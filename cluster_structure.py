import os, shutil, subprocess, glob, gzip
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from tqdm import tqdm
from proteinshake.utils import load, save, unzip_file, write_avro

import itertools, time
from joblib import Parallel, delayed
from collections import defaultdict
from main import get_dataset, replace_avro_files, transfer_dataset, transfer_file
from tqdm import tqdm
from functools import partialmethod

# slow down tqdm
tqdm.__init__ = partialmethod(tqdm.__init__, mininterval=60*60) # once per hour

THRESHOLDS = [0.5, 0.6, 0.7, 0.8, 0.9]
BATCH_SIZE = int(os.getenv('BATCHSIZE'))

def tmalign_wrapper(pdb1, pdb2):
    """Compute TM score with TMalign between two PDB structures.
    Parameters
    ----------
    pdb1: str
        Path to PDB.
    arg2 : str
        Path to PDB.
    Returns
    -------
    float
        TM score from `pdb1` to `pdb2`
    float
        TM score from `pdb2` to `pdb1`
    float
        RMSD between structures
    """
    assert shutil.which('TMalign') is not None,\
           "No TMalign installation found. Go here to install : https://zhanggroup.org/TM-align/TMalign.cpp"
    try:
        out = subprocess.run(['TMalign','-outfmt','2', pdb1, pdb2], stdout=subprocess.PIPE).stdout.decode()
        path1, path2, TM1, TM2, RMSD, ID1, ID2, IDali, L1, L2, Lali = out.split('\n')[1].split('\t')
    except Exception as e:
        print(e)
        return 0
    return float(TM1), float(TM2), float(RMSD)

def prepare(dataset):
    proteins = list(dataset.proteins()[0])
    path_dict = {dataset.get_id_from_filename(os.path.basename(f)):f for f in dataset.get_raw_files()}
    paths = [path_dict[p['protein']['ID']] for p in proteins]
    paths = [unzip_file(p, remove=False) if not os.path.exists(p.rstrip('.gz')) else p.rstrip('.gz') for p in tqdm(paths, desc='Unzipping')]
    pairs = itertools.combinations(range(len(paths)), 2)
    todo = (f'{paths[p1]} {paths[p2]}' for p1, p2 in pairs)
    job_id = 0
    while True:
        chunk_it = itertools.islice(todo, BATCH_SIZE)
        try:
            first_el = next(chunk_it)
        except StopIteration:
            break
        batch = itertools.chain((first_el,), chunk_it)
        with open(f'{dataset.root}/jobs/{job_id}.txt', 'w') as file:
            file.write('\n'.join(batch))
        job_id += 1
    print(f'Prepared {job_id} jobs.')


def compute_clusters_structure(dataset):
    job_path = os.path.expandvars(f'{dataset.root}/jobs/$SLURM_ARRAY_TASK_ID')
    #if os.path.exists(job_path+'.tmalign.json') and os.path.exists(job_path+'.tmalign.json'):
        #return
    with open(job_path+'.txt', 'r') as file:
        todo = list(map(lambda x: x.split(), file.readlines()))
    tmscore = defaultdict(lambda: {})
    rmsd = defaultdict(lambda: {})
    start = time.time()
    for pdb1,pdb2 in todo:
        name1 = dataset.get_id_from_filename(os.path.basename(pdb1))
        name2 = dataset.get_id_from_filename(os.path.basename(pdb2))
        tm1, tm2, _rmsd = tmalign_wrapper(pdb1, pdb2)
        tmscore[name1][name2], tmscore[name2][name1] = tm1, tm2
        rmsd[name1][name2], rmsd[name2][name1] = _rmsd, _rmsd
    save(dict(tmscore), job_path+'.tmalign.json')
    save(dict(rmsd), job_path+'.rmsd.json')
    job_id = os.path.basename(job_path)
    duration = int((time.time()-start)/60)
    print(f'Computed job {job_id} in {duration} min.')


def collect(dataset):

    proteins = list(dataset.proteins()[0])
    pdbids = [p['protein']['ID'] for p in proteins]
    path_dict = {dataset.get_id_from_filename(os.path.basename(f)):f for f in dataset.get_raw_files()}
    paths = [path_dict[id] for id in pdbids]
    num_proteins = len(pdbids)

    def get_matrix(file, default):
        d = defaultdict(dict)
        for path in glob.glob(f'{dataset.root}/jobs/*.{file}.json'):
            for k,v in load(path).items():
                d[k] = {**d[k], **v}
        mat = np.zeros((num_proteins, num_proteins))
        for i,id1 in enumerate(pdbids):
            for j,id2 in enumerate(pdbids):
                mat[i][j] = d[id1][id2] if i != j else default
        return mat

    tmscore = get_matrix('tmalign', default=1)
    with gzip.open(f'{dataset.root}/{dataset.name}.tmscore.npy.gz', 'w') as file:
        np.save(file=file, arr=tmscore)
    transfer_file(f'{dataset.root}/{dataset.name}.tmscore.npy.gz', 'structure')

    rmsd = get_matrix('rmsd', default=0)
    with gzip.open(f'{dataset.root}/{dataset.name}.rmsd.npy.gz', 'w') as file:
        np.save(file=file, arr=rmsd)
    transfer_file(f'{dataset.root}/{dataset.name}.rmsd.npy.gz', 'structure')

    DM = np.zeros((num_proteins, num_proteins))
    DM = []
    for i in range(num_proteins):
        for j in range(i+1, num_proteins):
            DM.append(1 - max(tmscore[i][j], tmscore[j][i]))
    DM = np.array(DM).reshape(-1, 1)

    for d in THRESHOLDS:
        clusterer = AgglomerativeClustering(n_clusters=None, distance_threshold=(1-d))
        clusterer.fit(DM)
        for i, p in enumerate(proteins):
            p['protein'][f'structure_cluster_{d}'] = int(clusterer.labels_[i])

    replace_avro_files(dataset, proteins)
    transfer_dataset(dataset, 'structure')
    print('All data collected.')



if __name__ == '__main__':
    ds, args = get_dataset()
    if args.prepare:
        os.makedirs(f'{ds.root}/jobs', exist_ok=True)
        prepare(ds)
    elif args.compute:
        compute_clusters_structure(ds)
    elif args.collect:
        collect(ds)
