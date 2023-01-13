import os, shutil, subprocess, glob
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from tqdm import tqdm
from proteinshake.utils import load, save, unzip_file, write_avro

import itertools
from joblib import Parallel, delayed
from collections import defaultdict
from main import get_dataset, replace_avro_files, transfer_dataset, transfer_file
from tqdm import tqdm
from functools import partialmethod

# slow down tqdm
tqdm.__init__ = partialmethod(tqdm.__init__, mininterval=60*60) # once per hour

THRESHOLDS = [0.5, 0.6, 0.7, 0.8, 0.9]

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
        return -1.
    return float(TM1), float(TM2), float(RMSD)

def prepare(dataset):
    proteins = list(dataset.proteins()[0])
    path_dict = {dataset.get_id_from_filename(os.path.basename(f)):f for f in dataset.get_raw_files()}
    paths = [path_dict[p['protein']['ID']] for p in proteins]
    paths = [unzip_file(p, remove=False) if not os.path.exists(p.rstrip('.gz')) else p.rstrip('.gz') for p in tqdm(paths, desc='Unzipping')]
    pdbids = [dataset.get_id_from_filename(p) for p in paths]
    pairs = list(itertools.combinations(range(len(paths)), 2))
    todo = [f'{paths[p1]} {paths[p2]}' for p1, p2 in pairs]
    BATCH_SIZE = 1000
    for job_id, i in enumerate(range(0,len(todo), BATCH_SIZE)):
        with open(f'{dataset.root}/jobs/{job_id}.txt', 'w') as file:
            file.write('\n'.join(todo[i:i+BATCH_SIZE]))


def compute_clusters_structure(dataset):
    job_path = os.path.expandvars(f'{dataset.root}/jobs/$SLURM_ARRAY_TASK_ID')
    with open(job_path+'.txt', 'r') as file:
        todo = list(map(lambda x: x.split(), file.readlines()))
    tmscore = defaultdict(lambda: {})
    rmsd = defaultdict(lambda: {})
    for pdb1,pdb2 in tqdm(todo, desc='TMalign'):
        name1 = dataset.get_id_from_filename(pdb1)
        name2 = dataset.get_id_from_filename(pdb2)
        tm1, tm2, _rmsd = tmalign_wrapper(pdb1, pdb2)
        tmscore[name1][name2], tmscore[name2][name1] = tm1, tm2
        rmsd[name1][name2], rmsd[name2][name1] = _rmsd, _rmsd
    save(dict(tmscore), job_path+'.tmalign.json')
    save(dict(rmsd), job_path+'.rmsd.json')


def collect(dataset):
    tmscore = {}
    for path in glob.glob(f'{dataset.root}/jobs/*.tmalign.json'):
        tmscore = {**tmscore, **load(path)}
    save(tmscore, f'{dataset.root}/{dataset.name}.tmalign.json')
    transfer_file(f'{dataset.root}/{dataset.name}.tmalign.json', 'structure')

    rmsd = {}
    for path in glob.glob(f'{dataset.root}/jobs/*.rmsd.json'):
        rmsd = {**rmsd, **load(path)}
    save(rmsd, f'{dataset.root}/{dataset.name}.rmsd.json')
    transfer_file(f'{dataset.root}/{dataset.name}.rmsd.json', 'structure')

    proteins = list(dataset.proteins()[0])
    path_dict = {dataset.get_id_from_filename(os.path.basename(f)):f for f in dataset.get_raw_files()}
    paths = [path_dict[p['protein']['ID']] for p in proteins]
    pdbids = [dataset.get_id_from_filename(p) for p in paths]

    num_proteins = len(pdbids)
    DM = np.zeros((num_proteins, num_proteins))
    DM = []
    for i in range(num_proteins):
        for j in range(i+1, num_proteins):
            DM.append(
                1 - max(
                    tmscore[pdbids[i]][pdbids[j]],
                    tmscore[pdbids[j]][pdbids[i]]
                )
            )
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
