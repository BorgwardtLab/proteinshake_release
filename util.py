import importlib, argparse, os, shutil, subprocess, random
from tqdm import tqdm
from functools import partialmethod
from proteinshake.utils import write_avro, zip_file, unzip_file

DATASETS = importlib.import_module('proteinshake.datasets')
TASKS = importlib.import_module('proteinshake.tasks')

def replace_avro_files(dataset, proteins):
    residue_proteins = list(dataset.proteins(resolution='residue'))
    atom_proteins = list(dataset.proteins(resolution='atom'))
    for r,a,p in zip(residue_proteins,atom_proteins,proteins):
        r['protein'] = p['protein']
        a['protein'] = p['protein']
    write_avro(list(residue_proteins), f'{dataset.root}/{dataset.name}.residue.avro')
    write_avro(list(atom_proteins), f'{dataset.root}/{dataset.name}.atom.avro')

def get_dataset(root, name, organism=None, n_jobs=1):
    Dataset = getattr(DATASETS, name)
    #Dataset.limit = 10
    if name == 'AlphaFoldDataset':
        return Dataset(root=f'{root}/{name}_{organism}', organism=organism, use_precomputed=False, n_jobs=n_jobs, skip_signature_check=True)
    else:
        return Dataset(root=f'{root}/{name}', use_precomputed=False, n_jobs=n_jobs, skip_signature_check=True)

def get_task(root, name, n_jobs=1):
    Task = getattr(TASKS, name)
    return Task(root=f'{root}/{Task.DatasetClass.__name__}', use_precomputed=False, n_jobs=n_jobs, skip_signature_check=True)

def get_paths(dataset):
    pdbids = [p['protein']['ID'] for p in dataset.proteins()]
    path_dict = {dataset.get_id_from_filename(os.path.basename(f)):f for f in dataset.get_raw_files()}
    paths = [path_dict[id] for id in pdbids]
    return pdbids, paths, path_dict

def split(wrapper, ds, pool, test_size, threshold, path_dict, pdbids, n=0.05, seed=42, verbose=False):
    random.seed(seed)
    test = set()
    n = int(test_size*n)
    with tqdm(total=test_size, desc='Sampling split') as pbar:
        while len(test) < test_size:
            query = random.choice(pool)
            cluster = wrapper(ds, query, threshold, path_dict)
            cluster = [c for c in cluster if c in pdbids]
            #if len(cluster) < n: continue
            pool = [p for p in pool if not p in cluster]
            if len(cluster) > n: cluster = random.sample(cluster, n)
            if len(cluster) > (test_size-len(test)): cluster = random.sample(cluster, (test_size-len(test)))
            test.update(cluster)
            pbar.update(len(cluster))
    return pool, test