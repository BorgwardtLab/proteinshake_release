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

def transfer_file(file, dest):
    if os.path.exists(file):
        if not file.endswith('.gz'):
            zip_file(file)
            file += '.gz'
        subprocess.call(['rsync', f'{file}', dest])

def transfer_dataset(ds, dest, remove=False):
    transfer_file(f'{ds.root}/{ds.name}.atom.avro', dest)
    transfer_file(f'{ds.root}/{ds.name}.residue.avro', dest)
    if remove:
        shutil.rmtree(ds.root)

def get_dataset(root, name, organism=None, n_jobs=1):
    if name.endswith('Task'):
        Task = getattr(TASKS, name)
        root = f'{root}/{Task.DatasetClass.__name__}'
        return Task(root=root, use_precomputed=False)
    else:
        Dataset = getattr(DATASETS, name)
        root = f'{root}/{name}'
        Dataset.limit = 100
        if name == 'AlphaFoldDataset':
            return Dataset(root=root, organism=organism, use_precomputed=False, n_jobs=n_jobs)
        else:
            return Dataset(root=root, use_precomputed=False, n_jobs=n_jobs)

def get_paths(dataset):
    pdbids = [p['protein']['ID'] for p in dataset.proteins()]
    path_dict = {dataset.get_id_from_filename(os.path.basename(f)):f for f in dataset.get_raw_files()}
    paths = [path_dict[id] for id in pdbids]
    return pdbids, paths, path_dict

def split(wrapper, ds, pool, testsize, threshold, n=1, seed=42):
    random.seed(seed)
    test = []
    with tqdm(total=testsize, desc='Sampling split') as pbar:
        while len(test) < testsize:
            query = random.choice(pool)
            cluster = wrapper(ds, query, threshold)
            if len(cluster) < n: continue
            cluster = random.sample(cluster, n)
            test.extend(cluster)
            pool = [p for p in pool if not p in cluster]
            pbar.update(len(cluster))
    return pool, test