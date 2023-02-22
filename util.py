import importlib, argparse, os, shutil, subprocess
from tqdm import tqdm
from functools import partialmethod
from proteinshake.utils import write_avro, zip_file

DATASETS = importlib.import_module('proteinshake.datasets')
TASKS = importlib.import_module('proteinshake.tasks')

def replace_avro_files(dataset, proteins):
    residue_proteins = list(dataset.proteins(resolution='residue')[0])
    atom_proteins = list(dataset.proteins(resolution='atom')[0])
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
        #Dataset.limit = 100
        if name == 'AlphaFoldDataset':
            return Dataset(root=root, organism=organism, use_precomputed=False, n_jobs=n_jobs)
        else:
            return Dataset(root=root, use_precomputed=False, n_jobs=n_jobs)