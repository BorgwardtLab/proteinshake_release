import os, tarfile, shutil, subprocess
import pandas as pd
from datetime import datetime
from proteinshake.datasets import __all__ as DATASETS
from proteinshake.tasks import __all__ as TASKS
from proteinshake.datasets.alphafold import AF_DATASET_NAMES
from proteinshake.utils import zip_file
from util import get_dataset
from random_split import compute_random_split
from sequence_split import compute_sequence_split
from structure_split import compute_structure_split

######## VARIABLES #########

TAG = '02JUN2023' #datetime.now().strftime('%d%b%Y').upper() # release tag, set to the current date
ROOT = os.path.expandvars(f'$LOCAL_SCRATCH/proteinshake/{TAG}') # root directory to save the datasets to
NJOBS = 50 # number of jobs

###########################

os.makedirs(f'{ROOT}/release', exist_ok=True)

# construct dataset iterator
TASK_DATASETS = [d for d in DATASETS if not d in ['Dataset','AlphaFoldDataset','RCSBDataset','ProteinLigandDecoysDataset']] # filter unlabeled datasets
DATASETS = [d for d in DATASETS if not d in ['Dataset','AlphaFoldDataset']] # filter parent class and AF
ALL_DATASETS = list(zip(DATASETS,[None]*len(DATASETS))) + list(zip(['AlphaFoldDataset']*len(AF_DATASET_NAMES), AF_DATASET_NAMES)) # zip with organism name

TASK_DATASETS = ['ProteinProteinInterfaceDataset']
DATASETS = ['ProteinProteinInterfaceDataset']
ALL_DATASETS = [('ProteinProteinInterfaceDataset',None)]

# download data
for name, organism in ALL_DATASETS:
    print(f'Downloading {name} {organism}')
    ds = get_dataset(ROOT, name, organism, NJOBS)
print('Downloaded all datasets.')

# random splitting
for name in TASK_DATASETS:
    ds = get_dataset(ROOT, name, None, NJOBS)
    compute_random_split(ds)
print('Random split ready.')

# sequence splitting
for name in TASK_DATASETS:
    ds = get_dataset(ROOT, name, None, NJOBS)
    compute_sequence_split(ds)
print('Sequence split ready.')

# structure splitting
for name in TASK_DATASETS:
    ds = get_dataset(ROOT, name, None, NJOBS)
    compute_structure_split(ds)
print('Structure split ready.')

# collecting release
print('Collecting...')
for name, organism in ALL_DATASETS:
    ds = get_dataset(ROOT, name, organism, NJOBS)
    for filename in [f'{ds.name}.atom.avro', f'{ds.name}.residue.avro'] + ds.additional_files:
        if not filename.endswith('.gz'):
            zip_file(f'{ds.root}/{filename}')
            filename += '.gz'
        shutil.copyfile(f'{ds.root}/{filename}', f'{ROOT}/release/{filename}')

    