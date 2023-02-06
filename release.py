import os
from proteinshake.datasets import __all__ as DATASETS
from proteinshake.tasks import __all__ as TASKS
from proteinshake.datasets.alphafold import AF_DATASET_NAMES
from util import get_dataset, transfer_dataset

# construct dataset iterator
TASK_DATASETS = [d for d in DATASETS if not d in ['Dataset','AlphaFoldDataset','RCSBDataset']] # filter unlabeled datasets
DATASETS = [d for d in DATASETS if not d in ['Dataset','AlphaFoldDataset']] # filter parent class and AF
ALL_DATASETS = list(zip(DATASETS,[None]*len(DATASETS))) + list(zip(['AlphaFoldDataset']*len(AF_DATASET_NAMES), AF_DATASET_NAMES)) # zip with organism name

# variables
SCRATCH = os.expandvars('')
RELEASE = os.expandvars('')
NJOBS = 130

# download data
for name, organism in ALL_DATASETS:
    ds = get_dataset(SCRATCH, name, organism, NJOBS)
print('Downloaded all datasets.')

# sequence clustering
for name in TASK_DATASETS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
    compute_clusters_sequence(ds, NJOBS)
print('Clustered all sequences.')

# structure clustering
for name in TASK_DATASETS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
    compute_clusters_structure(ds, NJOBS)
print('Clustered all structures.')

# compute tasks
for name in TASKS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
print('Computed all tasks.')

# transfer
for name, organism in ALL_DATASETS:
    ds = get_dataset(SCRATCH, name, organism, NJOBS)
    transfer_dataset(ds, RELEASE)

for name in TASK_DATASETS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
    transfer_file(f'{ds.root}/{ds.name}.cdhit.json', RELEASE)
    transfer_file(f'{ds.root}/{ds.name}.tm.npy', RELEASE)
    transfer_file(f'{ds.root}/{ds.name}.rmsd.npy', RELEASE)
    transfer_file(f'{ds.root}/{ds.name}.gdt.npy', RELEASE)
    