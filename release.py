import os, tarfile
from datetime import datetime
from proteinshake.datasets import __all__ as DATASETS
from proteinshake.tasks import __all__ as TASKS
from proteinshake.datasets.alphafold import AF_DATASET_NAMES
from util import get_dataset, transfer_dataset, transfer_file
from cluster_sequence import compute_clusters_sequence
from cluster_structure import compute_clusters_structure

AF_DATASET_NAMES = ['methanocaldococcus_jannaschii']

######## VARIABLES #########

TAG = '10FEB2023' #datetime.now().strftime('%d%b%Y').upper() # release tag, set to the current date
SCRATCH = os.path.expandvars(f'$LOCAL_SCRATCH/proteinshake/{TAG}') # local scratch directory to avoid IO bottlenecks during processing
RELEASE = os.path.expandvars(f'$GLOBAL_SCRATCH/Datasets/proteinshake/{TAG}/') # final release directory to save the result files
NJOBS = 130 # number of jobs

###########################

os.makedirs(RELEASE, exist_ok=True)

# construct dataset iterator
TASK_DATASETS = [d for d in DATASETS if not d in ['Dataset','AlphaFoldDataset','RCSBDataset']] # filter unlabeled datasets
DATASETS = [d for d in DATASETS if not d in ['Dataset','AlphaFoldDataset']] # filter parent class and AF
ALL_DATASETS = list(zip(DATASETS,[None]*len(DATASETS))) + list(zip(['AlphaFoldDataset']*len(AF_DATASET_NAMES), AF_DATASET_NAMES)) # zip with organism name
TASKS = [t for t in TASKS if t != 'Task']

TASK_DATASETS = ['EnzymeCommissionDataset', 'ProteinProteinInterfaceDataset', 'ProteinLigandInterfaceDataset', 'TMAlignDataset', 'SCOPDataset', 'PfamDataset', 'GeneOntologyDataset'] # fix order for now

# download data
for name, organism in ALL_DATASETS:
    ds = get_dataset(SCRATCH, name, organism, NJOBS)
print('Downloaded all datasets.')

# sequence clustering
for name in TASK_DATASETS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
    compute_clusters_sequence(ds)
print('Clustered all sequences.')

# structure clustering
for name in TASK_DATASETS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
    compute_clusters_structure(ds)
print('Clustered all structures.')

# compute tasks
for name in TASKS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
print('Computed all tasks.')

# transfer
print('Transferring...')
for name, organism in ALL_DATASETS:
    ds = get_dataset(SCRATCH, name, organism, NJOBS)
    transfer_dataset(ds, RELEASE)
for name in TASK_DATASETS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
    transfer_file(f'{ds.root}/{ds.name}.cdhit.json', RELEASE)
    transfer_file(f'{ds.root}/{ds.name}.tm.npy', RELEASE)
    transfer_file(f'{ds.root}/{ds.name}.rmsd.npy', RELEASE)
    transfer_file(f'{ds.root}/{ds.name}.gdt.npy', RELEASE)
for name in TASKS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
    transfer_file(f'{ds.root}/{ds.name}.json', RELEASE)

print('Archiving...')
with tarfile.open(f'{RELEASE}.tar', 'w') as tar:
    tar.add(RELEASE, arcname=os.path.basename(RELEASE))
    