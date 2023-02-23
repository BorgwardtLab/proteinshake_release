import os, tarfile
from datetime import datetime
from proteinshake.datasets import __all__ as DATASETS
from proteinshake.tasks import __all__ as TASKS
from proteinshake.datasets.alphafold import AF_DATASET_NAMES
from util import get_dataset, transfer_dataset, transfer_file
from surface import compute_surface
from sequence_split import compute_sequence_split
from structure_split import compute_structure_split

AF_DATASET_NAMES = ['methanocaldococcus_jannaschii']

######## VARIABLES #########

TAG = 'test' #datetime.now().strftime('%d%b%Y').upper() # release tag, set to the current date
SCRATCH = os.path.expandvars(f'$LOCAL_SCRATCH/proteinshake/{TAG}') # local scratch directory to avoid IO bottlenecks during processing
RELEASE = os.path.expandvars(f'$GLOBAL_SCRATCH/Datasets/proteinshake/{TAG}/') # final release directory to save the result files
NJOBS = 20 # number of jobs

###########################

os.makedirs(RELEASE, exist_ok=True)

# construct dataset iterator
TASK_DATASETS = [d for d in DATASETS if not d in ['Dataset','AlphaFoldDataset','RCSBDataset']] # filter unlabeled datasets
DATASETS = [d for d in DATASETS if not d in ['Dataset','AlphaFoldDataset']] # filter parent class and AF
ALL_DATASETS = list(zip(DATASETS,[None]*len(DATASETS))) + list(zip(['AlphaFoldDataset']*len(AF_DATASET_NAMES), AF_DATASET_NAMES)) # zip with organism name
TASKS = [t for t in TASKS if t != 'Task']

TASK_DATASETS = ['GeneOntologyDataset']
TASKS = []
ALL_DATASETS = [('GeneOntologyDataset',None)]

# download data
for name, organism in ALL_DATASETS:
    ds = get_dataset(SCRATCH, name, organism, NJOBS)
print('Downloaded all datasets.')

# sequence splitting
for name in TASK_DATASETS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
    pass#compute_sequence_split(ds)
print('Clustered all sequences.')

# structure splitting
for name in TASK_DATASETS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
    pass#compute_structure_split(ds)
print('Clustered all structures.')

# compute surface
for name, organism in ALL_DATASETS:
    ds = get_dataset(SCRATCH, name, organism, NJOBS)
    compute_surface(ds)
print('Clustered all sequences.')

# compute tasks
for name in TASKS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
print('Computed all tasks.')

# transfer
print('Transferring...')
for name, organism in ALL_DATASETS:
    ds = get_dataset(SCRATCH, name, organism, NJOBS)
    transfer_dataset(ds, RELEASE)
for name in TASKS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
    transfer_file(f'{ds.root}/{ds.name}.json', RELEASE)

# archiving
print('Archiving...')
with tarfile.open(f'{RELEASE}.tar', 'w') as tar:
    tar.add(RELEASE, arcname=os.path.basename(RELEASE))
    