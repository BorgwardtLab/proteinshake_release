import os, tarfile, shutil
from datetime import datetime
from proteinshake.datasets import __all__ as DATASETS
from proteinshake.tasks import __all__ as TASKS
from proteinshake.datasets.alphafold import AF_DATASET_NAMES
from util import get_dataset, transfer_dataset, transfer_file
from sequence_split import compute_sequence_split
from structure_split import compute_structure_split

AF_DATASET_NAMES = ['methanocaldococcus_jannaschii']

######## VARIABLES #########

TAG = 'test' #datetime.now().strftime('%d%b%Y').upper() # release tag, set to the current date
SCRATCH = os.path.expandvars(f'$LOCAL_SCRATCH/proteinshake/{TAG}') # local scratch directory to avoid IO bottlenecks during processing
DESTINATION = os.path.expandvars(f'$GLOBAL_SCRATCH/Datasets/proteinshake/{TAG}/') # final release directory to save the result files
NJOBS = 20 # number of jobs

###########################

os.makedirs(f'{SCRATCH}/release', exist_ok=True)
os.makedirs(DESTINATION, exist_ok=True)

# construct dataset iterator
TASK_DATASETS = [d for d in DATASETS if not d in ['Dataset','AlphaFoldDataset','RCSBDataset']] # filter unlabeled datasets
DATASETS = [d for d in DATASETS if not d in ['Dataset','AlphaFoldDataset']] # filter parent class and AF
ALL_DATASETS = list(zip(DATASETS,[None]*len(DATASETS))) + list(zip(['AlphaFoldDataset']*len(AF_DATASET_NAMES), AF_DATASET_NAMES)) # zip with organism name

TASK_DATASETS = ['GeneOntologyDataset']
ALL_DATASETS = [('GeneOntologyDataset',None)]

# download data
for name, organism in ALL_DATASETS:
    ds = get_dataset(SCRATCH, name, organism, NJOBS)
print('Downloaded all datasets.')

# random splitting
for name in TASK_DATASETS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
    compute_random_split(ds)
print('Random split ready.')

# sequence splitting
for name in TASK_DATASETS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
    compute_sequence_split(ds)
print('Sequence split ready.')

# structure splitting
for name in TASK_DATASETS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
    compute_structure_split(ds)
print('Structure split ready.')

# collecting release
print('Collecting...')
for name, organism in ALL_DATASETS:
    ds = get_dataset(SCRATCH, name, organism, NJOBS)
    shutil.copyfile(f'{ds.root}/{ds.name}.atom.avro', f'{SCRATCH}/release')
    shutil.copyfile(f'{ds.root}/{ds.name}.residue.avro', f'{SCRATCH}/release')
    for filename in ds.additional_files:
        if not filename.endswith('.gz'):
            zip_file(filename)
            filename += '.gz'
        shutil.copyfile(f'{ds.root}/{filename}', f'{SCRATCH}/release')
with tarfile.open(f'{SCRATCH}/release.tar', 'w') as tar:
    tar.add(f'{SCRATCH}/release', arcname=f'ProteinShake_Release_{TAG}')
subprocess.call(['rsync', f'{SCRATCH}/release.tar', DESTINATION])

    