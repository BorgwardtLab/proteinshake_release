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
from summary import summary

#AF_DATASET_NAMES = ['methanocaldococcus_jannaschii']

######## VARIABLES #########

TAG = '09MAR2023' #datetime.now().strftime('%d%b%Y').upper() # release tag, set to the current date
SCRATCH = os.path.expandvars(f'$LOCAL_SCRATCH/proteinshake/{TAG}') # local scratch directory to avoid IO bottlenecks during processing
DESTINATION = os.path.expandvars(f'$GLOBAL_SCRATCH/Datasets/proteinshake/{TAG}/') # final release directory to save the result files
NJOBS = 100 # number of jobs

###########################

os.makedirs(f'{SCRATCH}/release', exist_ok=True)
os.makedirs(DESTINATION, exist_ok=True)

# construct dataset iterator
TASK_DATASETS = [d for d in DATASETS if not d in ['Dataset','AlphaFoldDataset','RCSBDataset','ProteinLigandDecoysDataset']] # filter unlabeled datasets
DATASETS = [d for d in DATASETS if not d in ['Dataset','AlphaFoldDataset']] # filter parent class and AF
ALL_DATASETS = list(zip(DATASETS,[None]*len(DATASETS))) + list(zip(['AlphaFoldDataset']*len(AF_DATASET_NAMES), AF_DATASET_NAMES)) # zip with organism name

'''
# download data
for name, organism in ALL_DATASETS:
    print(f'Downloading {name} {organism}')
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
'''
# structure splitting
for name in TASK_DATASETS:
    ds = get_dataset(SCRATCH, name, None, NJOBS)
    compute_structure_split(ds)
print('Structure split ready.')

# summaries
df = []
for name, organism in ALL_DATASETS:
    print(f'Summarizing {name} {organism}')
    ds = get_dataset(SCRATCH, name, organism, NJOBS)
    df.append(summary(ds))
df = pd.DataFrame(df)
print(df)
df.to_csv(f'{SCRATCH}/release/summary.csv', index=False)
zip_file(f'{SCRATCH}/release/summary.csv')
os.remove(f'{SCRATCH}/release/summary.csv')

# collecting release
print('Collecting...')
for name, organism in ALL_DATASETS:
    ds = get_dataset(SCRATCH, name, organism, NJOBS)
    for filename in [f'{ds.name}.atom.avro', f'{ds.name}.residue.avro'] + ds.additional_files:
        if not filename.endswith('.gz'):
            zip_file(f'{ds.root}/{filename}')
            filename += '.gz'
        shutil.copyfile(f'{ds.root}/{filename}', f'{SCRATCH}/release/{filename}')
with tarfile.open(f'{SCRATCH}/release.tar', 'w') as tar:
    tar.add(f'{SCRATCH}/release', arcname=f'ProteinShake_Release_{TAG}')
subprocess.call(['rsync', f'{SCRATCH}/release.tar', DESTINATION])

    