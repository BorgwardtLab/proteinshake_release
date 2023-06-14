
from proteinshake.datasets import __all__ as DATASETS
from proteinshake.tasks import __all__ as TASKS
from proteinshake.datasets.alphafold import AF_DATASET_NAMES
from util import get_dataset

TAG = '02JUN2023'
ROOT = os.path.expandvars(f'$LOCAL_SCRATCH/proteinshake/{TAG}')
NJOBS = 1

TASKS = [t for t in TASKS if not t in ['Task']]
DATASETS = [d for d in DATASETS if not d in ['Dataset']]

# Overview Table with all datasets and tasks. Dataset: Name, Annotation, number of instances. Task: Name, type, level, metric.
dataset_table = []
for name in DATASETS:
    ds = get_dataset(ROOT, name, 'swissprot', NJOBS)
    dataset_table.append({
        'Name': name,
        'Size': len(ds.proteins()),
        'Description': ds.__doc__.split('.')[0]
    })
dataset_table = pd.DataFrame(dataset_table).to_html()
with open(f'{ROOT}/release/summary.csv','w') as file:
    file.write(dataset_table)


# For each task: label distribution in train/test/val
