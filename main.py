from proteinshake.utils import write_avro
import importlib, argparse, os, shutil, subprocess
from proteinshake.utils import zip_file
from tqdm import tqdm
from functools import partialmethod

# slow down tqdm
tqdm.__init__ = partialmethod(tqdm.__init__, mininterval=60*60) # once per hour

def replace_avro_files(dataset, proteins):
    residue_proteins = list(dataset.proteins(resolution='residue')[0])
    atom_proteins = list(dataset.proteins(resolution='atom')[0])
    for r,a,p in zip(residue_proteins,atom_proteins,proteins):
        r['protein'] = p['protein']
        a['protein'] = p['protein']
    write_avro(list(residue_proteins), f'{dataset.root}/{dataset.name}.residue.avro')
    write_avro(list(atom_proteins), f'{dataset.root}/{dataset.name}.atom.avro')

def transfer_file(file, folder):
    if os.path.exists(file):
        if not file.endswith('.gz'):
            zip_file(file)
            file += '.gz'
        subprocess.call(['rsync', f'{file}', os.path.expandvars(f'$SHAKE_STORE/{folder}/')])

def transfer_dataset(ds, folder):
    transfer_file(f'{ds.root}/{ds.name}.atom.avro', folder)
    transfer_file(f'{ds.root}/{ds.name}.residue.avro', folder)

def get_dataset():
    datasets = importlib.import_module('proteinshake.datasets')

    parser = argparse.ArgumentParser(description='Script to generate all datasets for release.')
    parser.add_argument('--njobs', type=int, help='Number of jobs.', default=10)
    parser.add_argument('--dataset', type=str, help='Name of the dataset class (case sensitive)', default='RCSBDataset')
    parser.add_argument('--organism', type=str, help='Organism (for AlphaFold datasets)', default='swissprot')
    parser.add_argument('--prepare', action='store_true')
    parser.add_argument('--compute', action='store_true')
    parser.add_argument('--collect', action='store_true')
    parser.add_argument('--task',  type=str, help='Task', default='none')
    args = parser.parse_args()

    DATASET, ORGANISM, TASK, n_jobs = args.dataset, args.organism, args.task, args.njobs
    NAME = f'{DATASET}_{ORGANISM}' if DATASET == 'AlphaFoldDataset' else DATASET
    ROOT = os.path.expandvars(f'$SHAKE_SCRATCH/{NAME}')

    Dataset = getattr(datasets, DATASET)
    #Dataset.limit = 100
    if DATASET == 'AlphaFoldDataset':
        ds = Dataset(root=ROOT, organism=ORGANISM, use_precomputed=False, n_jobs=n_jobs)
    elif DATASET in ['RCSBDataset','GeneOntologyDataset','EnzymeCommissionDataset','SCOPDataset','PfamDataset']:
        ds = Dataset(root=ROOT, use_precomputed=False, n_jobs=n_jobs, max_requests=5)
    else:
        ds = Dataset(root=ROOT, use_precomputed=False, n_jobs=n_jobs)

    if TASK != 'none':
        tasks = importlib.import_module('proteinshake.tasks')
        Task = getattr(tasks, TASK)
        task = Task(root=ROOT, use_precomputed=False)

    return ds, args

if __name__ == '__main__':
    ds, args = get_dataset()
    if args.task != 'none':
        transfer_file(f'{ds.root}/{args.task}.json.gz', 'tasks')
    else:
        transfer_dataset(ds, 'parsed')
