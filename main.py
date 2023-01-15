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
        zip_file(file)
        subprocess.call(['rsync', f'{file}.gz', os.path.expandvars(f'$SHAKE_STORE/{folder}/')])

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
    args = parser.parse_args()

    DATASET, ORGANISM, n_jobs = args.dataset, args.organism, args.njobs
    NAME = f'{DATASET}_{ORGANISM}' if DATASET == 'AlphaFoldDataset' else DATASET
    ROOT = os.path.expandvars(f'$SHAKE_SCRATCH/{NAME}')

    Dataset = getattr(datasets, DATASET)
    #Dataset.limit = 100
    if DATASET == 'AlphaFoldDataset':
        return Dataset(root=ROOT, organism=ORGANISM, use_precomputed=False, n_jobs=n_jobs), args
    elif DATASET in ['RCSBDataset','GeneOntologyDataset','EnzymeCommissionDataset','SCOPDataset','PfamDataset']:
        return Dataset(root=ROOT, use_precomputed=False, n_jobs=n_jobs, max_requests=5), args
    else:
        return Dataset(root=ROOT, use_precomputed=False, n_jobs=n_jobs), args

if __name__ == '__main__':
    ds, args = get_dataset()
    transfer_dataset(ds, 'parsed')
