from proteinshake.utils import write_avro
import importlib, argparse, os, shutil, subprocess
from proteinshake.utils import zip_file

def replace_avro_files(dataset, proteins):
    residue_proteins = list(dataset.proteins(resolution='residue')[0])
    atom_proteins = list(dataset.proteins(resolution='atom')[0])
    for r,a,p in zip(residue_proteins,atom_proteins,proteins):
        r['protein'] = p['protein']
        a['protein'] = p['protein']
    write_avro(list(residue_proteins), f'{dataset.root}/{dataset.name}.residue.avro')
    write_avro(list(atom_proteins), f'{dataset.root}/{dataset.name}.atom.avro')

def transfer_dataset(ds, folder):
    if os.path.exists(f'{ds.root}/{ds.name}.atom.avro') and os.path.exists(f'{ds.root}/{ds.name}.residue.avro'):
        zip_file(f'{ds.root}/{ds.name}.atom.avro')
        zip_file(f'{ds.root}/{ds.name}.residue.avro')
        subprocess.call(['rsync', f'{ds.root}/{ds.name}.atom.avro.gz', os.path.expandvars(f'$RELEASE_DIR/{folder}')])
        subprocess.call(['rsync', f'{ds.root}/{ds.name}.residue.avro.gz', os.path.expandvars(f'$RELEASE_DIR/{folder}')])

def get_dataset():
    datasets = importlib.import_module('proteinshake.datasets')

    parser = argparse.ArgumentParser(description='Script to generate all datasets for release.')
    parser.add_argument('--njobs', type=int, help='Number of jobs.', default=10)
    parser.add_argument('--dataset', type=str, help='Name of the dataset class (case sensitive)', default='RCSBDataset')
    parser.add_argument('--organism', type=str, help='Organism (for AlphaFold datasets)', default='swissprot')
    parser.add_argument('--download', action='store_true', help='Stops after download')
    args = parser.parse_args()

    DATASET, ORGANISM, n_jobs = args.dataset, args.organism, args.njobs
    NAME = f'{DATASET}_{ORGANISM}' if DATASET == 'AlphaFoldDataset' else DATASET
    ROOT = os.path.expandvars(f'$TMPDIR/proteinshake/$TAG/{NAME}')
    print(args)
    os.makedirs(ROOT, exist_ok=True)

    Dataset = getattr(datasets, DATASET)
    Dataset.limit = 10
    if args.download:
        ROOT = os.path.expandvars(f'$GLOBAL_SCRATCH/{NAME}')
        print(ROOT)
        Dataset.parse = lambda self: None
    else:
        if os.path.exists(ROOT):
            shutil.rmtree(ROOT)
        shutil.copytree(os.path.expandvars(f'$GLOBAL_SCRATCH/{NAME}'), ROOT)
    if DATASET == 'AlphaFoldDataset':
        return Dataset(root=ROOT, organism=ORGANISM, use_precomputed=False, n_jobs=n_jobs)
    else:
        return Dataset(root=ROOT, use_precomputed=False, n_jobs=n_jobs)

if __name__ == '__main__':
    ds = get_dataset()
    transfer_dataset(ds, 'parsed')
