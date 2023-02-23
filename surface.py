import os
import tempfile
import shutil
import subprocess
import pandas as pd
from joblib import Parallel, delayed
from proteinshake.utils import write_avro, unzip_file
from tqdm import tqdm
from util import get_paths

def dms_wrapper(pdb_path, d=0.5):
    """ Call DMS to compute a surface for the PDB.

    Usage: dms input_file [-a] [-d density] [-g file] [-i file] [-n] [-w radius] [-v] -o file
    -a	use all atoms, not just amino acids
    -d	change density of points
    -g	send messages to file
    -i	calculate only surface for specified atoms
    -n	calculate normals for surface points
    -w	change probe radius
    -v	verbose
    -o	specify output file name (required)

    See: https://www.cgl.ucsf.edu/Overview/software.html#dms
    """
    with tempfile.TemporaryDirectory() as tmp:
        dest = f'{tmp}/out.surf'
        assert shutil.which('dms') is not None, "DMS executable not in PATH go here to install https://www.cgl.ucsf.edu/Overview/software.html#dms."
        cmd = ['dms', pdb_path, '-n', '-d', str(d), '-o', dest]
        out = subprocess.run(cmd, capture_output=True, text=True)
        print(out.stdout)
        print(out.stderr)
        df = pd.read_csv(dest,
            delim_whitespace = True,
            header = None,
            names = [
                'residue_name',
                'residue_index',
                'atom_name',
                'x',
                'y',
                'z',
                'point_type',
                'area',
                'x_norm',
                'y_norm',
                'z_norm'
            ]
        )
        df = df.dropna(axis=0)
        return df['residue_index'].tolist()

    

def compute_surface(dataset):

    print(f'Computing surfaces {dataset.name}...')
    
    proteins_residue = list(dataset.proteins(resolution='residue'))
    proteins_atom = list(dataset.proteins(resolution='atom'))
    pdbids, paths, path_dict = get_paths(dataset)
    paths = [unzip_file(p, remove=False) if not os.path.exists(p.rstrip('.gz')) else p.rstrip('.gz') for p in tqdm(paths, desc='Unzipping')]
    #surfaces = Parallel(n_jobs=dataset.n_jobs)(delayed(dms_wrapper)(p) for p in paths)
    surfaces = [dms_wrapper(p) for p in paths]
    for pr, pa, s in zip(proteins_residue, proteins_atom, surfaces):
        pr['residue']['surface'] = [(i in s) for i in pr['residue']['residue_number']]
        pa['atom']['surface'] = [(i in s) for i in pa['atom']['residue_number']]
    write_avro(list(proteins_residue), f'{dataset.root}/{dataset.name}.residue.avro')
    write_avro(list(proteins_atom), f'{dataset.root}/{dataset.name}.atom.avro')