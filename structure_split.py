import os, itertools, tempfile, random, subprocess, shutil
from tqdm import tqdm
from util import replace_avro_files, get_paths, split

def foldseek_create_database(ds):
    pdb_path = f'{ds.root}/raw/files/'
    out_path = f'{ds.root}/raw/foldseek/'
    db_path = f'{out_path}/foldseekDB'
    os.makedirs(out_path, exist_ok=True)
    cmd = ['foldseek', 'createdb', pdb_path, db_path]
    out = subprocess.run(cmd, capture_output=True, text=True)
    cmd = ['foldseek', 'createindex', db_path, out_path]
    out = subprocess.run(cmd, capture_output=True, text=True)

def foldseek_wrapper(ds, query, threshold):

    assert shutil.which('cd-hit') is not None,\
    "FoldSeek installation not found. Go here https://github.com/steineggerlab/foldseek to install"

    db_path = f'{ds.root}/raw/foldseek/foldseekDB'
    out_path = f'{ds.root}/raw/foldseek/'
    out_file = f'{out_path}/output.m8'
    n_jobs = 0 if ds.n_jobs < 0 else ds.n_jobs
    try:
        cmd = ['foldseek', 'easy-search', query, db_path, out_file, out_path,
            '--threads', str(n_jobs),
            '--max-seqs', '1000000000',
            '--lddt-threshold', str(threshold),
            '--format-output', 'target'
        ]
        out = subprocess.run(cmd, capture_output=True, text=True)
        with open(out_file, 'r') as file:
            cluster = file.read().split()
            cluster = list(set([c.rstrip('_A') for c in cluster])) # remove chain ID
            return cluster
    except Exception as e:
        return []

def compute_structure_split(dataset, thresholds=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], testsize=0.1, valsize=0.1):

    print(f'Structure split {dataset.name}')
    
    proteins = list(dataset.proteins())
    pdbids, paths, path_dict = get_paths(dataset)
    foldseek_create_database(dataset)

    for threshold in thresholds:
        pool = [p for p in paths]
        n_test, n_val = int(len(pool)*testsize), int(len(pool)*valsize)
        pool, test = split(foldseek_wrapper, dataset, pool, n_test, threshold)
        train, val = split(foldseek_wrapper, dataset, pool, n_val, threshold)
        train, test, val = [dataset.get_id_from_filename(p) for p in train], [dataset.get_id_from_filename(p) for p in test], [dataset.get_id_from_filename(p) for p in val]
        for p in proteins:
            p['protein'][f'split_{threshold}'] = 'test' if p['protein']['ID'] in test else 'train'
            p['protein'][f'split_{threshold}'] = 'val' if p['protein']['ID'] in val else 'train'
    replace_avro_files(dataset, proteins)







    



