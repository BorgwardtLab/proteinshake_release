import os, itertools, tempfile, random, subprocess
from tqdm import tqdm
from util import replace_avro_files

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
    db_path = f'{ds.root}/raw/foldseek/foldseekDB'
    out_path = f'{ds.root}/raw/foldseek/'
    out_file = f'{out_path}/output.m8'
    try:
        cmd = ['foldseek', 'easy-search', query, db_path, out_file, out_path,
            '--threads', str(ds.n_jobs),
            '--max-seqs', '1000000000',
            '--lddt-threshold', str(threshold),
            '--format-output', 'target'
        ]
        out = subprocess.run(cmd, capture_output=True, text=True)
        with open(out_file, 'r') as file:
            return file.read().split()
    except Exception as e:
        return []

def split(ds, pool, testsize, threshold, n=10, seed=42):
    random.seed(seed)
    test = []
    i = 0
    while len(test) < testsize:
        i += 1
        print(f'\rSampling cluster\t{i}\t({len(test)}/{testsize})', end='')
        query = random.choice(pool)
        cluster = foldseek_wrapper(ds, query, threshold)
        cluster = list(set([c.rstrip('_A') for c in cluster])) # remove chain ID
        if len(cluster) < n: continue
        cluster = random.sample(cluster, n)
        test.extend(cluster)
        pool = [p for p in pool if not p in cluster]
    print()
    return pool, test

def get_paths(dataset):
    pdbids = [p['protein']['ID'] for p in dataset.proteins()]
    path_dict = {dataset.get_id_from_filename(os.path.basename(f)):f for f in dataset.get_raw_files()}
    paths = [path_dict[id] for id in pdbids]
    return paths

def compute_clusters_structure(dataset, thresholds=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], testsize=0.1, valsize=0.1):

    print(f'Structure clustering {dataset.name}')
    proteins = list(dataset.proteins())
    paths = get_paths(dataset)
    foldseek_create_database(dataset)

    for threshold in thresholds:
        pool = [p for p in paths]
        n_test, n_val = int(len(pool)*testsize), int(len(pool)*valsize)
        pool, test = split(dataset, pool, n_test, threshold)
        train, val = split(dataset, pool, n_val, threshold)
        train, test, val = [dataset.get_id_from_filename(p) for p in train], [dataset.get_id_from_filename(p) for p in test], [dataset.get_id_from_filename(p) for p in val]
        for p in proteins:
            p['protein'][f'split_{threshold}'] = 'test' if p['protein']['ID'] in test else 'train'
            p['protein'][f'split_{threshold}'] = 'val' if p['protein']['ID'] in val else 'train'
    replace_avro_files(dataset, proteins)







    



