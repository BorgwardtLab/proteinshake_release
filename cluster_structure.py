import os, shutil, subprocess
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from tqdm import tqdm
from proteinshake.utils import save, unzip_file, write_avro

import itertools
from joblib import Parallel, delayed
from collections import defaultdict
from main import get_dataset, replace_avro_files, transfer_dataset

def tmalign_wrapper(pdb1, pdb2):
    """Compute TM score with TMalign between two PDB structures.
    Parameters
    ----------
    pdb1: str
        Path to PDB.
    arg2 : str
        Path to PDB.
    Returns
    -------
    float
        TM score from `pdb1` to `pdb2`
    float
        TM score from `pdb2` to `pdb1`
    float
        RMSD between structures
    """
    assert shutil.which('TMalign') is not None,\
           "No TMalign installation found. Go here to install : https://zhanggroup.org/TM-align/TMalign.cpp"
    try:
        out = subprocess.run(['TMalign','-outfmt','2', pdb1, pdb2], stdout=subprocess.PIPE).stdout.decode()
        path1, path2, TM1, TM2, RMSD, ID1, ID2, IDali, L1, L2, Lali = out.split('\n')[1].split('\t')
    except Exception as e:
        print(e)
        return -1.
    return float(TM1), float(TM2), float(RMSD)

def compute_clusters_structure(dataset):
    """ Launch TMalign on all pairs of proteins in dataset.
    Assign a cluster ID to each protein at protein-level key 'structure_cluster'.
    Saves TMalign output to `dataset.root/{Dataset.__class__}.tmalign.json.gz`
    Parameters:
    -----------
    paths: list
        List of paths to original pdb files (after filtering).
    """
    proteins = list(dataset.proteins()[0])
    path_dict = {dataset.get_id_from_filename(os.path.basename(f)):f for f in dataset.get_raw_files()}
    paths = [path_dict[p['protein']['ID']] for p in proteins]

    dump_name = f'{dataset.name}.tmalign.json'
    dump_path = os.path.join(dataset.root, dump_name)

    if dataset.n_jobs == 1:
        print('Computing the TM scores with use_precompute = False is very slow. Consider increasing n_jobs.')

    paths = [unzip_file(p, remove=False) if not os.path.exists(p.rstrip('.gz')) else p.rstrip('.gz') for p in tqdm(paths, desc='Unzipping')]

    pdbids = [dataset.get_id_from_filename(p) for p in paths]
    pairs = list(itertools.combinations(range(len(paths)), 2))
    todo = [(paths[p1], paths[p2]) for p1, p2 in pairs]

    dist = defaultdict(lambda: {})

    output = Parallel(n_jobs=dataset.n_jobs)(
        delayed(tmalign_wrapper)(*pair) for pair in tqdm(todo, desc='TMalign')
    )

    for (pdb1, pdb2), d in zip(todo, output):
        name1 = dataset.get_id_from_filename(pdb1)
        name2 = dataset.get_id_from_filename(pdb2)
        # each value is a tuple (tm-core, RMSD)
        dist[name1][name2] = (d[0], d[2])
        dist[name2][name1] = (d[0], d[2])

    save(dist, dump_path)
    num_proteins = len(paths)
    DM = np.zeros((num_proteins, num_proteins))
    DM = []
    for i in range(num_proteins):
        for j in range(i+1, num_proteins):
            # take the largest TMscore (most similar) between both
            # directions and convert to a distance
            DM.append(
                1 - max(
                    dist[pdbids[i]][pdbids[j]][0],
                    dist[pdbids[j]][pdbids[i]][0]
                )
            )
    DM = np.array(DM).reshape(-1, 1)

    if isinstance(dataset.similarity_threshold_structure, float):
        thresholds = [dataset.similarity_threshold_structure]
    else:
        thresholds = dataset.similarity_threshold_structure

    for d in thresholds:
        clusterer = AgglomerativeClustering(n_clusters=None, distance_threshold=(1-d))
        clusterer.fit(DM)
        for i, p in enumerate(proteins):
            p['protein'][f'structure_cluster_{d}'] = int(clusterer.labels_[i])
    replace_avro_files(dataset, proteins)

if __name__ == '__main__':
    ds = get_dataset()
    compute_clusters_structure(ds)
    transfer_dataset(ds, 'structure')
