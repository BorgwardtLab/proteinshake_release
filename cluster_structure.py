import os
from proteinshake.utils import tmalign_wrapper, unzip_file
from tqdm import tqdm
import numpy as np
from joblib import Parallel, delayed

def gdt(path1, path2):
    return tmalign_wrapper(path1, path2)['GDT']

def compute_clusters_structure(dataset, thresholds=[0.5, 0.6, 0.7, 0.8, 0.9]):

    proteins = list(dataset.proteins()[0])[:10]
    pdbids = [p['protein']['ID'] for p in proteins]
    path_dict = {dataset.get_id_from_filename(os.path.basename(f)):f for f in dataset.get_raw_files()}
    paths = [path_dict[id] for id in pdbids]
    paths = [unzip_file(p, remove=False) if not os.path.exists(p.rstrip('.gz')) else p.rstrip('.gz') for p in tqdm(paths, desc='Unzipping')]
    num_proteins = len(proteins)

    C = 0 # current cluster center
    R = 0 # current cluster radius
    M = {} # storage for distances from C to all points
    D = {} # storage for radii

    # k-center algorithm
    while R*2 <= max(thresholds):
        id = proteins[C]['protein']['ID']
        M[id] = Parallel(n_jobs=dataset.n_jobs)(delayed(gdt)(paths[C], path_i) for path_i in paths) # compute GDT (similarity) from current center to all points
        m = np.max(np.vstack([_ for _ in M.values()]), axis=0) # get maximum similarity over clusters for all points
        C = np.argmin(m) # next cluster center
        R = np.min(m) # current cluster radius
        D[id] = R*2
        print(C,R)

    # assign cluster centers
    for threshold in thresholds:
        center_distances = np.stack([M[id] for id in M.keys() if D[id] <= threshold]) # select all centers with radius <= threshold
        print(center_distances.shape)
        center_ids = np.argmax(center_distances, axis=0) # get maximum similarity over clusters for all points
        print(center_ids.shape, threshold, center_ids)
        for p, c in zip(proteins, center_ids):
            p['protein'][f'sequence_cluster_{threshold}'] = c
    print(proteins[9]['protein'])


    



