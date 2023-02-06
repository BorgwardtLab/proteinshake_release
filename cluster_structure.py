import os
from proteinshake.utils import tmalign_wrapper, unzip_file
from tqdm import tqdm
import numpy as np
from joblib import Parallel, delayed

def gdt(path1, path2):
    return tmalign_wrapper(path1, path2)['GDT']

def compute_clusters_structure(dataset, thresholds=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]):

    proteins = list(dataset.proteins()[0])[:100]
    pdbids = [p['protein']['ID'] for p in proteins]
    path_dict = {dataset.get_id_from_filename(os.path.basename(f)):f for f in dataset.get_raw_files()}
    paths = [path_dict[id] for id in pdbids]
    paths = [unzip_file(p, remove=False) if not os.path.exists(p.rstrip('.gz')) else p.rstrip('.gz') for p in tqdm(paths, desc='Unzipping')]
    num_proteins = len(proteins)

    # cd-hit algorithm, compare Holm & Sander 1998
    for threshold in thresholds: # do for each threshold
        clusters, cluster_ids = [], []
        pool = np.argsort([len(p['protein']['sequence']) for p in proteins])[::-1] # sort by length (descending)
        print(f'Structure clustering threshold {threshold}')
        with tqdm(total=num_proteins) as pbar:
            while len(pool) > 0:
                clusters, pool = [*clusters, pool[0]], pool[1:] # choose largest protein as first representative
                cluster_ids.append(proteins[clusters[-1]]['protein']['ID'])
                proteins[clusters[-1]]['protein'][f'sequence_cluster_{threshold}'] = len(cluster_ids)-1
                M = np.zeros((len(clusters),len(pool)))
                for i,rep in enumerate(clusters):
                    M[i] = Parallel(n_jobs=dataset.n_jobs)(delayed(gdt)(paths[rep], paths[i]) for i in pool) # compute GDT (similarity) from current center to all points
                cluster_assignments = np.argmax(M, axis=0)[np.max(M, axis=0) >= threshold]
                protein_idx = pool[np.max(M, axis=0) >= threshold]
                for i,c in zip(protein_idx, cluster_assignments):
                    proteins[i]['protein'][f'sequence_cluster_{threshold}'] = c
                pool = pool[np.max(M, axis=0) < threshold]
                print(len(clusters),len(pool), len(clusters)*len(pool))
                pbar.set_description(f'{len(clusters)} Clusters')
                pbar.update((np.max(M, axis=0) >= threshold).sum()+1)
    print(proteins[8]['protein'])


    



