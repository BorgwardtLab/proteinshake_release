import os, itertools
from proteinshake.utils import tmalign_wrapper, unzip_file
from tqdm import tqdm
import numpy as np
from joblib import Parallel, delayed
from util import replace_avro_files
from sklearn.cluster import AgglomerativeClustering

def batched(iterable, size=10):
    iterator = iter(iterable)
    for first in iterator:
        yield itertools.chain([first], itertools.islice(iterator, size - 1))

def compute_clusters_structure(dataset, thresholds=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]):

    proteins = list(dataset.proteins()[0])
    pdbids = [p['protein']['ID'] for p in proteins]
    path_dict = {dataset.get_id_from_filename(os.path.basename(f)):f for f in dataset.get_raw_files()}
    paths = [path_dict[id] for id in pdbids]
    paths = [unzip_file(p, remove=False) if not os.path.exists(p.rstrip('.gz')) else p.rstrip('.gz') for p in tqdm(paths, desc='Unzipping')]
    num_proteins = len(proteins)
    savepath = f'{dataset.root}/{dataset.name}'

    TM, RMSD, GDT = [np.load(f'{savepath}.{m}.npy') if os.path.exists(f'{savepath}.{m}.npy') else np.ones((num_proteins,num_proteins), dtype=np.float16) * np.nan for m in ['tm','rmsd','gdt']]
    np.fill_diagonal(TM, 1.0), np.fill_diagonal(RMSD, 0.0), np.fill_diagonal(GDT, 1.0)

    combinations = itertools.combinations(range(num_proteins), 2)
    chunk_size = 10
    for chunk in tqdm(batched(combinations, chunk_size), desc='Structure clustering', total=int(np.ceil((num_proteins**2-num_proteins)/2/chunk_size))):
        chunk = np.array([(x,y) for x,y in chunk if np.isnan(TM[x,y]) or np.isnan(RMSD[x,y]) or np.isnan(GDT[x,y])])
        if len(chunk) == 0: continue
        d = Parallel(n_jobs=dataset.n_jobs)(delayed(tmalign_wrapper)(paths[i], paths[j]) for i,j in chunk)
        x,y = tuple(chunk[:,0]), tuple(chunk[:,1])
        TM[x,y] = [x['TM1'] for x in d]
        TM[y,x] = [x['TM2'] for x in d]
        RMSD[x,y] = [x['RMSD'] for x in d]
        RMSD[y,x] = [x['RMSD'] for x in d]
        GDT[x,y] = [x['GDT'] for x in d]
        GDT[y,x] = [x['GDT'] for x in d]
        # save periodically
        np.save(f'{savepath}.tm.npy', TM)
        np.save(f'{savepath}.rmsd.npy', RMSD)
        np.save(f'{savepath}.gdt.npy', GDT)
 
    # clustering
    for threshold in thresholds:
        clusterer = AgglomerativeClustering(n_clusters=None, distance_threshold=(1-threshold), metric='precomputed', linkage='average')
        clusterer.fit(1-GDT)
        for i, p in enumerate(proteins):
            p['protein'][f'structure_cluster_{threshold}'] = int(clusterer.labels_[i])
    replace_avro_files(dataset, proteins)
    print(proteins[0]['protein'])
    print('Sequence clustering done.')




    



