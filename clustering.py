import os
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from tqdm import tqdm
from proteinshake.utils import save, unzip_file, write_avro

from wrappers import cdhit_wrapper, madoka_wrapper

def replace_avro_files(dataset, proteins):
    residue_proteins = list(dataset.proteins(resolution='residue')[0])
    atom_proteins = list(dataset.proteins(resolution='atom')[0])
    for r,a,p in zip(residue_proteins,atom_proteins,proteins):
        r['protein'] = p['protein']
        a['protein'] = p['protein']
    write_avro(list(residue_proteins), f'{dataset.root}/{dataset.name}.residue.avro')
    write_avro(list(atom_proteins), f'{dataset.root}/{dataset.name}.atom.avro')

def compute_clusters_sequence(dataset):
    """ Use CDHit to cluster sequences. Assigns the field 'sequence_cluster' to an integer cluster ID for each protein.
    """
    proteins = list(dataset.proteins()[0])
    if isinstance(dataset.similarity_threshold_sequence, float):
        thresholds = [dataset.similarity_threshold_sequence]
    else:
        thresholds = dataset.similarity_threshold_sequence

    representatives = {}
    sequences = [p['protein']['sequence'] for p in proteins]
    ids = [p['protein']['ID'] for p in proteins]
    for threshold in thresholds:
        clusters, reps = cdhit_wrapper(ids, sequences, sim_thresh=threshold, n_jobs=dataset.n_jobs)
        representatives[threshold] = reps
        if clusters == -1:
            print("Sequence clustering failed.")
            return
        for p, c in zip(proteins, clusters):
            p['protein'][f'sequence_cluster_{threshold}'] = c
    save(representatives, f'{dataset.root}/{dataset.name}.cdhit.json')
    replace_avro_files(dataset, proteins)

def compute_clusters_structure(dataset):
    """ Launch TMalign on all pairs of proteins in dataset.
    Assign a cluster ID to each protein at protein-level key 'structure_cluster'.

    Saves TMalign output to `dataset.root/{Dataset.__class__}.tmalign.json.gz`
    """
    proteins = list(dataset.proteins()[0])
    path_dict = {dataset.get_id_from_filename(os.path.basename(f)):f for f in dataset.get_raw_files()}
    paths = [path_dict[p['protein']['ID']] for p in proteins]

    dump_name = f'{dataset.name}.tmalign.json'
    dump_path = os.path.join(dataset.root, dump_name)

    if dataset.n_jobs == 1:
        print('Computing the TM scores with use_precompute = False is very slow. Consider increasing n_jobs.')

    paths = [unzip_file(p, remove=False) if p.endswith('.gz') else p for p in tqdm(paths, desc='Unzipping')]
    pdbids = [dataset.get_id_from_filename(p) for p in paths]

    DM = madoka_wrapper(paths, dataset.n_jobs)

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
