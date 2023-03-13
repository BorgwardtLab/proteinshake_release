import shutil, subprocess, tempfile
from proteinshake.utils import save
import os.path as osp
from tqdm import tqdm
from util import replace_avro_files, get_paths, split
from collections import defaultdict

def cdhit_wrapper(ids, sequences, sim_thresh=0.6, n_jobs=1):
    """ Cluster sequences using CD-hit

    Choose of word size:
    -n 5 for thresholds 0.7 ~ 1.0
    -n 4 for thresholds 0.6 ~ 0.7
    -n 3 for thresholds 0.5 ~ 0.6
    -n 2 for thresholds 0.4 ~ 0.5

    Parameters
    -----------
    sequences: list
        List of protein sequences to cluster.

    Returns
    --------
    representatives: list
        List of sequence indices to preserve as representatives.
    """
    assert sim_thresh >= 0.4 and sim_thresh <= 1, "Similarity threshold not in [0.4, 1]"

    if sim_thresh >= 0.4 and sim_thresh < 0.5:
        word_size = 2
    elif sim_thresh >= 0.5 and sim_thresh < 0.6:
        word_size = 3
    elif sim_thresh >= 0.6 and sim_thresh < 0.7:
        word_size = 4
    else:
        word_size = 5

    assert shutil.which('cd-hit') is not None,\
    "CD-HIT installation not found. Go here https://github.com/weizhongli/cdhit to install"

    n_jobs = 0 if n_jobs < 0 else n_jobs

    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = osp.join(tmpdir, 'in.fasta')
        out_file = osp.join(tmpdir, 'out.fasta')
        with open(in_file, "w") as inp:
            for id, s in zip(ids,sequences):
                inp.write(f">{id}\n")
                inp.write(s + "\n")
        cmd = ['cd-hit',
            '-c', str(sim_thresh),
            '-i', in_file,
            '-n', str(word_size),
            '-o', out_file,
            '-T', str(n_jobs),
            '-M', "0" # unlimited memory
            ]
        subprocess.run(cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT
                    )
        # parse cluster assignments
        pdb2cluster = {}
        cluster2pdb = defaultdict(list)
        with open(out_file + ".clstr", "r") as out:
            for line in out:
                if line.startswith(">"):
                    clust_id = int(line.split()[1])
                    continue
                pdb_id = line.split(">")[1].split('.')[0]
                pdb2cluster[pdb_id] = clust_id
                cluster2pdb[clust_id].append(pdb_id)
        return pdb2cluster, cluster2pdb


def compute_sequence_split(dataset, thresholds=[0.5, 0.6, 0.7, 0.8, 0.9], test_ratio=0.1, val_ratio=0.1):
    """ Use CDHit to cluster sequences. Assigns the field 'sequence_cluster' to an integer cluster ID for each protein.
    """
    if osp.exists(f'{dataset.root}/{dataset.name}.cdhit.json'): return

    print(f'Sequence split {dataset.name}')

    proteins = list(dataset.proteins())
    sequences = [p['protein']['sequence'] for p in proteins]
    pdbids, paths, path_dict = get_paths(dataset)

    for threshold in thresholds:
        pdb2cluster, cluster2pdb = cdhit_wrapper(pdbids, sequences, sim_thresh=threshold, n_jobs=dataset.n_jobs)
        def split_wrapper(ds, query, threshold):
            if not query in pdb2cluster: return []
            return cluster2pdb[pdb2cluster[query]]
        pool = [p for p in pdbids]
        test_size, val_size = int(len(pool)*test_ratio), int(len(pool)*val_ratio)
        pool, test = split(split_wrapper, dataset, pool, test_size, threshold)
        train, val = split(split_wrapper, dataset, pool, val_size, threshold)
        train, test, val = [dataset.get_id_from_filename(p) for p in train], [dataset.get_id_from_filename(p) for p in test], [dataset.get_id_from_filename(p) for p in val]
        for p in proteins:
            if p['protein']['ID'] in test: p['protein'][f'sequence_split_{threshold}'] = 'test'
            elif p['protein']['ID'] in val: p['protein'][f'sequence_split_{threshold}'] = 'val'
            elif p['protein']['ID'] in train: p['protein'][f'sequence_split_{threshold}'] = 'train'
            else: p['protein'][f'sequence_split_{threshold}'] = 'none'
    replace_avro_files(dataset, proteins)

