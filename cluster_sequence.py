import shutil, subprocess, tempfile
from proteinshake.utils import save
import os.path as osp
from tqdm import tqdm

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
    assert sim_thresh >= 0.4 and sim_thresh <= 1, "Threshold not in [0.4, 1]"

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
        try:
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
        except Exception as e:
            print(traceback.format_exc())
            return -1
        else:
            # parse cluster assignments
            clusters = {}
            representatives = []
            with open(out_file + ".clstr", "r") as out:
                for line in out:
                    if line.startswith(">"):
                        clust_id = int(line.split()[1])
                        continue
                    pdb_id = line.split(">")[1].split('.')[0]
                    clusters[pdb_id] = clust_id
                    if line.endswith('*'):
                        representatives.append(pdb_id)
            clusters = [clusters[id] if id in clusters else -1 for id in ids]
            return clusters, representatives

def compute_clusters_sequence(dataset, thresholds=[0.5, 0.6, 0.7, 0.8, 0.9]):
    """ Use CDHit to cluster sequences. Assigns the field 'sequence_cluster' to an integer cluster ID for each protein.
    """
    print('Starting sequence clustering.')

    proteins = list(dataset.proteins()[0])
    representatives = {}
    sequences = [p['protein']['sequence'] for p in proteins]
    ids = [p['protein']['ID'] for p in proteins]

    for threshold in thresholds:
        print(f'Threshold {threshold}')
        clusters, reps = cdhit_wrapper(ids, sequences, sim_thresh=threshold, n_jobs=dataset.n_jobs)
        representatives[threshold] = reps
        if clusters == -1:
            print("Sequence clustering failed.")
            return
        for p, c in zip(proteins, clusters):
            p['protein'][f'sequence_cluster_{threshold}'] = c
    save(representatives, f'{dataset.root}/{dataset.name}.sequence_cluster_centers.json')
    replace_avro_files(dataset, proteins)
    print('Sequence clustering done.')

