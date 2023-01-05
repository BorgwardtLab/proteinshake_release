import os, itertools
import numpy as np
from collections import defaultdict
import os.path as osp
import tempfile
import shutil
import subprocess
import pandas as pd
import traceback
import re
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm

def madoka_wrapper(paths, n_jobs):
    assert shutil.which('MADOKA') is not None and shutil.which('dssp2') is not None,\
           "No MADOKA/dssp2 installation found. Go here to install : http://madoka.denglab.org/download.html"

    pairs = list(itertools.combinations(range(len(paths)), 2))
    todo = [(paths[p1], paths[p2]) for p1, p2 in pairs]
    dist = defaultdict(lambda: {})

    with tempfile.TemporaryDirectory() as tmpdir:

        def dssp(path):
            pdb = os.path.basename(path).rstrip('.pdb')
            subprocess.run(['dssp2','-i', path, '-o', f'{tmpdir}/{pdb}.sse'], stdout=subprocess.PIPE)

        def madoka(path1, path2):
            cwd = os.path.dirname(path1)
            pdb1 = os.path.basename(path1).rstrip('.pdb')
            pdb2 = os.path.basename(path2).rstrip('.pdb')
            subprocess.run(['MADOKA','-o', './', f'{pdb1}.sse', f'{pdb2}.sse'], cwd=tmpdir, stdout=subprocess.PIPE)
            with open(f'{tmpdir}/{pdb1}-{pdb2}-result.txt', 'r') as file:
                search = re.search('TM-score \d*.\d*', file.read())
                tm_score = float(search.group(0).split()[1]) if not search is None else 0
            return tm_score

        Parallel(n_jobs=n_jobs)(delayed(dssp)(path) for path in tqdm(paths, desc='DSSP'))
        output = Parallel(n_jobs=n_jobs)(delayed(madoka)(*pair) for pair in tqdm(todo, desc='MADOKA'))

        for (pdb1, pdb2), tm in zip(todo, output):
            dist[pdb1][pdb2] = tm
            dist[pdb2][pdb1] = tm

        num_proteins = len(paths)
        DM = np.zeros((num_proteins, num_proteins))
        DM = []
        for i in range(num_proteins):
            for j in range(i+1, num_proteins):
                DM.append(dist[paths[j]][paths[i]][0])
        DM = np.array(DM).reshape(-1, 1)
        return DM


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
