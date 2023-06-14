import os, itertools
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from collections import Counter
from plotly.subplots import make_subplots
from proteinshake.datasets import __all__ as DATASETS
from proteinshake.tasks import __all__ as TASKS
from proteinshake.datasets.alphafold import AF_DATASET_NAMES
from util import get_dataset, get_task

TAG = '09MAR2023'
ROOT = os.path.expandvars(f'$LOCAL_SCRATCH/proteinshake/{TAG}')
NJOBS = 20

TASKS = [t for t in TASKS if not t in ['Task']]
DATASETS = [d for d in DATASETS if not d in ['Dataset']]

html_pre = """
<html>
<head>
<style>
table {
  border-collapse: collapse;
  font-family: Monaco;
  font-size: 16px;
}
tr:nth-child(even){background-color: #f2f2f2;}
th {
  background-color: #333;
  color: #fff;
}
th, td {
  text-align: left;
  border: 3px solid #333;
  padding: 10px;
}
</style>
</head>
<body>
"""
html_post = """
</body></html>
"""


# Overview Table with all datasets and tasks. Dataset: Name, Annotation, number of instances. Task: Name, type, level, metric.
dataset_table = []
for name in DATASETS:
    print(name)
    ds = get_dataset(ROOT, name, 'swissprot', NJOBS)
    dataset_table.append({
        'Name': name,
        'Size': len(ds.proteins()),
        'Description': ds.__doc__.split('.')[0]
    })
dataset_table = pd.DataFrame(dataset_table).to_html(index=False)
with open(f'debug/datasets.html','w') as file:
    file.write(html_pre + dataset_table + html_post)

task_table = []
for name in TASKS:
    print(name)
    task = get_task(ROOT, name, NJOBS)
    task_table.append({
        'Name': name,
        'Type': task.type,
        'Input': task.input,
        'Output': task.output,
    })
task_table = pd.DataFrame(task_table).to_html(index=False)
with open(f'debug/tasks.html','w') as file:
    file.write(html_pre + task_table + html_post)


exit()
# For each task: label distribution in train/test/val
TASKS = ['GeneOntologyTask']
fig = make_subplots(rows=len(TASKS), cols=1)
for ax,name in enumerate(TASKS):
    task = get_task(ROOT, name, NJOBS)
    targets = np.array([task.target(p) for p in task.proteins], dtype=object)

    #if task.type == 'Multilabel Classification':
    targets = itertools.chain.from_iterable(targets)
    counts = dict(Counter(targets))
    plot = go.Bar(x=list(counts.keys()), y=list(counts.values()))


    #if task.type == 'Multiclass Classification':
    #    counts = dict(Counter(targets))
    #    plot = go.Bar(x=list(counts.keys()), y=list(counts.values()))
    #if task.type == 'Regression':
    #    plot = go.Histogram(x=targets, nbinsx=100)
    fig.add_trace(plot, row=ax+1, col=1)
fig.update_layout(height=400*len(TASKS), width=600, showlegend=False)
fig.write_html('debug/test.html')