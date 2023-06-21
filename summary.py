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
  font-size: 14px;
  border-radius: 5px;
  overflow:hidden;
}
tr {
  background-color: #fff;
}
td {
  text-align: left;
  border: 3px solid #e3e8ef;
  padding: 10px;
  color:#111729;
}
th {
  padding: 10px;
  text-align: left;
  background-color: #111729;
  border: 3px solid #111729;
  color:#fff;
}
tr:nth-child(even){
  background-color: #F2F5F9;
}
</style>
</head>
<body>
"""
html_post = """
</body></html>
"""

'''
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
with open(f'docs/datasets.html','w') as file:
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
with open(f'docs/tasks.html','w') as file:
    file.write(html_pre + task_table + html_post)


'''

# For each task: label distribution in train/test/val
add_js = 'cdn'
with open('docs/labels.html','w') as file:
    file.write("<html><head></head><body>" + "\n")
    for ax,name in enumerate(TASKS):
        print(name)
        task = get_task(ROOT, name, NJOBS)
        targets = list(task.train_targets) + list(task.test_targets) + list(task.val_targets)

        if task.type == 'Multiclass Classification' or task.type == 'Multilabel Classification':
            token_map = {v:k for k,v in task.token_map.items()}
            token_map = np.array([token_map[i] for i in range(len(token_map))])
            if task.type == 'Multilabel Classification':
                targets = [token_map[np.array(labels, dtype=bool)] for labels in targets]
                targets = itertools.chain.from_iterable(targets)
            else:
                targets = [token_map[label] for label in targets]
            labels,counts = zip(*Counter(targets).most_common())
            plot = px.bar(x=labels, y=counts, title=name, labels={'x':'Target Label', 'y':'Count'}, template='plotly_white')
            plot.update_xaxes(range=[-1, min(20, len(counts))])
        elif task.type == 'Binary Classification':
            if name == 'BindingSiteDetectionTask':
                targets = [task.target(i) for i in task.proteins]
            counts = [np.array(t).sum() for t in targets]
            plot = px.histogram(x=counts, nbins=100, title=name, labels={'x':'Number of Residue Contacts', 'count':'Count'}, template='plotly_white')
            plot.update_layout(yaxis_title="Count")
        elif task.type == 'Regression':
            plot = px.histogram(x=targets, nbins=100, title=name, labels={'x':'Target Value', 'count':'Count'}, template='plotly_white')
            plot.update_layout(yaxis_title="Count")
        elif task.type == 'Retrieval':
            lengths = [len(t) for t in targets]
            plot = px.histogram(x=lengths, nbins=100, title=name, labels={'x':'Number of similar instances', 'count':'Count'}, template='plotly_white')
            plot.update_layout(yaxis_title="Count")
        else:
            continue
        plot.update_layout(height=300)
        file.write(plot.to_html(full_html=False, include_plotlyjs=add_js))
        add_js = False
    file.write("</body></html>" + "\n")