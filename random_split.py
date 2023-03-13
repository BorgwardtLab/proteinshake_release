from sklearn.model_selection import train_test_split
from util import replace_avro_files, get_paths

def compute_random_split(dataset, test_ratio=0.1, val_ratio=0.1, seed=42):
    """ Compute a simple random split.
    """
    print(f'Random split {dataset.name}')

    proteins = list(dataset.proteins())
    indices = range(len(proteins))
    train, valtest = train_test_split(indices, test_size=test_ratio, random_state=seed)
    val, test = train_test_split(valtest, test_size=test_ratio/(test_ratio+val_ratio), random_state=seed)

    for i,p in enumerate(proteins):
        if p['protein']['ID'] in test: p['protein'][f'random_split'] = 'test'
        elif p['protein']['ID'] in val: p['protein'][f'random_split'] = 'val'
        else: p['protein'][f'random_split'] = 'train'
    replace_avro_files(dataset, proteins)