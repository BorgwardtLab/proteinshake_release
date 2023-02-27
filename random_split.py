from sklearn.model_selection import train_test_split

def compute_random_split(dataset, test_ratio=0.1, train_ratio=0.1, seed=42):
    """ Compute a simple random split.
    """
    print(f'Random split {dataset.name}')

    proteins = list(dataset.proteins())
    indices = range(len(proteins))
    train, valtest = train_test_split(indices, test_size=1-train_ratio, random_state=seed)
    val, test = train_test_split(valtest, test_size=test_ratio/(test_ratio+val_ratio), random_state=seed)

    for i,p in enumerate(proteins):
        p['protein'][f'random_split'] = 'test' if i in test else 'train'
        p['protein'][f'random_split'] = 'val' if i in val else 'train'
    replace_avro_files(dataset, proteins)