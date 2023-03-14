

def summary(ds):
    return {
        'Name': ds.name,
        'Size': len(ds.proteins())
    }