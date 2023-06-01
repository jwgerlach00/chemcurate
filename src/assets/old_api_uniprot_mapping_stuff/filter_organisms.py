import yaml
from yaml import SafeLoader


def filter_organisms(protein_num_threshold:int) -> dict:
    with open('uniprot_mapping.yaml') as f:
        uniprot_mapping = yaml.load(f, Loader=SafeLoader)

    keys = [x for x in uniprot_mapping.keys() if len(uniprot_mapping[x]['protein']) > protein_num_threshold]

    filtered = {}
    for key in keys:
        filtered[key] = uniprot_mapping[key]

    return filtered          


if __name__ == '__main__':
    filtered = filter_organisms(1000)

    with open('uniprot_mapping_filtered.yaml', 'w') as f:
        yaml.dump(filtered, f)
