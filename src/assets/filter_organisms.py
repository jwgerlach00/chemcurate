import yaml
from yaml import SafeLoader


if __name__ == '__main__':
    with open('uniprot_mapping.yaml') as f:
        uniprot_mapping = yaml.load(f, Loader=SafeLoader)

    threshold = 100
    keys = [x for x in uniprot_mapping.keys() if len(uniprot_mapping[x]['protein']) > threshold]

    out = {}
    for key in keys:
        out[key] = uniprot_mapping[key]

    with open('uniprot_mapping_filtered.yaml', 'w') as f:
        yaml.dump(out, f)
