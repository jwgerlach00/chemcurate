import yaml

with open('src/assets/uniprot_mapping_filtered.yaml') as f:
    uniprot_mapping = yaml.load(f, Loader=yaml.FullLoader)
    
targets = list(uniprot_mapping.keys())
print(len(targets))

with open('targets.txt', 'w') as f:
    for target in targets:
        f.write(f'{target}\n')
