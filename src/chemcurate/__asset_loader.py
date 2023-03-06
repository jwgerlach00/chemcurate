from importlib import resources
import yaml
from yaml.loader import SafeLoader


with resources.open_text('assets', 'uniprot_mapping.yaml') as f:
    uniprot_mapping = yaml.load(f, Loader=SafeLoader)
