import time
start_time = time.time()

import gemmi
import json
import pathlib
from rdkit import Chem
from rdkit.Chem import rdmolops

from rdkit import RDLogger
logger = RDLogger.logger()
logger.setLevel(RDLogger.CRITICAL)
import sys, logging, chemtempgen
logger = logging.getLogger('chemtempgen')
logger.setLevel(logging.WARNING)
handler = logging.StreamHandler(sys.stdout)

# Download modified_to_change_data.json
import urllib.request

url = "https://rna.bgsu.edu/modified/modified_to_change_data.json"
file_path = "modified_to_change_data.json"

try:
    urllib.request.urlretrieve(url, file_path)
    logging.info(f"File downloaded successfully: {file_path}")
except Exception as e:
    logging.error(f"Failed to download file. Error: {e}")

# Load JSON file
with open(file_path, 'r') as f:
    data = json.load(f)

# Filter by error
filter_key = 'error'
filter_value = 'not mappable'
count_error = sum(1 for value in data.values() if filter_key in value)
not_mappable = {key: value for key, value in data.items() if value.get(filter_key) == filter_value}
print(f"Total number of residues: {len(data)}")
print(f"Number of residues with error: {count_error}")
print(f"Number of residues not mappable: {len(not_mappable)}")

mappable = {key: value for key, value in data.items() if filter_key not in value}

# Make CC
from chemtempgen import * 
import copy

# Constants & Configuration
acidic_proton_loc_canonical = {
        # any carboxylic acid, sulfuric/sulfonic acid/ester, phosphoric/phosphinic acid/ester
        '[H][O]['+atom+'](=O)': 0 for atom in ('CX3', 'SX4', 'SX3', 'PX4', 'PX3')
    } | {
        # any thio carboxylic/sulfuric acid
        '[H][O]['+atom+'](=S)': 0 for atom in ('CX3', 'SX4')
    } | {
        '[H][SX2][a]': 0, # thiophenol
    } 
embed_allowed_smarts = "[O][PX4](=O)([O])[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2][H]"
cap_allowed_smarts = "[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2]"
pattern_to_label_mapping_standard = {'[PX4h1]': '5-prime', '[O+0X2h1:1]': '3-prime'}

variant_dict = {
        "":  ({"[O][PX4](=O)([O])[OX2][CX4]": {0} ,"[CX4]1[OX2][CX4][CX4][CX4]1[OX2][H]": {6}}, None), # embedded nucleotide 
        "3": ({"[O][PX4](=O)([O])[OX2][CX4]": {0}}, None), # 3' end nucleotide 
        "5p": ({"[CX4]1[OX2][CX4][CX4][CX4]1[OX2][H]": {6}}, None), # 5' end nucleotide (extra phosphate than canonical X5)
        "5": ({"[O][PX4](=O)([O])[OX2][CX4]": {0,1,2,3}, "[CX4]1[OX2][CX4][CX4][CX4]1[OX2][H]": {6}}, {"[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2]": {0}}), # 5' end nucleoside (canonical X5 in Amber)
    }

# Single Make Process
def make_variants(basename: str) -> list[ChemicalComponent]: 

    cc_from_cif = ChemicalComponent.from_cif(fetch_from_pdb(basename), basename)
    if cc_from_cif is None:
        return None

    editable = cc_from_cif.rdkit_mol.GetSubstructMatches(Chem.MolFromSmarts(embed_allowed_smarts))
    if not editable:
        logging.warning(f"Molecule doesn't contain embed_allowed_smarts: {embed_allowed_smarts} -> no templates will be made. ")
        return None

    cc_variants = []
    for suffix in variant_dict:
        cc = copy.deepcopy(cc_from_cif)
        cc.resname += suffix
        print(f"*** using CCD residue {basename} to construct {cc.resname} ***")

        cc = (
            cc
            .make_canonical(acidic_proton_loc = acidic_proton_loc_canonical) 
            .make_embedded(allowed_smarts = embed_allowed_smarts, 
                           leaving_smarts_loc = variant_dict[suffix][0])
            )
        if len(rdmolops.GetMolFrags(cc.rdkit_mol))>1:
            logging.warning(f"Molecule breaks into fragments during the deleterious editing of {cc.resname} -> skipping the vaiant... ")
            continue

        cc = (
            cc
            .make_capped(allowed_smarts = cap_allowed_smarts, 
                         capping_smarts_loc = variant_dict[suffix][1]) 
            .make_pretty_smiles()
            .make_link_labels_from_patterns(pattern_to_label_mapping = pattern_to_label_mapping_standard)
            )

        # Try Meeko check
        try:
            cc.ResidueTemplate_check()
        except Exception as e:
            logging.error(f"Template Failed to pass Meeko check. Error: {e}")
            continue
        
        # Check redundancy
        if any(cc == other_variant for other_variant in cc_variants):
            logging.error(f"Template Failed to pass redundancy check -> skipping the template... ")
            continue

        cc_variants.append(cc)
        print(f"*** finish making {cc.resname} ***")

    return cc_variants


# Check conflicts
default_json_fn = '/Users/amyhe/Desktop/0_forks/Meeko/meeko/data/residue_chem_templates.json'
with open(default_json_fn, 'r') as f:
    default_dict = json.load(f)
default_resnames = default_dict['ambiguous'].keys() | default_dict['residue_templates'].keys()


# Run Make Process
cc_byparent = {}
for cc_name in mappable: 

    if cc_name in default_resnames: 
        logging.error(f"Name conflict with existing residue: {cc_name}. Skipping... ")
        continue

    variant_list = make_variants(cc_name)
    if variant_list:
        cc_byparent[cc_name] = variant_list
    else:
        cc_byparent[cc_name] = None
        logging.warning(f"No template generated for {cc_name}. ")

# Write to JSON
write_to_single_json_fn = "NAKB_templates.json"
residue_counter = 0
variant_counter = 0

if write_to_single_json_fn:
    temp_list = [] 
    for cc_name in cc_byparent:
        if cc_byparent[cc_name] is not None:
            temp_list += cc_byparent[cc_name]
            residue_counter += 1
            variant_counter += len(cc_byparent[cc_name])
    export_chem_templates_to_json(temp_list, write_to_single_json_fn)
else:
    for cc_name in cc_byparent:
        if cc_byparent[cc_name] is not None:
            residue_counter += 1
            variant_counter += len(cc_byparent[cc_name])
            export_chem_templates_to_json(cc_byparent[cc_name], f"{cc_name}_templates.json")

end_time = time.time()
execution_time = end_time - start_time
print(f"Generated {variant_counter} templates for {residue_counter} residues")
print(f"Execution time: {execution_time} seconds")