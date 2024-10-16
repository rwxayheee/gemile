from chemtempgen import *

acidic_proton_loc_canonical = {
    # any carboxylic acid, sulfuric/sulfonic acid/ester, phosphoric/phosphinic acid/ester
    '[H][O]['+atom+'](=O)': 0 for atom in ('CX3', 'SX4', 'SX3', 'PX4', 'PX3')
} | {
    # any thio carboxylic/sulfuric acid
    '[H][O]['+atom+'](=S)': 0 for atom in ('CX3', 'SX4')
} | {
    '[H][SX2][a]': 0, # thiophenol
} 

source_cif = '/Users/amyhe/Desktop/0_forks/gemile/examples/1_embedded_PSU/PSU.cif'
basename = 'PSU'
cc = ChemicalComponent.from_cif(source_cif, basename)

cc.resname = 'PSU5'
cc = (
    cc
    .make_canonical(acidic_proton_loc = acidic_proton_loc_canonical) 
    .make_embedded(allowed_smarts = "[O][PX4](=O)([O])[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2][H]", 
                   leaving_smarts_loc = {"[O][PX4](=O)([O])[OX2][CX4]": {0,1,2,3}, "[CX4]1[OX2][CX4][CX4][CX4]1[OX2][H]": {6}}) 
    .make_capped(allowed_smarts = "[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2]", 
                 capping_smarts_loc = {"[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2]": {0}}) 
    .make_pretty_smiles()
    .make_link_labels_from_patterns(pattern_to_label_mapping = {'[PX4h1]': '5-prime', '[O+0X2h1:1]': '3-prime'})
    # equlavent to 
    # .make_link_labels_from_names(name_to_label_mapping = {'P': '5-prime', "O3'": '3-prime'})
)

export_chem_templates_to_json([cc], 'PSU_residue_templates.json')