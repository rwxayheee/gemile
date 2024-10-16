from chemtempgen import *

source_cif = '/Users/amyhe/Desktop/0_forks/gemile/examples/1_embedded_PSU/PSU.cif'
basename = 'PSU'

cc = ChemicalComponent.from_cif(source_cif, basename)
cc.make_canonical()
cc.make_embedded(allowed_smarts = "[O][PX4](=O)([O])[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2][H]", 
                 leaving_names = {"HO3'"}, 
                 leaving_smarts_loc = {"[O][PX4](=O)([O])[OX2][CX4]": {0,1,2,3}, "[CX4]1[OX2][CX4][CX4][CX4]1[OX2][H]": {6}})

cc.make_capped(allowed_smarts = "[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2]", 
               capping_smarts_loc = {"[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2]": {0}})

cc.smiles_exh, cc.atom_name = get_smiles_with_atom_names(cc.rdkit_mol)
cc.smiles_exh = make_pretty_smiles(cc.smiles_exh)
print(cc.smiles_exh, cc.atom_name)


cc.make_link_labels_from_names()
cc.resname = 'PSU5'
cc.parent = basename

export_chem_templates_to_json([cc], 'PSU_residue_templates.json')