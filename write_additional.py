import os
from chemtempgen import *
import traceback

leap_lib_path = '/Users/amyhe/Downloads/amber24_src/dat/leap/lib/'
AA_lib = (leap_lib_path+x for x in (
    # AA in ff19SB
    'amino19.lib',
    'amino19ipq_0.9.lib',
    'aminoct19ipq_0.9.lib',
    'aminont19ipq_0.9.lib',
    # mod AA
    'phosaa19SB.lib', 
    'mod_amino19.lib',
))

NA_lib = (leap_lib_path+x for x in (
    # RNA
    'RNA.lib',
    'terminalphos.LJbb-RNA.lib',
    # DNA
    'DNA.OL15.lib',
    'parmBSC1.lib',
    # mod NA
    'all_modrna08.lib',
))

class Amber_AA_recipe: 

    embed_allowed_smarts = "[NX3,NX4]([H])([H])[CX4][CX3](=O)[X]"
    cap_allowed_smarts = "[NH2,NH3][CX4][CX3](=O)"
    cap_protonate = True
    pattern_to_label_mapping_standard = {'[NX3h1]': 'N-term', '[CX3h1]': 'C-term'}

    variant_dict = {
            "":  ({"[NX3]([H])([H])[CX4][CX3](=O)[X]": {1, 6}}, None), # embedded amino acid
            "_N": ({"[NX3]([H])([H])[CX4][CX3](=O)[X]": {6}}, {"[NH2,NH3][CX4][CX3](=O)": {0}}), # N-term amino acid
            "_C": ({"[NX3]([H])([H])[CX4][CX3](=O)[X]": {1}}, None), # C-term amino acid
        }

class Amber_NA_recipe: 
        
    embed_allowed_smarts = "[X][PX4](=O)([O])[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2][H]"
    cap_allowed_smarts = "[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2]"
    cap_protonate = False
    pattern_to_label_mapping_standard = {'[PX4h1]': '5-prime', '[O+0X2h1]': '3-prime'}
    variant_dict = {
            "_":  ({"[PX4]([H])(=O)([O])[OX2][CX4]": {1}, "[OH][PX4](=O)([O])[OX2][CX4]": {0}, "[O-][PX4](=O)([O-])[OX2][CX4]": {0}, "[CX4]1[OX2][CX4][CX4][CX4]1[OX2][H]": {6}}, None), # embedded nucleotide 
            "_3": ({"[PX4]([H])(=O)([O])[OX2][CX4]": {1}, "[OH][PX4](=O)([O])[OX2][CX4]": {0}, "[O-][PX4](=O)([O-])[OX2][CX4]": {0}}, None), # 3' end nucleotide 
            "_5p": ({"[CX4]1[OX2][CX4][CX4][CX4]1[OX2][H]": {6}}, None), # 5' end nucleotide (extra phosphate than canonical X5)
            "_5": ({"[X][PX4](=O)([O])[OX2][CX4]": {0,1,2,3}, "[CX4]1[OX2][CX4][CX4][CX4]1[OX2][H]": {6}}, {"[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2]": {0}}), # 5' end nucleoside (canonical X5 in Amber)
        }

# Check conflicts
default_json_fn = '/Users/amyhe/Desktop/0_forks/gemile/residue_chem_templates.vanila.json'
with open(default_json_fn, 'r') as f:
    default_dict = json.load(f)
default_resnames = default_dict['ambiguous'].keys() | default_dict['residue_templates'].keys()


cc_list = []

for recipe, libfiles in [(Amber_AA_recipe, AA_lib), (Amber_NA_recipe, NA_lib)]:
    for lib in libfiles:

        json_fname = os.path.basename(lib)

        print("Lib Name:", lib)
        num_failures = 0

        residues = parmed.amber.offlib.AmberOFFLibrary.parse(lib)
        for resname in residues: 

            if resname in default_resnames or resname in set(cc.parent for cc in cc_list) or resname in set(cc.resname for cc in cc_list):
                print("Skipping existing residue: ", resname)
                continue
            
            try:
                # from lib
                cc_from_lib = ChemicalComponent.from_lib(lib, resname)
                if cc_from_lib is None:
                    print("Template can't be generated for: ", resname)
                    num_failures += 1
                
                else:
                    print("Lib to Template Success:", resname)
                    
                    cc_variants = add_variants(cc_orig = cc_from_lib, cc_list = [], 
                                embed_allowed_smarts = recipe.embed_allowed_smarts, 
                                cap_allowed_smarts = recipe.cap_allowed_smarts, cap_protonate = recipe.cap_protonate, 
                                pattern_to_label_mapping_standard = recipe.pattern_to_label_mapping_standard, 
                                variant_dict = recipe.variant_dict, 
                                acidic_proton_loc = {}
                                )
                    
                    cc_from_lib = cc_from_lib.make_pretty_smiles()
                    cc_from_lib.resname += "_fl"
                    amb_cc_list = [cc_from_lib]
                    
                    if not cc_variants: 
                        print("No template generated for: ", resname)
                    else:
                        for var in cc_variants: 
                            if any(var==cc for cc in amb_cc_list): 
                                continue
                            else:
                                amb_cc_list.append(var)

                    # from CCD
                    try: 
                        cc_from_ccd = build_noncovalent_CC(resname)

                        if not cc_from_ccd: 
                            print("No CCD template available for: ", resname)
                            continue

                        cc_from_ccd.resname += '_fl-ccd'
                        amb_cc_list += [cc_from_ccd]

                        cc_variants = build_linked_CCs(resname)

                        if not cc_variants: 
                            print("No CCD template generated for: ", resname)
                        else:
                            for var in cc_variants: 
                                var.resname += '-ccd'
                                if any(var==cc for cc in amb_cc_list): 
                                    continue
                                else:
                                    amb_cc_list.append(var)
                    except: 
                        continue

                    cc_list += amb_cc_list

                    
            except Exception as e:
                traceback.print_exc()
                print("Error during template generation for: ", resname)
                print(e)
                num_failures += 1
        
        if num_failures > 0:
            print("Num Failures:", num_failures, "Lib Name:", lib)

export_chem_templates_to_json(cc_list=cc_list, json_fname="all_additional.json")