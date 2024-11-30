import os
os.environ['PYTHONUNBUFFERED'] = '1'
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

    recipe_tag = "Amber_AA"
    embed_options = [("[NX4]([#1])([#1])([#1])[CX4][CX3](=O)[O,#1]", (1, 7)),
                     ("[NX3,NX4]([#1])([#1])[CX4][CX3](=O)[O,#1]", (1, 6)),
                     ("[NX3]([#1])[CX4][CX3](=O)[O,#1]", (1, 5))]
    cap_option = ("[NX2,NX3][CX4][CX3](=O)", {0})
    cap_allowed_smarts = cap_option[0]
    cap_protonate = True
    pattern_to_label_mapping_standard = {'[NX3h1]': 'N-term', '[CX3h1]': 'C-term'}

class Amber_NA_recipe: 
    
    recipe_tag = "Amber_NA"
    embed_options = [("[OH,O-,#1][PX4](=O)([O])[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2][H]", (0, 12)),
                     ("[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2][H]", (8,))]
    cap_option = ("[OX2][CX4][CX4]1[OX2][CX4][CX4][CX4]1[OX2]", {0})
    cap_allowed_smarts = cap_option[0]
    cap_protonate = False
    pattern_to_label_mapping_standard = {'[PX4h1]': '5-prime', '[O+0X2h1]': '3-prime'}

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
                    
                    for embed_option in recipe.embed_options: 
                        embed_allowed_smarts = embed_option[0]
                        matches = cc_from_lib.rdkit_mol.GetSubstructMatches(Chem.MolFromSmarts(embed_allowed_smarts))
                        print(len(matches), embed_allowed_smarts, Chem.MolToSmiles(cc_from_lib.rdkit_mol))
                        if len(matches)==1:
                            break
                    print("Lib to Template Success:", resname)
                    print(f"{embed_allowed_smarts=}")
                    print(f"{recipe.cap_allowed_smarts=}")
                    if recipe.recipe_tag == "Amber_AA": 
                        variant_dict = {
                            "":  ({embed_allowed_smarts: embed_option[1]}, None), # embedded amino acid
                            "_N": ({embed_allowed_smarts: {embed_option[1][1]}}, {recipe.cap_option[0]: recipe.cap_option[1]}), # N-term amino acid
                            "_C": ({embed_allowed_smarts: {embed_option[1][0]}}, None), # C-term amino acid
                        }
                    if recipe.recipe_tag == "Amber_NA": 
                        variant_dict = {
                            "_":  ({embed_allowed_smarts: embed_option[1]}, None), # embedded nucleotide / nucleoside
                        }
                        if len(embed_option[1]) > 1: 
                            variant_dict.update({"_3": ({embed_allowed_smarts: {embed_option[1][0]}}, None)}) # 3' end nucleotide 
                            variant_dict.update({"_5p": ({embed_allowed_smarts: {embed_option[1][1]}}, None)}) # 5' end nucleotide (extra phosphate than canonical X5)
                            variant_dict.update({"_5": ({embed_allowed_smarts: {0,1,2,3,embed_option[1][1]}}, {recipe.cap_option[0]: recipe.cap_option[1]})}) # 5' end nucleotide 

                    cc_variants = add_variants(cc_orig = cc_from_lib, cc_list = [], 
                                embed_allowed_smarts = embed_allowed_smarts, 
                                cap_allowed_smarts = recipe.cap_allowed_smarts, cap_protonate = recipe.cap_protonate, 
                                pattern_to_label_mapping_standard = recipe.pattern_to_label_mapping_standard, 
                                variant_dict = variant_dict, 
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