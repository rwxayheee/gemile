import os
from chemtempgen_for_amber import *
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

class Amber_AA_recipe(AA_recipe): 
    def __init__(self):
        super().__init__()

class Amber_NA_recipe(NA_recipe): 
    def __init__(self):
        super().__init__()

# Check conflicts
default_json_fn = '/Users/amyhe/Desktop/0_forks/Meeko/meeko/data/residue_chem_templates.json'
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
                    print("Success:", resname)
                    
                    cc_variants = add_variants(cc_orig = cc_from_lib, cc_list = [], 
                                embed_allowed_smarts = recipe.embed_allowed_smarts, 
                                cap_allowed_smarts = recipe.cap_allowed_smarts, cap_protonate = recipe.cap_protonate, 
                                pattern_to_label_mapping_standard = recipe.pattern_to_label_mapping_standard, 
                                variant_dict = recipe.variant_dict, 
                                acidic_proton_loc = {}
                                )
                    
                    if not cc_variants: 
                        print("No template generated for: ", resname)

                    cc_list += cc_variants

                    # from CCD
                    cc_from_ccd = build_noncovalent_CC(resname)

                    if not cc_from_ccd: 
                        print("No CCD template available for: ", resname)
                        continue

                    cc_from_ccd.resname += '-ccd'
                    cc_list += [cc_from_ccd]

                    cc_variants = build_linked_CCs(resname)

                    if not cc_variants: 
                        print("No CCD template generated for: ", resname)
                    else:
                        for cc in cc_variants:
                            cc.resname += '-ccd'
                        cc_list += cc_variants

                    
            except Exception as e:
                traceback.print_exc()
                print("Error during template generation for: ", resname)
                print(e)
                num_failures += 1
        
        if num_failures > 0:
            print("Num Failures:", num_failures, "Lib Name:", lib)

export_chem_templates_to_json(cc_list=cc_list, json_fname="all_additional.json")