import os
from chemtempgen_for_amber import *
import traceback

leap_lib_path = '/Users/amyhe/Downloads/amber24_src/dat/leap/lib/'
all_lib = (leap_lib_path+x for x in (
    # AA in ff19SB
    'amino19.lib',
    'amino19ipq_0.9.lib',
    'aminoct19ipq_0.9.lib',
    'aminont19ipq_0.9.lib',
    # mod AA
    'phosaa19SB.lib', 
    'mod_amino19.lib',
    # RNA
    'RNA.lib',
    'terminalphos.LJbb-RNA.lib',
    # DNA
    'DNA.OL15.lib',
    'parmBSC1.lib',
    # mod NA
    'all_modrna08.lib',
))

for lib in all_lib:
    json_fname = os.path.basename(lib)
    cc_list = []

    print("Lib Name:", lib)
    num_failures = 0

    residues = parmed.amber.offlib.AmberOFFLibrary.parse(lib)
    for resname in residues: 
        try:
            cc = ChemicalComponent.from_lib(lib, resname)
            if cc is None:
                print("We have a problem...", resname)
                num_failures += 1
            else:
                print("Success:", resname)
                cc_list.append(cc)
                
        except Exception as e:
            # traceback.print_exc()
            print("We have a problem...", resname)
            print(e)
            num_failures += 1
    
    if num_failures > 0:
        print("Num Failures:", num_failures, "Lib Name:", lib)
    export_chem_templates_to_json(cc_list=cc_list, json_fname=json_fname.replace(".lib", ".json"))