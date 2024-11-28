import json
import mols2grid
import pandas as pd
from rdkit import Chem

def write_to_html(json_fname, html_fname):
    with open(json_fname, 'r') as f:
        data = json.load(f)
    
    ambiguous, templates = (data["ambiguous"], data["residue_templates"])
    amb_dict = {v: k for k,v_set in ambiguous.items() for v in v_set}
    res_in_dict = {"resname": [res for res in templates]}

    df = pd.DataFrame.from_dict(res_in_dict)
    for colname in ("mol", "ID", "Smiles", "parent"): 
        df[colname] = ''

    for ind in df.index:
        resname = df["resname"][ind]
        df["parent"][ind] = amb_dict[resname] if resname in amb_dict else ''

        df[colname][ind] = int(ind+1)
        smiles = templates[resname]["smiles"]
        df["Smiles"][ind] = smiles
        ps = Chem.SmilesParserParams()
        ps.removeHs = False
        mol = Chem.MolFromSmiles(smiles, ps)
        for i, atom in enumerate(mol.GetAtoms()):
            atom.SetProp("atomNote", templates[resname]["atom_name"][i])
        df["mol"][ind] = mol

    mols2grid.save(df, mol_col="mol", #prerender = True, 
                   subset=['resname', 'parent', 'ID', 'img'], output=html_fname)

write_to_html("residue_chem_templates.json", "residue_chem_templates.html")