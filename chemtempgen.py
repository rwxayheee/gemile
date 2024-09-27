import gemmi

from rdkit import Chem
from rdkit import RDLogger
logger = RDLogger.logger()
logger.setLevel(RDLogger.CRITICAL)

def get_smiles_with_atom_names(mol: Chem.Mol) -> tuple[str, list[str]]:
    """Generate SMILES with atom names in the order of SMILES output."""
    # allHsExplicit may expose the implicit Hs of linker atoms to Smiles; the implicit Hs don't have names
    smiles_exh = Chem.MolToSmiles(mol, allHsExplicit=True)

    smiles_atom_output_order = mol.GetProp('_smilesAtomOutputOrder')
    delimiters = ['[', ']', ',']
    for delimiter in delimiters:
        smiles_atom_output_order = smiles_atom_output_order.replace(delimiter, ' ')
    smiles_output_order = [int(x) for x in smiles_atom_output_order.split()]

    atom_name = [mol.GetAtomWithIdx(atom_i).GetProp('atom_id') for atom_i in smiles_output_order]

    return smiles_exh, atom_name


def embed(mol: Chem.Mol, leaving_name: list[str]):
    """Remove atoms from the molecule based on the leaving_name set."""
    leaving_atoms = [atom for atom in mol.GetAtoms() if atom.GetProp('atom_id') in leaving_name]

    if not leaving_atoms: 
        return mol
    
    print(f"removing {leaving_name}...") 
    rwmol = Chem.RWMol(mol)
    # remove from the largest index to preserve smaller indexes
    for atom in sorted(leaving_atoms, key=lambda atom: atom.GetIdx(), reverse=True): 
        rwmol.RemoveAtom(atom.GetIdx())
    # linker atoms will have implicit Hs to compensate missing connections
    for atom in rwmol.GetAtoms():
        atom.UpdatePropertyCache()
    return rwmol.GetMol()


def deprotonate(mol, acidic_proton_loc: dict[str, int]) -> Chem.Mol:
    """Remove acidic protons from the molecule based on acidic_proton_loc"""
    # acidic_proton_loc is a mapping 
    # keys: smarts pattern of a fragment
    # value: the index (order in smarts) of the leaving proton

    # deprotonate all matched protons
    acidic_protons = []
    for smarts_pattern, idx in acidic_proton_loc.items():
        qmol = Chem.MolFromSmarts(smarts_pattern)
        acidic_protons.extend([mol.GetAtomWithIdx(match[idx]) for match in mol.GetSubstructMatches(qmol)])
    
    if not acidic_protons:
        return mol
    
    print(f"deprotonating {list(acidic_proton_loc.keys())}...") 
    rwmol = Chem.RWMol(mol)
    for atom in sorted(acidic_protons, key=lambda atom: atom.GetIdx(), reverse=True):
        rwmol.RemoveAtom(atom.GetIdx())
        neighbors = atom.GetNeighbors()
        if neighbors:
            neighbor_atom = rwmol.GetAtomWithIdx(neighbors[0].GetIdx())
            neighbor_atom.SetFormalCharge(neighbor_atom.GetFormalCharge() - 1)
    
    return rwmol.GetMol()


def is_chemical_element(symbol: str) -> bool:
    """Check if a string represents a valid chemical element."""
    try:
        return Chem.GetPeriodicTable().GetAtomicNumber(symbol) > 0
    # rdkit throws RuntimeError if invalid
    except RuntimeError:
        return False


def make_pretty_smiles(smi: str): 
    """Convert Smiles with allHsExplicit to pretty Smiles to be put on chem templates"""
    # collect the inside square brackets contents
    contents = set()
    inside_bracket = False
    for char in smi:
        if char == '[':
            inside_bracket = True
            content = ""
        elif char == ']' and inside_bracket:
            inside_bracket = False
            contents.add(content)
        elif inside_bracket:
            content += char

    for content in contents:
        # keep [H] for explicit Hs
        if content == 'H': 
            continue
        # drop H in the content to hide implicit Hs
        H_stripped = content.split('H')[0]
        # drop [ ] if the content is an uncharged element symbol
        if is_chemical_element(content) or is_chemical_element(H_stripped):
            smi = smi.replace(f"[{content}]", f"{H_stripped}" if 'H' in content else f"{content}")
    return smi


class ChemicalComponent:

    def __init__(self, rdkit_mol: Chem.Mol, resname: str, smiles_exh: str, atom_name: list[str], leaving_name: set[str]):
        self.rdkit_mol = rdkit_mol
        self.resname = resname
        self.smiles_exh = smiles_exh
        self.atom_name = atom_name
        self.leaving_name = leaving_name
        self.link_labels = {} # default to empty dict (free molecular form)
        self.parent = resname # default parent to itself
    
    @classmethod
    # requires gemmi
    def from_cif(cls, source_cif: str, resname: str):
        """Create ChemicalComponent from a chemical component dict file and a searchable residue name in file."""
        
        # Locate block by resname
        doc = gemmi.cif.read(source_cif)
        block = doc.find_block(resname)
        
        # Populate atom table
        atom_category = '_chem_comp_atom.'
        atom_attributes = ['atom_id', # atom names
                           'type_symbol', # atom elements
                           'charge', # (atomic) formal charges
                           'pdbx_leaving_atom_flag', # flags for leaving atoms after polymerization
                           ]
        atom_table = block.find(atom_category, atom_attributes)
        atom_cols = {attr: atom_table.find_column(f"{atom_category}{attr}") for attr in atom_attributes}

        # Summon rdkit atoms into empty RWMol
        rwmol = Chem.RWMol()
        atom_elements = atom_cols['type_symbol']

        for idx, element in enumerate(atom_elements):
            rdkit_atom = Chem.Atom(element)
            for attr in atom_attributes:
                rdkit_atom.SetProp(attr, atom_cols[attr][idx])
                # strip double quotes in names
                raw_name = rdkit_atom.GetProp('atom_id')
                rdkit_atom.SetProp('atom_id', raw_name.strip('"'))
            target_charge = atom_cols['charge'][idx]
            if target_charge!='0':
                rdkit_atom.SetFormalCharge(float(target_charge)) # this is usually int
            rwmol.AddAtom(rdkit_atom)

        # Map atom_id (atom names) with rdkit idx
        name_to_idx_mapping = {atom.GetProp('atom_id'): idx for (idx, atom) in enumerate(rwmol.GetAtoms())}

        # Populate bond table
        bond_category = '_chem_comp_bond.'
        bond_attributes = ['atom_id_1', # atom name 1
                           'atom_id_2', # atom name 2
                           'value_order', # bond order
                           ]
        bond_table = block.find(bond_category, bond_attributes)
        bond_cols = {attr: bond_table.find_column(f"{bond_category}{attr}") for attr in bond_attributes}

        # Connect atoms by bonds
        bond_type_mapping = {
            'SING': Chem.BondType.SINGLE,
            'DOUB': Chem.BondType.DOUBLE,
            'TRIP': Chem.BondType.TRIPLE,
            'AROM': Chem.BondType.AROMATIC
        }
        bond_types = bond_cols['value_order']

        for bond_i, bond_type in enumerate(bond_types):
            rwmol.AddBond(name_to_idx_mapping[bond_cols['atom_id_1'][bond_i].strip('"')], 
                          name_to_idx_mapping[bond_cols['atom_id_2'][bond_i].strip('"')], 
                          bond_type_mapping.get(bond_type, Chem.BondType.UNSPECIFIED))

        # Finish eidting mol    
        for atom in rwmol.GetAtoms():
            atom.UpdatePropertyCache()
        rdkit_mol = rwmol.GetMol()
            
        # Get Smiles with explicit Hs and ordered atom names
        smiles_exh, atom_name = get_smiles_with_atom_names(rdkit_mol)

        # Get set of leaving atom names
        leaving_name = {atom_cols['atom_id'][atom_i].strip('"')
                        for atom_i, flag in enumerate(atom_cols['pdbx_leaving_atom_flag']) 
                        if flag == 'Y'}
        print(f"leaving atoms from pdbx annotations: {leaving_name}")
        
        return cls(rdkit_mol, resname, smiles_exh, atom_name, leaving_name)


    def make_canonical(self):
        """Deprotonate acidic groups til the canonical (most deprotonated) state."""
        self.rdkit_mol = deprotonate(self.rdkit_mol, acidic_proton_loc = {'[H][O][PX4](=O)([O])[OX2]': 0})
        self.smiles_exh, self.atom_name = get_smiles_with_atom_names(self.rdkit_mol)

    def make_embedded(self):
        """Remove leaving atoms from the molecule."""
        self.rdkit_mol = embed(self.rdkit_mol, self.leaving_name)
        self.smiles_exh, self.atom_name = get_smiles_with_atom_names(self.rdkit_mol)

    def make_link_labels_from_names(self, name_to_label_mapping = {'P': '5-prime', "O3'": '3-prime'}):
        """Map atom names to link labels based on a given mapping."""
        link_labbels = {}
        for atom in self.rdkit_mol.GetAtoms():
            if atom.GetProp('atom_id') in name_to_label_mapping:
                if atom.GetNumImplicitHs() > 0:
                    name = atom.GetProp('atom_id')
                    link_labbels[str(self.atom_name.index(name))] = name_to_label_mapping[name]
        self.link_labels = link_labbels


# XXX make embedding variants by PATTERNS not names, assign/check link_labels and update smiles/idx
# XXX read from prepared? enumerate in stepwise? all protonation state variants, alter charge and update smiles/idx
# XXX group by parent: variants for "ambiguous"; tag synonym Amber and CCD resnames


if __name__ == '__main__':

    import json
    
    # """Download components.cif"""
    # import subprocess, sys
    # result = subprocess.run(["curl", "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif"], capture_output=True, text=True)
    # if result.returncode != 0:
    #    print(f"Unable to download components.cif from files.wwpdb.org")
    #    sys.exit(2)

    """Make chemical templates"""
    source_cif = 'components.cif'
    basenames = ['A', 'U', 'C', 'G', 'DA', 'DT', 'DC', 'DG']
    NA_ccs = []

    # nucleotide
    variant_dict = {
        "_": {}, # free nucleotide monophosphate
        "":  {'OP3', 'HOP3', "HO3'"}, # embedded nucleotide 
        "3": {'OP3', 'HOP3'}, # 3' end nucleotide 
        "5p":  {"HO3'"}, # 5' end nucleotide (extra phosphate than canonical X5)
        }

    for basename in basenames:
        for suffix in variant_dict:
            cc = ChemicalComponent.from_cif(source_cif, basename)
            cc.resname += suffix
            print(f"using CCD residue {basename} to construct {cc.resname}")

            cc.leaving_name = variant_dict[suffix]
            print(f"setting leaving atoms to {cc.leaving_name}")

            cc.make_canonical()
            cc.make_embedded()

            cc.smiles_exh = make_pretty_smiles(cc.smiles_exh)
            cc.make_link_labels_from_names()
            print(f"setting link labels to {cc.link_labels}")
            cc.parent = basename

            print(f"*** finish making {cc.resname} ***")
            NA_ccs.append(cc)

    # nucleoside (no phosphate)
    resname_mapping ={
        'A': 'ADN',
        'U': 'URI',
        'C': 'CTN',
        'G': 'GMP',
        'DA': '3D1',
        'DT': 'THM',
        'DC': 'DCZ',
        'DG': 'GNG'
        }
    variant_dict = {
        "N": {}, # free nucleoside 
        "5":  {"HO3'"}, # 5' end nucleoside (canonical X5 in Amber)
        # there shouldn't be a 3' end nucleoside variant (always has a phosphate)
        }

    for basename in basenames:
        for suffix in variant_dict:
            cc = ChemicalComponent.from_cif(source_cif, resname_mapping[basename])
            cc.resname = basename+suffix
            print(f"using CCD residue {resname_mapping[basename]} to construct {cc.resname}")

            cc.leaving_name = variant_dict[suffix]
            print(f"setting leaving atoms to {cc.leaving_name}")

            cc.make_canonical()
            cc.make_embedded()

            cc.smiles_exh = make_pretty_smiles(cc.smiles_exh)
            cc.make_link_labels_from_names()
            print(f"setting link labels to {cc.link_labels}")
            cc.parent = basename

            print(f"*** finish making {cc.resname} ***")
            NA_ccs.append(cc)


    """Export chemical templates"""
    residue_templates = {}

    for cc in NA_ccs:
        residue_templates[cc.resname] = {
            "smiles": cc.smiles_exh,
            "atom_name": cc.resname+".atom_names",
        }
        if cc.link_labels:
            residue_templates[cc.resname]["link_labels"] = cc.resname+".link_labels"
        else:
            residue_templates[cc.resname]["link_labels"] = {}

    data_to_export = {"residue_templates": residue_templates}
    json_str = json.dumps(data_to_export, indent = 4)

    # format link_labels and atom_name to one line
    for cc in NA_ccs:
        single_line_atom_name = json.dumps(cc.atom_name, separators=(', ', ': '))
        json_str = json_str.replace(json.dumps(data_to_export["residue_templates"][cc.resname]["atom_name"], indent = 4), single_line_atom_name)
        if cc.link_labels:
            single_line_link_labels = json.dumps(cc.link_labels, separators=(', ', ': '))
            json_str = json_str.replace(json.dumps(data_to_export["residue_templates"][cc.resname]["link_labels"], indent = 4), single_line_link_labels)

    new_template_json = 'NA_residue_templates.json'
    with open(new_template_json, 'w') as f:
        f.write(json_str)


    """Export ambiguous mapping"""
    ambiguous_dict = {basename:[] for basename in basenames}
    for cc in NA_ccs:
        ambiguous_dict[cc.parent].append(cc.resname)

    data_to_export = {"ambiguous": {basename:basename+'.resnames' for basename in basenames}}
    json_str = json.dumps(data_to_export, indent = 4)

    # format ambiguous resnames to one line
    for basename in ambiguous_dict:
        single_line_resnames = json.dumps(ambiguous_dict[basename], separators=(', ', ': '))
        json_str = json_str.replace(json.dumps(data_to_export["ambiguous"][basename], indent = 4), single_line_resnames)

    new_ambiguous_json = 'NA_ambiguous.json'
    with open(new_ambiguous_json, 'w') as f:
        f.write(json_str)

