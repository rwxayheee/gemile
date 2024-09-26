from rdkit import Chem
import gemmi

def get_smiles_with_atom_names(mol):
    # Generate SMILES with explicit hydrogens
    smiles_exh = Chem.MolToSmiles(mol, allHsExplicit=True)

    # Get SMILES output order from molecule properties
    smiles_atom_output_order = mol.GetProp('_smilesAtomOutputOrder')
        
    # Clean the order string to get a list of indices
    delimiters = ['[', ']', ',']
    for delimiter in delimiters:
        smiles_atom_output_order = smiles_atom_output_order.replace(delimiter, ' ')
    smiles_output_order = [int(x) for x in smiles_atom_output_order.split()]

    # Re-index atom names based on the SMILES output order
    atom_name = []
    for atom_i in smiles_output_order:
        atom_name.append(mol.GetAtomWithIdx(atom_i).GetProp('atom_id'))  # Get atom name using 'atom_id' property

    return smiles_exh, atom_name


def embed(mol, leaving_name: list[str]):

    init_atoms = mol.GetAtoms()
    leaving_atoms = [atom for atom in init_atoms if atom.GetProp('atom_id') in leaving_name]

    if len(leaving_atoms) == 0: 
        print(f"no leaving atoms for {leaving_name}...")
        return mol
    else:
        sorted_leaving_atoms = sorted(leaving_atoms, key=lambda atom: atom.GetIdx(), reverse=True)

    rwmol = Chem.RWMol(mol)
    for atom in sorted_leaving_atoms:
        rwmol.RemoveAtom(atom.GetIdx())
    
    return rwmol.GetMol()


def deprotonate(mol, acidic_proton_loc = {'[H][O][PX4](=O)([O])[OX2]': 0}):

    acidic_protons = []
    for smarts_pattern in acidic_proton_loc:
        qmol = Chem.MolFromSmarts(smarts_pattern)
        matches = mol.GetSubstructMatches(qmol)
        for match in matches:
            acidic_protons.append(mol.GetAtomWithIdx(match[acidic_proton_loc[smarts_pattern]]))
    
    if len(acidic_protons) == 0:
        print(f"no acidic protons for {acidic_proton_loc.keys()}")
        return mol
    else:
        sorted_acidic_protons = sorted(acidic_protons, key=lambda atom: atom.GetIdx(), reverse=True)
    
    rwmol = Chem.RWMol(mol)
    for atom in sorted_acidic_protons:
        rwmol.RemoveAtom(atom.GetIdx())
        neighbors = atom.GetNeighbors()
        if neighbors:
            neighbor_idx = neighbors[0].GetIdx()
            neighbor_atom = rwmol.GetAtomWithIdx(neighbor_idx)
            neighbor_atom_charge = neighbor_atom.GetFormalCharge()
            neighbor_atom.SetFormalCharge(neighbor_atom_charge-1)
    
    return rwmol.GetMol()


def is_chemical_element(symbol):
    # Get the periodic table from RDKit
    periodic_table = Chem.GetPeriodicTable()
    
    # Try to get the atomic number of the element using the symbol
    try:
        atomic_number = periodic_table.GetAtomicNumber(symbol)
        # If atomic_number is greater than 0, it's a valid element
        return atomic_number > 0
    except:
        # If symbol is not a valid element, GetAtomicNumber raises a ValueError
        return False


def make_pretty_smiles(smi: str): 
    # remove square brackets if atom does not have a charge;
    # make H implicit

    contents = set()

    for char in smi:
        if char == '[':
            inside_bracket = True
            content = ""  # Start a new element
        elif char == ']' and inside_bracket:
            inside_bracket = False
            contents.add(content)  # Add the completed element
        elif inside_bracket:
            content += char  # Collect characters inside brackets

    for content in contents:
        if content=='H':
            pass
        elif is_chemical_element(content):
            smi = smi.replace(f"[{content}]", f"{content}")
        elif 'H' in content:
            H_stripped = content.split('H')[0]
            if is_chemical_element(H_stripped):
                smi = smi.replace(f"[{content}]", f"{H_stripped}")

    return smi


def make_link_labels_from_names(names: list[str], name_to_label_mapping = {'P': '5-prime', "O3'": '3-prime'}):
    # find link labels by names
    link_labels = {}
    for name in name_to_label_mapping: 
        if name in names:
            link_labels[str(names.index(name))] = name_to_label_mapping[name]

    return link_labels


class ChemicalComponent:

    def __init__(self, rdkit_mol: Chem.Mol, \
                 resname: str, smiles_exh: str, atom_name: list[str], leaving_name: set[str]):
        self.rdkit_mol = rdkit_mol
        self.resname = resname
        self.smiles_exh = smiles_exh
        self.atom_name = atom_name
        self.leaving_name = leaving_name
        self.link_labels = {} # default to empty dict (free molecular form)
        self.parent = resname # default parent to itself
    
    @classmethod # requires gemmi
    def from_cif(cls, source_cif: str, CCD_resname: str):
        resname = CCD_resname # default resname to CCD_resname

        # Locate block
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
                raw_name = rdkit_atom.GetProp('atom_id')
                rdkit_atom.SetProp('atom_id', raw_name.strip('"'))
            
            target_charge = rdkit_atom.GetProp('charge')
            if target_charge!='0':
                rdkit_atom.SetFormalCharge(int(target_charge))
            
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
            if bond_type not in bond_type_mapping:
                raise ValueError(f"value_order {bond_type} from input is not in {bond_type_mapping.keys()}")
            else:
                rwmol.AddBond(name_to_idx_mapping[bond_cols['atom_id_1'][bond_i].strip('"')], \
                            name_to_idx_mapping[bond_cols['atom_id_2'][bond_i].strip('"')], \
                            bond_type_mapping[bond_type])

        smiles_exh, atom_name = get_smiles_with_atom_names(rwmol)
        for atom in rwmol.GetAtoms():
            atom.UpdatePropertyCache()
        rdkit_mol = rwmol.GetMol()

        leaving_name = {atom_cols['atom_id'][atom_i].strip('"')
                        for atom_i, atom_leaving_flag in enumerate(atom_cols['pdbx_leaving_atom_flag']) 
                        if atom_leaving_flag == 'Y'}
        print(f"leaving atoms from pdbx annotations: {leaving_name}")
            

        return cls(rdkit_mol, resname, smiles_exh, atom_name, leaving_name)


    def make_embedded(self): 
        newmol = embed(self.rdkit_mol, self.leaving_name)
        print(f"Current leaving atoms: {self.leaving_name}")
        
        self.smiles_exh, self.atom_name = get_smiles_with_atom_names(newmol)
        self.rdkit_mol = newmol
    

    def make_canonical(self): 
        newmol = deprotonate(self.rdkit_mol)

        self.smiles_exh, self.atom_name = get_smiles_with_atom_names(newmol)
        self.rdkit_mol = newmol


# XXX make embedding variants, assign/check link_labels and update smiles/idx
# XXX make? read? protonation variants, alter charge and update smiles/idx
# XXX group by parent: variants for "ambiguous"


if __name__ == '__main__':

    source_cif = '/Users/amyhe/Desktop/7_Mk_for_NA/0_ccd/components.cif'
    basenames = ['A', 'U', 'C', 'G', 'DA', 'DT', 'DC', 'DG']
    NA_ccs = []

    # nucleotide
    variant_dict = {
    "_": {}, # free nucleotide
    "":  {'OP3', 'HOP3', "HO3"}, # embedded nucleotide
    "5": {"HO3'"} # 5' end nucleotide
    }

    for basename in basenames:
        for suffix in variant_dict:
            cc = ChemicalComponent.from_cif(source_cif, basename)
            cc.leaving_name = variant_dict[suffix]

            cc.make_canonical()
            cc.make_embedded()

            cc.smiles_exh = make_pretty_smiles(cc.smiles_exh)
            cc.link_labels = make_link_labels_from_names(cc.atom_name)
            cc.resname += suffix
            cc.parent = basename

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
    ''
    }
    variant_dict = {
    "N": {}, # free nucleoside
    "3":  {'OP3', 'HOP3'}, # 3' end nucleoside
    }

    for basename in basenames:
        for suffix in variant_dict:
            cc = ChemicalComponent.from_cif(source_cif, resname_mapping[basename])
            cc.leaving_name = variant_dict[suffix]

            cc.make_canonical()
            cc.make_embedded()

            cc.smiles_exh = make_pretty_smiles(cc.smiles_exh)
            cc.link_labels = make_link_labels_from_names(cc.atom_name)
            cc.resname += suffix
            cc.parent = basename

            NA_ccs.append(cc)
