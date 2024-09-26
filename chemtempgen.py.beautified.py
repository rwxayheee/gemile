from rdkit import Chem
import gemmi

def get_smiles_with_atom_names(mol: Chem.Mol) -> tuple[str, list[str]]:
    """Generate SMILES with atom names in the order of SMILES output."""
    smiles_exh = Chem.MolToSmiles(mol, allHsExplicit=True)
    smiles_atom_output_order = mol.GetProp('_smilesAtomOutputOrder')
    # Clean and convert output order string to list of indices
    delimiters = ['[', ']', ',']
    for delimiter in delimiters:
        smiles_atom_output_order = smiles_atom_output_order.replace(delimiter, ' ')
    smiles_output_order = [int(x) for x in smiles_atom_output_order.split()]
    
    # Re-index atom names based on the SMILES output order
    atom_name = [mol.GetAtomWithIdx(atom_i).GetProp('atom_id') for atom_i in smiles_output_order]
    return smiles_exh, atom_name

def embed(mol: Chem.Mol, leaving_name: list[str]) -> Chem.Mol:
    """Remove atoms by leaving_name list from the molecule."""
    leaving_atoms = [atom for atom in mol.GetAtoms() if atom.GetProp('atom_id') in leaving_name]
    if not leaving_atoms:
        print(f"No leaving atoms for {leaving_name}...")
        return mol

    rwmol = Chem.RWMol(mol)
    for atom in sorted(leaving_atoms, key=lambda atom: atom.GetIdx(), reverse=True):
        rwmol.RemoveAtom(atom.GetIdx())
    return rwmol.GetMol()

def deprotonate(mol: Chem.Mol, acidic_proton_loc: dict[str, int]) -> Chem.Mol:
    """Remove acidic protons from the molecule based on patterns."""
    acidic_protons = []
    for smarts_pattern, idx in acidic_proton_loc.items():
        qmol = Chem.MolFromSmarts(smarts_pattern)
        acidic_protons.extend([mol.GetAtomWithIdx(match[idx]) for match in mol.GetSubstructMatches(qmol)])
    
    if not acidic_protons:
        print(f"No acidic protons for {list(acidic_proton_loc.keys())}")
        return mol

    rwmol = Chem.RWMol(mol)
    for atom in sorted(acidic_protons, key=lambda atom: atom.GetIdx(), reverse=True):
        rwmol.RemoveAtom(atom.GetIdx())
        neighbors = atom.GetNeighbors()
        if neighbors:
            neighbor_atom = rwmol.GetAtomWithIdx(neighbors[0].GetIdx())
            neighbor_atom.SetFormalCharge(neighbor_atom.GetFormalCharge() - 1)
    return rwmol.GetMol()

def is_chemical_element(symbol: str) -> bool:
    """Check if the symbol represents a valid chemical element."""
    try:
        return Chem.GetPeriodicTable().GetAtomicNumber(symbol) > 0
    except RuntimeError:
        return False

def make_pretty_smiles(smi: str) -> str:
    """Remove square brackets from SMILES representation where not needed."""
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
        if content == 'H':
            continue
        H_stripped = content.split('H')[0]
        if is_chemical_element(content) or is_chemical_element(H_stripped):
            smi = smi.replace(f"[{content}]", f"{H_stripped}" if 'H' in content else f"{content}")
    return smi

class ChemicalComponent:
    """Class to represent and manipulate a chemical component."""
    def __init__(self, rdkit_mol: Chem.Mol, resname: str, smiles_exh: str, atom_name: list[str], leaving_name: set[str]):
        self.rdkit_mol = rdkit_mol
        self.resname = resname
        self.smiles_exh = smiles_exh
        self.atom_name = atom_name
        self.leaving_name = leaving_name
        self.link_labels = {}
        self.parent = resname
    
    @classmethod
    def from_cif(cls, source_cif: str, CCD_resname: str):
        """Create ChemicalComponent from a CIF file and residue name."""
        doc = gemmi.cif.read(source_cif)
        block = doc.find_block(CCD_resname)
        
        atom_category = '_chem_comp_atom.'
        atom_attributes = ['atom_id', 'type_symbol', 'charge', 'pdbx_leaving_atom_flag']
        atom_table = block.find(atom_category, atom_attributes)
        atom_cols = {attr: atom_table.find_column(f"{atom_category}{attr}") for attr in atom_attributes}

        rwmol = Chem.RWMol()
        for idx, element in enumerate(atom_cols['type_symbol']):
            rdkit_atom = Chem.Atom(element)
            for attr in atom_attributes:
                rdkit_atom.SetProp(attr, atom_cols[attr][idx].strip('"'))
            if atom_cols['charge'][idx] != '0':
                rdkit_atom.SetFormalCharge(int(atom_cols['charge'][idx]))
            rwmol.AddAtom(rdkit_atom)

        name_to_idx_mapping = {atom.GetProp('atom_id'): idx for idx, atom in enumerate(rwmol.GetAtoms())}
        bond_category = '_chem_comp_bond.'
        bond_attributes = ['atom_id_1', 'atom_id_2', 'value_order']
        bond_table = block.find(bond_category, bond_attributes)
        bond_cols = {attr: bond_table.find_column(f"{bond_category}{attr}") for attr in bond_attributes}

        bond_type_mapping = {'SING': Chem.BondType.SINGLE, 'DOUB': Chem.BondType.DOUBLE, 'TRIP': Chem.BondType.TRIPLE, 'AROM': Chem.BondType.AROMATIC}
        for bond_i, bond_type in enumerate(bond_cols['value_order']):
            rwmol.AddBond(name_to_idx_mapping[bond_cols['atom_id_1'][bond_i].strip('"')],
                          name_to_idx_mapping[bond_cols['atom_id_2'][bond_i].strip('"')],
                          bond_type_mapping.get(bond_type, Chem.BondType.UNSPECIFIED))
            
        for atom in rwmol.GetAtoms():
            atom.UpdatePropertyCache()

        smiles_exh, atom_name = get_smiles_with_atom_names(rwmol)
        leaving_name = {atom_cols['atom_id'][i].strip('"') for i, flag in enumerate(atom_cols['pdbx_leaving_atom_flag']) if flag == 'Y'}
        
        return cls(rwmol.GetMol(), CCD_resname, smiles_exh, atom_name, leaving_name)
    
    def make_canonical(self):
        """Deprotonate the molecule to make it canonical."""
        self.rdkit_mol = deprotonate(self.rdkit_mol, acidic_proton_loc = {'[H][O][PX4](=O)([O])[OX2]': 0})
        self.smiles_exh, self.atom_name = get_smiles_with_atom_names(self.rdkit_mol)

    def make_embedded(self):
        """Embed leaving atoms into the molecule."""
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
    

if __name__ == '__main__':
    source_cif = '/Users/amyhe/Desktop/7_Mk_for_NA/0_ccd/components.cif'
    basenames = ['A', 'U', 'C', 'G', 'DA', 'DT', 'DC', 'DG']
    NA_ccs = []

    variant_dict = {
        "_": {},  # free nucleotide
        "":  {'OP3', 'HOP3', "HO3'"},  # embedded nucleotide
        "5": {"HO3'"}  # 5' end nucleotide
    }

    for basename in basenames:
        for suffix, leaves in variant_dict.items():
            cc = ChemicalComponent.from_cif(source_cif, basename)
            cc.leaving_name = leaves

            cc.make_canonical()
            cc.make_embedded()

            cc.smiles_exh = make_pretty_smiles(cc.smiles_exh)
            cc.make_link_labels_from_names()
            cc.resname += suffix
            cc.parent = basename

            NA_ccs.append(cc)

    # Nucleoside (no phosphate)
    resname_mapping = {
        'A': 'ADN', 'U': 'URI', 'C': 'CTN', 'G': 'GMP',
        'DA': '3D1', 'DT': 'THM', 'DC': 'DCZ', 'DG': 'GNG'
    }
    variant_dict = {
        "N": {},  # free nucleoside
        "3": {'OP3', 'HOP3'},  # 3' end nucleoside
    }

    for basename in basenames:
        for suffix, leaves in variant_dict.items():
            cc = ChemicalComponent.from_cif(source_cif, resname_mapping[basename])
            cc.leaving_name = leaves

            cc.make_canonical()
            cc.make_embedded()

            cc.smiles_exh = make_pretty_smiles(cc.smiles_exh)
            cc.make_link_labels_from_names()
            cc.resname += suffix
            cc.parent = basename

            NA_ccs.append(cc)
