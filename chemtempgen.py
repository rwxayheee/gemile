from rdkit import Chem
import gemmi

class ChemicalComponent:

    def __init__(self, mol_from_gemmi: Chem.Mol, \
                 resname: str, smiles_exh: str, atom_name: list[str], leaving_atom_idx: list[int]):
        self.mol_from_gemmi = mol_from_gemmi
        self.resname = resname
        self.smiles_exh = smiles_exh
        self.atom_name = atom_name
        self.leaving_atom_idx = leaving_atom_idx
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
                           'charge', # formal? charges
                           'pdbx_leaving_atom_flag', # flags for leaving atoms after polymerization
                           ]
        atom_table = block.find(atom_category, atom_attributes)


        # Summon rdkit atoms into empty RWMol
        mol = Chem.RWMol()
        atom_elements = atom_table.find_column(f"{atom_category}type_symbol")

        for idx, element in enumerate(atom_elements):
            rdkit_atom = Chem.Atom(element)
            
            for attribute in atom_attributes:
                rdkit_atom.SetProp(attribute, atom_table.find_column(attribute)[idx])
            
            target_charge = rdkit_atom.GetProp('charge')
            if target_charge!='0':
                rdkit_atom.SetFormalCharge(int(target_charge))
            
            mol.AddAtom(rdkit_atom)


        # Map atom_id (atom names) with rdkit idx
        name_to_idx_mapping = {atom.GetProp('atom_id'): idx for (idx, atom) in enumerate(mol.GetAtoms())}
                
        # Populate bond table
        bond_category = '_chem_comp_bond.'
        bond_attributes = ['atom_id_1', # atom name 1
                           'atom_id_2', # atom name 2
                           'value_order', # bond order
                           ]
        bond_table = block.find(bond_category, bond_attributes)


        # Connect atoms by bonds
        bond_type_mapping = {
            'SING': Chem.BondType.SINGLE,
            'DOUB': Chem.BondType.DOUBLE,
            'TRIP': Chem.BondType.TRIPLE,
            'AROM': Chem.BondType.AROMATIC
        }
        bond_types = bond_table.find_column(f"{bond_category}value_order")

        for bond_i, bond_type in enumerate(bond_types):
            if bond_type not in bond_type_mapping:
                raise ValueError(f"value_order {bond_type} from input is not in {bond_type_mapping.keys()}")
            else:
                mol.AddBond(name_to_idx_mapping[bond_table.find_column(f"{bond_category}atom_id_1")[bond_i]], \
                            name_to_idx_mapping[bond_table.find_column(f"{bond_category}atom_id_2")[bond_i]], \
                            bond_type_mapping[bond_type])


        # Export Smiles with explicit Hs and mol    
        smiles_exh = Chem.MolToSmiles(mol, allHsExplicit=True)
        mol_from_gemmi = mol.GetMol()
        

        # Get Smiles output order
        smilesAtomOutputOrder = mol.GetProp('_smilesAtomOutputOrder')
        delimiters = ['[', ']', ',']
        for delimiter in delimiters:
            smilesAtomOutputOrder = smilesAtomOutputOrder.replace(delimiter, ' ')
        smiles_output_order = [int(x) for x in smilesAtomOutputOrder.split()]

        # Re-index atom names and leaving atom idx
        atom_name = []
        leaving_atom_idx = []
        for idx, atom_i in enumerate(smiles_output_order):
            atom_name_raw = atom_table.find_column(f"{atom_category}atom_id")[atom_i] 
            # could be like '"H5\'"', '"H5\'\'"', which have a pair of excess double quotes
            atom_name.append(atom_name_raw.strip('"'))

            leaving_atom_flag = atom_table.find_column(f"{atom_category}pdbx_leaving_atom_flag")[atom_i]
            if leaving_atom_flag == 'Y':
                leaving_atom_idx.append(idx)


        return cls(mol_from_gemmi, resname, smiles_exh, atom_name, leaving_atom_idx)


# XXX make embedding variants, assign/check link_labels and update smiles/idx
# XXX make? read? protonation variants, alter charge and update smiles/idx
# XXX group by parent: variants for "ambiguous"