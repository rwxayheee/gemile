import re
from collections import namedtuple

def parse_bld_file(bld_filename): 
    with open(bld_filename, 'r') as f:
        bld_lines = f.readlines()
    bld_blocks = []
    bld_block = []
    for line in bld_lines:
        if line.startswith("f_m_ct"):
            if len(bld_block) > 0: 
                bld_blocks.append(bld_block)
            bld_block = [line]
        else:
            bld_block.append(line)
    if len(bld_block) > 0: 
        bld_blocks.append(bld_block)
    return bld_blocks

def parse_bld_block(bld_block): 
    # assumes first subsection is title
    title_block = []
    line_number = 1
    max_line_number = len(bld_block) - 2
    line = bld_block[line_number]
    while "{" not in line and line_number<max_line_number: 
        title_block.append(line.strip())
        line_number += 1
        line = bld_block[line_number]
    title_dict = parse_mapping(title_block)

    # assumes atom section starts with m_atom[<number of atoms>]
    # followed by bond section that starts with m_bond[<number of bonds>]
    atom_block = []
    bond_block = []
    while "m_atom" not in line and line_number<max_line_number:
        line_number += 1
        line = bld_block[line_number]
    line_number += 1
    line = bld_block[line_number]
    while "}" not in line and line_number<max_line_number:
        atom_block.append(line.strip())
        line_number += 1
        line = bld_block[line_number]
    while "m_bond" not in line and line_number<max_line_number:
        line_number += 1
        line = bld_block[line_number]
    line_number += 1
    line = bld_block[line_number]
    while "}" not in line and line_number<max_line_number:
        bond_block.append(line.strip())
        line_number += 1
        line = bld_block[line_number]

    if atom_block: 
        atom_list = parse_table(atom_block)
    else:
        return {"title": title_dict, "atom_list": [], "bond_list": []}

    if bond_block: 
        bond_list = parse_table(bond_block)
    else:
        return {"title": title_dict, "atom_list": atom_list, "bond_list": []}

    return {"title": title_dict, "atom_list": atom_list, "bond_list": bond_list}

def parse_mapping(input_block): 
    parsed_dict = {}
    keys = []
    is_val_row = False
    for line in input_block: 
        line = line.strip()
        if ":::" in line:
            is_val_row = True
        elif not is_val_row: 
            keys.append(line)
        elif keys:
            key = keys.pop(0)
            parsed_dict[key] = line
    return parsed_dict

def parse_table(input_block): 
    obj_fields = []
    obj_list = []
    is_obj_row = False
    bld_obj = None
    for line in input_block:
        line = line.strip()
        if ":::" in line: 
            bld_obj = namedtuple('bld_obj', obj_fields)
            is_obj_row = True
        elif not is_obj_row: 
            if " " not in line: 
                obj_fields.append(line)
            else: 
                obj_fields.append("index")
        elif obj_fields and is_obj_row: 
            obj_list.append(bld_obj._make(parse_line_to_tuple(line)))
    return obj_list

def parse_line_to_tuple(line):
    """
    parse line
    '1 18 -0.095400 1.152600 1.546600 1 a A 70 " DA " " OP3" "  c2" 8 -1'
    into tuple
    ("1", "18", "-0.095400", "1.152600", "1.546600", "1", "a", "A", "70", " DA ", " OP3", "  c2", "8", "-1")
    """
    matches = re.findall(r'("[^"]*"|\S+)', line)
    return tuple(match.strip('"') if match.startswith('"') and match.endswith('"') else match for match in matches)

bld_blocks = parse_bld_file("/Users/amyhe/Desktop/0_forks/template_mod/res/dna2.bld")
bld_block_1 = bld_blocks[1]
parse_bld_block(bld_block_1)


