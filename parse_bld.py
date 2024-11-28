import re
from collections import namedtuple

def parse_bld_file(bld_filename): 
    bld_blocks = []
    with open(bld_filename, 'r') as f:
        bld_lines = (line.strip() for line in f) 
        current_block = []
        for line in bld_lines:
            if line.startswith("f_m_ct"):
                if current_block: 
                    bld_blocks.append(current_block)
                current_block = [line] 
            else:
                current_block.append(line)
        if current_block:
            bld_blocks.append(current_block)
    return bld_blocks

def extract_block(start_marker, end_marker, bld_block, start_index):
    """Extracts lines between a start marker and end marker."""
    block = []
    line_number = start_index
    max_line_number = len(bld_block) - 1
    # Find the start marker
    while start_marker not in bld_block[line_number] and line_number < max_line_number:
        line_number += 1
    # Skip the line with the start marker
    line_number += 1
    # Extract lines until the end marker
    while line_number < max_line_number and end_marker not in bld_block[line_number]:
        block.append(bld_block[line_number])
        line_number += 1
    return block, line_number

def extract_title(bld_block):
    title_block = []
    line_number = 1  # Start after the first line
    max_line_number = len(bld_block) - 1
    while line_number < max_line_number and "{" not in bld_block[line_number]:
        title_block.append(bld_block[line_number])
        line_number += 1
    return title_block, line_number

def parse_bld_block(bld_block):
    title_block, line_number = extract_title(bld_block)
    atom_block, line_number = extract_block("m_atom", "}", bld_block, line_number)
    bond_block, _ = extract_block("m_bond", "}", bld_block, line_number)

    return {
        "title_block": parse_mapping(title_block),
        "atom_block": parse_table(atom_block),
        "bond_block": parse_table(bond_block),
    }

def parse_mapping(input_block): 
    parsed_dict = {}
    keys = []
    is_val_row = False
    for line in input_block: 
        line = line
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
        line = line
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

bld_blocks = parse_bld_file("/Users/amyhe/Desktop/0_forks/template_mod/res/dna.bld")
bld_block_1 = bld_blocks[1]
parse_bld_block(bld_block_1)

