
# gemmi import
import time
start_time = time.time()
import gemmi
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit import RDLogger
logger = RDLogger.logger()
logger.setLevel(RDLogger.CRITICAL)
import sys, logging, chemtempgen
logger = logging.getLogger('chemtempgen')
logger.setLevel(logging.WARNING)
handler = logging.StreamHandler(sys.stdout)

# prody import
import argparse
import json
import math
from os import linesep as os_linesep
import pathlib
import sys
import numpy as np
from meeko.reactive import atom_name_to_molsetup_index, assign_reactive_types_by_index
from meeko import PDBQTMolecule
from meeko import RDKitMolCreate
from meeko import MoleculePreparation
from meeko import MoleculeSetup
from meeko import ResidueChemTemplates
from meeko import PDBQTWriterLegacy
from meeko import LinkedRDKitChorizo
from meeko import LinkedRDKitChorizoEncoder
from meeko import ChorizoCreationError
from meeko import reactive_typer
from meeko import get_reactive_config
from meeko import gridbox
from meeko import __file__ as pkg_init_path
from rdkit import Chem
import prody

SUPPORTED_PRODY_FORMATS = {"pdb": prody.parsePDB, "cif": prody.parseMMCIF}

read_with_prody = "/Users/amyhe/Downloads/1d3g.pdb"
ext = pathlib.Path(read_with_prody).suffix[1:].lower()
parser = SUPPORTED_PRODY_FORMATS[ext]
input_obj = parser(read_with_prody, altloc="all")

templates_fn = "/Users/amyhe/Desktop/0_forks/Meeko/meeko/data/residue_chem_templates.json"
with open(templates_fn) as f:
    res_chem_templates = json.load(f)
templates = ResidueChemTemplates.from_dict(res_chem_templates)

mk_prep = MoleculePreparation()

chorizo = LinkedRDKitChorizo.from_prody(
    input_obj,
    templates,
    mk_prep,
    set_template=None,
    residues_to_delete=None,
    allow_bad_res=False,
    bonds_to_delete=None,
    blunt_ends=None,
    wanted_altloc=None,
    default_altloc="A",
)

