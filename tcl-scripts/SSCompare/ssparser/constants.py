import pathlib
import os

"""
Contains constants for parsing and plotting each ABeta monomer
"""

file_cwd = pathlib.Path(__file__).parent.resolve()
DATASET_DIR = os.path.join(file_cwd, "../datasets/abeta-fibril-csv/")

ABETA_LONG_SEQUENCE = ['ASP', 'ALA', 'GLU', 'PHE', 'ARG', 'HIS', 'ASP', 'SER', 'GLY', 'TYR', 'GLU', 'VAL', 'HIS', 'HIS', 'GLN', 'LYS', 'LEU', 'VAL', 'PHE', 'PHE', 'ALA', 'GLU', 'ASP', 'VAL', 'GLY', 'SER', 'ASN', 'LYS', 'GLY', 'ALA', 'ILE', 'ILE', 'GLY', 'LEU', 'MET', 'VAL', 'GLY', 'GLY', 'VAL', 'VAL', 'ILE', 'ALA']
ABETA_LONG_STRING = "".join(ABETA_LONG_SEQUENCE)
ABETA_LONG_RESIDUE_COUNT = 42

ABETA_SHORT_SEQUENCE = ['ASP', 'ALA', 'GLU', 'PHE', 'ARG', 'HIS', 'ASP', 'SER', 'GLY', 'TYR', 'GLU', 'VAL', 'HIS', 'HIS', 'GLN', 'LYS', 'LEU', 'VAL', 'PHE', 'PHE', 'ALA', 'GLU', 'ASP', 'VAL', 'GLY', 'SER', 'ASN', 'LYS', 'GLY', 'ALA', 'ILE', 'ILE', 'GLY', 'LEU', 'MET', 'VAL', 'GLY', 'GLY', 'VAL', 'VAL']
ABETA_SHORT_STRING = "".join(ABETA_SHORT_SEQUENCE)
ABETA_SHORT_RESIDUE_COUNT = 40

