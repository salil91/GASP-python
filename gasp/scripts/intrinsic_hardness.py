#!/usr/bin/env python
# coding: utf-8

"""
Intrinsic Harness from GASP results:

This script reads the structures from the POSCAR files output by GASP
and calculates the intrinsic hardness of the structures based on
the modified electonegativity (EN) and bond strength (BS) models.
Output file: ./hard_data

Usage: python intrinsic_hardness.py /path/to/run_data/file
POSCAR files must be in the same directory as the run_data file
"""

# Import libraries and read files
from pymatgen.core import Structure
from pymatgen.analysis.local_env import CrystalNN
import pandas as pd
import numpy as np
import math
import os
import sys
import warnings

warnings.filterwarnings('ignore')


# Read run_data
run_data_path = sys.argv[1]
with open(run_data_path, 'r') as f_gasp:
    lines = f_gasp.readlines()

# Prepare output file
run_data_dir = os.path.dirname(run_data_path)
with open(os.path.join(run_data_dir,"hard_data"), 'w') as f_hard:
    for line_number, line in enumerate(lines[:4]):
        if line_number == 2:
            line = line.strip()
            f_hard.writelines(line)
            f_hard.writelines("\t\t H_EN\t\t H_BS\t\t H_C")
        else:
            f_hard.writelines(line)
    f_hard.writelines("\n")


# Look up electronegativities
#df_EN = pd.read_excel("EN_lookup.xlsx", header = None, names = range(1,13))
#df_EN.index = np.arange(1, len(df_EN)+1)
d_EN = {
1: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
2: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
3: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
4: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
5: {1: float('nan'),  2: 1.641,  3: 1.641,  4: 1.641,  5: float('nan'),  6: 1.415,  7: 1.415,  8: float('nan'),  9: 1.415,  10: float('nan'),  11: float('nan'),  12: float('nan')},
6: {1: float('nan'),  2: 2.619,  3: 2.619,  4: 2.5,  5: float('nan'),  6: 2.347,  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
7: {1: float('nan'),  2: float('nan'),  3: 3.565,  4: 3.437,  5: float('nan'),  6: 3.085,  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
8: {1: float('nan'),  2: float('nan'),  3: 4.561,  4: 4.375,  5: float('nan'),  6: 3.75,  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
9: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
10: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
11: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
12: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
13: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
14: {1: float('nan'),  2: 1.645,  3: 1.645,  4: 1.645,  5: float('nan'),  6: 1.565,  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
15: {1: float('nan'),  2: float('nan'),  3: 2.187,  4: 2.187,  5: float('nan'),  6: 2.074,  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
16: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
17: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
18: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
19: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
20: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
21: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
22: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
23: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
24: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: 0.659,  8: float('nan'),  9: 0.659,  10: float('nan'),  11: float('nan'),  12: 0.659},
25: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
26: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
27: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
28: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
29: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
30: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
31: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
32: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
33: {1: float('nan'),  2: float('nan'),  3: 2.039,  4: 2.039,  5: float('nan'),  6: 1.925,  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
34: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
35: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
36: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
37: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
38: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
39: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
40: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
41: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
42: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: 2.171,  11: float('nan'),  12: float('nan')},
43: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
44: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
45: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
46: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
47: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
48: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
49: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
50: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
51: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
52: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
53: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
54: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
55: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
56: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
57: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
58: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
59: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
60: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
61: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
62: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
63: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
64: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
65: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
66: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
67: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
68: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
69: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
70: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
71: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
72: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
73: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
74: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: 2.171,  11: float('nan'),  12: float('nan')},
75: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: 0.71,  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
76: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: 1.504,  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')}
}

# Perform calculations on each POSCAR file
# Each line of the run_data file corresponds to a particular organism, and has it's own POSCAR file.
for line_number, line in enumerate(lines[4:]):
    id = line.split()[0]

    # Read POSCAR corresponding to id
    poscar = os.path.join(run_data_dir,str("POSCAR." + id))
    with open(poscar, "r") as f_poscar:
        formula = f_poscar.readlines()[0]
    #print("\nOrganism", id, "\nFormula =", formula.strip())

    # Read the structure and calculate unit cell volume
    structure = Structure.from_file(poscar)
    vol = structure.volume

    # Find bonds and calculate bond length
    bonds = []
    for site_index, atom in enumerate(structure):
        nn_object = CrystalNN()
        neighbors = nn_object.get_nn_info(structure, site_index)
        CN1 = nn_object.get_cn(structure, site_index)
        for neighbor in neighbors:
            if neighbor['site_index'] < site_index: continue
            CN2 = nn_object.get_cn(structure, neighbor['site_index'])
            #bl = math.dist(structure[site_index].coords, neighbor['site'].coords)
            bl = np.linalg.norm(structure[site_index].coords - neighbor['site'].coords)            #print(atom.specie, neighbor['site'].specie, "\t", site_index, neighbor['site_index'], "\t", CN1, CN2, "\t", bl)
            bonds.append({"site_1": site_index, "atom_1": atom.specie,
                      "Z_1": atom.specie.Z, "CN_1": CN1, "EN_1": d_EN[atom.specie.Z][CN1],
                      "site_2": neighbor['site_index'], "atom_2": neighbor['site'].specie,
                      "Z_2": neighbor['site'].specie.Z, "CN_2": CN2, "EN_2": d_EN[neighbor['site'].specie.Z][CN2],
                      "bond_length": bl})
    N = len(bonds)

    # Calculate Intrinsic hardness
    prod_EN = 1
    prod_BS = 1
    prod_C = 1

    for bond in bonds:
        # EN model
        fi_EN = 0.25*abs(bond["EN_1"] - bond["EN_2"]) / (bond["EN_1"]*bond["EN_2"])**0.5
        X_ij = ((bond["EN_1"]/bond["CN_1"]) * (bond["EN_2"]/bond["CN_2"]))**0.5
        prod_EN = prod_EN * (X_ij * math.exp(-2.7*fi_EN))

        # BS model
        fi_BS = ((bond["EN_1"] - bond["EN_2"]) / (bond["EN_1"] + bond["EN_2"]))**2
        S_ij = (1/bond["CN_1"]/bond["CN_2"]) * (bond["EN_1"]/0.481*bond["EN_2"]/0.481)**0.5 / bond["bond_length"]
        prod_BS = prod_BS * S_ij * math.exp(-2.8*fi_BS)

        # Cheenady model
        fi_C = ((bond["EN_1"] - bond["EN_2"]) / (bond["EN_1"] + bond["EN_2"]))**2
        Z_ij = (bond["EN_1"]/bond["CN_1"]) * (bond["EN_2"]/bond["CN_2"])
        prod_C = prod_C * (Z_ij**0.35) * (bond["bond_length"]**-1.48) * math.exp(-2.55*fi_C)

    H_EN = (469.27*N/vol) * prod_EN**(1/N) - 7.24
    H_BS = (1450*N/vol) * prod_BS**(1/N)
    H_C = (865*(N/vol)**1.13) * prod_C**(1/N)

    #print("Vickers hardness, as per EN model = {:.2f}".format(H_EN))
    #print("Vickers hardness, as per BS model = {:.2f}".format(H_BS))
    #print("Average Vickers hardness = {:.2f}".format(H_C))

    # Write output to file
    with open(os.path.join(run_data_dir,"hard_data"), 'a') as f_hard:
        f_hard.writelines(line.strip())
        f_hard.writelines("\t\t {:.2f}\t\t {:.2f}\t\t {:.2f}\n".format(H_EN, H_BS, H_C))

# Print the ouput file
with open(os.path.join(run_data_dir,"hard_data"), 'r') as f:
    print(f.read())
