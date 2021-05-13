# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

"""
Intrinsic Harness from GASP results:

This script reads the structures from the POSCAR files output by GASP
and calculates the intrinsic hardness of the structures based on
the modified electonegativity (EN) and bond strength (BS) models.
Output file: hardness_data

Usage: python intrinsic_hardness.py /path/to/run_data/file
POSCAR files must be in the same directory as the run_data file
"""

# Import libraries and read files
from gasp import electronegativities
from pymatgen.core import Structure
from pymatgen.analysis.local_env import CrystalNN
import numpy as np
import math
import os
import sys
import warnings

warnings.filterwarnings('ignore')

def main():
    # Read run_data
    run_data_path = sys.argv[1]
    with open(run_data_path, 'r') as f_gasp:
        lines = f_gasp.readlines()

    # Prepare output file
    run_data_dir = os.path.dirname(run_data_path)
    with open(os.path.join(run_data_dir,"hardness_data"), 'w') as f_hard:
        for line_number, line in enumerate(lines[:4]):
            if line_number == 2:
                line = line.strip()
                f_hard.writelines(line)
                f_hard.writelines("\t\t H_EN\t\t H_BS\t\t H_C")
            else:
                f_hard.writelines(line)
        f_hard.writelines("\n")


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
            try:
                neighbors = nn_object.get_nn_info(structure, site_index)
            except:
                continue
            # if not neighbors: continue
            CN1 = nn_object.get_cn(structure, site_index)
            for neighbor in neighbors:
                if neighbor['site_index'] < site_index: continue
                CN2 = nn_object.get_cn(structure, neighbor['site_index'])
                if CN1==0 or CN2==0: continue
                #bl = math.dist(structure[site_index].coords, neighbor['site'].coords)
                bl = np.linalg.norm(structure[site_index].coords - neighbor['site'].coords)            #print(atom.specie, neighbor['site'].specie, "\t", site_index, neighbor['site_index'], "\t", CN1, CN2, "\t", bl)
                bonds.append({"site_1": site_index, "atom_1": atom.specie,
                          "Z_1": atom.specie.Z, "CN_1": CN1,
                          "EN_1": electronegativities.get_EN(atom.specie.Z, CN1),
                          "site_2": neighbor['site_index'], "atom_2": neighbor['site'].specie,
                          "Z_2": neighbor['site'].specie.Z, "CN_2": CN2,
                          "EN_2": electronegativities.get_EN(neighbor['site'].specie.Z, CN2),
                          "bond_length": bl})
        N = len(bonds)

        #Temporary fix for when no bonds are found
        if N == 0:
            with open(os.path.join(run_data_dir,"hardness_data"), 'a+') as f_hard:
                f_hard.writelines(line.strip())
                f_hard.writelines("\t\t nan\t\t nan\t\t nan\n".format(H_EN, H_BS, H_C))
        else:
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


        # Write output to file
        with open(os.path.join(run_data_dir,"hardness_data"), 'a') as f_hard:
            f_hard.writelines(line.strip())
            f_hard.writelines("\t\t {:.2f}\t\t {:.2f}\t\t {:.2f}\n".format(H_EN, H_BS, H_C))

    # Print the ouput file
    with open(os.path.join(run_data_dir,"hardness_data"), 'r') as f:
        print(f.read())


if __name__ == "__main__":
    main()
