import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

import os

csv_dir = "abeta-csv"
residue_count = 42
total_chain_count = 0

alpha = np.zeros(residue_count, dtype=int)
beta  = np.zeros(residue_count, dtype=int)
turn  = np.zeros(residue_count, dtype=int)
other = np.zeros(residue_count, dtype=int)

#alpha, beta, turn, others
for fname in os.listdir(csv_dir):
    infile = os.path.join(csv_dir, fname)
    df = pd.read_csv(infile)

    cur_chain = None
    prev_res_index = 0

    # write secondary strcuture at position to cumulative array
    for index, row in df.iterrows():
        resid, chain, structure = row["PDB ResID"], row["Chain"], row["SS Code"]
        res_index = resid - 1 # resID is 1-indexed

        # check if this is a new chain
        if cur_chain != chain:
            # consider all unmarked values as "other"
            for i in range(0, res_index):
                other[i] += 1
        
            total_chain_count += 1
            cur_chain = chain

        # alpha-helix
        if structure in ["H", "G", "I"]:
            alpha[res_index] += 1
        
        # beta strand or sheet
        elif structure in ["B", "b", "E"]:
            beta[res_index] += 1
        
        # turn
        elif structure in ["T"]:
            turn[res_index] += 1

        # other, listed as a coil
        elif structure in ["C"]:
            other[res_index] += 1

        prev_res_index = res_index


legend_dict = {"Alpha": alpha, "Beta": beta, "Turn": turn, "Other": other}

x = np.arange(residue_count)  
width = 0.25  # the width of the bars
multiplier = 3

fig, ax = plt.subplots(layout='constrained')

for attribute, measurement in legend_dict.items():
    # normalize over total numeber of chains
    measurement = measurement / total_chain_count

    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    multiplier += 1

ax.set_ylabel("Frequency of Secondary Structure")
ax.set_xlabel("Residue ID")
ax.set_title("Frequency of Secondary Structure per Residue of ABeta-Monomer")
ax.legend(loc='upper right', ncols=3)

plt.show()
