import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

import os

abeta_sequence = ['ASP', 'ALA', 'GLU', 'PHE', 'ARG', 'HIS', 'ASP', 'SER', 'GLY', 'TYR', 'GLU', 'VAL', 'HIS', 'HIS', 'GLN', 'LYS', 'LEU', 'VAL', 'PHE', 'PHE', 'ALA', 'GLU', 'ASP', 'VAL', 'GLY', 'SER', 'ASN', 'LYS', 'GLY', 'ALA', 'ILE', 'ILE', 'GLY', 'LEU', 'MET', 'VAL', 'GLY', 'GLY', 'VAL', 'VAL', 'ILE', 'ALA']
abeta_str = "".join(abeta_sequence)
csv_dir = "abeta-fibril-csv"
residue_count = 42
total_chain_count = 0
long_chain_chount = 0

alpha = np.zeros(residue_count, dtype=int)
beta  = np.zeros(residue_count, dtype=int)
turn  = np.zeros(residue_count, dtype=int)
other = np.zeros(residue_count, dtype=int)

#alpha, beta, turn, others
for fname in os.listdir(csv_dir):
    infile = os.path.join(csv_dir, fname)
    df = pd.read_csv(infile)        

    # iterate though each chain in the complex
    chain_names = set(df["Chain"])
    for cname in chain_names:
        chain = df.loc[df['Chain'] == cname]
        
        # check if the current chain is an ABeta monomer, sequence its structure if it is
        aa_sequence = list(chain["Residue Name"])
        if "".join(aa_sequence) in abeta_str:
            print(f"Sequencing {infile} Chain {cname}")
            
            total_chain_count += 1
            if len(aa_sequence) == 42:
                long_chain_chount += 1

            chain_start = int(min(chain["PDB ResID"]))
            chain_stop  = int(max(chain["PDB ResID"]))

            # iterate through each residue in the chain. If it is in the structure,
            # add it to out sequence
            for res_index in range(residue_count):
                # select data by residueID. Note that this is 1-indexed
                resid = res_index + 1
                residue_data = chain.loc[chain["PDB ResID"] == resid]

                if not residue_data.empty:
                    structure = list(residue_data["SS Code"])[0]
                    
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
                    else:
                        other[res_index] += 1
                
                else:
                    other[res_index] += 1


legend_dict = {"Alpha": alpha, "Beta": beta, "Turn": turn}

x = np.arange(residue_count)  
width = 0.25  # the width of the bars
multiplier = 3

fig, ax = plt.subplots(layout='constrained')
plt.xticks(np.arange(1,43, dtype=int), abeta_sequence, rotation=-60)

# add line at 50%
plt.axhline(y=0.2, color='r', linestyle='-')

for index, (attribute, measurement) in enumerate(legend_dict.items()):
    # normalize over total numeber of chains
    if index < 40:
        measurement = measurement / total_chain_count
    else:
        measurement = measurement / long_chain_chount

    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    multiplier += 1

ax.set_ylabel("Frequency of Secondary Structure")
ax.set_xlabel("Residue ID")
ax.set_title("Frequency of Secondary Structure per Residue of ABeta-Monomer")
ax.legend(loc='upper right')

plt.show()

