import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

import os

abeta_sequence = ['ASP', 'ALA', 'GLU', 'PHE', 'ARG', 'HIS', 'ASP', 'SER', 'GLY', 'TYR', 'GLU', 'VAL', 'HIS', 'HIS', 'GLN', 'LYS', 'LEU', 'VAL', 'PHE', 'PHE', 'ALA', 'GLU', 'ASP', 'VAL', 'GLY', 'SER', 'ASN', 'LYS', 'GLY', 'ALA', 'ILE', 'ILE', 'GLY', 'LEU', 'MET', 'VAL', 'GLY', 'GLY', 'VAL', 'VAL', 'ILE', 'ALA']
abeta_sequence = "".join(abeta_sequence)
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

    # iterate though each chain in the complex
    chain_names = set(df["Chain"])
    for cname in chain_names:
        chain = df.loc[df['Chain'] == cname]
        
        # check if the current chain is an ABeta monomer, sequence its structure if it is
        aa_sequence = list(chain["Residue Name"])
        if "".join(aa_sequence) in abeta_sequence:
            print(f"Sequencing {infile} Chain {cname}")
            
            total_chain_count += 1

            # TODO: POSSIBLE BUG here
            chain_start = int(min(chain["PDB ResID"]))
            chain_stop  = int(max(chain["PDB ResID"]))

            # iterate through each residue in the chain. If it is in the structure,
            # add it to out sequence
            for resid in range(chain_start, chain_stop):
                # select residue by resid
                residue_data = chain.loc[chain["PDB ResID"] == resid]

                if not residue_data.empty:
                    structure = list(residue_data["SS Code"])[0]
                    
                    # alpha-helix
                    if structure in ["H", "G", "I"]:
                        alpha[resid-1] += 1
                    
                    # beta strand or sheet
                    elif structure in ["B", "b", "E"]:
                        beta[resid-1] += 1
                    
                    # turn
                    elif structure in ["T"]:
                        turn[resid-1] += 1

                    # other, listed as a coil
                    else:
                        other[resid-1] += 1
                
                else:
                    other[resid-1] += 1


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
ax.legend(loc='upper right', ncols=2)

plt.show()

