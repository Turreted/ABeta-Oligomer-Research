import pandas as pd
import numpy as np
import os
import re
import pickle

from stride import parse_with_stride
from ssparser import parse_dataset_freq
from constants import *

"""
Scoring algorithm for ABeta-chains generated by AlphaFold

1. Generates a strcuture template based on ABeta monomers in fibrils or complexes
that are taken from the PDB. The files in the database are specified by DATABASE_DIR.
This template contains the frequency of each secondary structure at each position 
(ie the number of times the ith residue was observed to have that secondary strcuture
divided by the total number of chain)

2. Selects 'target residues' which are residues which frequentley have secondary 
strcutures of a specific type

3. For each result generated by AlphaFold, iterate over each target residue and
see if it matches the template. If it does, add the frequency of the secondary
structure in the template to the overall score.

4. Scale the final score by dividing it by the sum of the frequencies at each
position for each chain, thus scaling each residue by the frequency in the 
template.
"""


# returns the .pkl and .pdb output of an AlphaFold attempt
def get_full_runs(data_dir: str) -> tuple:
    # r = re.compile(r"relaxed_model_.*_multimer_.*_pred_.*.pdb")
    r = re.compile(r"relaxed_model_.*_pred_.*.pdb")
    pairs = []

    for model in list(filter(r.match, os.listdir(data_dir))):
        run_output = model.strip(".pdb").strip("relaxed")
        pairs.append(
            (
                model,
                [
                    f
                    for f in os.listdir(data_dir)
                    if run_output in f and f.endswith(".pkl")
                ][0],
            )
        )
    return pairs


# data_dir = "/Users/gideon/Documents/UChicago/Research/box/monomer/2M4J-exclude-templates"
data_dir = (
    "/Users/gideon/Documents/UChicago/Research/box/multimer/2M4J-dimer-unmodified"
)

get_full_runs(data_dir)


beta_threshold = 0.3
turn_threshold = 0.2
template_dict = parse_dataset_freq()
beta_template = template_dict["Beta"]
turn_template = template_dict["Turn"]

# select all residue indices where the frequency of secondary strcuture is above a threshold
target_beta_residues = [
    ri for ri, freq in enumerate(beta_template) 
    if freq > beta_threshold
]
target_turn_residues = [
    ri for ri, freq in enumerate(turn_template) 
    if freq > turn_threshold
]

print(
    f"Beta-Strand target residues (threshold of {beta_threshold*100}%): {target_beta_residues}"
)
print(
    f"Turn target residues (threshold of {turn_threshold*100}%):        {target_turn_residues}"
)

res = pd.DataFrame(columns=["Filename", "Beta Score", "Turn Score", "AlphaFold Score"])

for sname, mname in get_full_runs(data_dir):
    print(f"Processing {sname} {mname}")
    full_struct_path = os.path.join(data_dir, sname)
    full_metadata_path = os.path.join(data_dir, mname)

    struct = parse_with_stride(full_struct_path)

    chain_names = set(struct["Chain"])
    chain_count = len(chain_names)
    beta_score = 0
    turn_score = 0

    for cname in chain_names:
        chain = struct.loc[struct["Chain"] == cname]

        for res_index in range(ABETA_LONG_RESIDUE_COUNT):
            resid = res_index + 1
            residue_data = chain.loc[chain["PDB ResID"] == resid]

            if not residue_data.empty:
                res_struct_code = list(residue_data["SS Code"])[0]

                if (
                    res_struct_code in ["B", "b", "E"]
                    and res_index in target_beta_residues
                ):
                    beta_score += beta_template[res_index]

                elif (
                    res_struct_code in ["T"] 
                    and res_index in target_turn_residues
                ):
                    turn_score += turn_template[res_index]

    # extract AlphaFold score
    mf = open(full_metadata_path, "rb")
    pf = pickle.load(mf)
    plddt_scores = pf["plddt"]
    target_residues = set(target_beta_residues).union(target_turn_residues)
    af_score = np.mean([plddt_scores[i] for i in target_residues])

    # take sum of all frequencies so we can divide by then and get scaled result
    beta_sum = sum([beta_template[ri] for ri in target_beta_residues])
    turn_sum = sum([turn_template[ri] for ri in target_turn_residues])

    res.loc[len(res.index)] = [
        os.path.basename(sname),
        round(beta_score / (beta_sum * chain_count), 3),
        round(turn_score / (turn_sum * chain_count), 3),
        round(af_score, 3),
    ]

res = res.sort_values(by=['AlphaFold Score'], ascending=False)
print(res.head())
res.to_csv("testout.csv")
