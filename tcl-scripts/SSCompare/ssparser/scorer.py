import pandas as pd
import numpy as np
import os
import re
import pickle

from stride import parse_with_stride
from ssparser import parse_dataset_freq
from constants import *


# get file that are rankedn_.pdn
def get_structure_files(data_dir: str):
    return sorted(
        [
            os.path.join(data_dir, fname)
            for fname in os.listdir(data_dir)
            if "ranked" in fname
        ]
    )


# returns the .pkl and .pdb output of an AlphaFold attempt
def get_full_runs(data_dir: str) -> tuple:
    #r = re.compile(r"relaxed_model_.*_multimer_.*_pred_.*.pdb")
    r = re.compile(r"relaxed_model_.*_pred_.*.pdb")
    pairs = []

    for model in list(filter(r.match, os.listdir(data_dir))):
        run_output = model.strip(".pdb").strip("relaxed")
        pairs.append(
            (
                model,
                [f for f in os.listdir(data_dir)
                 if run_output in f and f.endswith(".pkl")][0],
            )
        )
    return pairs


data_dir = (
    "/Users/gideon/Documents/UChicago/Research/box/monomer/2M4J-exclude-templates"
)
get_full_runs(data_dir)


beta_threshold = 0.3
turn_threshold = 0.2
freq_dict = parse_dataset_freq()
beta_freq = freq_dict["Beta"]
turn_freq = freq_dict["Turn"]

# select all residue indices where the frequency of secondary strcuture is above a threshold
target_beta_residues = [
    ri for ri, freq in enumerate(beta_freq) if freq > beta_threshold
]
target_turn_residues = [
    ri for ri, freq in enumerate(turn_freq) if freq > turn_threshold
]
# scale each residue by its frequency then multiply by the sum
beta_residue_scale = []
turn_residue_scale = []

print(target_beta_residues, target_turn_residues)
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
        # print(chain["SS Code"])

        for res_index in target_beta_residues:
            resid = res_index + 1
            residue_data = chain.loc[chain["PDB ResID"] == resid]

            if not residue_data.empty:
                res_struct_code = list(residue_data["SS Code"])[0]

                if res_struct_code in ["B", "b", "E"]:
                    beta_score += 1

        for res_index in target_turn_residues:
            resid = res_index + 1
            residue_data = chain.loc[chain["PDB ResID"] == resid]

            if not residue_data.empty:
                res_struct_code = list(residue_data["SS Code"])[0]

                if res_struct_code in ["T"]:
                    turn_score += 1

    # extract AlphaFold score
    mf = open(full_metadata_path, "rb")
    pf = pickle.load(mf)
    
    plddt_scores = pf["plddt"]
    target_residues = set(target_beta_residues).union(target_turn_residues)
    af_score = np.mean([plddt_scores[i] for i in target_residues])

    # Scale each value by the frequency in the tempate so that stronger template
    # residues lead to better results

    # take sum of all frequencies so we can divide by then and get scaled result
    beta_sum = sum([beta_freq[ri] for ri in target_beta_residues])
    turn_sum = sum([turn_freq[ri] for ri in target_turn_residues])

    res.loc[len(res.index)] = [
        os.path.basename(sname),
        round(beta_score / (chain_count * len(target_beta_residues)), 3),
        round(turn_score / (chain_count * len(target_turn_residues)), 3),
        round(af_score, 3),
    ]

print(res.head())
res.to_csv("testout.csv")
