import re
import os

# returns the .pkl and .pdb output of an AlphaFold attempt
def get_full_runs(data_dir: str, multimer=False) -> tuple:
    # r = re.compile(r"relaxed_model_.*_multimer_.*_pred_.*.pdb")
    re_str = r"relaxed_model_.*_multimer_.*_pred_.*.pdb" if multimer else r"relaxed_model_.*_pred_.*.pdb"
    r = re.compile(re_str) 
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

# get the .pdb files of each output
def get_run_pdb(data_dir: str, multimer=False) -> list:
    re_str = r"relaxed_model_.*_multimer_.*_pred_.*.pdb" if multimer else r"relaxed_model_.*_pred_.*.pdb"
    r = re.compile(re_str) 
    return list(filter(r.match, os.listdir(data_dir)))
