import pickle
import sys
import json
import numpy as np

"""
Converts an AF-2.X run to an AF-2.3.X run
"""

def main():
    try:
        filename = sys.argv[1]
        fp = open(filename, "rb")
    except:
        print("Error: script must be run with a valid .pkl and .pdb input file")
        print("Example: ./AF-convert.py result_model_1_pred.pkl")
        return 1

    outfile = f"{ filename.strip('.pkl') }_vmd.csv"
    p = pickle.load(fp)
    print(p.keys())

    plddt = {"confidenceScore": p["plddt"].tolist()}
    pae = [
        {"predicted_aligned_error": p["predicted_aligned_error"].tolist(), 
        "max_predicted_aligned_error": p["max_predicted_aligned_error"].tolist()}]

    model = filename.strip("result_").strip(".pkl")
    plddt_name = f"confidence_{model}.json"
    pae_name = f"pae_{model}.json"

    with open(plddt_name, "w") as outfile:
        outfile.write(json.dumps(plddt, indent=4))

    with open(pae_name, "w") as outfile:
        outfile.write(json.dumps(pae, indent=4))

if __name__ == "__main__":
    main()
