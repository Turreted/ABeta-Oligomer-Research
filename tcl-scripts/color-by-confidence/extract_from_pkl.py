import pickle
import sys
import numpy as np

def main():
    try:
        filename = sys.argv[1]
        fp = open(filename, "rb")
    except:
        print("Error: script must be run with a valid .pkl and .pdb input file")
        print("Example: ./extract_from_pkl.py result_model_1_pred.pkl")
        return 1

    outfile = f"{ filename.strip('.pkl') }_vmd.csv"
    p = pickle.load(fp)
    confidence_array = p['plddt']
    np.savetxt(outfile, confidence_array, delimiter="\n", fmt='%f')

if __name__ == "__main__":
    main()

