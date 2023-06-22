from ssparser.plotter import plot_monomer_hits
from ssparser.ssparser import parse_dataset_freq
from ssparser.scorer import generate_score_table
from ssparser.utils import get_run_pdb
from ssparser.stride import parse_with_stride
from ssparser.constants import *


alphafold_dir = "inputs/2M4J-monomer-exclude-all"
output_dir = "outputs/2M4J-monomer-exclude-all"
multimer = False

print("Generating scoring table...")
generate_score_table(alphafold_dir, output_dir)

print("Generating frequency table from fibril structures...")
structure_freq = parse_dataset_freq(dataset_dir=DATASET_DIR)

print("Generating monomer secondary strcuture charts...")
for pdb in get_run_pdb(alphafold_dir):
    fullpath = os.path.join(alphafold_dir, pdb)
    df = parse_with_stride(fullpath)

    # output for each chain, pad with coils
    if multimer:
        chain_iter = 0
        for chain_code in set(df["Chain"]):
            chain = df.loc[df["Chain"] == chain_code]
            ss = list(chain["SS Code"]) + ["C", "C"]

            outfile = (
                f"{os.path.join(output_dir, pdb.strip('.pdb'))}-chain-{chain_iter}.png"
            )
            plot_monomer_hits(ss, structure_freq, output_file=outfile)
            chain_iter += 1

            print(f"Output {outfile}")
    else:
        # extract monomer SS from datafile, pad with coils
        ss = list(df["SS Code"]) + ["C", "C"]

        outfile = os.path.join(output_dir, pdb.replace("pdb", "png"))
        plot_monomer_hits(ss, structure_freq, output_file=outfile)
        print(f"Output {outfile}")
