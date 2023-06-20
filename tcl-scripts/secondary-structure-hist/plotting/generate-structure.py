import subprocess
import tempfile
import re
import os

"""
Generates the structure of a given PDB file using the STRIDE algorithm, which must
be installed to the user's $PATH. Instructions for installing the STRIDE executable
can be found at https://webclu.bio.wzw.tum.de/stride/install.html.
"""

# overwrite with path to your own binary, if necessary
STRIDE_EXE = "stride"


def parse_with_stride(infile: str, outfile: str):
    """
    Output Format:
        6-8  Residue name
        10-10 Protein chain identifier
        12-15 PDB	residue	number
        17-20 Ordinal residue number
        25-25 One	letter secondary structure code	**)
        27-39 Full secondary structure name
        43-49 Phi	angle
        53-59 Psi	angle
        65-69 Residue solvent accessible area
    """

    print(f"Running for {infile}...")
    tfile = tempfile.NamedTemporaryFile(mode="w", delete=False)
    tfilename = tfile.name

    # pipe output to tmpfile, wait for it to complete
    subprocess.check_call(
        [STRIDE_EXE, infile], stdout=tfile, stderr=subprocess.PIPE
    )

    csv_contents = [
        "Residue Name,Chain,PDB ResID,Ordinal ResID,SS Code,SS Full,Phi Angle,Psi Angle,Solvent Area,PDB Code"
    ]
    with open(tfilename, "r") as f:
        for line in f.readlines():
            if line.startswith("ASG"):
                csv_contents.append(re.sub(r"\s+", ",", line[5:-1]))

    with open(outfile, "w") as f:
        for line in csv_contents:
            f.write(f"{line}\n")


in_dir = "abeta-ent"
out_dir = "abeta-csv"

for fname in os.listdir(in_dir):
    infile = os.path.join(in_dir, fname)
    outfile = os.path.join(out_dir, fname.replace(".ent", ".csv"))
    parse_with_stride(infile, outfile)
