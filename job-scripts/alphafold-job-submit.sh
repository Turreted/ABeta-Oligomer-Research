#!/bin/bash
#SBATCH --job-name=alphafold2
#SBATCH --partition=gpu
#SBATCH --nodes=4
#SBATCH --time=03:00:00
#SBATCH --account=pi-haddadian
#SBATCH --partition=caslake
#SBATCH --mem=64G

module load alphafold/2.2.0 cuda/11.3

# This file is intened to be used for the default run_alphafold script installed
# on Midway3. It runs v2.2.0 by default.

# Input parameters so we can easily configure input/output, all other parameters
# should be modified from within the script
# shamelessley ripped from https://github.com/kalininalab/alphafold_non_docker/blob/main/run_alphafold.sh
usage() {
        echo ""
        echo "Please make sure all required parameters are given"
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        echo "-o <output_dir>       Path to a directory that will store the results."
        echo "-f <fasta_paths>      Path to FASTA files containing sequences. If a FASTA file contains multiple sequences, then it will be folded as a multimer. To fold more sequences one after another, write the files separated by a comma"
        echo ""
        exit 1
}

while getopts ":d:o:f:t:g:r:e:n:a:m:c:p:l:b:" i; do
        case "${i}" in
        o)
                output_dir=$(realpath $OPTARG)
        ;;
        f)
                fasta_path=$(realpath $OPTARG)
        ;;
        t)
        ;;
        esac
done

if [[ "$output_dir" == "" || "$fasta_path" == ""  ]] ; then
    usage
fi

download_dir=/software/alphafold-data-2.2/

echo "Alphafold executable: $(which run_alphafold)"
echo "Data directory: $data_dir"
echo "Output directory: $output_dir"
echo "Fasta files: $fasta_path"
echo "Running with GPU"
echo "Started job at $date"
echo ""

run_alphafold \
	-f $fasta_path \
	-t 2020-05-14 \
	-m monomer \
	-p full_dbs \
	-g true \
	-d $download_dir \
	-o $output_dir
