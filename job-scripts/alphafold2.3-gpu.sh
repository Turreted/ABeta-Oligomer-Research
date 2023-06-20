#!/bin/bash
#SBATCH --job-name=alphafold-doris
#SBATCH --account=pi-haddadian
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:2
#SBATCH --constraint=v100
#SBATCH --mem=64G

# This file is intened to be used for a custom build of AlphaFold v2.3.0.
# The path to the executable, database, and run-alphafold directory must 
# be specified within the script.

# load so we can load conda locally
module load alphafold/2.2.0 cuda/11.3
conda activate alphafold-2.3

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

ALPHAFOLD_DIR=~/alphafold-2.3.1
ALPHAFOLD_EXE=$ALPHAFOLD_DIR/run_alphafold.sh
data_dir=/scratch/midway3/gideonm/alphafold-database/

echo "Alphafold executable: $ALPHAFOLD_EXE"
echo "Data directory: $data_dir"
echo "Output directory: $output_dir"
echo "Fasta files: $fasta_path"
echo "Running with GPU"
echo "Started job at $date"
echo ""

# This is the only way to get the script to work
cd $ALPHAFOLD_DIR

$ALPHAFOLD_EXE \
        -f $fasta_path \
        -t 2020-05-14 \
        -m monomer \
        -g true \
        -c reduced_dbs \
        -d $data_dir \
        -o $output_dir
