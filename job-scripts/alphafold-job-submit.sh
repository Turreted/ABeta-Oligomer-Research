#!/bin/bash
#SBATCH --job-name=alphafold2
#SBATCH --partition=gpu
#SBATCH --nodes=4
#SBATCH --time=03:00:00
#SBATCH --account=pi-haddadian
#SBATCH --partition=caslake
#SBATCH --mem=64G

module load alphafold/2.2.0 cuda/11.3

echo "GPUs available: $CUDA_VISIBLE_DEVICES"
echo "CPU cores: $SLURM_CPUS_PER_TASK"

echo "Using `which run_alphafold`"

OUTPUT_DIR=/home/gideonm/scratch/gideonm/doris-alphafold/output
DOWNLOAD_DIR=/software/alphafold-data-2.2/

run_alphafold \
	-f 24MJ.fasta \
	-t 2020-05-14 \
	-m monomer \
	-p full_dbs \
	-g true \
	-o $OUTPUT_DIR
