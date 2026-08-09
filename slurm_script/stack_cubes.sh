#!/bin/bash
#SBATCH --job-name=stack_cubes
#SBATCH --partition=obs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --output=slurm_script/process_cubes_%A.out
#SBATCH --error=slurm_script/process_cubes_%A.err

module load conda
conda activate /carnegie/scidata/groups/muse_stacking/conda_envs/cassi

python stack_cubes.py --/carnegie/scidata/groups/muse_stacking/CASSI_26/cassi_outputs/masked -o stacked.fits
