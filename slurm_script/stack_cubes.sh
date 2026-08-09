#!/bin/bash
#SBATCH --job-name=stack_cubes
#SBATCH --partition=obs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --output=slurm_script/stack_cubes_%A.out
#SBATCH --error=slurm_script/stack_cubes_%A.err

module load conda
conda activate /carnegie/scidata/groups/muse_stacking/conda_envs/cassi

cd /carnegie/scidata/groups/muse_stacking/CASSI_26

python stack_cubes.py --dir cassi_outputs/masked -o cassi_outputs/stacked/stacked.fits
