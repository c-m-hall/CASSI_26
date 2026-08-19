#!/bin/bash
#SBATCH --job-name=stack_cubes
#SBATCH --partition=obs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=10:00:00
#SBATCH --output=/carnegie/scidata/groups/muse_stacking/CASSI_26/slurm_script/stack_cubes_%A.out
#SBATCH --error=/carnegie/scidata/groups/muse_stacking/CASSI_26/slurm_script/stack_cubes_%A.err

module load conda
conda activate /carnegie/scidata/groups/muse_stacking/conda_envs/cassi

cd /carnegie/scidata/groups/muse_stacking/CASSI_26

python stack_cubes.py --dir cassi_outputs/masked \
    -o /carnegie/scidata/groups/muse_stacking/CASSI_26/cassi_output/results/150fields_final_stack.fits \
    --method median
