#!/bin/bash
#SBATCH --job-name=center_spaxel_check
#SBATCH --output=center_spaxel_check_%j.out
#SBATCH --error=center_spaxel_check_%j.err
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=slurm_script/check_spaxel_%A.out
#SBATCH --error=slurm_script/check_spaxel_%A.err





module load conda
conda activate /carnegie/scidata/groups/muse_stacking/conda_envs/cassi

cd /carnegie/scidata/groups/muse_stacking/CASSI_26

python check_center_spaxel_mask.py /carnegie/scidata/groups/muse_stacking/CASSI_26/cassi_outputs/masked --check-nan --glob "*_masked.fits" --mode plane -o results.csv
