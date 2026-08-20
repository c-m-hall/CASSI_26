#!/bin/bash
#SBATCH --job-name=center_spaxel_check
#SBATCH --output=slurm_script/check_spaxel_%A.out
#SBATCH --error=slurm_script/check_spaxel_%A.err
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# ------------------------------------------------------------------------
# Runs check_center_spaxel_mask.py against a directory of *_masked.fits
# science cubes. In these files, apply_mask.py has already converted
# masked (mask==0) voxels to NaN -- so we use:
#   --mode cube        : check DATA ext0 for NaN (this IS the default,
#                         included explicitly here for clarity)
#   --wave-check mid    : only check the middle wavelength channel at the
#                         center spaxel (fast; equivalent to the old
#                         "--mode plane" behavior). Use --wave-check all
#                         to require every wavelength channel be NaN.
# ------------------------------------------------------------------------

set -euo pipefail

module load conda
conda activate /carnegie/scidata/groups/muse_stacking/conda_envs/cassi
cd /carnegie/scidata/groups/muse_stacking/CASSI_26

python check_center_spaxel_mask.py \
    /carnegie/scidata/groups/muse_stacking/CASSI_26/cassi_outputs/masked \
    --mode cube \
    --wave-check mid \
    --glob "*_masked.fits" \
    -o results.csv
