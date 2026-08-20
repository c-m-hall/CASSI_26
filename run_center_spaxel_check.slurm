#!/bin/bash
#SBATCH --job-name=center_spaxel_check
#SBATCH --output=center_spaxel_check_%j.out
#SBATCH --error=center_spaxel_check_%j.err
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# ------------------------------------------------------------------------
# Runs check_center_spaxel_mask.py with --check-nan against a directory
# of *_masked.fits science cubes (NaN = masked center spaxel, per
# apply_mask.py's convention -- see script docstring for details).
#
# EDIT THESE BEFORE SUBMITTING:
#   SCRIPT_PATH  - path to check_center_spaxel_mask.py
#   DATA_DIR     - directory containing the *_masked.fits cubes
#   OUT_CSV      - where to write the results CSV
#   FILE_GLOB    - pattern matching your masked cube files
#   MODE         - "all" (default), "any", or "middle"
# ------------------------------------------------------------------------

set -euo pipefail

SCRIPT_PATH="/carnegie/scidata/groups/muse_stacking/CASSI_26/check_center_spaxel_mask.py"
DATA_DIR="/carnegie/scidata/groups/muse_stacking/CASSI_26/cassi_outputs/masked"
OUT_CSV="/carnegie/scidata/groups/muse_stacking/CASSI_26/cassi_outputs/results/center_spaxel_nan_results.csv"
FILE_GLOB="*_masked.fits"
MODE="middle"

# Activate your environment here if needed, e.g.:
# module load python/3.11
# source /path/to/venv/bin/activate

mkdir -p "$(dirname "$OUT_CSV")"

python3 "$SCRIPT_PATH" "$DATA_DIR" \
    --check-nan \
    --glob "$FILE_GLOB" \
    --mode "$MODE" \
    -o "$OUT_CSV"
