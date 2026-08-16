#!/bin/bash
#SBATCH --job-name=nodownload_process
#SBATCH --output=slurm_script/nodownload_process_%j.out
#SBATCH --error=slurm_script/nodownload_process_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4



module load conda
conda activate /carnegie/scidata/groups/muse_stacking/conda_envs/cassi
cd /carnegie/scidata/groups/muse_stacking/CASSI_26
 

 
# 1. build the local pairs table from whatever cubes are ALREADY on disk
# in ./scratch/ (e.g. if a previous download job was stopped/killed
# partway through and never reached its own final pairs_table.write()).
python build_pairs_from_existing.py \
    --n-fields 150 \
    --sample muse_x_milliquas_sample.fits \
    --dest ./scratch/ \
    --pairs-out muse_x_milliquas_pairs_local.fits
 
# 2. run the rest of the pipeline (convert_wl -> psf_sub -> mask -> crop ->
# resample -> apply mask) on exactly the cubes that pairs table covers
python3 <<'PYEOF'
import mainscript
 
mainscript.main(
    sample_path='muse_x_milliquas_sample.fits',
    pairs_path='muse_x_milliquas_pairs_local.fits',
    cube_dir='./scratch/',
    sex_config='/carnegie/scidata/groups/muse_stacking/CASSI_26/wl_eso.sex',
    sex_binary='sex',
    resume_from='download',
)
PYEOF
