#!/bin/bash
#SBATCH --job-name=stack_cubes
#SBATCH --partition=obs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --output=/carnegie/scidata/groups/muse_stacking/CASSI_26/slurm_script/stack_cubes_%A.out
#SBATCH --error=/carnegie/scidata/groups/muse_stacking/CASSI_26/slurm_script/stack_cubes_%A.err

module load conda
conda activate /carnegie/scidata/groups/muse_stacking/conda_envs/cassi

cd /carnegie/scidata/groups/muse_stacking/CASSI_26

python stack_cubes.py ADP.2016-06-21T00:31:05.324_vac.fits ADP.2016-07-12T14:10:36.607_vac.fits ADP.2016-08-08T11:40:19.878_vac.fits ADP.2016-11-16T15:11:46.454_vac.fits -o four_stacked.fits

