#!/bin/bash
# Job name: jds3_test
#SBATCH --job-name=jds2
#
# Project (change to your :
#SBATCH --account=uio
#
# Wall clock limit:
#SBATCH --time=100:0:0
#
# Max memory usage per core (MB):
#SBATCH --mem-per-cpu=61G
#
## Number of tasks (nodes):
##SBATCH --ntasks-per-node=12
##

##Setup job environment
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors


#Load R module

module load R

#Set variable with number of processors and load mpi

# export OMP_NUM_<span id="cke_bm_104E" style="display: none" data-cke-bookmark="1"> </span>THREADS=1
module load R # openmpi.gnu
mpirun -np 1 R CMD BATCH --no-save --no-restore /usit/abel/u1/alexaei/scripts/projects/02_hundred_lakes/archaea/dada2_abel.R

#End of script
