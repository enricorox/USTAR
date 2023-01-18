#!/bin/bash
#SBATCH --job-name compile-bcalm
#SBATCH --output out_compile-bcalm.txt
#SBATCH --error err_compile-bcalm.txt
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH
#SBATCH --time 10:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --partition allgroups
#SBATCH --mem 20G

set -e

cd "/nfsd/bcb/bcbg/prog_9/bin/bcalm"
srun cmake .
srun make -j 16
