#!/bin/bash

#SBATCH --job-name fasterq
#SBATCH --output out_fasterq.txt
#SBATCH --error err_fasterq.txt
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH
#SBATCH --time 60:00
#SBATCH --cpus-per-task 8
#SBATCH --partition allgroups
#SBATCH --mem 8G

source ../configure

set -e

cd $SRR_DIR

for sra in */*".sra"; do
  s=$(basename "${sra%".sra"}")
  echo
  echo "*** Converting $sra..."
  [[ -f "${s}.fasta" || -f "${s}_1.fasta" || -f "${s}_2.fasta" ]] \
  || srun fasterq-dump --threads $NTHREAD --progress --fasta --split-files --outdir "$s" "$sra"
done

wait

echo "*** All converted!"
