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

cd /nfsd/bcb/bcbg/prog_9/SRR


set -e
NTHREAD=8
PATH="/nfsd/bcb/bcbg/prog_9/bin/sratoolkit.3.0.1-ubuntu64/bin":$PATH

for sra in */*.sra; do
  s=${sra%".sra"}
  echo
  echo "*** Converting $s..."
  [[ -f "${s}.fasta" || -f "${s}_1.fasta" || -f "${s}_2.fasta" ]] \
  || srun fasterq-dump --threads $NTHREAD --progress --fasta --split-files "$s" --outdir "$s"
done

wait

echo "*** All converted!"
