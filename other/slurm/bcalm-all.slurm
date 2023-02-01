#!/bin/bash
#SBATCH --job-name bcalm
#SBATCH --output bcalm.txt
#SBATCH --error bcalm.txt
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH
#SBATCH --time 10:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --partition allgroups
#SBATCH --mem 20G

set -e
NTHREAD=8
abundance_min=1
K=(15 21 41)
PATH="/nfsd/bcb/bcbg/prog_9/bin/bcalm":$PATH

cd SRR

for k in "${K[@]}"; do
  for s in */*.fasta; do
    outfile_base="${s:0:-6}.a${abundance_min}.k${k}"
    outfile="${outfile_base}.unitigs.fa"
    echo "*** Running bcalm with NTHREAD=$NTHREAD, abundance_min=$abundance_min, s=$s..."
    [ ! -f "$outfile" ] || continue
    srun bcalm -max-memory 16000 -nb-cores $NTHREAD -kmer-size "$k" -abundance-min $abundance_min \
    -in "$s" -all-abundance-counts -out "$outfile"
  done
done

cd ..

echo "Computed unitigs for all sequences!"
