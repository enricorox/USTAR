#!/bin/bash

#SBATCH --job-name kmers-stats
#SBATCH --output out_kmers-stats.txt
#SBATCH --error err_kmers-stats.txt
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH
#SBATCH --time 60:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --partition allgroups
#SBATCH --mem 16G

set -e

kmer_sizes=(15 21 41)
a=1

# ustars path
PATH="/nfsd/bcb/bcbg/prog_9/USTAR":$PATH

cd "/nfsd/bcb/bcbg/prog_9"

for k in "${kmer_sizes[@]}"; do
  for s in SRR/*/*.fasta; do
    s_name=$(basename "$s" ".fasta")
    s_base=${s%".fasta"}
    kmers_list="${s_base}-a$a-k$k-kmers.txt"

    echo "Computing stats of $s_name ($kmers_list)"
    srun ustars -i "$kmers_list" -k "$k" > "${kmers_list%".txt"}-stats.txt"
  done
done
cat SRR/*/*-stats.txt > SRR/stats.txt
echo "Done. Kmers stats are in SRR/stats.txt"
