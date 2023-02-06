#!/bin/bash

#SBATCH --job-name jellyfish
#SBATCH --output out_jellyfish.txt
#SBATCH --error err_jellyfish.txt
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
# jellyfish path
PATH="/nfsd/bcb/bcbg/prog_9/bin":$PATH

cd "/nfsd/bcb/bcbg/prog_9/USTAR"

for k in "${kmer_sizes[@]}"; do
  for s in ../SRR/*/*.fasta; do
    s_name=$(basename "$s" ".fasta")
    s_base=${s%".fasta"}

    kmer_list="${s_base}-a${a}-k${k}-kmers.txt"
    echo
    echo "*** Extracting kmers of $s_name..."

    if [[ -f "${kmer_list}" ]]; then
      echo "Already processed!"
      continue
    fi

    srun jellyfish-linux count -t 4 -m "$k" -C -s 100M -L "$a" "$s"
    srun jellyfish-linux dump -c mer_counts.jf > "${s_base}-kmers-unsorted.txt"

    srun sort "${s_base}-kmers-unsorted.txt" -o "$kmer_list"

    rm -f "${s_base}-kmers-unsorted.txt" "mer_counts.jf"

    echo "*** kmers of $s_name extracted!"
  done
done
