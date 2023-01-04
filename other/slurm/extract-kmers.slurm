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

K=31

# jellyfish path
PATH="/nfsd/bcb/bcbg/prog_9/bin":$PATH

cd "/nfsd/bcb/bcbg/prog_9/USTAR"

for s in ../SRR/*/*.fasta; do
  s_name=$(basename "$s" ".fasta")
  s_base=${s%".fasta"}

  echo
  echo "*** Extracting kmers of $s_name..."

  [[ -f "${s_base}-kmers.txt" ]] \
  || srun jellyfish-linux count -t 4 -m "$K" -C -s 100M -L 2 "$s" \
  && srun jellyfish-linux dump -c mer_counts.jf > "${s_base}-kmers-unsorted.txt"

  srun sort "${s_base}-kmers-unsorted.txt" -o "${s_base}-kmers.txt"

  rm "${s_base}-kmers-unsorted.txt" "mer_counts.jf"

  echo "*** kmers of $s_name extracted!"
done