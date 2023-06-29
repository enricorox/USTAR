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

source ../configure

abundance_min=1
KMER_SIZES=(13 17)

cd $SRR_DIR

for k in "${KMER_SIZES[@]}"; do
  for s in */*.fasta; do
    outfile_base="${s:0:-6}.a${abundance_min}.k${k}"
    outfile="${outfile_base}.unitigs.fa"
    # SRR/SRR001665/SRR001665_1.fasta
    seq=$(basename "${s:0:-8}")

    # skip edge case (no 41-mers)
    # shellcheck disable=SC2128
    if [[ $k == 41 && $seq == "SRR001665" ]]; then
      echo "Skipping SRR001665..."
      continue;
    fi

    echo "*** Running bcalm with NTHREAD=$NTHREAD, abundance_min=$abundance_min, s=$s..."
    [ ! -f "$outfile" ] || continue
    srun bcalm -max-memory 16000 -nb-cores $NTHREAD -kmer-size "$k" -abundance-min $abundance_min \
    -in "$s" -all-abundance-counts -out "$outfile_base"

    rm -f "${outfile}"*glue* "${outfile}"*h5
  done
done

echo "Computed unitigs for all sequences!"
