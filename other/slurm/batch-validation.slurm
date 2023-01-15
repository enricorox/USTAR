#!/bin/bash

#SBATCH --job-name validation
#SBATCH --output out_validation.txt
#SBATCH --error err_validation.txt
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH
#SBATCH --time 24:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --partition allgroups
#SBATCH --mem 16G

set -e

K=31

ROOT_DIR="/nfsd/bcb/bcbg/prog_9"
SRR="${ROOT_DIR}/SRR"

# ustarx
PATH="${ROOT_DIR}/USTAR:$PATH"

cd "${ROOT_DIR}/results"

for s_name in SRR*; do
  kmers_list=("${SRR}/"*"/${s_name}-kmers.txt")
  cd "$s_name"
  params_dir=("m"*)
  params_dir+=("ust")

  for param_dir in "${params_dir[@]}"; do
    echo "Doing validation for $s_name, $param_dir"
    cd "$param_dir"
    fasta=(*.fa)
    counts=(*.counts)

    if [[ "$param_dir" == *"rle"* ]]; then
      encoding="rle"
    else
      encoding="plain"
    fi

    # shellcheck disable=SC2128
    [[ -f "ustar-kmers.txt" ]] || srun ustarx -k "$K" -i "$fasta" -c "$counts" -e "$encoding" -s

    # shellcheck disable=SC2128
    if cmp "ustar-kmers.txt" "$kmers_list"; then
      echo "SUCCESS: (${s_name}, ${param_dir}) is correct!"
      echo "Correct!" > validity.txt
    else
      echo "FAILURE: (${s_name}, ${param_dir}) is NOT correct!"
      echo "Incorrect!" > validity.txt
    fi

    cd ..
  done

  cd ..
done

echo "All validation done!"
echo "Check incorrect files with 'grep $ROOT_DIR/results/*/*/validity.txt Incorrect'"
