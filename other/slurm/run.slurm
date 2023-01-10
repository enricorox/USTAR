#!/bin/bash
#SBATCH --job-name run-ustar
#SBATCH --output out_run-ustar.txt
#SBATCH --error err_run-ust.txt
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH
#SBATCH --time 48:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --partition allgroups
#SBATCH --mem 16G

# compute file size
dim(){
  stat -c %s "${1}"
}

# exit if errors occur
set -e
# expand globs
shopt -s nullglob

ROOT_DIR="/nfsd/bcb/bcbg/prog_9"

cd "${ROOT_DIR}/USTAR"

# ustar & MFCompressC
PATH=${ROOT_DIR}/USTAR:${ROOT_DIR}/bin:$PATH

params=("-sf -xf" "-sf -xf -e flip_rle"\
        "-s-aa -x=a" "-s-aa -x=a -e flip_rle"\
        ""
)

mkdir -p "results" && cd "results"

for bcalm_unitigs_abs in "$ROOT_DIR"/SRR/SRR*/SRR*".a2.unitigs.fa"; do
  # derive strings
  bcalm_unitigs=$(basename "$bcalm_unitigs_abs")
  seq_name=${bcalm_unitigs:0:-14} # remove ".ax.unitigs.fa" from "SRR001665_2.unitigs.fa"
  seq_dir=$(echo "$seq_name" | awk -F '_' '{printf $1};') # remove _1 or _2 if any

  results_csv="${ROOT_DIR}/${seq_name}.csv"
  rm -f ./*.csv

  echo "*** Processing $bcalm_unitigs ($seq_dir:$seq_name) >> $results_csv"

  mkdir -p "$seq_name" && cd "$seq_name"

  ## reset outputs vars
  # data
  headers="${seq_name},"
  n_kmers="n_kmers,"
  ust_c_lb="CL_LB,"
  ust_c="CL,"
  outfile_fa_sizes="fasta size,"
  outfile_fa_zipped_sizes="fasta zipped size,"
  outfile_counts_sizes="counts size,"
  outfile_counts_zipped_sizes="counts zipped size,"
  # ratios (to compute in worksheet)
  gap="gap (lower is better),"
  ratio_counts_ust_out="ratio counts,"

  ## file names
  outfile_stat="${bcalm_unitigs}.stats.txt"
  outfile_fa="${bcalm_unitigs}.ustar.fa"
  outfile_fa_zipped="${outfile_fa}.mfc"
  outfile_counts="${bcalm_unitigs}.ustar.counts"
  outfile_counts_zipped="${outfile_counts}.bz3"

  for param in "${params[@]}"; do
    echo "** params: $param"

    mkdir -p "m${param}" && cd "m${param}"

    # need to remove .ustar.fa in SRR001665_1.a2.unitigs.fa.ustar.fa.ustar.fa
    # shellcheck disable=SC2086
    [ -f "$outfile_fa" ] || srun ustar -k 31 -i "$bcalm_unitigs_abs" -o "${outfile_fa:0:-9}" \
    -c "$outfile_counts" $param > $outfile_stat
    echo "* Compressing ${outfile_fa}..."
    [ -f "$outfile_fa_zipped" ] || srun MFCompressC -o "$outfile_fa_zipped" -3 "$outfile_fa"
    echo "* Compressing ${outfile_counts}..."
    [ -f "$outfile_counts_zipped" ] || srun bzip3 -k "$outfile_counts"

    headers+="${param},"
    n_kmers+="$(grep "number of kmers" "$outfile_stat" | awk -F ':' '{print $2;}'),"
    ust_c_lb+="," # to add in worksheet
    ust_c+="$(grep "cumulative length (CL)" "$outfile_stat" | awk -F '=' '{print $2;}'),"
    outfile_fa_sizes+="$(dim "$outfile_fa"),"
    outfile_fa_zipped_sizes+="$(dim "$outfile_fa_zipped"),"
    outfile_counts_sizes+="$(dim "$outfile_counts"),"
    outfile_counts_zipped_sizes+="$(dim "$outfile_counts_zipped"),"

    cd ..
  done
  # save results_csv
  printf "%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n,\n%s\n%s\n,\n,\n" \
  "$headers" "$n_kmers" "$ust_c_lb" "$ust_c" \
  "$outfile_fa_sizes" "$outfile_fa_zipped_sizes" "$outfile_counts_sizes" "$outfile_counts_zipped_sizes" \
  "$gap" "$ratio_counts_ust_out" >> "$results_csv"
  cd ..
done

cd ..

# all results in one file
cat ./*.csv > results_all.csv
