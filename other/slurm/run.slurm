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

# run original UST
run_ust(){
  mkdir -p "ust" && cd "ust"
  ## file names
  outfile_stat="${bcalm_unitigs}.stats.txt"
  outfile_fa="${bcalm_unitigs}.ust.fa"
  outfile_fa_zipped="${outfile_fa}.mfc"
  outfile_counts="${bcalm_unitigs}.ust.counts"
  outfile_counts_zipped="${outfile_counts}.bz3"

  echo "* Running UST..."
  [[ -f "$outfile_fa" ]] || srun ust -k 31 -a 1 -i "$bcalm_unitigs_abs"
  echo "* Compressing ${outfile_fa}..."
  [[ -f "$outfile_fa_zipped" ]] || srun MFCompressC -3 -o "$outfile_fa_zipped" "$outfile_fa"
  echo "* Compressing ${outfile_counts}..."
  [[ -f "$outfile_counts_zipped" ]] || srun bzip3 "$outfile_counts"

  # harvest graph properties
  headers+="ust,"
  n_kmers+="$(grep "N_KMER" "$outfile_stat" | awk -F '=' '{print $2;}'),"
  n_unitigs+="$(grep "V_BCALM" "$outfile_stat" | awk -F '=' '{print $2;}'),"
  n_iso+="$(grep -m 1 "N_ISOLATED" "$outfile_stat" | awk -F '=' '{print $2;}')," # match only first occurrence
  n_arcs+="$(grep "E_BCALM" "$outfile_stat" | awk -F '=' '{print $2;}'),"
  ust_c_lb+="$(grep "C_LB" "$outfile_stat" | awk -F '=' '{print $2;}')," # to add in worksheet


  # harvest file sizes
  ust_c+="$(grep "C_USTITCH_twoway" "$outfile_stat" | awk -F '=' '{print $2;}'),"
  outfile_fa_sizes+="$(dim "$outfile_fa"),"
  outfile_fa_zipped_sizes+="$(dim "$outfile_fa_zipped"),"
  outfile_counts_sizes+="$(dim "$outfile_counts"),"
  outfile_counts_zipped_sizes+="$(dim "$outfile_counts_zipped"),"

  cd ..
}

# exit if errors occur
set -e
# expand globs
shopt -s nullglob

ROOT_DIR="/nfsd/bcb/bcbg/prog_9"

cd "${ROOT_DIR}/"

# ust & ustar & MFCompressC
PATH="${ROOT_DIR}/UST:${ROOT_DIR}/USTAR:${ROOT_DIR}/bin:$PATH"

params=("-sf -xf" "-sf -xf -e flip_rle"\
        "-s-aa -x=a" "-s-aa -x=a -e flip_rle"\
        "-s-l x+l" "-s-l x+l -e flip_rle"
)

mkdir -p "results" && cd "results"

for bcalm_unitigs_abs in "$ROOT_DIR"/SRR/SRR*/SRR*".a2.unitigs.fa"; do
  # derive strings
  bcalm_unitigs=$(basename "$bcalm_unitigs_abs")
  seq_name=${bcalm_unitigs:0:-14} # remove ".ax.unitigs.fa" from "SRR001665_2.unitigs.fa"
  seq_dir=$(echo "$seq_name" | awk -F '_' '{printf $1};') # remove _1 or _2 if any

  results_csv="${ROOT_DIR}/results/${seq_name}/${seq_name}.csv"

  echo "*** Processing $bcalm_unitigs ($seq_dir:$seq_name) >>> $results_csv"

  mkdir -p "$seq_name" && cd "$seq_name"

  ## reset outputs vars
  # data
  headers="${seq_name},"
  n_kmers="n_kmers,"
  n_unitigs="n_unitigs,"
  n_iso="n_iso,"
  n_arcs="n_arcs,"
  ust_c_lb="CL_LB,"
  ust_c="CL,"
  outfile_fa_sizes="fasta size,"
  outfile_fa_zipped_sizes="fasta zipped size,"
  outfile_counts_sizes="counts size,"
  outfile_counts_zipped_sizes="counts zipped size,"
  # ratios (to compute in worksheet)
  gap="gap (lower is better),"
  ratio_counts_ust_out="ratio counts,"

  # run ust here
  run_ust

  ## file names
  outfile_stat="${bcalm_unitigs}.stats.txt"
  outfile_fa="${bcalm_unitigs}.ustar.fa"
  outfile_fa_zipped="${outfile_fa}.mfc"
  outfile_counts="${bcalm_unitigs}.ustar.counts"
  outfile_counts_zipped="${outfile_counts}.bz3"

  for param in "${params[@]}"; do
    echo "** params: $param"

    mkdir -p "m${param}" && cd "m${param}"

    echo "* Running USTAR..."
    # shellcheck disable=SC2086
    [ -f "$outfile_fa" ] || srun ustar -k 31 -i "$bcalm_unitigs_abs" -o "${outfile_fa}" \
      $param -c "$outfile_counts" > $outfile_stat
    echo "* Compressing ${outfile_fa}..."
    [ -f "$outfile_fa_zipped" ] || srun MFCompressC -o "$outfile_fa_zipped" -3 "$outfile_fa"
    echo "* Compressing ${outfile_counts}..."
    [ -f "$outfile_counts_zipped" ] || srun bzip3 "$outfile_counts"

    headers+="${param},"
    ust_c+="$(grep "cumulative length (CL)" "$outfile_stat" | awk -F ':' '{print $2;}' | tr -d " "),"
    outfile_fa_sizes+="$(dim "$outfile_fa"),"
    outfile_fa_zipped_sizes+="$(dim "$outfile_fa_zipped"),"
    outfile_counts_sizes+="$(dim "$outfile_counts"),"
    outfile_counts_zipped_sizes+="$(dim "$outfile_counts_zipped"),"

    cd ..
  done

  # save results_csv
  printf "%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n,\n%s\n%s\n,\n,\n" \
    "$headers" "$n_kmers" "$n_unitigs" "$n_iso" "$n_arcs" "$ust_c_lb" "$ust_c" \
    "$outfile_fa_sizes" "$outfile_fa_zipped_sizes" "$outfile_counts_sizes" "$outfile_counts_zipped_sizes" \
    "$gap" "$ratio_counts_ust_out" > "$results_csv"
  cd ..
done

cd "${ROOT_DIR}/results"

# all results in one file
cat "${ROOT_DIR}/results"/*/*.csv > results_all.csv

echo "Done. Results are saved in $(realpath results_all.csv)"
