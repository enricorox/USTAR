#!/bin/bash
#SBATCH --job-name run-ustar-nofilter
#SBATCH --output out_run-ustar-nofilter.txt
#SBATCH --error err_run-ustar-nofilter.txt
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH
#SBATCH --time 48:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --partition allgroups
#SBATCH --mem 200G

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
  [[ -f "$outfile_fa" ]] || srun ust -k "$k" -a 1 -i "$bcalm_unitigs_abs"
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

  # print DBG info
  printf "%s\n%s\n%s\n%s\n%s\n%s\n" "$seq_name" "$n_kmers" "$n_unitigs" "$n_iso" "$n_arcs" "$ust_c_lb" > "../profile-k${k}.csv"
  cd ..
}

# exit if errors occur
set -e
# expand globs
shopt -s nullglob

source ../configure

cd "${ROOT_DIR}/"

# default params
params=("-s+aa -x=a"\
        "-s+aa -x-c"\
        "-s+aa -x=a -e avg_flip_rle"\
        "-s+aa -x-c -e avg_flip_rle"\
        "-s+aa -x=ma" "-s+aa -x=ma -e avg_flip_rle"
)

# only the best
params=("-s+aa -x-c")

RESULTS_DIR="$ROOT_DIR/results-nofilter"
mkdir -p $RESULTS_DIR && cd $RESULTS_DIR

kmer_sizes=(13 17)
for k in "${kmer_sizes[@]}"; do
  for bcalm_unitigs_abs in "$ROOT_DIR"/../SRR/SRR*/SRR*".a1.k${k}.unitigs.fa"; do
    # derive strings
    bcalm_unitigs=$(basename "$bcalm_unitigs_abs")
    seq_name=${bcalm_unitigs:0:-18} # remove ".ax.kyz.unitigs.fa" from "SRR001665_2.ax.kyz.unitigs.fa"
    seq_dir=$(echo "$seq_name" | awk -F '_' '{printf $1};') # remove _1 or _2 if any

    results_csv="$RESULTS_DIR/${seq_name}/${seq_name}-k${k}.csv"

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
    outfile_fa_sizes="fasta_size,"
    outfile_fa_zipped_sizes="fasta_zipped_size,"
    outfile_counts_sizes="counts size,"
    outfile_counts_zipped_sizes="counts_zipped_size,"
    # ratios (to compute in worksheet)
    gap="gap (lower is better),"
    ratio_counts_ust_out="ratio_counts,"

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
      [ -f "$outfile_fa" ] || srun ustar -k "$k" -i "$bcalm_unitigs_abs" -o "${outfile_fa}" \
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
    printf "%s\n%s\n%s\n%s\n%s\n%s\n,\n%s\n%s\n,\n,\n" \
      "$headers" "$ust_c" \
      "$outfile_fa_sizes" "$outfile_fa_zipped_sizes" "$outfile_counts_sizes" "$outfile_counts_zipped_sizes" \
      "$gap" "$ratio_counts_ust_out" > "$results_csv"
    cd ..
  done

  # all results in one file
  cat ${RESULTS_DIR}/*/*"-k${k}.csv" > "${RESULTS_DIR}/results_k${k}.csv"
done

echo "Done. Results are saved in $(realpath results_all.csv)"
