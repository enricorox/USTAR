#!/bin/bash
#SBATCH --job-name run-counters
#SBATCH --output out_run-counters.txt
#SBATCH --error err_run-counters.txt
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH
#SBATCH --time 48:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --partition allgroups
#SBATCH --mem 200G

# compute file size
dim(){
  stat -c %s "${1}"
}

# exit if errors occur
set -e
# expand globs
shopt -s nullglob

source configure

# activate conda for dsk & kmc & seqtk
. /nfsd/bcb/bcbg/rossigno/anaconda3/bin/activate

cd "${ROOT_DIR}/"

RESULTS_DIR="$ROOT_DIR/results-nofilter"
mkdir -p $RESULTS_DIR && cd $RESULTS_DIR

kmer_sizes=(21)

headers=",bcalm2,kmc,dsk,squeakr"
for k in "${kmer_sizes[@]}"; do
  for fasta_abs in "$ROOT_DIR"/../SRR/SRR*/SRR*".fasta"; do
    # derive strings
    fasta=$(basename "$fasta_abs")
    seq_name=${fasta:0:-6} # remove ".fasta" from "SRR001665_2.fasta"
    # seq_dir=$(echo "$seq_name" | awk -F '_' '{printf $1};') # remove _1 or _2 if any

    results_csv="$RESULTS_DIR/${seq_name}/${seq_name}-k${k}.csv"

    echo "*** Processing $fasta >>> $results_csv"

    mkdir -p "$seq_name" && cd "$seq_name"

    ## reset outputs vars
    # data
    row="$seq_name,"

    ## bcalm2
    echo "** Running BCALM2..."
    outfile="${fasta}.unitigs.fa"
    outfile_zipped="${outfile}.bz3"
    [ -f "$outfile" ] || srun bcalm -max-memory "$MAX_MEMORY" -nb-cores "$NTHREAD" \
      -kmer-size "$k" -abundance-min 1 -in "$fasta_abs" -all-abundance-counts -out "${outfile:0:-11}"
    echo "* Compressing ${outfile}..."
    [ -f "$outfile_zipped" ] || srun bzip3 -j "$NTHREAD" "$outfile"
    row+=$(dim "$outfile_zipped"),

    ## kmc
    echo "** Running KMC..."
    outfile="${fasta}.kmc"
    outfile_zipped="${outfile}.kmc_pre.bz3"
    [ -f "$outfile".kmc_pre ] || srun kmc -k"$k" -t"$NTHREAD" -m"$MAX_MEMORY" -fm -ci1 "$fasta_abs" "$outfile" .
    echo "* Compressing ${outfile}.kmc_{pre,suf}..."
    [ -f "$outfile_zipped" ] || srun bzip3 -j "$NTHREAD" -B "$outfile".kmc_pre "$outfile".kmc_suf
    row+=$(($(dim "$outfile_zipped") + $(dim "$outfile".kmc_suf.bz3))),

    ## dsk
    echo "** Running DSK..."
    outfile="${fasta}.dsk"
    outfile_zipped="${outfile}.h5.bz3"
    [ -f "$outfile".h5 ] || srun dsk  -nb-cores "$NTHREAD" -file "$fasta_abs" -kmer-size "$k" -max-memory "$MAX_MEMORY" -abundance-min 1 -out "$outfile"
    echo "* Compressing ${outfile}..."
    [ -f "$outfile_zipped" ] || srun bzip3 -j "$NTHREAD" "$outfile".h5
    row+=$(dim "$outfile_zipped"),

    ## squeakr
    echo "** Converting to fastq"
    fastq_abs=${fasta}.fq
    [ -f "$fastq_abs" ] || srun seqtk seq -F Q -U "$fasta_abs" > "$fastq_abs"
    outfile="$(basename "$fastq_abs").squeakr"
    outfile_zipped="${outfile}.bz3"
    echo "** Running Squeakr..."
    [ -f "$outfile" ] || [ -f "$outfile_zipped" ] || srun squeakr count -e -k "$k" -c 1 -t "$NTHREAD" --log-slots 38 -o ./"$outfile" "$fastq_abs"
    echo "* Compressing ${outfile}..."
    [ -f "$outfile_zipped" ] || srun bzip3 -j "$NTHREAD" "$outfile"
    row+=$(dim "$outfile_zipped"),

    # save results_csv
    printf "%s\n" \
      "$row" > "$results_csv"
    cd ..
  done

  # all results in one file
  printf "%s\n" \
        $headers > "${RESULTS_DIR}/results_k${k}.csv"
  cat "${RESULTS_DIR}"/*/*"-k${k}.csv" >> "${RESULTS_DIR}/results_k${k}.csv"
done

echo "Done. Results are saved in $(realpath results_all.csv)"
