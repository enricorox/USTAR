#!/bin/bash
set -e
ROOT_DIR="/nfsd/bcb/bcbg/prog_9/USTAR"
PATH="/nfsd/bcb/bcbg/prog_9/bin/sratoolkit.3.0.1-ubuntu64/bin":$PATH

sequences=$(grep -vE "#" $ROOT_DIR/other/datasets/sequences.txt)

cd $ROOT_DIR/..
mkdir -p SRR && cd SRR

for s in $sequences; do
  echo
  echo "*** Downloading $s..."
  prefetch --progress "$s"
  echo "*** $s downloaded!"
done