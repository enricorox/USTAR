# USTAR (Unitig STitch Advanced constRuction)
## Overview
USTAR is a kmers set compressor, with counting.

It is based on the ideas of [UST](https://github.com/medvedevgroup/UST) 
and [prophAsm](https://github.com/prophyle/prophasm) 
for computing an SPSS representation (aka simplitigs) for the given kmers set.

You will find four executable and one bash script:
* `ustar`: the main program
* `ustarx`: a kmers extractor
* `ustars`: compute kmers statistics
* `ustar-test`: used for debug
* `validate`: a validation script

## Dependencies
There are no dependencies. 
However, you'll need [BCALM2](https://github.com/GATB/bcalm) 
in order to compute a compacted de Brujin graph (cdBG) of your multi-fasta file.

## How to download and compile
* `git clone https://github.com/enricorox/USTAR`.
* `cd USTAR`
* `cmake . && make -j 4`.

## How to run USTAR
Run BCALM2 first: 
* `./bcalm -kmer-size <kmer-size> -in <your-multi-fasta> -all-abundance-counts`

Then ustar:
* `./ustar -k <kmer-size> -i <bcalm-output>`

To use the best heuristic, add `-s+aa -x-c`

See the help `./ustar -h` for details and advanced options.

## How to validate the output
You can check that the output file contains the same kmers of
your bcalm file with your preferred kmer counter.

If you want to check that __kmers and counts__ are correct,
* `./ustarx -k <kmer-size> -i <ustar-fasta> -c <ustar-counts> -s`
* `./validate <kmer-size> <your-multi-fasta> <ustar-kmers-counts>` 

Note that you'll need to install [Jellyfish-2](https://github.com/zippav/Jellyfish-2) in order to use `validate`.
