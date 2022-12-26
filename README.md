# USTAR (Unitig STitch STar)
## What is that?
Another kmers set compressor, with counting.

It is based on the ideas of [UST](https://github.com/medvedevgroup/UST) 
and [prophAsm](https://github.com/prophyle/prophasm) 
for computing an SPSS representation (aka simplitigs) for the given kmers set.

## Dependencies
There are no dependencies. 
However, you'll need [BCALM2](https://github.com/GATB/bcalm) 
in order to compute a compacted de Brujin graph (cdBG) of your multi-fasta file.

## How to download and compile
* `git clone https://github.com/enricorox/USTAR`.
* `cd USTAR`
* `cmake . && make`.

## How to run USTAR
* Run BCALM2: `./bcalm -kmer-size 31 -in <your-multi-fasta> -all-abundance-counts` 
* `./ustar -k 31 -i <bcalm_file> -e avg_flip_rle`

See the help `./ustar -h` for details.

## Cool! But it is correct?
Sure! 
You can check that the output file contains the same kmers of
your multi-fasta file with your preferred kmer counter (just ignore counts because they are in a different file).

For your convenience, you can just run `./validate <kmer-size> <your-multi-fasta> <ustar-output>`,
but you'll need to install [DSK](https://github.com/GATB/dsk) and change the path in the script.
