# USTAR (Unitig STitch STar)
## What is that?
Another kmers set compressor, with counting.

It is based on the ideas of [UST](https://github.com/medvedevgroup/UST) 
and [prophAsm](https://github.com/prophyle/prophasm) 
for computing an SPSS representation (aka simplitigs) for the given kmers set.

## Dependencies
There are no dependencies. 
However, you'll need [BCALM2](https://github.com/GATB/bcalm) 
in order to compute a cDBG of your multi-fasta file.

## How to compile
Just run `cmake .` and `make`.

## How to run
See the help `./ustar -h`.

## Cool! But it is correct?
Sure! 
You can check that the output file contains the same kmers of
your multi-fasta file with your preferred kmer counter (just ignore the counts!).

For your convenience, you can just run `./validate <kmer-size> <ustar-output> <your-multi-fasta>`,
but you'll need to install [DSK](https://github.com/GATB/dsk) and change the path in the script.
