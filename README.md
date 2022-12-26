# USTAR (Unitig STitch STar)
## What is that?
Another kmers set compressor, with counting.

It is based on the ideas of [UST](https://github.com/medvedevgroup/UST) and [prophAsm](https://github.com/prophyle/prophasm) 
for computing an SPSS representation (aka simplitigs) for the given kmers set.

## Dependencies
There are no dependencies. However, you'll need [BCALM2](https://github.com/GATB/bcalm) 
in order to compute a cDBG of your multi-fasta file.

## How to compile
Just run `cmake .` and `make`.

## How to run
See the help `./ustar -h`.