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

## References

Please cite [our paper](https://link.springer.com/chapter/10.1007/978-981-99-7074-2_16) in your research:
```
@InProceedings{10.1007/978-981-99-7074-2_16,
author="Rossignolo, Enrico and Comin, Matteo",
editor="Guo, Xuan and Mangul, Serghei and Patterson, Murray and Zelikovsky, Alexander",
title="USTAR: Improved Compression of k-mer Sets with Counters Using de Bruijn Graphs",
booktitle="Bioinformatics Research and Applications",
year="2023",
publisher="Springer Nature Singapore",
address="Singapore",
pages="202--213",
abstract="A fundamental operation in computational genomics is to reduce the input sequences to their constituent k-mers. Finding a space-efficient way to represent a set of k-mers is important for improving the scalability of bioinformatics analyses. One popular approach is to convert the set of k-mers into a de Bruijn graph and then find a compact representation of the graph through the smallest path cover.",
isbn="978-981-99-7074-2"
}
```

