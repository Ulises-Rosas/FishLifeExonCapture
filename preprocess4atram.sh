#!/bin/sh

# this only needs to be run once in the entire pipeline


# this assumes mixed end data. change if you're using single-end or two paired-end files. see aTRAM docs.

for f in *.fastq;
do
atram_preprocessor.py -b ${f%.*} --mixed-ends $f;
done
