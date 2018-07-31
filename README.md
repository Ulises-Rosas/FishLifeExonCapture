## FishLifeExonCapture

# Tutorial for assembling exon capture data across the diversity of fishes

Software required:

aTRAM

samtools (most versions should work)

BLAST+ (command line)

bwa

Velvet

Biopython

Exonerate

Mafft

TranslatorX

# Step 1: Organize your fastq files into separate directories for each species

You should establish a main working directory for all of your samples. Starting with a list of demultiplexed fastq files from the sequencer, we need to make separate directories to process our samples through the pipeline. These files are commonly gzip compressed (so you might see .fastq.gz). If you need to expand the files you can run the following:

```
for f in *gz;
do
gunzip $f;
done
```

Now make directories for each of the fastq files, and move them to those directories:

```
for f in *fastq;
do
mkdir ${f%.*};
mv $f ${f%.*}/;
done
```

# Step 2: Run Trimmomatic to quality trim the sequences

Run the trimmomatic-loop.sh script in the main working directory. 

```
../FishLifeExonCapture/trimmomatic-loop.sh
```

Note: This script will look for a file called 'adapters.fa' in the main working directory, to trim out adapter contamination. Different sets of adapters are used for different library preparations, so you should check which are appropriate. Trimmomatic has all of these sequences packaged with it, so you can just move this to the 'adapters.fa' file. The adapters are proprietary Illumina sequences, so I have not included them here.

When the script is finished, each .fastq file will have an associated .trimmed.fastq file.

# Step 3: Map raw reads back to representative bait sequences

aTRAM2.0 currently only takes a single reference sequence as input. To avoid having to choose that sequence for each individual project, and to increase the over all number of loci assembled, we do a initial de novo assembly partitioned by locus for each taxon. The product of that de novo assembly is used as the reference sequence for aTRAM in the following steps. This might be redundant in some cases, but tends to produce longer, more accurate assemblies, without having to choose a reference sequence for every taxon.

This script maps the raw reads with bwa against the reference sequences that all exon baits were designed on, as well as coding mitochondrial genes. These sequences are included in the all_Master.fasta file.
