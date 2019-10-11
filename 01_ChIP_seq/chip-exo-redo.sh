#!/usr/bin/env bash

#
# Pipeline to reproduce the data analysis for the paper:
#
# Genome-Wide Mapping of Binding Sites Reveals Multiple Biological Functions of the Transcription Factor Cst6p in Saccharomyces cerevisiae
#
# http://mbio.asm.org/content/7/3/e00559-16.long#F1
#

#
# This script requires Entrez-Direct, bwa, bamUtils and GEM installed.
#

# Bail on errors.
set -ueo pipefail

# All shortcuts to programs and data are defined here.
GEM=~/src/gem/gem.jar
DIST=~/src/gem/Read_Distribution_ChIP-exo.txt

# Reference genome.
REF=refs/saccer3.fa

# Glucose samples.
GLU1=bam/trimmed-SRR3033154.bam
GLU2=bam/trimmed-SRR3033155.bam

# Ethanol samples.
ETH1=bam/trimmed-SRR3033156.bam
ETH2=bam/trimmed-SRR3033157.bam

# Remote urls.
URL=http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz

# Directories to organize data in.
mkdir -p refs
mkdir -p data
mkdir -p bam

# Download the chromosomes and build a genome.
# The sequences are by chromosome.
curl $URL | tar zxv
curl http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes > refs/sacCer3.chrom.sizes

# Move the files
mv *.fa refs

# Create genome.
cat refs/chr*.fa > $REF

# Index the reference.
bwa index $REF
samtools faidx $REF

# Download the sequencing data in FASTQ format.

# Get run info from SRA.
esearch -db sra -query PRJNA306490  | efetch -format runinfo > runinfo.csv

# Isolate just the run ids for the run.
cat runinfo.csv | cut -f 1 -d , | grep SRR > runids.txt

# Download the fastq files.
cat runids.txt | parallel --eta --verbose  "fastq-dump -O data --split-files -F {}"

# Run the alignments in paired end mode.
#cat runids.txt | parallel --eta --verbose "bwa mem -t 8 $REF data/{}_1.fastq data/{}_2.fastq | samtools sort -@ 8 > bam/{}.bam"

# The data is aligned in single end mode and only file 1 is used.
cat runids.txt | parallel --eta --verbose "bwa mem -t 8 $REF data/{}_1.fastq | samtools sort -@ 8 > bam/{}.bam"

# Index the results in BAM format.
cat runids.txt | parallel --eta --verbose "samtools index bam/{}.bam"

# Trim the BAM file to keep left/right borders.
cat runids.txt | parallel --eta --verbose "bam trimBam bam/{}.bam bam/temp-{}.bam -R 70 --clip"

# Sort trimmed BAM file.
cat runids.txt | parallel --eta --verbose "samtools sort -@ 8 bam/temp-{}.bam > bam/trimmed-{}.bam"

# Get rid of temporary BAM files.
rm -f bam/temp*

# Reindex trimmed bam files.
cat runids.txt | parallel --eta --verbose "samtools index bam/trimmed-{}.bam"

# Create the coverage files for all BAM files.
ls bam/*.bam | parallel --eta --verbose "bedtools genomecov -ibam {} -g $REF.fai -bg | sort -k1,1 -k2,2n > {.}.bedgraph"

# Generate all bigwig coverages from bedgraphs.
ls bam/*.bedgraph | parallel --eta --verbose "bedGraphToBigWig {} $REF.fai {.}.bw"

# Run GEM on glucose samples.
java -Xmx3G -jar $GEM --d $DIST --g refs/sacCer3.chrom.sizes --genome refs --exptCond1 $GLU1  --exptCond1 $GLU2 --f SAM --out glucose --outBED --k 5  --smooth 3 --mrc 20 --k_win 30

# Run GEM on ethanol samples.
java -Xmx3G -jar $GEM --d $DIST --g refs/sacCer3.chrom.sizes --genome refs --exptCond1 $ETH1  --exptCond1 $ETH2 --f SAM --out ethanol --outBED --k 5  --smooth 3 --mrc 20 --k_win 30
