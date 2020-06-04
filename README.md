# PipeRiboseq
A comprehensive pipeline for Ribo-seq data analysis
## Software prerequisites
This pipeline is designed to run on Linux servers, and requires the following softwares:
```
Python2
STAR
bowtie2
bedtools
samtools
salmon
FeatureCount (from Subread)
fastqc (optional)
cufflinks (optional)
ribotish (optional, it requires Python2)
```
Besides the pipeline script PipeRiboseq.sh, dependencies are in ./bin folder

Two UCSC tools (from http://hgdownload.soe.ucsc.edu/admin/exe/) are used: bedGraphToBigWig and bigWigToBedGraph. Other scripts were generated from this project.

To save time, you can directly use STAR and featureCounts program in the ./bin folder (just add it to $PATH), without installing them again.

## Pipeline setup

Here is an example of mm10 genome setup.

1, Download scripts from github to Linux server:

```git clone https://github.com/sunyumail93/PipeRiboseq.git
mv PipeRiboseq PipelineHomeDir```

2, Set up index files for genome mapping

2a, Download whole genome fasta sequence from UCSC goldenpath:

`cd PipelineHomeDir/mm10/Sequence
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip *`

2b, Set up index files:
`cd PipelineHomeDir/mm10/Index`
