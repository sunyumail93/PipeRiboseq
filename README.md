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
Besides the pipeline script PipeRiboseq.sh, dependencies are in /bin folder

Two UCSC tools (from http://hgdownload.soe.ucsc.edu/admin/exe/) are used: bedGraphToBigWig and bigWigToBedGraph. Other scripts were generated from this project.

To safe time, you can directly use STAR and featureCounts program in /bin folder (just add it to $PATH), rathan than installing them again.
