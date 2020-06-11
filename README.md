# PipeRiboseq
A comprehensive pipeline for Ribo-seq data analysis
## Software prerequisites
This pipeline is designed to run on Linux servers, and requires the following softwares:
```
R
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

```
git clone https://github.com/sunyumail93/PipeRiboseq.git
mv PipeRiboseq PipelineHomeDir
```

2, Set up index files for genome mapping

2a, Download whole genome fasta sequence from UCSC goldenpath:

```
cd PipelineHomeDir/mm10/Sequence
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip *
```

2b, Set up index files:
```
cd PipelineHomeDir/mm10/Index

#STAR index:
mkdir STARIndex
STAR --runMode genomeGenerate --genomeDir STARIndex --genomeFastaFiles ../Sequence/mm10.fa --sjdbGTFfile ../Annotation/mm10.RefSeq.reduced.bed12.geneid.gtf --sjdbOverhang 100

#salmon index:
salmon index -t ../Sequence/mm10.RefSeq.reduced.bed12.fa -i SalmonIndex --type quasi -k 31

#miRNA and rRNA bowtie2 index:
mkdir miRNAIndex
mkdir rRNAIndex
bowtie2-build ../../Sequence/mm10.rRNA.fa ./rRNAIndex/rRNAIndex
bowtie2-build ../../Sequence/mm10.miRNA.fa ./miRNAIndex/miRNAIndex
```

## Pipeline components
```
PipelineHomeDir/
    ├── PipeRiboseq.sh
    ├── bin/
    └── mm10/
      └── Annotation/
        ├── mm10.RefSeq.reduced.bed12
        ├── mm10.RefSeq.reduced.mRNA.bed12
        ├── mm10.RefSeq.reduced.bed12.geneid.gtf
        └── mm10.uniqMatching.txt
      └── Index/
        ├── miRNAIndex/
        ├── rRNAIndex/
        ├── SalmonIndex/
        └── STARIndex/
      └── Sequence/
        ├── mm10.fa
        ├── mm10.fai
        ├── mm10.ChromInfo.txt
        ├── mm10.miRNA.fa
        ├── mm10.rRNA.fa
        └── mm10.RefSeq.reduced.bed12.fa
```

Notes: 
1, For Annotation folder, download GTF file from UCSC table browser. `reduced`: Only one location was chosen when one gene duplicates at multiple genomic loci.

2, `uniqMatching.txt` file contains one-to-one matching from transcript to gene name.

3, For Index folder, indexes are not included in this github directory, but need to be created during set up.

4, For Sequence folder, `RefSeq.reduced.bed12.fa` converts from `RefSeq.reduced.bed12` file using bedtools. genome.fa file also needs to be downloaded from UCSC goldenpath.

## Usage

Type the pipeline name, then you will see the manual page:

`PipeRiboseq.sh`

Manual page:

![](images/Usages.png)

The trimming pipeline also has manual page:

`PipeSETrimmer.sh`

## Examples

A regular run using mostly default parameters:

`PipeRiboseq.sh -i Data.fastq.gz -g mm10 -normCDS`

More parameters used, and plot given genes in list file (mRNAs in the list file must be Refseq mRNA ID: NM_xxx):

`PipeRiboseq.sh -i Data.fastq.gz -g mm10 -noqc -noriboqc -p 4 -normCDS -m 3 -plotRNA list`

## Run a real data to test the pipeline

1, Download data

Use a public dataset: [GEO SRA: SRR989509](https://www.ncbi.nlm.nih.gov/sra/?term=SRR989509&utm_source=gquery&utm_medium=search)

`fastq-dump --split-3 SRR989509`

2, Trim adaptor using the trimming script from this repository: PipeSETrimmer.sh

`PipeSETrimmer.sh`

3, Run Ribo-seq pipeline:

`PipeRiboseq.sh -i Data.trimmed.fastq.gz -g mm10 -p 4 -normCDS`

## Outputs

1, Length distribution of Ribo-seq mapped reads:

![](images/Lendis.png)

2, Ribo-seq corrected 5´-ends around top expressed mRNAs' AUG:

![](images/Meta.AUG.png)

3, Ribo-seq corrected 5´-ends around top expressed mRNAs' STOP codon:

![](images/Meta.STOP.png)

4, Ribo-seq QC by [RiboTISH](https://github.com/zhpn1024/ribotish):

![](images/RiboTISHQC.png)
