##RNASequel

RNAsequel runs as a post-processing step on top of an RNA-seq aligner and systematically detects and corrects many common alignment artifacts. Its key innovations are a two-pass splice junction detection system that combines de novo splice junction prediction with a splice junction database, and the use of an empirically determined estimate of the fragment size distribution for use in resolving read pairs. 

##RNASequel Dependencies
- Boost (www.boost.org)
- Samtools (https://github.com/samtools/htslib)

##Installation

The easiest way to obtain the RNASequel source code is to download the latest release LINK.

RNASequel can be built by typing 

```bash
cd <RNASequel Directory>
make
cp src/rnasequel <install directory>
```

The dependencies for RNASequel can be specified using the following variables:

- ```BOOST_ROOT  -- The path to the boost library root (default: /usr)```
- ```LIBBAM_ROOT -- The path to the libbam.a and the samtools header files (default: /usr)```
- ```BOOST_SUFFIX -- The version / compiler suffix used on boost library includes (Not usually necessary)```
- ```LDADD -- Extra libraries, for example, on some systems -lrt needs to be included```

For Example:

```bash
make BOOST_ROOT=/usr/local LIBBAM_ROOT=/usr/local
```


##Usage:   rnasequel [command] options

###Commands:
- index            Reference genome fasta file indexing
- transcriptome    Transcriptome index generation
- merge            Reference / Transcriptome alignment merging

Additional command line options can be viewed by using the -h flag for example:

```bash
rnasequel merge -h
```

##Notes:
For merging the alignments RNASequel requires that the bam files be sorted lexicographically using the same sorted scheme as samtools sort -n
This can be most easily accomplished by renaming the reads 1..N prior to the alignment bamfiles can also be sorting using
samtools sort -n prior to running the rnasequel merge command.

##Example Usage with BWA-mem
```bash
#Index the genome fasta file this only has to be done once
rnasequel index genome.fa
bwa index genome.fa

#Generate a transcriptome using a genes.gtf file and a denovo_alignment by STAR or another spliced read aligner
rnasequel transcriptome -g genes.gtf -r genome.fa –n 76 -b denovo_alignment.bam -o tx

# Index the transcriptome using BWA
bwa index tx.fa

# Map read 1 and 2 individually to the reference genome
bwa mem –L 2,2 -c 20000 -M -k 15 -a -t 8 -B 2 genome.fa {reads1 or 2} | samtools view -bS -F 4 - > {ref 1 or 2.bam}

# Map read 1 and 2 individually to the transcriptome
bwa mem –L 2,2 -c 20000 -M -k 15 -a -t 8 -B 2 tx.fa {read 1 or 2} | samtools view -bS -F 4 - > {juncs 1 or 2.bam}
```
