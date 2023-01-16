#!/bin/bash

#You need a reference file Bos_taurus.fa and the sequence reads
#Download a vcf file of all known variants
#gunzip the file
#gunzip bos_taurus.vcf.gz
#Run the quality control (multiqc) for the sequence reads
Multiqc . {This is a command line tool that can be used to generate a single report summarizing the
results of multiple sequences.}

#Index the reference genome
hisat2-build Bos_taurus.fa Bos_taurus.idx {This command is used to build a HISAT2 index of a
reference genome in FASTA format (Bos_taurus.fa) and save it as Bos_taurus.idx. The index is used
to align reads to the reference genome with HISAT2.}

samtools faidx Bos_taurus.fa {The command &quot;samtools faidx&quot; is used to index a FASTA file, in this
case &quot;Bos_taurus.fa&quot;. The indexed file can then be used with other tools, such as &quot;samtools&quot; to more
efficiently retrieve sequences from the genome.}

#Create sequence dictionary
gatk CreateSequenceDictionary --REFERENCE Bos_taurus.fa - -OUTPUT Bos_taurus.dict {This is
used to create a sequence dictionary file for a reference genome in FASTA format. The command
takes two arguments:
REFERENCE: The input FASTA file (Bos_taurus.fa)
OUTPUT: The name of the output sequence dictionary file (Bos_taurus.dict).
It creates a dictionary file, which contains metadata about the reference genome, such as the names
and lengths of the contigs. This file is required by some GATK tools and other bioinformatics tools
that work with the reference genome.}

#Create an index file for a VCF
gatk IndexFeatureFile -I bos_taurus.vcf {This command is used to create an index file for a VCF
(variant call format) file. The &quot;-I&quot; flag is used to specify the input VCF file, in this case
&quot;bos_taurus.vcf&quot;. This command is typically used to improve the performance of GATK when
working with large VCF files. The index file allows GATK to quickly locate and process specific regions
of the VCF file rather than having to scan the entire file every time.}

#Create a txt file called bos_accession.txt and copy all your accession numbers into that file

#Run the following for loop.

for sample in `cat bos_accession.txt`;
do
#Define your variables
R1=&quot;${sample}&quot;_1.fastq
R2=&quot;${sample}&quot;_2.fastq
#Align and get a sorted bam file output. Output it to HISAT2 dir
hisat2 -x Bos_taurus.idx -1 $R1 -2 $R2\
--rg-id RG1 --rg &quot;PL:ILLUMINA&quot;\
--rg SM:&quot;${sample}&quot; | samtools sort -\
-o HISAT2/&quot;${sample}&quot;_sorted.bam

#Index the sorted.bam files
samtools index HISAT2/&quot;${sample}&quot;_sorted.bam
#Mark duplicates
gatk MarkDuplicates\
--INPUT HISAT2/&quot;${sample}&quot;_sorted.bam\
--OUTPUT HISAT2/&quot;${sample}&quot;_sorted_dedup.bam\
--METRICS_FILE HISAT2/&quot;${sample}&quot;_dup_metrics.txt\
--REMOVE_DUPLICATES true --CREATE_INDEX true

#SplitNCigar Reads
gatk SplitNCigarReads\

-I HISAT2/&quot;${sample}&quot;_sorted_dedup.bam\
-O HISAT2/&quot;${sample}&quot;_cigarReads.bam -R Bos_taurus.fa

#Base recalibration
gatk BaseRecalibrator\
--input HISAT2/&quot;${sample}&quot;_cigarReads.bam\
--output HISAT2/&quot;${sample}&quot;_baseRecal.table\
--reference Bos_taurus.fa\
--known-sites bos_taurus.vcf


#Apply BQSR
gatk ApplyBQSR\

--bqsr-recal-file HISAT2/&quot;${sample}&quot;_baseRecal.table\
--input HISAT2/&quot;${sample}&quot;_sorted_dedup.bam\
--output HISAT2/&quot;${sample}&quot;_BQSR.bam\
--reference Bos_taurus.fa

#Variant calling
gatk HaplotypeCaller\
--reference Bos_taurus.fa\
--input HISAT2/&quot;${sample}&quot;_BQSR.bam\
--output HISAT2/&quot;${sample}&quot;_haplotype.g.vcf.gz\
--ERC GVCF

done
