#!/bin/sh

# input arguments
n=$1
dataset=$2
bn=$(basename $n _read1_trimmed.fastq.gz)

#directory paths
cwd="/users/lserrano/xhernandez/tRNA_mapping"

#change here your directory paths
dataset="hydroseq_PE_test"
ngsDir="${cwd}/Data/NGS/${dataset}"
scriptDir="${cwd}"
genomeDir="${cwd}/Data/Genomes/rRNA"
workDir="${cwd}/Analysis/${dataset}_rRNA"

#program paths
#change here your program paths
samtools="samtools"				
#1.3.1
segemehl="segemehl.x"				
#0.3.1
bedtools="bedtools"				
#2.27.1
picard="picard"					
#2.18.17

# variables
#change here the names of the variable
genomeName="rRNA"


$segemehl --evalue 500 --differences 3 --maxinterval 1000 --accuracy 85 --index ${genomeDir}/${genomeName}.idx --database ${genomeDir}/${genomeName}.fa --nomatchfilename ${bn}_unmatched.fastq --query $n --mate ${ngsDir}/${bn}/${bn}_read2_trimmed.fastq.gz | $samtools view -bS | $samtools sort -T ${bn} -o ${bn}.bam
gzip ${bn}_unmatched.fastq
$samtools index ${bn}.bam
$picard BuildBamIndex I=${bn}.bam O=${bn}.bai
$bedtools multicov -bams ${bn}.bam -bed ${genomeDir}/${genomeName}.bed > ${bn}.rRNA.txt
