#!/bin/sh

#directory paths
cwd=$(pwd)

#specify multimapper handling
# set to uniq for only uniq mapped reads
# set to phased for only phased mapped reads# set to all for all reads
multimapperHandling="uniq"


#change here your directory paths
ngsDir="${cwd}/TestData/NGS/H.sapiens"
scriptDir="${cwd}/"
genomeDir="${cwd}/TestData/Genomes/H.sapiens"
workDir="${cwd}/TestAnalysis/"
adapterFile="${cwd}/TestData/NGS/adapter.fa"

#program paths
#change here your program paths
bbduk="bbduk.sh"				
#BBMap version 38.22
fastqc="fastqc"					
#0.11.4
samtools="samtools"				
#1.3.1
tRNAscanSE="tRNAscan-SE"			
#2.0
bedtools="bedtools"				
#2.27.1
segemehl="segemehl.x"				
#0.3.1
picard="picard"					
#2.18.17
gatk="gatk3"					
#3.8

# variables
#change here the names of the variable
genomeName="hg38.genomic"
tRNAName="hg38.tRNAscan"

###pre-mapping against artificial genome
mkdir -p ${workDir}/mapping
cd ${workDir}/mapping

for n in $(ls ${ngsDir}/*trimmed.fastq.gz)
do
    bn=$(basename $n _trimmed.fastq.gz)

    $segemehl --evalue 500 --differences 3 --maxinterval 1000 --accuracy 80 --index ${genomeDir}/${genomeName}_artificial.idx --database ${genomeDir}/${genomeName}_artificial.fa --nomatchfilename ${bn}_unmatched.fastq --query $n -o ${bn}.sam
    gzip ${bn}_unmatched.fastq

    ##remove all reads mapping at least once to the genome
    perl ${scriptDir}/removeGenomeMapper.pl ${genomeDir}/${tRNAName}_pre-tRNAs.fa ${bn}.sam ${bn}_filtered.sam

    ##remove pre-tRNA reads, keep only mature tRNA reads
    perl ${scriptDir}/removePrecursor.pl  ${genomeDir}/${tRNAName}_pre-tRNAs.bed12 ${bn}_filtered.sam $n > ${bn}_filtered.fastq
    gzip ${bn}_filtered.fastq
done
