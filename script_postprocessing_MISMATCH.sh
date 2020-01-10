#!/bin/sh

# input arguments
n=$1

# directory paths and names
cwd="/users/lserrano/xhernandez/tRNA_mapping"
scriptDir="${cwd}/"
genomeDir="${cwd}/Data/Genomes/H.sapiens"
genomeName="hg38.genomic"
tRNAName="hg38.tRNAscan"

#program paths
#change here your program paths
samtools="samtools"				
#1.3.1
picard="picard"					
#2.18.17


bn=$(basename $n _filtered.fastq.gz)

## Separate reads wrt mismatches
$samtools view -h ${bn}.mmHandled.bam | grep -P '\t\d{2}=\t'\|'^@' | $samtools view -bS | $samtools sort -T ${bn} -o ${bn}.0mismatches.bam
$samtools view -h ${bn}.mmHandled.bam | grep -v -P '\t\d{2}=\t' | $samtools view -bS | $samtools sort -T ${bn} -o ${bn}.mismatches.bam
$samtools view -h ${bn}.mismatches.bam | grep -P '\t\d{0,2}={0,1}1X\d{0,2}={0,1}\t'\|'^@' | $samtools view -bS | $samtools sort -T ${bn} -o ${bn}.1mismatches.bam
$samtools view -h ${bn}.mismatches.bam | grep -v -P '\t\d{0,2}={0,1}1X\d{0,2}={0,1}\t' | $samtools view -bS | $samtools sort -T ${bn} -o ${bn}.moremismatches.bam

#indexing
$samtools index ${bn}.0mismatches.bam
$samtools index ${bn}.mismatches.bam
$samtools index ${bn}.1mismatches.bam
$samtools index ${bn}.moremismatches.bam

$picard BuildBamIndex I=${bn}.0mismatches.bam O=${bn}.0mismatches.bai
$picard BuildBamIndex I=${bn}.mismatches.bam O=${bn}.mismatches.bai
$picard BuildBamIndex I=${bn}.1mismatches.bam O=${bn}.1mismatches.bai
$picard BuildBamIndex I=${bn}.moremismatches.bam O=${bn}.moremismatches.bai

## Expression analysis
$samtools idxstats ${bn}.mmHandled.bam > ${bn}.mmHandled.expression.txt
$samtools idxstats ${bn}.0mismatches.bam > ${bn}.0mismatches.expression.txt
$samtools idxstats ${bn}.mismatches.bam > ${bn}.mismatches.expression.txt
$samtools idxstats ${bn}.1mismatches.bam > ${bn}.1mismatches.expression.txt
$samtools idxstats ${bn}.moremismatches.bam > ${bn}.moremismatches.expression.txt

## Compute expression
python ${scriptDir}/compute_expression.py ${bn}.mmHandled.expression.txt ${bn}.ASE.csv ${genomeDir}/${tRNAName}_clusterInfo.fa ${genomeDir}/tRNAs.txt ${bn}.mmHandled.RPM.csv
python ${scriptDir}/compute_expression.py ${bn}.0mismatches.expression.txt ${bn}.ASE.csv ${genomeDir}/${tRNAName}_clusterInfo.fa ${genomeDir}/tRNAs.txt ${bn}.0mismatches.RPM.csv
python ${scriptDir}/compute_expression.py ${bn}.mismatches.expression.txt ${bn}.ASE.csv ${genomeDir}/${tRNAName}_clusterInfo.fa ${genomeDir}/tRNAs.txt ${bn}.mismatches.RPM.csv
python ${scriptDir}/compute_expression.py ${bn}.1mismatches.expression.txt ${bn}.ASE.csv ${genomeDir}/${tRNAName}_clusterInfo.fa ${genomeDir}/tRNAs.txt ${bn}.1mismatches.RPM.csv
python ${scriptDir}/compute_expression.py ${bn}.moremismatches.expression.txt ${bn}.ASE.csv ${genomeDir}/${tRNAName}_clusterInfo.fa ${genomeDir}/tRNAs.txt ${bn}.moremismatches.RPM.csv

