#!/bin/sh

# input arguments
n=$1
dataset=$2

#specify multimapper handling
# set to uniq for only uniq mapped reads
# set to phased for only phased mapped reads
# set to all for all reads
multimapperHandling=$3

# directory paths and names
cwd="/users/lserrano/xhernandez/tRNA_mapping"
scriptDir="${cwd}/"
genomeDir="${cwd}/Data/Genomes/H.sapiens"
workDir="${cwd}/Analysis/${dataset}"
genomeName="hg38.genomic"
tRNAName="hg38.tRNAscan"

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


bn=$(basename $n _read1_filtered.fastq.gz)

#post-mapping against cluster
$segemehl --evalue 500 --differences 3 --maxinterval 1000 --accuracy 85 --index ${genomeDir}/${tRNAName}_cluster.idx --database ${genomeDir}/${tRNAName}_cluster.fa --nomatchfilename ${bn}_unmatched.fastq --query $n --mate ${workDir}/mapping/${bn}/${bn}_read2_filtered.fastq.gz | $samtools view -bS | $samtools sort -T ${bn} -o ${bn}.bam
gzip ${bn}_unmatched.fastq

##preparing bam file for indel realignment
#indexing
$samtools index ${bn}.bam
$picard BuildBamIndex I=${bn}.bam O=${bn}.bai

#add read groups to bam file
$picard AddOrReplaceReadGroups I=${bn}.bam O=${bn}.mod.bam RGPL=RNASeqReadSimulator RGLB=Simlib RGPU=unit1 RGSM=36bam

#indexing
$samtools index ${bn}.mod.bam
$picard BuildBamIndex I=${bn}.mod.bam O=${bn}.mod.bai

#modify mapping quality to 60 (otherwise all were removed)
$gatk -T PrintReads -R ${genomeDir}/${tRNAName}_cluster.fa -I ${bn}.mod.bam -o ${bn}.temp.bam -rf UnmappedRead -rf ReassignMappingQuality --default_mapping_quality 60
mv -f ${bn}.temp.bam ${bn}.mod.bam
rm -f ${bn}.mod.bai ${bn}.mod.bam.bai

#indexing
$samtools index ${bn}.mod.bam
$picard BuildBamIndex I=${bn}.mod.bam O=${bn}.mod.bai

##realignment
$gatk -T RealignerTargetCreator -R ${genomeDir}/${tRNAName}_cluster.fa -I ${bn}.mod.bam -o ${bn}.temp.intervals
$gatk -T IndelRealigner -R ${genomeDir}/${tRNAName}_cluster.fa -I ${bn}.mod.bam -targetIntervals ${bn}.temp.intervals -o ${bn}.realigned.bam --maxReadsForRealignment 500000
rm -f ${bn}.temp.intervals

##filter multimapped reads
if [ "$multimapperHandling" == "uniq" ]; then
  $samtools view -h ${bn}.realigned.bam | grep -P 'NH:i:1\D'\|'^@' | $samtools view -bS | $samtools sort -T ${bn} -o ${bn}.mmHandled.bam
elif [ "$multimapperHandling" == "quant" ]; then
  $samtools sort -n -T ${bn} -O sam ${bn}.realigned.bam  > ${bn}.nSorted.sam
  python ${scriptDir}/distribute_reads.py ${bn}.nSorted.sam ${genomeDir}/${tRNAName}_clusterInfo.fa ${bn}.nSorted.quant.sam ${bn}.nSorted.multimappers.sam
  $samtools view -bS ${bn}.nSorted.quant.sam | samtools sort -T ${bn} -o ${bn}.mmHandled.quant.bam
  $samtools view -h ${bn}.realigned.bam | grep -P 'NH:i:1\D'\|'^@' | $samtools view -bS | $samtools sort -T ${bn} -o ${bn}.mmHandled.bam
elif [ "$multimapperHandling" == "phased" ]; then
  $samtools sort -n -T ${bn} -O sam ${bn}.realigned.bam  > ${bn}.nSorted.sam
  perl ${scriptDir}/multimapperPhasing.pl -ed 0 -id 0 -verbose 0 -sam ${bn}.nSorted.sam -out ${bn}.nSorted.phased.sam
  $samtools view -bS ${bn}.nSorted.phased.sam | samtools sort -T ${bn} -o ${bn}.mmHandled.bam
elif [ "$multimapperHandling" == "all" ]; then
  cp ${bn}.realigned.bam ${bn}.mmHandled.bam
else
  echo "Unkown parameter for multimapperHandling; set to 'uniq', 'all', or 'phased'";
  exit;
fi

#indexing
$samtools index ${bn}.mmHandled.bam
$picard BuildBamIndex I=${bn}.mmHandled.bam O=${bn}.mmHandled.bai

##modification site calling
$gatk -R ${genomeDir}/${tRNAName}_cluster.fa -T UnifiedGenotyper -I ${bn}.mmHandled.bam -o ${bn}.GATK.vcf -stand_call_conf 50.0
grep -i -v lowqual ${bn}.GATK.vcf > ${bn}.GATK_filtered.vcf
    
## Allele-specific expression analysis
$samtools idxstats ${bn}.mmHandled.quant.bam > ${bn}.expression.txt
$gatk -R ${genomeDir}/${tRNAName}_cluster.fa -T ASEReadCounter -I ${bn}.mmHandled.bam -o ${bn}.ASE.csv -sites ${bn}.GATK_filtered.vcf

## Compute expression
python ${scriptDir}/compute_expression.py ${bn}.expression.txt ${bn}.ASE.csv ${genomeDir}/${tRNAName}_clusterInfo.fa ${genomeDir}/tRNAs.txt ${bn}.RPM.csv

## Stats
$samtools stats ${bn}.mmHandled.quant.bam | grep ^RL | cut -f 2- > ${bn}.readlengths.txt
$samtools depth ${bn}.mmHandled.quant.bam > ${bn}.depth.txt
