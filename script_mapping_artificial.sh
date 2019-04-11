#!/bin/sh

#program paths
segemehl="segemehl.x"				
#0.3.1
# input arguments
n=$1
dataset=$2

# directory paths and names
cwd="/users/lserrano/xhernandez/tRNA_mapping"
scriptDir="${cwd}/"
genomeDir="${cwd}/Data/Genomes/H.sapiens"
workDir="${cwd}/Analysis/${dataset}"
genomeName="hg38.genomic"
tRNAName="hg38.tRNAscan"

# COpy needed files on local node
cp ${genomeDir}/{${genomeName}_artificial.idx,${genomeName}_artificial.fa,${tRNAName}_pre-tRNAs.fa,${tRNAName}_pre-tRNAs.bed12} $TMPDIR
cp $n $TMPDIR
cp ${scriptDir}/{removeGenomeMapper.pl,removePrecursor.pl} $TMPDIR

cd $TMPDIR

###pre-mapping against artificial genome
bn=$(basename $n _trimmed.fastq.gz)
$segemehl --evalue 500 --differences 3 --maxinterval 1000 --accuracy 80 --index $TMPDIR/${genomeName}_artificial.idx --database $TMPDIR/${genomeName}_artificial.fa --nomatchfilename $TMPDIR/${bn}_unmatched.fastq --query $TMPDIR/${bn}_trimmed.fastq.gz -o $TMPDIR/${bn}.sam
gzip $TMPDIR/${bn}_unmatched.fastq

##remove all reads mapping at least once to the genome
perl $TMPDIR/removeGenomeMapper.pl $TMPDIR/${tRNAName}_pre-tRNAs.fa $TMPDIR/${bn}.sam $TMPDIR/${bn}_filtered.sam

##remove pre-tRNA reads, keep only mature tRNA reads
perl $TMPDIR/removePrecursor.pl  $TMPDIR/${tRNAName}_pre-tRNAs.bed12 $TMPDIR/${bn}_filtered.sam $TMPDIR/${bn}_trimmed.fastq.gz > $TMPDIR/${bn}_filtered.fastq
gzip $TMPDIR/${bn}_filtered.fastq

# Copy results from local node to cwd
cp $TMPDIR/{${bn}.sam,${bn}_filtered.sam,${bn}_unmatched.fastq.gz,${bn}_filtered.fastq.gz} ${workDir}/mapping/${bn}
