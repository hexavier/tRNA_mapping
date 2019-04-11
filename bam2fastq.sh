#!/bin/sh

n=$1
bn=$(basename $n .bam)
dn=$(dirname $n)
cwd="/users/lserrano/xhernandez/tRNA_mapping"

samtools="samtools"				
#1.3.1

# Transform bam to fastq
$samtools fastq $n > ${dn}/${bn}.fastq
gzip ${dn}/${bn}.fastq

