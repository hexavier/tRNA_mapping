#!/bin/sh

i=$1
dataset=$2
cwd="/users/lserrano/xhernandez/tRNA_mapping"
ngsDir="${cwd}/Data/NGS/${dataset}"

fastqc="fastqc"					
#0.11.4

## pre- and post-quality control
bi=$(basename $i _trimmed.fastq)
gzip $i
$fastqc -q $i.gz
