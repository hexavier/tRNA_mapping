#!/bin/sh

i=$1
dataset=$2
adapter=$3
cwd="/users/lserrano/xhernandez/tRNA_mapping"
ngsDir="${cwd}/Data/NGS/${dataset}"

bbduk="bbduk.sh"				
#BBMap version 38.22
fastqc="fastqc"					
#0.11.4

# change initial heap size of Java VM
export _JAVA_OPTIONS="-Xms2000M -Xmx3500M"

## adapter and quality trimming using BBduk
bi=$(basename $i .fastq.gz)
$fastqc -q $i
$bbduk in=$i  out=${ngsDir}/${bi}/${bi}_trimmed.fastq ref=${adapter} mink=8 ktrim=r k=10 rcomp=t hdist=1 qtrim=rl trimq=25 minlength=10 maxlength=50 2> ${ngsDir}/${bi}/${bi}.bbduk.log
