#directory paths
cwd=$(pwd)

#specify multimapper handling
# set to uniq for only uniq mapped reads
# set to phased for only phased mapped reads
# set to all for all reads
multimapperHandling="uniq"


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

###pre-mapping against artificial genome
mkdir -p ${workDir}/mapping
cd ${workDir}/mapping

for n in $(ls ${ngsDir}/*/*_read1_trimmed.fastq.gz)
do
    bn=$(basename $n _read1_trimmed.fastq.gz)
    mkdir ${workDir}/mapping/${bn}
    cd ${workDir}/mapping/${bn}
    qsub -V -cwd -l h_rt=48:00:00 -q long-sl7 ${scriptDir}/script_mapping_PE_rRNA.sh $n ${dataset}
done


