#directory paths
cwd=$(pwd)

#specify multimapper handling
# set to uniq for only uniq mapped reads
# set to phased for only phased mapped reads
# set to all for all reads
multimapperHandling="uniq"


#change here your directory paths
ngsDir="${cwd}/Data/NGS/HEK"
scriptDir="${cwd}/"
genomeDir="${cwd}/Data/Genomes/H.sapiens"
workDir="${cwd}/Analysis/HEK"
adapterFile="${cwd}/Data/NGS/adapter.fa"

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



###genome preparation
cd ${genomeDir}

### pretRNA clustering
##only identical pretRNAs were clustered
perl ${scriptDir}/clustering.pl ${tRNAName}_pre-tRNAs.fa ${tRNAName}_cluster_pre-tRNA.fa ${tRNAName}_clusterInfo_pre-tRNA.fa

##indexing pretRNA cluster
$samtools faidx ${tRNAName}_cluster_pre-tRNA.fa
$segemehl -x ${tRNAName}_cluster_pre-tRNA.idx -d ${tRNAName}_cluster_pre-tRNA.fa
$picard CreateSequenceDictionary R=${tRNAName}_cluster_pre-tRNA.fa O=${tRNAName}_cluster_pre-tRNA.dict



###pre-mapping against artificial genome
cd ${workDir}/mapping

for n in $(ls ${ngsDir}/*trimmed.fastq.gz)
do
    ##keep only pre-tRNA reads
    perl ${scriptDir}/keepPrecursor.pl  ${genomeDir}/${tRNAName}_pre-tRNAs.bed12 ${bn}_filtered.sam $n > ${bn}_pretRNA.fastq
    gzip ${bn}_filtered.fastq
done



###post-processing
mkdir -p ${workDir}/postprocessing_precursors
cd ${workDir}/postprocessing_precursors

for n in $(ls ${workDir}/mapping/*pretRNA.fastq.gz)
do
    bn=$(basename $n _pretRNA.fastq.gz)
    mkdir ${workDir}/postprocessing_precursors/${bn}
    cd ${workDir}/postprocessing_precursors/${bn}
    qsub -V -cwd -l virtual_free=80G -q long-sl7 ${scriptDir}/script_postprocessing_pretRNA.sh $n
done

