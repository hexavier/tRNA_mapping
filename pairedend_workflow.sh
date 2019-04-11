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
genomeDir="${cwd}/Data/Genomes/H.sapiens"
workDir="${cwd}/Analysis/${dataset}"
adapterFile="${cwd}/Data/NGS/adapters_Jochen.fa"

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




## pre- and post-quality control, adapter and quality trimming using BBduk
cd ${ngsDir}

## adapter and quality trimming using BBduk
for i in $(ls ${ngsDir}/*/*_read1.fastq.gz)
do
  bi=$(basename $i _read1.fastq.gz)
  cd ${ngsDir}/${bi}
  qsub -V -cwd ${scriptDir}/script_trimming_PE.sh $i ${dataset} ${adapterFile}
done


## pre- and post-quality control
for i in $(ls ${ngsDir}/*/*trimmed.fastq)
do
  s=${i##*/}
  bi=${s%_read*_trimmed.fastq}
  cd ${ngsDir}/${bi}
  qsub -V -cwd ${scriptDir}/script_quality.sh $i ${dataset}
done



###genome preparation
#cd ${genomeDir}

##genome indexing using samtools
#$samtools faidx ${genomeName}.fa

## scan for tRNA nuclear
#$tRNAscanSE -q -b ${tRNAName}.nuc.bed12  ${genomeName}.fa

## scan for mitochondrial tRNA, consider: tRNAscanSE finds only 21 mt tRNA
#cat ${genomeName}.fa | perl -lane 'BEGIN{$c=0;}if(m/^>NC_012920.1 Homo sapiens mitochondrion, complete genome$/){$c=1}elsif(m/^>/){$c=0;}#print if $c' > ${genomeName}.chrM.fa
#$tRNAscanSE -q -O -b ${tRNAName}.chrM.bed12 ${genomeName}.chrM.fa

#grep -v chrM ${tRNAName}.nuc.bed12 > ${tRNAName}.nuc_mod.bed12
#cat ${tRNAName}.nuc_mod.bed12 ${tRNAName}.chrM.bed12 > ${tRNAName}.bed12

##mask found tRNAs genomic
#$bedtools maskfasta -fi ${genomeName}.fa -fo ${genomeName}.masked.fa -mc N -bed ${tRNAName}.bed12


###create pre-tRNA library
##add 50 nt 5' and 3' flanking regions
#perl ${scriptDir}/modBed12.pl ${tRNAName}.bed12 ${tRNAName}_pre-tRNAs.bed12

##remove introns, make fasta from bed12
#$bedtools getfasta -name -split -s -fi ${genomeName}.fa -bed ${tRNAName}_pre-tRNAs.bed12 -fo ${tRNAName}_pre-tRNAs.fa

##add pre-tRNAs as extra chromosoms to the genome (get the artificial genome)
#cat ${genomeName}.masked.fa ${tRNAName}_pre-tRNAs.fa > ${genomeName}_artificial.fa

##indexing artificial genome
#$samtools faidx ${genomeName}_artificial.fa
#$segemehl -x ${genomeName}_artificial.idx -d ${genomeName}_artificial.fa


###create mature tRNA library
##remove introns, make fasta from bed12
#$bedtools getfasta -name -split -s -fi ${genomeName}.fa -bed ${tRNAName}.bed12 -fo ${tRNAName}.fa

##add CCA tail to tRNA chromosomes
##remove pseudogenes
#perl ${scriptDir}/addCCA.pl ${tRNAName}.fa ${tRNAName}_mature.fa



###mature tRNA clustering
##only identical tRNAs were clustered
#perl ${scriptDir}/clustering.pl ${tRNAName}_mature.fa ${tRNAName}_cluster.fa ${tRNAName}_clusterInfo.fa

##indexing tRNA cluster
#$samtools faidx ${tRNAName}_cluster.fa
#$segemehl -x ${tRNAName}_cluster.idx -d ${tRNAName}_cluster.fa
#$picard CreateSequenceDictionary R=${tRNAName}_cluster.fa O=${tRNAName}_cluster.dict



###pre-mapping against artificial genome
mkdir -p ${workDir}/mapping
cd ${workDir}/mapping

for n in $(ls ${ngsDir}/*/*_read1_trimmed.fastq.gz)
do
    bn=$(basename $n _read1_trimmed.fastq.gz)
    mkdir ${workDir}/mapping/${bn}
    cd ${workDir}/mapping/${bn}
    qsub -V -cwd -l virtual_free=60G,h_rt=48:00:00 -q long-sl7 ${scriptDir}/script_mapping_artificial_PE.sh $n ${dataset}
done


###post-processing
mkdir -p ${workDir}/postprocessing
cd ${workDir}/postprocessing

for n in $(ls ${workDir}/mapping/*/*filtered.fastq.gz)
do
    bn=$(basename $n _filtered.fastq.gz)
    mkdir ${workDir}/postprocessing/${bn}
    cd ${workDir}/postprocessing/${bn}

    qsub -V -cwd ${scriptDir}/script_postprocessing.sh $n
done

cd ${workDir}/postprocessing
python ${scriptDir}/extract_data.py ${ngsDir}/metadata.txt
#python ${scriptDir}/reads_stats.py ${genomeDir}/${tRNAName}_clusterInfo.fa ${genomeDir}/tRNAs.txt results_nomod.csv ${ngsDir}/metadata.txt
qsub -V -cwd -l h_rt=96:00:00 -q long-sl7 /users/lserrano/xhernandez/tRNA_mapping/script_read_stats.sh ${dataset}
