cwd="/users/lserrano/xhernandez/tRNA_mapping"
scriptDir="${cwd}"
genomeDir="${cwd}/Data/Genomes/H.sapiens"
tRNAName="hg38.tRNAscan"
dataset=$1
ngsDir="${cwd}/Data/NGS/${dataset}"

python ${scriptDir}/reads_stats.py ${genomeDir}/${tRNAName}_clusterInfo.fa ${genomeDir}/tRNAs.txt results_nomod.csv ${ngsDir}/metadata.txt
