#!/bin/bash
# Download from NCBI genomes, annotations and protein sequences

# Paths
export file_path=$1
export genome_path=$2


nb_chevron=$(grep -c "^>" ${genome_path})
nb_by=$(echo "$((${nb_chevron}/20))")
nb_by=${nb_chevron}
echo ${nb_by}

awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%3==0){file=sprintf("'${file_path}'split_%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < ${genome_path}

list_file=$(ls ${file_path})
echo ${list_file}
for file in ${list_file}
do
echo ${file_path}/${file}
tRNAscan-SE -E -o ${file_path}/${file}_tRNA.tab -S off ${file_path}/${file}
#conda run -n R_env tRNAscan-SE -E -o ${file_path}/${file}_tRNA.tab -S off ${file_path}/${file}

find ${file_path} -type f -name *_tRNA.tab -exec cat {} + > $3
done
