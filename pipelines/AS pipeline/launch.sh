#!/usr/bin/env bash
echo '     ################################'
echo "############ RNA_seq Analysis ############"
echo "######  Credit: Florian Bénitière   ######"
echo '     ################################'

echo
echo -n "Enter the PROJECT FILE and press [ENTER]: "
read projet
echo
pathData=/beegfs/data/XXXXXX/Projet-${projet}/
echo The directory path is: ${pathData}


mkdir -p ${pathData}Alignements/
mkdir -p ${pathData}Analyses/
mkdir -p ${pathData}Analyses-RNAseq/
mkdir -p ${pathData}Annotations/
mkdir -p ${pathData}Genomes/
mkdir -p ${pathData}RNAseq_table/
mkdir -p ${pathData}Busco/

mkdir -p ${pathData}Output
mkdir -p ${pathData}Output/Annotations_treatment/
mkdir -p ${pathData}Output/Splicing_analysis/
mkdir -p ${pathData}Output/Retention_analysis/
mkdir -p ${pathData}Output/Expression/
mkdir -p ${pathData}Output/Compilation/
mkdir -p ${pathData}Output/Sequencing_depth/
mkdir -p ${pathData}Output/Hisat/Indexer/
mkdir -p ${pathData}Output/Hisat/Alignement/


#echo -n "Enter the SNAKEMAKE DIRECTORY and press [ENTER]: "
#read dir_version
dir_version=$(basename "$PWD")
echo
pathScriptsPl=/beegfs/home/XXXXXX/Scripts/Projet-SplicedVariants/SVR_estimation/${dir_version}/config
echo The species config file is: ${pathScriptsPl}

echo
echo -n "Enter the EXCEL NAME and press [ENTER]: "
read excel_file
echo
export pathHome=/beegfs/home/XXXXXX/  # (2) Path to change
export Excel_dataset=${pathData}Fichiers-data/${excel_file}.xls  # (2) Path to change

echo The excel file is: ${Excel_dataset}




echo
read -p "Do you want to DETECT RNAseq data (y/n) ? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    Rscript ${pathHome}/Scripts/Projet-SplicedVariants/SVR_estimation/data_source_generator.R ${pathData}RNAseq_table/ ${pathHome} ${Excel_dataset} # Generate RNAseq_list files for each species
fi


echo
read -p "Do you want to download NCBI data (y/n) ? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    ncbi_data="Yes"
else
    ncbi_data="No"
fi


echo
echo -n "Enter the BUSCO DATASET and press [ENTER]: "
read BUSCO_dataset


echo
echo -n "Enter the number of JOB and press [ENTER]: "
read job_number

echo
read -p "Do you want to LAUNCH (y/n) ? " -n 1 -r
echo    # (optional) move to a new line
echo 'The command is: snakemake --cluster "sbatch -J {params.name} -p {params.partition} -N 1 --ntasks={params.ntasks} --mem={params.mem} -t {params.time} -o {params.out} -e {params.err}" --config BUSCO_dataset='${BUSCO_dataset}' ncbi_data='${ncbi_data} 'excel_file='${excel_file}' projet='${projet}' dir_version='${dir_version}' -j '${job_number}' --rerun-incomplete -n'

if [[ $REPLY =~ ^[Yy]$ ]]
then
    snakemake --cluster "sbatch -J {params.name} -p {params.partition} -N 1 --ntasks={params.ntasks} --mem={params.mem} -t {params.time} -o {params.out} -e {params.err}" --config BUSCO_dataset=${BUSCO_dataset} ncbi_data=${ncbi_data} excel_file=${excel_file} projet=${projet} dir_version=${dir_version} -j ${job_number} --rerun-incomplete
else
    snakemake --cluster "sbatch -J {params.name} -p {params.partition} -N 1 --ntasks={params.ntasks} --mem={params.mem} -t {params.time} -o {params.out} -e {params.err}" --config BUSCO_dataset=${BUSCO_dataset} ncbi_data=${ncbi_data} excel_file=${excel_file} projet=${projet} dir_version=${dir_version} -j ${job_number} -n --rerun-incomplete
    echo 'The command is: snakemake --cluster "sbatch -J {params.name} -p {params.partition} -N 1 --ntasks={params.ntasks} --mem={params.mem} -t {params.time} -o {params.out} -e {params.err}" --config BUSCO_dataset='${BUSCO_dataset}' ncbi_data='${ncbi_data} 'excel_file='${excel_file}' projet='${projet}' dir_version='${dir_version}' -j '${job_number}' --rerun-incomplete -n'
fi

