#!/bin/bash
echo '    #################################'
echo "############ dN/dS calculation ############"
echo "######  Credit: Florian Bénitière  ######"
echo '    #################################'



export pathData=/beegfs/XXXXX # (1) Path to change
export pathHome=/beegfs/XXXXX  # (2) Path to change
export pathdNdS=${pathData}DnDs// # (3) Path to change

# Paths
export pathFolder=${pathdNdS}
export pathScripts=${pathHome}Scripts/Projet-SplicedVariants/phylogenie_dNdS_pipeline/
export data_sequence=${pathFolder}GC_content_and_dNdS.tab
export tree_input=${pathFolder}/RAxML/concatenatAAS.aln.raxml.support

Rscript --vanilla ${pathScripts}/GC_content.R ${pathFolder}CDS/ ${pathFolder}PRANK_CDS ${tree_input} ${data_sequence}

Rscript --vanilla ${pathScripts}/partitioned_to_subset_along_GC3.R ${pathFolder}PRANK_CDS_at_least_85percent/ ${pathFolder}subset_200_ksites_GC3/ ${data_sequence} 200000


