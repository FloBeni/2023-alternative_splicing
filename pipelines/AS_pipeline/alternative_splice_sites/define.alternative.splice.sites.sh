#!/bin/bash

export sp=$1
export sample=$2

##############################################################

export pathResults=/beegfs/data/fbenitiere/Projet-SplicedVariants/Analyses/Species_${sp}
export pathAnnot=/beegfs/data/fbenitiere/Projet-SplicedVariants/Annotations/Species_${sp}
export pathScripts=/beegfs/home/fbenitiere/Scripts/Projet-SplicedVariants/alternative_splice_sites

##############################################################

perl ${pathScripts}/define.alternative.splice.sites.pl --pathExonBlocks=${pathAnnot}/exon_union.txt --pathIntronLibrary=${pathResults}/IntronLibrary.txt --pathIntronQuantification=${pathResults}/${sample}/junctions.txt --pathOutput=alternative_splice_sites.txt

##############################################################
