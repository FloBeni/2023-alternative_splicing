# Calcule du dN/dS par branche et par espèce
library(stringr)
library(ggplot2)
library(ape)
library(seqinr)

# DNDS INFORMATION PER SEQUENCE
get_newick_value = function(arbre_phylo){
  nodes <- sapply(arbre_phylo$tip.label,function(x,y) which(y==x),y=arbre_phylo$tip.label)
  edge.length <- setNames(arbre_phylo$edge.length[sapply(nodes,function(x,y) which(y==x),y=arbre_phylo$edge[,2])],names(nodes))
  return(edge.length)
}


args = (commandArgs(TRUE))
path = args[1]
output_path = args[2]

list_busco <- list.dirs(path,full.names = F,recursive = F)


data_sequence = data.frame()
tree = read.tree(paste(path,"/concatenatAAS.aln.raxml.support",sep=""))
for ( busco_gene in list_busco ){ print(busco_gene)# through fasta file
  for (seq_id in tree$tip.label ){ # through sequences
    data_sequence = rbind(data_sequence,data.frame(
      species = seq_id,
      id = busco_gene
    ))
  }
}

data_sequence[,paste("alignment_nucl_length",sep="")] = NA
data_sequence[,paste("sequence_nucl_length",sep="")] = NA
data_sequence[,paste("branch_length",sep="")] = NA
data_sequence[,paste("GC_content",sep="")] = NA

for ( busco_id in unique(data_sequence$id)){print(busco_id)
  if ( grepl("dNdS_per_gene",path) ){
    fasta = read.fasta(paste(path,busco_id,"/",busco_id,".fa-prank.aln",sep=""))
  } else {
    fasta = read.fasta(paste(path,busco_id,"/concatenatCDS.aln",sep=""))
  }

  if ( grepl("dNdS_per_gene",path) ){
    tree = get_newick_value(read.tree(paste(path,busco_id,"/concatenatAAS.aln.raxml.support.filter.dnd_1",sep="")))
  } else {
    tree = get_newick_value(read.tree(paste(path,busco_id,"/concatenatAAS.aln.raxml.support.dnd_1",sep="")))
  }


  data_sequence[data_sequence$id == busco_id,paste("alignment_nucl_length",sep="")] = length(fasta[[1]])
  data_sequence[data_sequence$id == busco_id,paste("sequence_nucl_length",sep="")] =
    unlist(lapply(fasta,function(x) sum(table(x[!grepl("-",x)]))))[data_sequence[data_sequence$id == busco_id,]$species]

  data_sequence[data_sequence$id == busco_id,paste("GC_content",sep="")] = unlist(lapply(fasta,function(x) sum(table(x[!grepl("-",x) & !grepl("t",x) & !grepl("a",x)])) / sum(table(x[!grepl("-",x) ]))
  ))[data_sequence[data_sequence$id == busco_id,]$species]

  data_sequence[data_sequence$id == busco_id,paste("branch_length",sep="")] = tree[data_sequence[data_sequence$id == busco_id,]$species]

}


protocol_subst = list("dNdS_SW" = c("_W->W","_S->S","_S->W","_W->S"),
                      "dNdS"="")

for (protocol in names(protocol_subst)){print(protocol)
  for (subst in protocol_subst[protocol][[1]]){print(subst)

    if (protocol == "dNdS_SW"){
      dS_name = paste("counts_dS_X",subst,sep="")
      dN_name = paste("counts_dN_X",subst,sep="")
    } else if (protocol == "dNdS"){
      dS_name = paste("counts_dS",subst,sep="")
      dN_name = paste("counts_dN",subst,sep="")
    }


    data_sequence[,paste("num_dS",subst,sep="")] = NA
    data_sequence[,paste("den_dS",subst,sep="")] = NA
    data_sequence[,paste("num_dN",subst,sep="")] = NA
    data_sequence[,paste("den_dN",subst,sep="")] = NA

    for ( busco_id in unique(data_sequence$id)){print(busco_id)
      dSBioPP = get_newick_value(read.tree(paste(path,busco_id,"/",protocol,"/",dS_name,".dnd",sep="")))
      dSBioPP_norm = get_newick_value(read.tree(paste(path,busco_id,"/",protocol,"/",dS_name,"_norm.dnd",sep="")))
      dNBioPP = get_newick_value(read.tree(paste(path,busco_id,"/",protocol,"/",dN_name,".dnd",sep="")))
      dNBioPP_norm = get_newick_value(read.tree(paste(path,busco_id,"/",protocol,"/",dN_name,"_norm.dnd",sep="")))

      data_sequence[data_sequence$id == busco_id,paste("num_dS",subst,sep="")] = dSBioPP[data_sequence[data_sequence$id == busco_id,]$species]
      data_sequence[data_sequence$id == busco_id,paste("den_dS",subst,sep="")] = dSBioPP_norm[data_sequence[data_sequence$id == busco_id,]$species]
      data_sequence[data_sequence$id == busco_id,paste("num_dN",subst,sep="")] = dNBioPP[data_sequence[data_sequence$id == busco_id,]$species]
      data_sequence[data_sequence$id == busco_id,paste("den_dN",subst,sep="")] = dNBioPP_norm[data_sequence[data_sequence$id == busco_id,]$species]

    } 
  }
}
colnames(data_sequence$sum_den_dN)


write.table(data_sequence,output_path, row.names=F, col.names=T, sep="\t", quote=F)


