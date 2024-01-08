{
  source("figure/figure_main_generator/library_path.R")
  
  ## Abstract ##
  print("Abstract")
  
  data_1 = read.delim("data/Data1_supp.tab",comment.char = "#")
  data_1 = data_1[data_1$nb_busco > .8 * 978,]
  data_1 = data_1[data_1$CoverageBuscoExon > 196 ,]
  print(paste("To test this hypothesis, we analyzed ",
              sum(data_1$nb_rnaseq)," transcriptome sequencing samples to quantify AS in ",nrow(data_1),
              " metazoan species spanning a wide range of Ne values",sep=""))
  
  
  ## Introducton ##
  print("Introducton")
  print(paste("we quantified AS rates in ",nrow(data_1)," metazoan species",sep=""))
  
  
  ## Results ##
  print("Results")
  print("Genomic and transcriptomic data collection")
  data_1 = read.delim("data/Data1_supp.tab",comment.char = "#")
  print(paste("we examined a collection of ",nrow(data_1)," species ",sep=""))
  print(paste("single-copy orthologous genes shared across metazoan (n=",
              978," genes)",sep=""))
  
  print(paste(
    "We retained for further analyses those species for which at least 80% of the BUSCO metazoan gene set could be unambiguously identified (N=",
    sum(data_1$nb_busco > .8*978)," species).",sep=""))
  
  
  
  data_1 = read.delim("data/Data1_supp.tab",comment.char = "#")
  data_1 = data_1[data_1$nb_busco > .8*978,]
  data_1 = data_1[data_1$CoverageBuscoExon > 196 ,]
  
  
  print(paste("Our final dataset thus consisted of ",nrow(data_1)," species (",nrow(data_1[data_1$clade %in% c("Mammalia", "Crocodylia","Aves", "Chondrichthyes" ),]),
              " vertebrates and ",
              nrow(data_1[!data_1$clade %in% c("Mammalia", "Crocodylia","Aves", "Chondrichthyes" ),])," insects;",sep=""))
  
  
  print(paste(", and of ",sum(data_1$nb_rnaseq)," RNA-seq samples (",
              round(sum(data_1$nb_rnaseq)/nrow(data_1))," per species on average).",sep=""))
  
  print(paste("among BUSCO genes ranges from ",round(min(data_1$analyzable_intron_busco),1)," to ",
              round(max(data_1$analyzable_intron_busco),1)," (which represents "
              ,round(min(data_1$analyzable_intron_busco/data_1$annotated_intron_busco*100),1),"% to ",round(max(data_1$analyzable_intron_busco/data_1$annotated_intron_busco*100),1),
              "% of their introns;",sep=""))
  
  
  data_10 = read.delim("data/Data10_supp.tab",comment.char = "#")
  data_10 = data_10[data_10$species %in% data_1$species,]
  print(paste("very homogenous in terms of sequencing technology (",round(sum(data_10$Platform == "ILLUMINA") / nrow(data_10)*100,0),"% of RNA-seq samples sequenced with Illumina platforms)",sep=""))
  
  ######
  print("Proxies for the effective population size (Ne)")
  print(paste("To evaluate this ratio, we aligned 922 BUSCO genes, reconstructed the phylogenetic tree of the ",nrow(data_1)," species",sep=""))
  
  ylabel="body_size"
  xlabel="longevity"
  shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)
  pval = summary(pgls(log10(ylabel)~log10(xlabel),shorebird))$coefficients[2,4]
  ylabel="dNdS_200k"
  xlabel="longevity"
  shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)
  pval=c(pval,summary(pgls((ylabel)~log10(xlabel),shorebird))$coefficients[2,4])
  print(paste("these different proxies of Ne are positively correlated with each other (p < ",
              formatC(max(pval), format = "e", digits = 1),", Fig. 1B,C)",sep=""))
  
  ######
  print("Alternative splicing rates are negatively correlated with Ne proxies")
  print(paste("In vertebrates, each BUSCO gene contains on average ",round(mean(data_1[data_1$clade %in% c("Mammalia", "Crocodylia","Aves", "Chondrichthyes" ),]$mean_intron_per_gene_busco),1),
              " introns (Supplementary Tab. 1). The intron density is more variable among insect clades, ranging from ",round(mean(data_1[data_1$clade %in% c("Diptera" ),]$mean_intron_per_gene_busco),1),
              " introns per BUSCO gene in Diptera to ",round(mean(data_1[data_1$clade %in% c("Blattodea" ),]$mean_intron_per_gene_busco),1)," in Blattodea.",sep=""))
  
  
  print(paste("As expected, most major introns have GT/AG splice sites (",round(mean(data_1$splsite_gtag_major_busco)*100,1),"% on average across species)",sep=""))
  print(paste("and only a small fraction have non-canonical boundaries (",
              round(mean(data_1$splsite_gcag_major_busco)*100,1), "% GC/AG and ",
              round(mean(data_1$splsite_atac_major_busco)*100,1), "% AT/AC).",sep=""))
  
  print(paste("The fraction of non-canonical splice sites is slightly higher among minor introns (",
              round(mean(data_1$splsite_gcag_minor_busco)*100,1), "% GC/AG and ", 
              round(mean(data_1$splsite_atac_minor_busco)*100,1), "% AT/AC).",sep=""))
  
  
  print(paste("The proportion of major introns for which AS has been detected (i.e. with Na > 0) ranges from ",
              round(min(data_1$prop_major_sv_busco)*100,1),"% to ",round(max(data_1$prop_major_sv_busco)*100,1),"%",sep=""))
  
  
  svr_busco = tapply(data_1$mean_as_busco,data_1$species,sum)
  
  print(paste("The average AS rate for BUSCO genes varies by a factor of ",round(max(svr_busco) / min(svr_busco)) ," among species, from ",
              round(svr_busco["Drosophila_grimshawi"],1),"% in Drosophila grimshawi (Diptera) to ",round(svr_busco["Megachile_rotundata"],1),"% in Megachile rotundata (Hymenoptera) (",
              round(svr_busco["Homo_sapiens"],1),"% in humans).",sep=""))
  
  data_4 = read.delim("data/Data4_supp.tab",comment.char = "#")
  data_4 = data_4[data_4$species != "Dendroctonus_ponderosae",]
  data_4 = data_4[which(data_4$organs %in% c("Brain", "Cerebellum", "Heart", "Kidney","Head", "Liver", "Ovary", "Testis")),] ## pour enlèver KidneyTestis et WholeBrain, pour lesquels on n'avait pas toutes les espèces.
  
  a = anova(lm(data_4$SVR~(data_4$species+data_4$organs)))
  output_anova = 100*a[["Sum Sq"]]/sum(a[["Sum Sq"]])
  
  print(paste( "Specifically, in an ANOVA analysis performed on the average AS rate across BUSCO gene introns, with the species and the organ of origin as explanatory variables, the species factor explained ",
               round(output_anova[1]),"% of the total variance, while the organ factor explained only ",round(output_anova[2]),"%.",sep=""))
  
  
  ########
  print("Functional vs. non-functional alternative splicing")
  
  print(paste("Of note, the subset of rare SVs represents the vast majority of the SV repertoire (from ",round(min(data_1$prop_rare_sv_protein_coding),1),
              "% to ",round(max(data_1$prop_rare_sv_protein_coding),1),"% depending on the species; Supplementary Tab. 1).",sep=""))
  
  
  #########
  print("Investigating selective pressures on minor splice sites")
  
  data_5 = read.delim(paste("data/Data5_supp.tab",sep=""),comment.char = "#")
  data_5$color_group = factor(data_5$color_group,levels = c("red","green","blue"))
  data_5$occurences = data_5$mean_polymorphism * 2 * data_5$Nb_introns_minor
  data_5$length = 2 * data_5$Nb_introns_minor
  data_5 = data_5[order(data_5$group),]
  dt = data_5[data_5$filtering == "Drosophila_melanogaster_abundant_sv",]
  
  print(paste("the SNP density at minor splice signals of abundant SVs is much lower than in corresponding controls (from ",round(max(dt$percent_difference_vs_control,na.rm = T),0),
        "% to ",round(min(dt$percent_difference_vs_control,na.rm = T),0),"%,",sep=""))
  
  print(paste("and than in minor splice signals of rare SVs (from ",round(max(dt[!is.na(dt$pvalue_vs_control),]$difference_vs_rare,na.rm = T),0),
        "% to ",round(min(dt[!is.na(dt$pvalue_vs_control),]$difference_vs_rare,na.rm = T),0),"%,",sep=""))
  
  dt = data_5[grepl("Homo_sapiens",data_5$filtering) & grepl("AG unshared minor splice site into major intron",data_5$group),]
  print(paste("(AG) located within the major intron show a weak but significant SNP deficit relative to corresponding control sites (p-value < ",
        formatC(max(dt$pvalue_vs_control,na.rm = T), format = "e", digits = 1), ")",sep=""))
  
  #########
  print("The splicing rate of rare SVs is negatively correlated with gene expression levels")
  
  print(paste("The selected subset represents ",round(min(data_1$prop_major_nt_sup100 * 100),1),"% to ",
              round(max(data_1$prop_major_nt_sup100 * 100),1),"% of introns of each species (median=",round(mean(data_1$prop_major_nt_sup100 * 100),1),"%).",sep=""))
  
  print(paste("AS rate and gene expression level in ",sum(data_1$cor_as_fpkm_low_as < 0)," out of the ",nrow(data_1 ),
              " species (significant with p < 0.05, in ", sum(data_1$cor_as_fpkm_low_as < 0 & data_1$pval_cor_as_fpkm_low_as < 0.05),"/",nrow(data_1 )," species;",sep=""))
  
  print(paste("Interestingly, when we performed this analysis on all introns (including those with abundant SVs, which are enriched in functional variants), then most species (",
              sum(data_1$cor_as_fpkm_all_as < 0 & data_1$pval_cor_as_fpkm_all_as < 0.05),"/53) still showed a negative correlation between AS rate and gene expression level",sep=""))
  
  
  
  ## Discussion ##
  print("Discussion")
  
  print(paste("publicly available RNA-seq data across a large set of ",nrow(data_1)," species",sep=""))
  
  print(paste("With this setting, on average ",round(mean(data_1$analyzable_intron_busco/data_1$annotated_intron_busco*100),1),"% of BUSCO annotated introns could be analyzed in each species.",sep=""))
  
  print(paste("We observed a ",round(max(svr_busco) / min(svr_busco)) ,"-fold variation in the average AS rate of BUSCO introns across species (Fig. 3A).",sep=""))
  
  print(paste("Rare SVs represent the vast majority of the repertoire of splicing isoforms detected in a given transcriptome (from ",
              round(min(data_1$prop_rare_sv_protein_coding),1),
              "% to ",
              round(max(data_1$prop_rare_sv_protein_coding),1),
              "% according to the species; Supplementary Tab. 1).",sep=""))
  
  
  
  
  print(paste("In agreement with previous work, we observed that AS rates tend to be high in vertebrates (average=",
              round( mean(data_1[data_1$clade %in% c("Mammalia", "Crocodylia","Aves", "Chondrichthyes" ),]$mean_as_busco),1),
              "%). and notably in primates (average=",
              round( mean(data_1[data_1$species %in% c("Homo_sapiens","Macaca_mulatta" ),]$mean_as_busco),1)
              ,"%).",sep=""))
  
  print(paste("Given that intron density varies widely across clades (from ",
              round(mean(data_1[data_1$clade %in% c("Diptera" ),]$mean_intron_per_gene_busco),1),
              " introns per gene in diptera to ",round(mean(data_1[data_1$clade %in% c("Mammalia", "Crocodylia","Aves", "Chondrichthyes" ),]$mean_intron_per_gene_busco),1),
              " introns per gene in vertebrates; Supplementary Tab. 1),",sep=""))
}
