# Generate Data 5
options(stringsAsFactors = F, scipen = 999)
pathData="/home/XXXXX/data/Projet-SplicedVariants/"
# pathData="/beegfs/data/XXXXX/Projet-SplicedVariants/"

freq=0

data_5 = data.frame()


### DROSO FREQUENT
species = "Drosophila_melanogaster"
polymorphisme = read.delim(file=paste(pathData , "Annotations/" , species , "/polymorphism/by_minor_intron.tab",sep=""))

polymorphisme = polymorphisme[polymorphisme$which_shared_site != "both",]
polymorphisme = polymorphisme[polymorphisme$criptic_intron == "False",]
polymorphisme = polymorphisme[polymorphisme$into_cds == "True",]
polymorphisme = polymorphisme[polymorphisme$mira > 0.05 ,]
{
  # Pour Dmel on ne souhaite pas dissocier les CpG des non CpG car il n'y a pas d'imapct du CpG, ainsi si pour le CpG controle on n'a rien en CpG alors on prend les non CpG et on analysera le CpG controle qui reviendra a tout analyser
  polymorphisme[is.na(polymorphisme$GT_CpG.splice3.20_60bp),"GT_CpG.splice3.20_60bp"] = polymorphisme[is.na(polymorphisme$GT_CpG.splice3.20_60bp),"GT_noCpG.splice3.20_60bp"]
  polymorphisme[is.na(polymorphisme$GT_CpG.splice3.20_60bp_intron),"GT_CpG.splice3.20_60bp_intron"] = polymorphisme[is.na(polymorphisme$GT_CpG.splice3.20_60bp_intron),"GT_noCpG.splice3.20_60bp_intron"]
  
  ########### AG analyses
  
  polymorphisme_data = data.frame()
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$AG.splice5.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG control 20-60bp into major intron 5' side",
    mean_polymorphism = mean(polymorphisme_select$AG.splice5.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$AG.splice5.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG control 20-60bp into major exon 5' side",
    mean_polymorphism = mean(polymorphisme_select$AG.splice5.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$specific_splsite_into_shared_major == "True",]
  
  control_dt = polymorphisme_data[polymorphisme_data$group == "AG control 20-60bp into major intron 5' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$specific_splsite_into_shared_major == "False",]
  
  control_dt = polymorphisme_data[polymorphisme_data$group == "AG control 20-60bp into major exon 5' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  # GT
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$GT_CpG.splice3.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT control 20-60bp into major intron 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_CpG.splice3.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$GT_CpG.splice3.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT control 20-60bp into major exon 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_CpG.splice3.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3"& polymorphisme$specific_splsite_into_shared_major == "True",]
  
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT control 20-60bp into major intron 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3"  & polymorphisme$specific_splsite_into_shared_major == "False",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT control 20-60bp into major exon 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  
  rownames(polymorphisme_data) = polymorphisme_data$group
  df = polymorphisme_data[c("AG unshared major splice site","AG unshared minor splice site into major intron","AG unshared minor splice site into major exon",
                            "AG control 20-60bp into major intron 5' side","AG control 20-60bp into major exon 5' side","GT unshared major splice site",
                            "GT unshared minor splice site into major intron","GT unshared minor splice site into major exon","GT control 20-60bp into major intron 3' side",
                            "GT control 20-60bp into major exon 3' side"),]
  df$color_group = c("green","red","red","blue","blue","green","red","red","blue","blue")
}

df$filtering = paste(species,"abundant_sv",sep="_")
data_5 = rbind(data_5,df)






#### HOMO SAPIENS FREQUENT

species = "Homo_sapiens"
{
  list_files = list.files(paste(pathData,"Annotations/Homo_sapiens/SNP_chr_v2/polymorphism_",freq,"_v2",sep=""))
  polymorphisme = data.frame()
  for (file in list_files){print(file)
    polymorphisme = rbind(polymorphisme,read.delim(file=paste(pathData,"Annotations/Homo_sapiens/SNP_chr_v2/polymorphism_",freq,"_v2/",
                                                              file,sep="")))
  }
  all_polymorphisme = polymorphisme
  
  polymorphisme = polymorphisme[polymorphisme$which_shared_site != "both",]
  polymorphisme = polymorphisme[polymorphisme$criptic_intron == "False",]
  polymorphisme = polymorphisme[polymorphisme$into_cds == "True" ,]
  polymorphisme$mira = polymorphisme$n1 / (polymorphisme$n1_major + polymorphisme$n2_major)
  polymorphisme = polymorphisme[polymorphisme$mira > 0.05 ,]  
  
}
{
  ########### AG analyses
  
  polymorphisme_data = data.frame()
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$AG.splice5.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG control 20-60bp into major intron 5' side",
    mean_polymorphism = mean(polymorphisme_select$AG.splice5.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$AG.splice5.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG control 20-60bp into major exon 5' side",
    mean_polymorphism = mean(polymorphisme_select$AG.splice5.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$specific_splsite_into_shared_major == "True",]
  
  control_dt = polymorphisme_data[polymorphisme_data$group == "AG control 20-60bp into major intron 5' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$specific_splsite_into_shared_major == "False",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "AG control 20-60bp into major exon 5' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  # GT
  
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$GT_noCpG.splice3.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG control 20-60bp into major intron 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_noCpG.splice3.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$GT_noCpG.splice3.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG control 20-60bp into major exon 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_noCpG.splice3.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$major_splice5_CpG == "No",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$major_splice5_CpG == "No",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "No",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "No" & polymorphisme$specific_splsite_into_shared_major == "True",]
  
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT no CpG control 20-60bp into major intron 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "No" & polymorphisme$specific_splsite_into_shared_major == "False",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT no CpG control 20-60bp into major exon 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$GT_CpG.splice3.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG control 20-60bp into major intron 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_CpG.splice3.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$GT_CpG.splice3.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG control 20-60bp into major exon 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_CpG.splice3.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$major_splice5_CpG == "Yes",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$major_splice5_CpG == "Yes",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "Yes",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "Yes" & polymorphisme$specific_splsite_into_shared_major == "True",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT CpG control 20-60bp into major intron 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "Yes" & polymorphisme$specific_splsite_into_shared_major == "False",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT CpG control 20-60bp into major exon 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  
  rownames(polymorphisme_data) = polymorphisme_data$group
  df = polymorphisme_data[c("AG unshared major splice site","AG unshared minor splice site into major intron","AG unshared minor splice site into major exon",
                              "AG control 20-60bp into major intron 5' side","AG control 20-60bp into major exon 5' side","GT no CpG unshared major splice site",
                           "GT no CpG unshared minor splice site into major intron" ,"GT no CpG unshared minor splice site into major exon",
                            "GT no CpG control 20-60bp into major intron 3' side","GT no CpG control 20-60bp into major exon 3' side"),]
  
  df$color_group = c("green","red","red","blue","blue","green","red","red","blue","blue")
}

df$filtering = paste(species,"abundant_sv",sep="_")
data_5 = rbind(data_5,df)




######## RARE VARIANTS

### DROSO FREQUENT


species = "Drosophila_melanogaster"
polymorphisme = read.delim(file=paste(pathData,"Annotations/",species,"/polymorphism/by_minor_intron.tab",sep=""))

polymorphisme = polymorphisme[polymorphisme$which_shared_site != "both",]
polymorphisme = polymorphisme[polymorphisme$criptic_intron == "False",]
polymorphisme = polymorphisme[polymorphisme$into_cds == "True",]
polymorphisme$mira = polymorphisme$n1 / (polymorphisme$n1_major + polymorphisme$n2_major)
polymorphisme = polymorphisme[polymorphisme$mira <= 0.05 ,]

{
  # Pour Dmel on ne souhaite pas dissocier les CpG des non CpG car il n'y a pas d'imapct du CpG, ainsi si pour le CpG controle on n'a rien en CpG alors on prend les non CpG et on analysera le CpG controle qui reviendra a tout analyser
  polymorphisme[is.na(polymorphisme$GT_CpG.splice3.20_60bp),"GT_CpG.splice3.20_60bp"] = polymorphisme[is.na(polymorphisme$GT_CpG.splice3.20_60bp),"GT_noCpG.splice3.20_60bp"]
  polymorphisme[is.na(polymorphisme$GT_CpG.splice3.20_60bp_intron),"GT_CpG.splice3.20_60bp_intron"] = polymorphisme[is.na(polymorphisme$GT_CpG.splice3.20_60bp_intron),"GT_noCpG.splice3.20_60bp_intron"]
  
  ########### AG analyses
  
  polymorphisme_data = data.frame()
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$AG.splice5.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG control 20-60bp into major intron 5' side",
    mean_polymorphism = mean(polymorphisme_select$AG.splice5.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$AG.splice5.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG control 20-60bp into major exon 5' side",
    mean_polymorphism = mean(polymorphisme_select$AG.splice5.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$specific_splsite_into_shared_major == "True",]
  
  control_dt = polymorphisme_data[polymorphisme_data$group == "AG control 20-60bp into major intron 5' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$specific_splsite_into_shared_major == "False",]
  
  control_dt = polymorphisme_data[polymorphisme_data$group == "AG control 20-60bp into major exon 5' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  # GT
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$GT_CpG.splice3.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT control 20-60bp into major intron 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_CpG.splice3.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$GT_CpG.splice3.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT control 20-60bp into major exon 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_CpG.splice3.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3"& polymorphisme$specific_splsite_into_shared_major == "True",]
  
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT control 20-60bp into major intron 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3"  & polymorphisme$specific_splsite_into_shared_major == "False",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT control 20-60bp into major exon 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  
  rownames(polymorphisme_data) = polymorphisme_data$group
  df = polymorphisme_data[c("AG unshared major splice site","AG unshared minor splice site into major intron","AG unshared minor splice site into major exon",
                            "AG control 20-60bp into major intron 5' side","AG control 20-60bp into major exon 5' side","GT unshared major splice site",
                            "GT unshared minor splice site into major intron","GT unshared minor splice site into major exon","GT control 20-60bp into major intron 3' side",
                            "GT control 20-60bp into major exon 3' side"),]
  df$color_group = c("green","red","red","blue","blue","green","red","red","blue","blue")
}
df$filtering = paste(species,"rare_sv",sep="_")
data_5 = rbind(data_5,df)


### HOMO SAPIENS FREQUENT

species = "Homo_sapiens"

{
  
  polymorphisme = all_polymorphisme
  
  polymorphisme = polymorphisme[polymorphisme$which_shared_site != "both",]
  polymorphisme = polymorphisme[polymorphisme$criptic_intron == "False",]
  polymorphisme = polymorphisme[polymorphisme$into_cds == "True" ,]
  polymorphisme$mira = polymorphisme$n1 / (polymorphisme$n1_major + polymorphisme$n2_major)
  polymorphisme = polymorphisme[polymorphisme$mira  <= 0.05  ,]  
  
}
{
  ########### AG analyses
  
  polymorphisme_data = data.frame()
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$AG.splice5.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG control 20-60bp into major intron 5' side",
    mean_polymorphism = mean(polymorphisme_select$AG.splice5.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$AG.splice5.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG control 20-60bp into major exon 5' side",
    mean_polymorphism = mean(polymorphisme_select$AG.splice5.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$specific_splsite_into_shared_major == "True",]
  
  control_dt = polymorphisme_data[polymorphisme_data$group == "AG control 20-60bp into major intron 5' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$specific_splsite_into_shared_major == "False",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "AG control 20-60bp into major exon 5' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  # GT
  
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$GT_noCpG.splice3.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG control 20-60bp into major intron 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_noCpG.splice3.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$GT_noCpG.splice3.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG control 20-60bp into major exon 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_noCpG.splice3.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$major_splice5_CpG == "No",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$major_splice5_CpG == "No",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "No",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "No" & polymorphisme$specific_splsite_into_shared_major == "True",]
  
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT no CpG control 20-60bp into major intron 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "No" & polymorphisme$specific_splsite_into_shared_major == "False",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT no CpG control 20-60bp into major exon 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$GT_CpG.splice3.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG control 20-60bp into major intron 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_CpG.splice3.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$GT_CpG.splice3.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG control 20-60bp into major exon 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_CpG.splice3.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$major_splice5_CpG == "Yes",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$major_splice5_CpG == "Yes",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "Yes",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "Yes" & polymorphisme$specific_splsite_into_shared_major == "True",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT CpG control 20-60bp into major intron 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "Yes" & polymorphisme$specific_splsite_into_shared_major == "False",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT CpG control 20-60bp into major exon 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  
  rownames(polymorphisme_data) = polymorphisme_data$group
  df = polymorphisme_data[c("AG unshared major splice site","AG unshared minor splice site into major intron","AG unshared minor splice site into major exon",
                            "AG control 20-60bp into major intron 5' side","AG control 20-60bp into major exon 5' side","GT no CpG unshared major splice site",
                            "GT no CpG unshared minor splice site into major intron" ,"GT no CpG unshared minor splice site into major exon",
                            "GT no CpG control 20-60bp into major intron 3' side","GT no CpG control 20-60bp into major exon 3' side"),]
  
  df$color_group = c("green","red","red","blue","blue","green","red","red","blue","blue")
}

df$filtering = paste(species,"rare_sv",sep="_")
data_5 = rbind(data_5,df)


### CpG Homo sapiens


species = "Homo_sapiens"
{
  
  polymorphisme = all_polymorphisme
  
  
  polymorphisme = polymorphisme[polymorphisme$which_shared_site != "both",]
  polymorphisme = polymorphisme[polymorphisme$criptic_intron == "False",]
  polymorphisme = polymorphisme[polymorphisme$into_cds == "True" ,]
  polymorphisme = polymorphisme[polymorphisme$mira > 0.05 ,]  
  
}
{
  ########### AG analyses
  
  polymorphisme_data = data.frame()
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$AG.splice5.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG control 20-60bp into major intron 5' side",
    mean_polymorphism = mean(polymorphisme_select$AG.splice5.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$AG.splice5.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG control 20-60bp into major exon 5' side",
    mean_polymorphism = mean(polymorphisme_select$AG.splice5.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$specific_splsite_into_shared_major == "True",]
  
  control_dt = polymorphisme_data[polymorphisme_data$group == "AG control 20-60bp into major intron 5' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$specific_splsite_into_shared_major == "False",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "AG control 20-60bp into major exon 5' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  # GT
  
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$GT_noCpG.splice3.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG control 20-60bp into major intron 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_noCpG.splice3.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$GT_noCpG.splice3.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG control 20-60bp into major exon 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_noCpG.splice3.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$major_splice5_CpG == "No",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$major_splice5_CpG == "No",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "No",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "No" & polymorphisme$specific_splsite_into_shared_major == "True",]
  
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT no CpG control 20-60bp into major intron 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "No" & polymorphisme$specific_splsite_into_shared_major == "False",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT no CpG control 20-60bp into major exon 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$GT_CpG.splice3.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG control 20-60bp into major intron 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_CpG.splice3.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$GT_CpG.splice3.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG control 20-60bp into major exon 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_CpG.splice3.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$major_splice5_CpG == "Yes",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$major_splice5_CpG == "Yes",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "Yes",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "Yes" & polymorphisme$specific_splsite_into_shared_major == "True",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT CpG control 20-60bp into major intron 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "Yes" & polymorphisme$specific_splsite_into_shared_major == "False",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT CpG control 20-60bp into major exon 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  
  rownames(polymorphisme_data) = polymorphisme_data$group
  df = polymorphisme_data[c("GT CpG unshared major splice site","GT CpG unshared minor splice site into major intron" ,"GT CpG unshared minor splice site into major exon",
                           "GT CpG control 20-60bp into major intron 3' side","GT CpG control 20-60bp into major exon 3' side"    ),]
  
  df$color_group = c("green","red","red","blue","blue")
}

df$filtering = paste(species,"CpG_abundant_sv",sep="_")
data_5 = rbind(data_5,df)
### CpG Homo sapiens


species = "Homo_sapiens"
{
  
  polymorphisme = all_polymorphisme
  
  
  polymorphisme = polymorphisme[polymorphisme$which_shared_site != "both",]
  polymorphisme = polymorphisme[polymorphisme$criptic_intron == "False",]
  polymorphisme = polymorphisme[polymorphisme$into_cds == "True" ,]
  polymorphisme = polymorphisme[polymorphisme$mira <= 0.05 ,]  
  
}
{
  ########### AG analyses
  
  polymorphisme_data = data.frame()
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$AG.splice5.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG control 20-60bp into major intron 5' side",
    mean_polymorphism = mean(polymorphisme_select$AG.splice5.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$AG.splice5.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG control 20-60bp into major exon 5' side",
    mean_polymorphism = mean(polymorphisme_select$AG.splice5.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$AG.splice5.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$specific_splsite_into_shared_major == "True",]
  
  control_dt = polymorphisme_data[polymorphisme_data$group == "AG control 20-60bp into major intron 5' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$specific_splsite_into_shared_major == "False",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "AG control 20-60bp into major exon 5' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "AG unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  # GT
  
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$GT_noCpG.splice3.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG control 20-60bp into major intron 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_noCpG.splice3.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$GT_noCpG.splice3.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG control 20-60bp into major exon 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_noCpG.splice3.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_noCpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$major_splice5_CpG == "No",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$major_splice5_CpG == "No",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "No",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "No" & polymorphisme$specific_splsite_into_shared_major == "True",]
  
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT no CpG control 20-60bp into major intron 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "No" & polymorphisme$specific_splsite_into_shared_major == "False",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT no CpG control 20-60bp into major exon 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT no CpG unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  
  
  polymorphisme_select =  polymorphisme[!is.na(polymorphisme$GT_CpG.splice3.20_60bp_intron),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG control 20-60bp into major intron 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_CpG.splice3.20_60bp_intron)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp_intron),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[!is.na(polymorphisme$GT_CpG.splice3.20_60bp),]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG control 20-60bp into major exon 3' side",
    mean_polymorphism = mean(polymorphisme_select$GT_CpG.splice3.20_60bp)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$GT_CpG.splice3.20_60bp),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice5" & polymorphisme$major_splice5_CpG == "Yes",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG shared splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_shared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_shared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$major_splice5_CpG == "Yes",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared major splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_major_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_major_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "Yes",]
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared minor splice site",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = NA,
    percent_difference_vs_control = NA,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "Yes" & polymorphisme$specific_splsite_into_shared_major == "True",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT CpG control 20-60bp into major intron 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared minor splice site into major intron",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  polymorphisme_select = polymorphisme[polymorphisme$which_shared_site == "splice3" & polymorphisme$splice5_CpG == "Yes" & polymorphisme$specific_splsite_into_shared_major == "False",]
  control_dt = polymorphisme_data[polymorphisme_data$group == "GT CpG control 20-60bp into major exon 3' side",]
  prop_test = prop.test(x = sum(polymorphisme_select$snp_notshared_site), n = nrow(polymorphisme_select)*2 , p = control_dt$error_bar, alternative = "less")
  polymorphisme_data = rbind(polymorphisme_data,data.frame(
    group = "GT CpG unshared minor splice site into major exon",
    mean_polymorphism = mean(polymorphisme_select$snp_notshared_site)/2,
    Nb_introns_minor = nrow(polymorphisme_select),
    svr_mean = mean(polymorphisme_select$mira)/2,
    
    pvalue_vs_control = prop_test$p.value,
    percent_difference_vs_control = (mean(polymorphisme_select$snp_notshared_site)/2 - control_dt$mean_polymorphism) / (control_dt$mean_polymorphism) * 100 ,
    error_bar = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[1],
    error_bar_2 = prop.test(sum(polymorphisme_select$snp_notshared_site),n=2*nrow(polymorphisme_select))$conf.int[2]
  ))
  
  
  
  
  rownames(polymorphisme_data) = polymorphisme_data$group
  df = polymorphisme_data[c("GT CpG unshared major splice site","GT CpG unshared minor splice site into major intron" ,"GT CpG unshared minor splice site into major exon",
                            "GT CpG control 20-60bp into major intron 3' side","GT CpG control 20-60bp into major exon 3' side"    ),]
  
  df$color_group = c("green","red","red","blue","blue")
}

df$filtering = paste(species,"CpG_rare_sv",sep="_")
data_5 = rbind(data_5,df)

data_5[!is.na(data_5$pvalue_vs_control) & data_5$pvalue_vs_control < 0.05,"significance"] = "*"
data_5[!is.na(data_5$pvalue_vs_control) & data_5$pvalue_vs_control < 0.005,"significance"] = "**"
data_5[!is.na(data_5$pvalue_vs_control) & data_5$pvalue_vs_control < 0.0005,"significance"] = "***"

data_5[grepl("_abundant",data_5$filtering),"difference_vs_rare"] = (data_5[grepl("_abundant",data_5$filtering),]$mean_polymorphism - data_5[grepl("_rare",data_5$filtering),]$mean_polymorphism) /
  data_5[grepl("_rare",data_5$filtering),]$mean_polymorphism * 100

write.table(data_5,paste("data/Data5_supp.tab",sep=""), row.names=F, col.names=T, sep="\t", quote=F)

