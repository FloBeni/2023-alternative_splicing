source("figure/figure_supp_generator/library_path.R")

############## Supplementary Figure 11

imgA = load.image(paste(path_require,"pipeline.png",sep=""))

{
  pdf(file=paste(path_pannel,"Figure11_supp.pdf",sep=""), width=7, height=4.5)
  
  m=matrix(rep(1,10*17), nrow=17)
  layout(m)
  
  par(mar=c(0, 0, 1, 0))
  plot(imgA, axes = F)
  dev.off()
}
