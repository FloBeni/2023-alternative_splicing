source("figure/figure_supp generator/library_path.R")

#### Figure 10

imgA = load.image(paste(path_require,"pipeline.png",sep=""))

{
  pdf(file=paste(path_pannel,"Figure10_supp.pdf",sep=""), width=7, height=4.5)
  
  m=matrix(rep(1,10*17), nrow=17)
  layout(m)
  
  par(mar=c(0, 0, 1, 0))
  plot(imgA, axes = F)
  mtext("A", side=2,at=5,adj=-2, line=1, font=2, cex=1.2,las=2)
  dev.off()
}
