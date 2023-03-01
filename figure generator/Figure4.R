source("figure generator/library_path.R")

#### Figure 4

imgA = load.image(paste(path_require,"mira.png",sep=""))

{
  pdf(file=paste(path_pannel,"Figure4.pdf",sep=""), width=6.75, height=2)
  
  m=matrix(rep(1,10*7), nrow=7)
  m
  layout(m)
  
  par(mar=c(0, 3, 0, 1))
  plot(imgA, axes = F)
  mtext("A",at=80,adj=0, side=2, line=1, font=2, cex=1.7,las=2)
  dev.off()
}
