source("figure/figure_refPCI_generator/library_path.R")


############## PCI referee Pannel 1 A
data_11 = read.delim(paste("data/Data11_supp.tab",sep=""),comment.char = "#")
data_11$major_introns.Var1 = as.character(data_11$major_introns.Var1 )
data_11$major_introns_busco.Var1 = as.character(data_11$major_introns_busco.Var1 )
p1A = ggplot(data_11,aes(y=major_introns.Freq * 100,x=major_introns.Var1 ))  + geom_boxplot(fill=set_color[2]) + 
  theme_bw()+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) +  scale_y_continuous(breaks=seq(0,100,20), labels=paste(seq(0,100,20),"%"),limits = c(0,70)) + 
  ggtitle("(all protein-coding genes)") + theme(legend.position = "none") + 
  ylab("Fraction of introns") + xlab("Phase") 
p1A

jpeg(paste(path_figure,"refPCI_p1A.jpg",sep=""), width = 3600/resolution, height = 2500/resolution,res=350/resolution)
print(p1A)
dev.off()


############## PCI referee Pannel 1 B
p1B = ggplot(data_11,aes(y=major_introns_busco.Freq*100,x=major_introns_busco.Var1 ))  + geom_boxplot(fill=set_color[1]) + 
  theme_bw()+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) + ggtitle("(BUSCO genes)") + theme(legend.position = "none") + ylab("Fraction of introns")+ xlab("Phase") + 
  scale_y_continuous(breaks=seq(0,100,20), labels=paste(seq(0,100,20),"%"),limits = c(0,70)) 
p1B

jpeg(paste(path_figure,"refPCI_p1B.jpg",sep=""), width = 3600/resolution, height = 2500/resolution,res=350/resolution)
print(p1B)
dev.off()


############## PCI referee Pannel 1 C
dt = data.frame(
  values = data_11[grepl("0",data_11$major_introns.Var1),]$major_introns_svr - data_11[grepl("1",data_11$major_introns.Var1),]$major_introns_svr,
  group = "P0 - P1"
)
dt = rbind(dt,data.frame(
  values = data_11[grepl("0",data_11$major_introns.Var1),]$major_introns_svr - data_11[grepl("2",data_11$major_introns.Var1),]$major_introns_svr,
  group = "P0 - P2"
))
dt = rbind(dt,data.frame(
  values = data_11[grepl("2",data_11$major_introns.Var1),]$major_introns_svr - data_11[grepl("1",data_11$major_introns.Var1),]$major_introns_svr,
  group = "P2 - P1"
))

p1C = ggplot(dt,aes(y=values ,x=group ))  + geom_boxplot(fill=set_color[2]) + 
  theme_bw()+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) + ggtitle("(all protein-coding genes)") + theme(legend.position = "none") + ylab("Average AS rate ")+ xlab("") + 
  scale_y_continuous(breaks=seq(-1,1,.2), labels=paste(seq(-1,1,.2),"%")) +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron differences")))
p1C

jpeg(paste(path_figure,"refPCI_p1C.jpg",sep=""), width = 3600/resolution, height = 2500/resolution,res=350/resolution)
print(p1C)
dev.off()


############## PCI referee Pannel 1 D
dt = data.frame(
  values = data_11[grepl("0",data_11$major_introns_busco.Var1),]$major_introns_busco_svr - data_11[grepl("1",data_11$major_introns_busco.Var1),]$major_introns_busco_svr,
  group = "P0 - P1"
)
dt = rbind(dt,data.frame(
  values = data_11[grepl("0",data_11$major_introns_busco.Var1),]$major_introns_busco_svr - data_11[grepl("2",data_11$major_introns_busco.Var1),]$major_introns_busco_svr,
  group = "P0 - P2"
))
dt = rbind(dt,data.frame(
  values = data_11[grepl("2",data_11$major_introns_busco.Var1),]$major_introns_busco_svr - data_11[grepl("1",data_11$major_introns_busco.Var1),]$major_introns_busco_svr,
  group = "P2 - P1"
))

p1D = ggplot(dt,aes(y=values ,x=group ))  + geom_boxplot(fill=set_color[1]) + 
  theme_bw()+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) + ggtitle("(BUSCO genes)") + theme(legend.position = "none") + ylab("Average AS rate ")+ xlab("") + 
  scale_y_continuous(breaks=seq(-1,1,.2), labels=paste(seq(-1,1,.2),"%")) +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron differences")))
p1D

jpeg(paste(path_figure,"refPCI_p1D.jpg",sep=""), width = 3600/resolution, height = 2500/resolution,res=350/resolution)
print(p1D)
dev.off()

############## PCI referee Figure 1

imgA = load.image(paste(path_figure,"p1A.jpg",sep=""))
imgB = load.image(paste(path_figure,"p1B.jpg",sep=""))
imgC = load.image(paste(path_figure,"p1C.jpg",sep=""))
imgD = load.image(paste(path_figure,"p1D.jpg",sep=""))

{
  pdf(file=paste(path_pannel,"Figure1.pdf",sep=""), width=5.4, height=4)
  
  m=matrix(c(1,3,2,4), nrow=2)
  
  
  m
  layout(m)
  
  par(mar=c(1, 0, 1, 0))
  plot(imgA, axes = F)
  mtext("A", side=2,at=0,adj=-1.5, line=1, font=2, cex=1.1,las=2)
  par(mar=c(1, 0, 1, 0))
  plot(imgB, axes = F)
  mtext("B", side=2,at=0,adj=-1.5, line=1, font=2, cex=1.1,las=2)
  par(mar=c(1, 0, 1, 0))
  plot(imgC, axes = F)
  mtext("C", side=2,at=50,adj=-1.5, line=1, font=2, cex=1.1,las=2)
  par(mar=c(1, 0, 1, 0))
  plot(imgD, axes = F)
  mtext("D", side=2,at=60,adj=-1.5, line=1, font=2, cex=1.1,las=2)
  dev.off()
}









