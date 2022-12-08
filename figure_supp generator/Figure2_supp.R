source("figure_supp generator/library_path.R")


############## Supplementary Pannel 2 A
data_3 = read.delim(paste("data/Data3_supp.tab",sep=""),comment.char = "#")

p1 =  ggplot(data_3,aes(x=rate*100,group=species,y=Freq))  + theme_bw() + ylab("Proportion of introns")+ 
  geom_line(size=0.5,col="grey")+
  geom_point(pch=21,col="grey",fill="grey",size=.5) + ggtitle("All introns (all protein-coding genes)")+ 
  geom_line(data=data_3[data_3$species=="Drosophila_melanogaster",],size=.7,col="red") +
  geom_line(data=data_3[data_3$species=="Homo_sapiens",],size=.7,col="#66281A") +
  geom_point(data=data_3[data_3$species=="Drosophila_melanogaster",],alpha=0.7,pch=21,size=1,fill="red") +
  geom_point(data=data_3[data_3$species=="Homo_sapiens",],alpha=0.7,pch=21,size=1,fill="#66281A") +theme_bw()+
  theme(
    axis.title.x = element_text(color="black", size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=31, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) + scale_x_continuous(breaks=seq(0,100,25), labels=paste(seq(0,100,25),"%"))+ scale_y_continuous(breaks=seq(0,100,25), labels=paste(seq(0,100,25),"%"),limits=c(0,75)) +
  labs(x=expression(paste("RAS ",italic("per")," intron")))
p1


p = ggdraw() + draw_plot(p1, 0, 0, 1, 1)  
p

jpeg(paste(path_figure,"p5_hist_ras.png",sep=""), width = 6000/resolution, height = 2500/resolution,res=350/resolution)
print(p)
dev.off()







############## Supplementary Pannel 2 B

p1 =  ggplot(data_3,aes(x=rate*100,group=species,y=Freq.1))  + theme_bw() + ylab("Proportion of introns")+ 
  geom_line(size=0.5,col="grey")+
  geom_point(pch=21,col="grey",fill="grey") +ggtitle("All introns (all protein-coding genes)")+
  geom_line(data=data_3[data_3$species=="Drosophila_melanogaster",],size=.7,col="red") +
  geom_line(data=data_3[data_3$species=="Homo_sapiens",],size=.7,col="#66281A") +
  geom_point(data=data_3[data_3$species=="Drosophila_melanogaster",],alpha=0.7,pch=21,size=1,fill="red") +
  geom_point(data=data_3[data_3$species=="Homo_sapiens",],alpha=0.7,pch=21,size=1,fill="#66281A") +theme_bw()+
  theme(
    axis.title.x = element_text(color="black", size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=31, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) + scale_x_continuous(breaks=seq(0,100,25), labels=paste(seq(0,100,25),"%"))+ scale_y_continuous(breaks=seq(0,100,25), labels=paste(seq(0,100,25),"%"),limits=c(0,75))+
  labs(x=expression(paste("RANS ",italic("per")," intron")))
p1


p = ggdraw() + draw_plot(p1, 0, 0, 1, 1)  
p

jpeg(paste(path_figure,"p6_hist_rans.png",sep=""), width = 6000/resolution, height = 2500/resolution,res=350/resolution)
print(p)
dev.off()




############## Supplementary Figure 2 
imgA = load.image(paste(path_figure,"p5_hist_ras.png",sep=""))
imgB = load.image(paste(path_figure,"p6_hist_rans.png",sep=""))

{
  pdf(file=paste(path_pannel,"Figure2_supp.pdf",sep=""), width=4, height=4)
  
  m=matrix(c(1,2), nrow=2)
  
  
  m
  layout(m)
  
  par(mar=c(1, 0, 1, 0))
  plot(imgA, axes = F)
  mtext("A", side=2,at=0,adj=-2, line=1, font=2, cex=1.1,las=2)
  par(mar=c(1, 0, 1, 0))
  plot(imgB, axes = F)
  mtext("B", side=2,at=0,adj=-2, line=1, font=2, cex=1.1,las=2)
  dev.off()
}


