source("figure/figure_main generator/library_path.R")

############## Pannel 2 A
data_3 = read.delim(paste("data/Data3_supp.tab",sep=""),comment.char = "#")

p1 =  ggplot(data_3,aes(x=rate*100,group=species,y=Freq))  + theme_bw() + ylab("Proportion of introns")+ 
  geom_line(size=0.2,col="grey")+
  geom_point(pch=21,col="grey",fill="grey") + ggtitle("All introns (all protein-coding genes)")+ 
  geom_line(data=data_3[data_3$species=="Drosophila_melanogaster",],size=2,col="red") +
  geom_line(data=data_3[data_3$species=="Homo_sapiens",],size=2,col="#66281A") +
  geom_point(data=data_3[data_3$species=="Drosophila_melanogaster",],size=4,pch=21,fill="red") +
  geom_point(data=data_3[data_3$species=="Homo_sapiens",],size=4,pch=21,fill="#66281A") + theme_bw()+
  theme(
    axis.title.x = element_text(color="black", size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=31, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) + scale_x_continuous(breaks=seq(0,100,25), labels=paste(seq(0,100,25),"%"))+ scale_y_continuous(breaks=seq(0,100,25), labels=paste(seq(0,100,25),"%"),limits=c(0,80)) +
  labs(x=expression(paste("RAS ",italic("per")," intron")))
p1


p = ggdraw() + draw_plot(p1, 0, 0, 1, 1)  
p



{
  p = ggdraw()+ draw_plot(p,0,0,1,1)  + 
    draw_image(paste(path_require,"Drosophila_melanogaster_red.png",sep=""),.26,.775,.19,.065)+
    draw_image(paste(path_require,"human_brown.png",sep=""),.245,.73,.08,.17) 
  p
  jpeg(paste(path_figure,"p5_hist_ras.jpg",sep=""), width = 4000/resolution, height = 2500/resolution,res=250/resolution)
  print(p)
  dev.off()
}

############## Pannel 2 B

p1 =  ggplot(data_3,aes(x=rate*100,group=species,y=Freq.1))  + theme_bw() + ylab("Proportion of introns")+ 
  geom_line(size=0.2,col="grey")+
  geom_point(pch=21,col="grey",fill="grey") +ggtitle("All introns (all protein-coding genes)")+
  geom_line(data=data_3[data_3$species=="Drosophila_melanogaster",],size=2,col="red") +
  geom_line(data=data_3[data_3$species=="Homo_sapiens",],size=2,col="#66281A") +
  geom_point(data=data_3[data_3$species=="Drosophila_melanogaster",],size=4,pch=21,fill="red") +
  geom_point(data=data_3[data_3$species=="Homo_sapiens",],size=4,pch=21,fill="#66281A") +theme_bw()+
  theme(
    axis.title.x = element_text(color="black", size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=31, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) + scale_x_continuous(breaks=seq(0,100,25), labels=paste(seq(0,100,25),"%"))+ scale_y_continuous(breaks=seq(0,100,25), labels=paste(seq(0,100,25),"%"),limits=c(0,80))+
  labs(x=expression(paste("RANS ",italic("per")," intron")))
p1


p = ggdraw() + draw_plot(p1, 0, 0, 1, 1)  
p


{
  p = ggdraw()+ draw_plot(p,0,0,1,1)  + 
    draw_image(paste(path_require,"Drosophila_melanogaster_red.png",sep=""),.26,.775,.19,.065)+
    draw_image(paste(path_require,"human_brown.png",sep=""),.245,.73,.08,.17) 
  p
  jpeg(paste(path_figure,"p6_hist_rans.jpg",sep=""), width = 4000/resolution, height = 2500/resolution,res=250/resolution)
  print(p)
  dev.off()
}

#### Figure 2

imgA = load.image(paste(path_require,"ns_na_nu.png",sep=""))
imgB = load.image(paste(path_figure,"p5_hist_ras.jpg",sep=""))
imgC = load.image(paste(path_figure,"p6_hist_rans.jpg",sep=""))
imgD = load.image(paste(path_require,"mira.png",sep=""))

{
  pdf(file=paste(path_pannel,"Figure2.pdf",sep=""), width=6.75, height=6/1.1)
  
  m=matrix(rep(1,10*17), nrow=17)
  
  for(i in 6:11){
    m[i,]=c(rep(2,5),rep(3,5))
  }
  
  for(i in 12:17){
    m[i,]=c(rep(4,10))
  }
  m
  layout(m)
  
  par(mar=c(0, 0, 1, 0))
  plot(imgA, axes = F)
  mtext("A", side=2,at=10,adj=-5, line=1, font=2, cex=1.2,las=2)
  par(mar=c(0, 0, 1, 0))
  plot(imgB, axes = F)
  mtext("B",side=2,at=30,adj=-1.5,  line=1, font=2, cex=1.2,las=2)
  par(mar=c(0, 0, 1, 0))
  plot(imgC, axes = F)
  mtext("C", side=2,at=40,adj=-1.5, line=1, font=2, cex=1.2,las=2)
  par(mar=c(0, 0, 1, 0))
  plot(imgD, axes = F)
  mtext("D", side=2,at=10,adj=-5, line=1, font=2, cex=1.2,las=2)
  dev.off()
}
