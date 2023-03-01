source("figure generator/library_path.R")


############## Pannel 5 A
data_2 = read.delim(paste("data/Data2_supp.tab",sep=""),comment.char = "#")
data_2 = data_2[data_2$species %in% arbrePhylo$tip.label,]
data_2$average_mira=data_2$average_mira * 100
data_2$framepreserving_proportion=data_2$framepreserving_proportion * 100

p1 = ggplot(data_2,aes(x=average_mira,y=framepreserving_proportion,fill=species,label=species))  + theme_bw()+
  scale_x_log10(breaks=c(0.001,0.01,0.1,1,10, 100), labels=c("0.001 %","0.01 %","0.1 %","1 %","10 %", "100 %"))+
  coord_cartesian(xlim=c(0.001,100)) +  ylab("Proportion of\nframe-preserving variants")+ 
  scale_y_continuous(breaks=seq(0,100,10), labels=paste(seq(0,100,10),"%")) +
  geom_line(size=0.5,col="grey")+ xlab("Minor intron relative abundance (MIRA)")+
  geom_point(size=1,pch=21,col="grey",fill="grey") +
  geom_vline(xintercept=5, linetype="dashed", color = "black", size=1,alpha=1)+
  geom_hline(yintercept=33, linetype="dashed", color = "black", size=1,alpha=0.5)+
  geom_text(label="5%",x=.5,y=70,cex=13, family="serif")+
  geom_text(label="33%",x=2,y=37,cex=13, family="serif",alpha=0.5)+
  
  geom_line(data=data_2[data_2$species=="Drosophila_melanogaster",],size=2,aes(col=species),col="red")+
  geom_point(data=data_2[data_2$species=="Drosophila_melanogaster",],pch=21,size=4,aes(fill=species),fill="red") +
  
  geom_line(data=data_2[data_2$species=="Apis_mellifera",],size=2,aes(col=species),col="#ba8e18")+
  geom_point(data=data_2[data_2$species=="Apis_mellifera",],pch=21,size=4,aes(fill=species),fill="#ba8e18") +
  
  geom_line(data=data_2[data_2$species=="Homo_sapiens",],size=2,aes(col=species),col="#66281A")+
  geom_point(data=data_2[data_2$species=="Homo_sapiens",],pch=21,size=4,aes(fill=species),fill="#66281A") +
  theme(
    axis.title.x = element_text(color="black", size=37,margin = margin(t = 15, r = 0, b = 0, l = 0),family="serif"),
    axis.title.y = element_text(color="black", size=37, family="serif"),
    axis.text.y =  element_text(color="black", size=33, family="serif"),
    axis.text.x =  element_text(color="black", size=33, family="serif"),
    title =  element_text(color="black", size=31, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  )

p1

jpeg(paste(path_figure,"prop_md3_53sp.jpg",sep=""), width = 9900/resolution, height = 6000/resolution,res=700/resolution)
print(p1)
dev.off()




############## Pannel 5 B

data_8 = read.delim(paste("data/Data8_supp.tab",sep=""),comment.char = "#")
data_8$clade = data_1[data_8$species,]$clade
data_8$longevity = data_1[data_8$species,]$longevity
data_8$ratio =  data_8$prop_fp_sv_abundant * 100

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_8[,"species"],
                                                     xlabel=data_8[,"longevity"],
                                                     ylabel=data_8[,"ratio"]), species, vcv=TRUE)




p4 = ggplot(data_8, aes(x=longevity,y=ratio,fill=clade))  + theme_bw() + 
  xlab("Longevity (days, log scale)") + ylab("Proportion of frame-preserving SVs")+
  geom_point(pch=21,alpha=0.7, size=7) + 
  scale_y_continuous(breaks=seq(0,100,10), labels=paste(seq(0,100,10),"%")) +
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,50,100,1000,10000,50000), limits=c(7,50000))+ 
  scale_fill_manual("Clades",values=vectorColor) + labs(fill="Clades") +
  scale_color_manual("Clades",values=vectorColor) + ggtitle("Abundant SVs (all protein-coding genes)")+
  theme(
    axis.title.x = element_text(color="black", size=31,margin = margin(t = 15, r = 0, b = 0, l = 0),family="serif"),
    axis.title.y = element_text(color="black", size=31,margin = margin(t = 0, r = 20, b = 0, l = 0), family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.395, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  )  + annotation_logticks(sides="b")

p4

jpeg(paste(path_figure,"md3_10_percent.jpg",sep=""), width = 8300/resolution, height = 5300/resolution,res=700/resolution)
print(p4)
dev.off()




############## Pannel 5 C

ylabel = "mean_as_busco_low_as"
xlabel = "longevity"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p6=ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual(values=vectorColor)+ ggtitle("Low-AS major introns (BUSCO genes)")+
  theme_bw() +labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Longevity (days, log scale)") +
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,50,100,1000,10000,50000), limits=c(7,53000)) +
  scale_y_continuous(breaks=seq(0,10,0.2), labels=paste("",seq(0,10,0.2),"%")) +
  theme(
    axis.title.x = element_text(color="black", size=31,margin = margin(t = 15, r = 0, b = 0, l = 0),family="serif"),
    axis.title.y = element_text(color="black", size=31,margin = margin(t = 0, r = 20, b = 0, l = 0), family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif"),
    plot.caption = element_text(hjust = 0.735, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  ) + theme(legend.position = "none") + annotation_logticks(sides="b")
p6


jpeg(paste(path_figure,"busco_notfunctional_svr.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p6)
dev.off()





############## Pannel 5 D

ylabel="mean_as_busco_high_as"
xlabel="longevity"


shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p5 = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual(values=vectorColor)+ ggtitle("High-AS major introns (BUSCO genes)") +
  theme_bw() +
  ylab("Average AS rate per intron")+labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Longevity (days, log scale)")  +
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,50,100,1000,10000,50000), limits=c(7,53000)) +
  scale_y_continuous(breaks=seq(10,25,2), labels=paste(seq(10,25,2),"%")) +
  theme(
    axis.title.x = element_text(color="black", size=31,margin = margin(t = 15, r = 0, b = 0, l = 0),family="serif"),
    axis.title.y = element_text(color="black", size=31,margin = margin(t = 0, r = 20, b = 0, l = 0), family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif"),
    plot.caption = element_text(hjust = 0.71, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  ) + theme(legend.position = "none") + annotation_logticks(sides="b")

p5

jpeg(paste(path_figure,"busco_functional_svr.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p5)
dev.off()





######### Figure 5

{
  p = ggdraw()+ draw_plot(p1,0,0,1,1)  + 
    draw_image(paste(path_require,"Drosophila_melanogaster_red.png",sep=""),.29,.845,.19,.065)+
    draw_image(paste(path_require,"bee_yellow.png",sep=""),.275,.78,.085,.2) +
    draw_image(paste(path_require,"human_brown.png",sep=""),.2,.80,.08,.17) 
  p
  jpeg(paste(path_figure,"prop_md3_53sp.jpg",sep="") , width = 9900/resolution, height = 6000/resolution,res=700/resolution)
  print(p)
  dev.off()
}

imgA = load.image(paste(path_figure,"prop_md3_53sp.jpg",sep=""))
imgB = load.image(paste(path_figure,"md3_10_percent.jpg",sep=""))
imgC = load.image(paste(path_figure,"busco_notfunctional_svr.jpg",sep=""))
imgD = load.image(paste(path_figure,"busco_functional_svr.jpg",sep=""))



{ 
  pdf(file=   paste(path_pannel,"Figure5.pdf",sep="")  , width=10, height=7)
  
  m = matrix(rep(NA,10*10), nrow=10)
  
  for(i in 1:5){
    m[i,]=c(rep(1,5),rep(2,5))
  }
  
  for(i in 6:10){
    m[i,]=c(rep(3,5),rep(4,5))
  }
  
  m
  layout(m)
  
  showAxes = F
  
  par(mar=c(0, 1,0, 0))
  plot(imgA, axes=showAxes)
  mtext("A",at=0, adj=-1, side=2, line=1, font=2, cex=1.7,las=2)
  
  par(mar=c(0,2, 0.5, 0))
  plot(imgB, axes=showAxes)
  mtext("B",at=0,adj=0.5, side=2, line=1, font=2, cex=1.7,las=2)
  
  par(mar=c(1, 5.5,2, 5))
  plot(imgC, axes=showAxes)
  mtext("C",adj=1,at=20, side=2, line=1, font=2, cex=1.7,las=2)
  
  par(mar=c(1, 2,2, 8.5))
  plot(imgD, axes=showAxes)
  mtext("D",adj=0.5,at=0, side=2, line=1, font=2, cex=1.7,las=2)
  dev.off()
}