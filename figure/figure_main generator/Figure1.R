source("figure/figure_main generator/library_path.R")


############## Pannel 1 B
ylabel="body_size"
xlabel="longevity"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p1=ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,50,100,1000,10000,50000))+ 
  scale_y_log10(breaks=c(0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000),limits = c(0.01,1000))+ theme_bw() +
  ylab("Body length (cm, log scale)")+
  xlab("Longevity (days, log scale)")+ theme(
    axis.title.x = element_text(color="black", size=31,family="serif"),
    axis.title.y = element_text(color="black", size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=31, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.4, face= "italic", size=23, family="serif"),
    plot.caption.position =  "plot"
  )+
  labs(
    caption = substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls(log10(ylabel)~log10(xlabel),shorebird))))
  )+ annotation_logticks()
p1

jpeg(paste(path_figure,"history_traits_cor.jpg",sep=""), width = 8500/resolution, height = 6000/resolution,res=700/resolution)
print(p1)
dev.off()






############## Pannel 1 C
ylabel="dNdS_200k"
xlabel="longevity"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p2=ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ 
  scale_x_log10(breaks=c(0.01,0.05,0.1,0.5,1,5,10,50,100,1000,10000,50000)) + theme_bw() +
  ylab("dN/dS")+
  xlab("Longevity (days, log scale)")+ theme(
    axis.title.x = element_text(color="black", size=31,family="serif"),
    axis.title.y = element_text(color="black", size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=31 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.38, face= "italic", size=23, family="serif"),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  ) + scale_y_continuous(breaks=c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12), labels =c(0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12)) + annotation_logticks(sides="b")
p2

jpeg(paste(path_figure,"dNdS_longevity_cor.jpg",sep=""),width = 8500/resolution, height = 6000/resolution,res=700/resolution)
print(p2)
dev.off()

############## Pannel 1 A
nodes<-sapply(arbrePhylo$tip.label,function(x,y) which(y==x),y=arbrePhylo$tip.label)
colorSelected = vectorColor[data_1[arbrePhylo$tip.label,"clade"]]
arbrePhylo$tip.label = paste("  ",data_1[arbrePhylo$tip.label,"species"])
arbrePhylo$tip.label = str_replace_all(arbrePhylo$tip.label,"_"," ")

jpeg(paste(path_figure,"tree.jpg",sep=""), width = 6000/resolution, height = 7000/resolution,res=500/resolution)
plot(arbrePhylo, 
     use.edge.length = T, 
     tip.color=colorSelected, 
     align.tip.label = TRUE,
     family="mono",cex=1.3,lwd=2)
add.scale.bar(lwd=3)
dev.off()





#### Figure 1

imgB = load.image(paste(path_figure,"history_traits_cor.jpg",sep="") )
imgC = load.image(paste(path_figure,"dNdS_longevity_cor.jpg",sep="") )



monkey<-readPNG(paste(path_require,"monkey.png",sep=""))
fly<-readPNG(paste(path_require,"fly.png",sep=""))
ants<-readPNG(paste(path_require,"ants.png",sep=""))
wasp<-readPNG(paste(path_require,"wasp.png",sep=""))
bee<-readPNG(paste(path_require,"bee.png",sep=""))
hemiptera<-readPNG(paste(path_require,"hemiptera.png",sep=""))
crocodile<-readPNG(paste(path_require,"crocodile.png",sep=""))
termite<-readPNG(paste(path_require,"termite.png",sep=""))
Bombyx_mori<-readPNG(paste(path_require,"Bombyx_mori.png",sep=""))
Drosophila_melanogaster<-readPNG(paste(path_require,"Drosophila_melanogaster.png",sep=""))


{
  pdf(file= paste(path_pannel,"Figure1.pdf",sep=""), width=6.75, height=5.5)
  
  m=matrix(rep(NA,100*10), nrow=10)
  
  for(i in 1:100){
    m[,i]=rep(1)
  }
  
  for(i in 50:100){
    m[,i]=c(rep(2,5),rep(3,5))
  }
  layout(m)
  
  par(mar=c(0, 2, 0, 0))
  plot(arbrePhylo,
       use.edge.length = T,
       tip.color=colorSelected,
       align.tip.label = TRUE,
       family="mono",
       cex=0.6,
       edge.width=0.7,
       x.lim =c(0,2.2))
  
  xfly=1.8
  yfly=30
  rasterImage(fly,xleft=0+xfly, ybottom=0+yfly, xright=.2/1.5+xfly, ytop=2.5/1.5+yfly)
  
  xants=1.9
  yants=50
  rasterImage(ants,xleft=0+xants, ybottom=0+yants, xright=.2/1.5+xants, ytop=3/1.5+yants)
  
  xwasp=1.9
  ywasp=40.6
  rasterImage(wasp,xleft=0+xwasp, ybottom=0+ywasp, xright=.3/1.5+xwasp, ytop=3.1/1.5+ywasp)
  
  xbee=1.7
  ybee=45
  rasterImage(bee,xleft=0+xbee, ybottom=0+ybee, xright=.3/1.6+xbee, ytop=2.5/1.5+ybee)
  
  xBombyx_mori=1.7
  yBombyx_mori=26.4
  rasterImage(Bombyx_mori,xleft=0+xBombyx_mori, ybottom=0+yBombyx_mori, xright=.3/1.5+xBombyx_mori, ytop=2/1.5+yBombyx_mori)
  
  xhemiptera=2
  yhemiptera=23
  rasterImage(hemiptera,xleft=0+xhemiptera, ybottom=0+yhemiptera, xright=.2/1.5+xhemiptera, ytop=2.5/1.5+yhemiptera)
  
  xtermite=1.9
  ytermite=18.5
  rasterImage(termite,xleft=0+xtermite, ybottom=0+ytermite, xright=.2/1.2+xtermite, ytop=2.5/1.2+ytermite)
  
  xmonkey=1.8
  ymonkey=9
  rasterImage(monkey,xleft=0+xmonkey, ybottom=0+ymonkey, xright=.3/1.5+xmonkey, ytop=4/1.5+ymonkey)
  
  xcrocodile=1.9
  ycrocodile=2
  rasterImage(crocodile,xleft=0+xcrocodile, ybottom=0+ycrocodile, xright=.35/1.5+xcrocodile, ytop=3/1.5+ycrocodile)
  
  xDrosophila_melanogaster=1.9
  yDrosophila_melanogaster=34.5
  rasterImage(Drosophila_melanogaster,xleft=0+xDrosophila_melanogaster, ybottom=0+yDrosophila_melanogaster, xright=.15/1.6+xDrosophila_melanogaster, ytop=1.5/1.5+yDrosophila_melanogaster)
  
  mtext("A",at=49.4,adj=0.1, side=2, line=1, font=2, cex=1.4,las=2)
  plot(imgB, axes=FALSE)
  mtext("B", side=2,at=1, line=1, font=2, cex=1.4,las=2)
  plot(imgC, axes=FALSE)
  mtext("C", side=2,at=1, line=1, font=2, cex=1.4,las=2)
  dev.off()
}

