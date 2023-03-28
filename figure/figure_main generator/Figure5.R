source("figure/figure_main generator/library_path.R")


############## Pannel 5 A
data_5 = read.delim(paste("data/Data5_supp.tab",sep=""),comment.char = "#")
data_5$color_group = factor(data_5$color_group,levels = c("red","green","blue"))

df = data_5[data_5$filtering == "Drosophila_melanogaster_abundant_sv",]
{
  df$pos=c(3.73,3.13,4.48,2.46,1.75
           ,2.19,2.78,1.45,3.46,4.16)
  #, 3963396 SNPs
  p1  = ggplot(df,aes(x=pos,y=mean_polymorphism,fill=color_group)) + geom_col(width=0.1,col="black") + theme_bw()+
    geom_errorbar(aes(ymin=error_bar, ymax=error_bar_2),width=0.03,show.legend=FALSE)+ggtitle("Abundant SVs (all protein-coding genes)")+ 
    geom_text(data=df,aes(x=pos-0.07,y=mean_polymorphism+0.004, family="serif",label=paste(round(Nb_introns_minor,3))),angle=90,vjust=0,size=6)+
    theme(legend.position = "none") + xlim(c(1,5))  + ylab("SNP density (per bp)") +labs(y=expression(paste("SNP density (",italic("per")," bp)")))+
    theme( 
      axis.title.x = element_text(color=NA, size=NA,family="serif"),
      axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=28, family="serif"),
      axis.text.y =  element_text(color="black", size=22, family="serif"),
      axis.text.x =  element_text(color=NA, size=NA, family="serif"),
      title =  element_text(color="black", size=20, family="serif"),
      text =  element_text(color="black", size=31, family="serif"),
      legend.text =  element_text(color="black", size=26, family="serif"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank() ,
      panel.grid.major.y = element_line( size=.1, color="black" ) ,
    )  + ylim(c(0,0.055))
  p1
  
}



############## Pannel 5 C

df = data_5[data_5$filtering == "Homo_sapiens_abundant_sv",]

{
  df$pos=c(3.73,3.13,4.48,2.46,1.75
           ,2.19,2.78,1.45,3.46,4.16)
  p2  = ggplot(df,aes(x=pos,y=mean_polymorphism,fill=color_group)) + geom_col(width=0.1,col="black") + theme_bw()+
    geom_errorbar(aes(ymin=error_bar, ymax=error_bar_2),width=00.03,show.legend=FALSE)+ggtitle("  ")+
    geom_text(data=df,aes(x=pos-0.08,y=mean_polymorphism-0.002, family="serif",label=paste(round(Nb_introns_minor,3))),angle=90,vjust=0,size=6)+
    theme(legend.position = "none") + xlim(c(1,5))  + ylab("SNP density (per bp)") +labs(y=expression(paste("SNP density (",italic("per")," bp)")))+
    theme( 
      axis.title.x = element_text(color=NA, size=NA,family="serif"),
      axis.title.y = element_text(color="black", size=28,margin = margin(t = 0, r = 20, b = 0, l = 0), family="serif"),
      axis.text.y =  element_text(color="black", size=22, family="serif"),
      axis.text.x =  element_text(color=NA, size=NA, family="serif"),
      title =  element_text(color="black", size=20, family="serif"),
      text =  element_text(color="black", size=31, family="serif"),
      legend.text =  element_text(color="black", size=26, family="serif"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank() ,
      panel.grid.major.y = element_line( size=.1, color="black" ) 
    ) 
  
  
  p2
  
}
{
  p = ggdraw() + draw_plot(p1, 0, 0.55, 1, 0.42) + draw_plot(p2, 0, 0.13, 1, 0.42) + 
    draw_image(paste(path_require,"polymorphism_position.png",sep=""),0.08,-0.41,0.923,1) + 
    draw_image(paste(path_require,"Drosophila_melanogaster.png",sep=""),.86,.85,0.12,.06) + 
    draw_image(paste(path_require,"human.png",sep=""),.85,.39,0.15,.11)
  p
  
  resolution=1
  jpeg(paste(path_figure,"compile_polymorphism_frequent.jpg",sep=""), width = 3600/resolution, height = 4000/resolution,res=350/resolution)
  print(p)
  dev.off()
}



############## Pannel 5 B

df = data_5[data_5$filtering == "Drosophila_melanogaster_rare_sv",]

{
  df$pos=c(3.73,3.13,4.48,2.46,1.75
           ,2.19,2.78,1.45,3.46,4.16)
  
  p1  = ggplot(df,aes(x=pos,y=mean_polymorphism,fill=color_group)) + geom_col(width=0.1,col="black") + theme_bw()+
    geom_errorbar(aes(ymin=error_bar, ymax=error_bar_2),width=0.03,show.legend=FALSE)+ggtitle("Rare SVs (all protein-coding genes)")+ 
    geom_text(data=df,aes(x=pos-0.07,y=mean_polymorphism+0.004, family="serif",label=paste(round(Nb_introns_minor,3))),angle=90,vjust=0,size=6)+
    theme(legend.position = "none") + xlim(c(1,5))  + ylab("SNP density (per bp)") +labs(y=expression(paste("SNP density (",italic("per")," bp)")))+
    theme( 
      axis.title.x = element_text(color=NA, size=NA,family="serif"),
      axis.title.y = element_text(color="black", size=28,margin = margin(t = 0, r = 20, b = 0, l = 0), family="serif"),
      axis.text.y =  element_text(color="black", size=22, family="serif"),
      axis.text.x =  element_text(color=NA, size=NA, family="serif"),
      title =  element_text(color="black", size=20, family="serif"),
      text =  element_text(color="black", size=31, family="serif"),
      legend.text =  element_text(color="black", size=26, family="serif"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank() ,
      panel.grid.major.y = element_line( size=.1, color="black" ) 
    )  + ylim(c(0,0.055))
  
  
  p1
  
  p = ggdraw() + draw_plot(p1, 0, 0.1, 1, 0.8) + 
    draw_image(paste(path_require,"polymorphism_position.png",sep=""),0.055,-0.43,0.95,1)
  p
  
  resolution=1
}



############## Pannel 5 D
df = data_5[data_5$filtering == "Homo_sapiens_rare_sv",]
{
  df$pos=c(3.73,3.13,4.48,2.46,1.75
           ,2.19,2.78,1.45,3.46,4.16)
  
  p2  = ggplot(df,aes(x=pos,y=mean_polymorphism,fill=color_group)) + geom_col(width=0.1,col="black") + theme_bw()+
    geom_errorbar(aes(ymin=error_bar, ymax=error_bar_2),width=00.03,show.legend=FALSE)+ggtitle(" ")+ 
    geom_text(data=df,aes(x=pos-0.08,y=mean_polymorphism-0.002, family="serif",label=paste(round(Nb_introns_minor,3))),angle=90,vjust=0,size=6)+
    theme(legend.position = "none") + xlim(c(1,5))  + ylab("SNP density (per bp)") +labs(y=expression(paste("SNP density (",italic("per")," bp)")))+
    theme( 
      axis.title.x = element_text(color=NA, size=NA,family="serif"),
      axis.title.y = element_text(color="black", size=28,margin = margin(t = 0, r = 20, b = 0, l = 0), family="serif"),
      axis.text.y =  element_text(color="black", size=22, family="serif"),
      axis.text.x =  element_text(color=NA, size=NA, family="serif"),
      title =  element_text(color="black", size=20, family="serif"),
      text =  element_text(color="black", size=31, family="serif"),
      legend.text =  element_text(color="black", size=26, family="serif"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank() ,
      panel.grid.major.y = element_line( size=.1, color="black" ) 
    ) 
  
  
  p2
  
}


{
  p = ggdraw() + draw_plot(p1, 0, 0.55, 1, 0.42) + draw_plot(p2, 0, 0.13, 1, 0.42) + 
    draw_image(paste(path_require,"polymorphism_position.png",sep=""),0.08,-0.41,0.923,1) + 
    draw_image(paste(path_require,"Drosophila_melanogaster.png",sep=""),.86,.85,0.12,.06) + 
    draw_image(paste(path_require,"human.png",sep=""),.85,.39,0.15,.11)
  p
  
  resolution=1
  jpeg(paste(path_figure,"compile_polymorphism_rare.jpg",sep=""), width = 3600/resolution, height = 4000/resolution,res=350/resolution)
  print(p)
  dev.off()
}







#### Figure 5


imgA = load.image(paste(path_figure,"compile_polymorphism_frequent.jpg",sep=""))
imgB = load.image(paste(path_figure,"compile_polymorphism_rare.jpg",sep=""))

{
  pdf(file=  paste(path_pannel,"Figure5.pdf",sep=""), width=9, height=5)
  
  m=matrix(rep(NA,90*60), nrow=90)
  
  for(i in 1:90){
    m[i,]=c(rep(1,30),rep(2,30))
  }
  
  m
  layout(m)
  
  par(mar=c(0, 1.5, 0, 0))
  plot(imgA, axes=F)
  mtext("A",at=100,adj=0, side=2, line=1, font=2, cex=1.3,las=2)
  mtext("C",at=1900,adj=0, side=2, line=1, font=2, cex=1.3,las=2)
  
  par(mar=c(0, 1.5, 0, 0))
  plot(imgB, axes=F)
  mtext("B",at=100,adj=0, side=2, line=1, font=2, cex=1.3,las=2)
  mtext("D",at=1900,adj=0, side=2, line=1, font=2, cex=1.3,las=2)
  dev.off()
}
