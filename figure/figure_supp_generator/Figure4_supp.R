source("figure/figure_supp generator/library_path.R")



############## Supplementary Pannel 4 A

data_5 = read.delim(paste("data/Data5_supp.tab",sep=""),comment.char = "#")
data_5[!is.na(data_5$pvalue_vs_control) & data_5$pvalue_vs_control < 0.05,"significance"] = "*"
data_5[!is.na(data_5$pvalue_vs_control) & data_5$pvalue_vs_control < 0.005,"significance"] = "**"
data_5[!is.na(data_5$pvalue_vs_control) & data_5$pvalue_vs_control < 0.0005,"significance"] = "***"

data_5$color_group = factor(data_5$color_group,levels = c("red","green","blue"))

df = data_5[data_5$filtering == "Homo_sapiens_CpG_abundant_sv",]

df$pos=c(2.19,2.78,1.45,3.46,4.16)

p4A  = ggplot(df,aes(x=pos,y=mean_polymorphism,fill=color_group)) + geom_col(width=0.1,col="black") + theme_bw()+
  geom_errorbar(aes(ymin=error_bar, ymax=error_bar_2),width=00.03,show.legend=FALSE)+ggtitle("Abundant SVs (all protein-coding genes)")+ 
  geom_text(data=df,aes(x=pos-0.07,y=mean_polymorphism+0.004, family="serif",label=paste(round(Nb_introns_minor,3))),angle=90,vjust=0,size=6)+
  # geom_text(data=df,aes(x=pos+0.07,y=mean_polymorphism+0.004, family="serif",label=significance),angle=90,vjust=0,size=6)+
  theme(legend.position = "none") + xlim(c(1,5))  +labs(y=expression(paste("SNP density (",italic("per")," bp)")))+
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
  ) +scale_y_continuous(limits=c(0,0.185))


p4A

resolution=1


p4A = ggdraw() + draw_plot(p4A, 0, 0.25, 1, .7) + draw_image(paste(path_require,"polymorphism_position_CpG.png",sep=""),0.09,-0.31,0.91,1)+ 
  draw_image(paste(path_require,"human.png",sep=""),.85,.65,0.15,.17)
p4A

resolution=1
jpeg(paste(path_figure,"supp_p4A.jpg",sep=""), width = 3600/resolution, height = 2500/resolution,res=350/resolution)
print(p4A)
dev.off()



############## Supplementary Pannel 4 B
df = data_5[data_5$filtering == "Homo_sapiens_CpG_rare_sv",]

df$pos=c(2.19,2.78,1.45,3.46,4.16)



p4B  = ggplot(df,aes(x=pos,y=mean_polymorphism,fill=color_group)) + geom_col(width=0.1,col="black") + theme_bw()+
  geom_errorbar(aes(ymin=error_bar, ymax=error_bar_2),width=00.03,show.legend=FALSE)+ggtitle("Rare SVs (all protein-coding genes)")+ 
  geom_text(data=df,aes(x=pos-0.07,y=mean_polymorphism+0.004, family="serif",label=paste(round(Nb_introns_minor,3))),angle=90,vjust=0,size=6)+
  # geom_text(data=df,aes(x=pos+0.07,y=mean_polymorphism+0.004, family="serif",label=significance),angle=90,vjust=0,size=6)+
  theme(legend.position = "none") + xlim(c(1,5)) +labs(y=expression(paste("SNP density (",italic("per")," bp)")))+
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
  ) +scale_y_continuous(limits=c(0,0.185))


p4B

resolution=1


p4B = ggdraw() + draw_plot(p4B, 0, 0.25, 1, .7) + draw_image(paste(path_require,"polymorphism_position_CpG.png",sep=""),0.09,-0.31,0.91,1)+ 
  draw_image(paste(path_require,"human.png",sep=""),.85,.65,0.15,.17)
p4B

resolution=1
jpeg(paste(path_figure,"supp_p4B.jpg",sep=""), width = 3600/resolution, height = 2500/resolution,res=350/resolution)
print(p4B)
dev.off()



############## Supplementary Figure 4

imgA = load.image(paste(path_figure,"supp_p4A.jpg",sep=""))
imgB = load.image(paste(path_figure,"supp_p4B.jpg",sep=""))

{
  pdf(file=paste(path_pannel,"Figure4_supp.pdf",sep=""), width=4, height=6)
  
  m=matrix(c(1,2), nrow=2)
  m
  layout(m)
  
  par(mar=c(1, 0, 1, 0))
  plot(imgA, axes = F)
  mtext("A", side=2,at=0,adj=-3, line=1, font=2, cex=1.3,las=2)
  par(mar=c(1, 0, 1, 0))
  plot(imgB, axes = F)
  mtext("B", side=2,at=0,adj=-3, line=1, font=2, cex=1.3,las=2)
  dev.off()
}



