source("figure/figure_supp generator/library_path.R")
color="red"

############## Supplementary Pannel 1 A
species="Homo_sapiens"
data_7 = read.delim(paste("data/Data7_supp.tab",sep=""),comment.char = "#")

df = data_7[data_7$species==species,]

p3 = ggplot(df[df$echantillon == "all introns",],aes(x=sequencing_depth,y=annot_N1N2_sup10/annotated_intron*100)) + 
  geom_point(size=5,shape=21,alpha=0.8,aes(fill=no_compiled_rna_seq)) +
  xlab("Median read coverage on BUSCO genes (reads/bp)") +
  geom_vline(xintercept=200, linetype="dashed", color = color, size=1,alpha=0.5) + theme_bw()+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) + 
  scale_y_continuous(breaks=c(0,25,50,75,100),limits= c(0,100), labels = paste(c(0,25,50,75,100),"%",sep="")) +
  scale_x_continuous(limits = c(0,800))+ggtitle("Annotated introns (BUSCO genes)")+
  ylab(paste("Proportion with N\u226510" ) )+ 
  scale_fill_gradient("Number of\nRNA-seq\nsamples", low = "grey", high = color)
p3

p = ggdraw() + draw_plot(p3, 0, 0, 1, 1)+
  draw_image(paste(path_require,"human.png",sep=""),.61,.2,0.15,.2)

p

jpeg(paste(path_figure,"p1_prop_seqdepth_hsap.jpg",sep=""), width = 3600/resolution, height = 2500/resolution,res=270/resolution)
print(p)
dev.off()


############## Supplementary Pannel 1 B
species="Drosophila_melanogaster"
df = data_7[data_7$species==species,]


p4 =ggplot(df[df$echantillon == "all introns",],aes(x=sequencing_depth,y=annot_N1N2_sup10/annotated_intron*100)) + 
  geom_point(size=5,shape=21,alpha=0.8,aes(fill=no_compiled_rna_seq)) +
  xlab("Median read coverage on BUSCO genes (reads/bp)") +ggtitle("Annotated introns (BUSCO genes)")+
  geom_vline(xintercept=200, linetype="dashed", color = color, size=1,alpha=0.5) + theme_bw()+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  )+ scale_x_continuous(breaks=c(0,200,400,600,800,1000,1200))+ 
  scale_y_continuous(breaks=c(0,25,50,75,100),limits= c(0,100), labels = paste(c(0,25,50,75,100),"%",sep=""))+
  ylab(paste( "Proportion with N\u226510" ) )+
  scale_fill_gradient("Number of\nRNA-seq\nsamples",
                      low = "grey", high = color)
p4


p = ggdraw() + draw_plot(p4, 0, 0, 1, 1) +
  draw_image(paste(path_require,"Drosophila_melanogaster.png",sep=""),.61,.2,0.15,.11)
p

jpeg(paste(path_figure,"p2_prop_seqdepth_dmel.jpg",sep=""), width = 3600/resolution, height = 2500/resolution,res=270/resolution)
print(p)
dev.off()



############## Supplementary Pannel 1 C

species="Homo_sapiens"
df = data_7[data_7$species==species,]

p3 = ggplot(df[df$echantillon == "all introns",],aes(x=sequencing_depth,y=average_svr_busco*100)) + 
  geom_point(size=5,shape=21,alpha=0.8,aes(fill=no_compiled_rna_seq)) +
  xlab("Median read coverage on BUSCO genes (reads/bp)") +
  geom_vline(xintercept=200, linetype="dashed", color = color, size=1,alpha=0.5) + theme_bw()+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  ) + 
  scale_y_continuous(breaks=seq(0.5,5,0.5), labels=paste(seq(0.5,5,0.5),"%")) +
  scale_x_continuous(limits = c(0,800))+
   ggtitle("Major introns (BUSCO genes)")+labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  scale_fill_gradient("Number of\nRNA-seq\nsamples", low = "grey", high = color)
p3

p = ggdraw() + draw_plot(p3, 0, 0, 1, 1)+
  draw_image(paste(path_require,"human.png",sep=""),.61,.2,0.15,.2)

p

jpeg(paste(path_figure,"p3_as_seqdepth_dmel.jpg",sep=""), width = 3600/resolution, height = 2500/resolution,res=270/resolution)
print(p)
dev.off()


############## Supplementary Pannel 1 D


species="Drosophila_melanogaster"
df = data_7[data_7$species==species,]

p4 =ggplot(df[df$echantillon == "all introns",],aes(x=sequencing_depth,y=average_svr_busco*100)) + 
  geom_point(size=5,shape=21,alpha=0.8,aes(fill=no_compiled_rna_seq)) +
  xlab("Median read coverage on BUSCO genes (reads/bp)") +
  geom_vline(xintercept=200, linetype="dashed", color = color, size=1,alpha=0.5) + theme_bw()+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 20, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif")
  )+ ggtitle("Major introns (BUSCO genes)")+
  scale_x_continuous(breaks=c(0,200,400,600,800,1000,1200))+ 
  scale_y_continuous(breaks=seq(0,5,0.5), labels=paste(seq(0,5,0.5),"%")) +
   labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  scale_fill_gradient("Number of\nRNA-seq\nsamples",
                      low = "grey", high = color)
p4


p = ggdraw() + draw_plot(p4, 0, 0, 1, 1) +
  draw_image(paste(path_require,"Drosophila_melanogaster.png",sep=""),.61,.2,0.15,.11)
p

jpeg(paste(path_figure,"p4_as_seqdepth_dmel.jpg",sep=""), width = 3600/resolution, height = 2500/resolution,res=270/resolution)
print(p)
dev.off()



############## Supplementary Figure 1

imgA = load.image(paste(path_figure,"p1_prop_seqdepth_hsap.jpg",sep=""))
imgB = load.image(paste(path_figure,"p2_prop_seqdepth_dmel.jpg",sep=""))
imgC = load.image(paste(path_figure,"p3_as_seqdepth_dmel.jpg",sep=""))
imgD = load.image(paste(path_figure,"p4_as_seqdepth_dmel.jpg",sep=""))

{
  pdf(file=paste(path_pannel,"Figure1_supp.pdf",sep=""), width=5.4, height=4)
  
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










