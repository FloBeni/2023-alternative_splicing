source("figure generator/library_path.R")


############## Pannel 2 B
ylabel="mean_as_busco"
xlabel="longevity"

shorebird <- comparative.data(arbrePhylo, data.frame(species=data_1[,"species"],xlabel=data_1[,xlabel],ylabel=data_1[,ylabel]), species, vcv=TRUE)

p1 = ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade,text=species) )+ geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual("Clades",values=vectorColor)+ ggtitle("Major introns (BUSCO genes)")+ 
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,50,100,1000,10000,50000), limits=c(7,50000)) + theme_bw() +
  scale_y_continuous(breaks=seq(0.5,4.5,0.5), labels=paste(seq(0.5,4.5,0.5),"%"),limits=c(.5,4)) +
  labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Longevity (days, log scale)")+ theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26 ,family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 10)),
    plot.caption = element_text(hjust = 0.4, face= "italic", size=23),
    plot.caption.position =  "plot"
  )+
  labs(
    caption =substitute(paste("PGLS model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(pgls((ylabel)~log10(xlabel),shorebird))))
  )  + annotation_logticks(sides="b")
p1



jpeg(paste(path_figure,"svr_busco_longevity.jpg",sep=""), width = 8500/resolution, height = 6000/resolution,res=700/resolution)
print(p1)
dev.off()




############## Pannel 2 C
point_color = c("Brain"="#5c89b2","Cerebellum"="#7bc7eb","Heart"="#8c1b1f","Liver"="#54833b","Kidney"="#c19e40",
                "Testis"="#e16527","Ovary"="#ac3495","Head"="#A6CEE3")
 # "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#fdfd99" "#e2cc1a"
point_shape=c("Gallus gallus"=5,
              "Homo sapiens"=0,
              "Macaca mulatta"=8,
              "Monodelphis domestica"=2,
              "Mus musculus"=1,
              "Oryctolagus cuniculus"=4,
              "Rattus norvegicus"=6,
              "Dendroctonus ponderosae"=9
              )

data_4 = read.delim("data/Data4_supp.tab",comment.char = "#")
data_4$species_name = str_replace_all(data_4$species,"_"," ")
data_4$longevity = data_1[data_4$species,]$longevity

data_4 = data_4[data_4$organs %in% names(point_color),]

life_span_order = tapply(data_4$longevity,data_4$species_name,mean)
life_span_order = life_span_order[order(life_span_order,decreasing = T)]

data_4$species_name = factor( data_4$species_name, levels = names(life_span_order))



p2 = ggplot(font.label = c(50, "plain"),font.legend= c(20, "plain"),font.x= c(20, "plain"),font.y= c(20, "plain"),
            data_4, aes(x=longevity, y=SVR,group=species_name,shape=species_name)) +   ggtitle("Major introns (BUSCO genes)")+ 
  geom_point(size=5,alpha=1,aes(color=organs),stroke=1.5)+ scale_color_manual("Organs",values=point_color) +
  scale_y_continuous(breaks=seq(0.5,4.5,0.5),labels=paste(seq(0.5,4.5,0.5),"%"),limits=c(0.5,4)) +
  scale_shape_manual("Species",values=point_shape) + labs(fill="Species") +
  xlab("Longevity (days, log scale)") + theme_bw() + theme(
    axis.title.x = element_text(color="black",margin = margin(t = 15, r = 0, b = 0, l = 0), size=31,family="serif"),
    axis.title.y = element_text(color="black",margin = margin(t = 0, r = 20, b = 0, l = 0), size=31, family="serif"),
    axis.text.y =  element_text(color="black", size=26, family="serif"),
    axis.text.x =  element_text(color="black", size=26, family="serif"),
    title =  element_text(color="black", size=26, family="serif"),
    text =  element_text(color="black", size=31, family="serif"),
    legend.text =  element_text(color="black", size=26, family="serif"),
    plot.caption = element_text(hjust = 0.365, face= "italic", size=23),
    plot.caption.position =  "plot"
  ) + labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  scale_x_log10(breaks=c(0.05,0.1,0.5,1,5,10,50,100,1000,10000,50000), limits=c(7,50000)) +
  guides(shape = guide_legend(override.aes = list(stroke=1.1), 
                              label.theme = element_text(color="black", size=26,face="italic", family="serif",vjust = 1.5,margin = margin(t = 5)))
         ,color = guide_legend(override.aes = list(shape = 16, size = 6), 
                               label.theme = element_text(color="black", size=26, family="serif",vjust = 1.5,margin = margin(t = 5))))+
  labs(
    caption = substitute(paste("LM model:"," R"^2,pgls_eq), list(pgls_eq=lm_eqn(lm((data_4$SVR)~log10(data_4$longevity)))))
  ) + annotation_logticks(sides="b")
p2

jpeg(paste(path_figure,"C_Moreira_svr.jpg",sep=""), width = 8500/resolution, height = 5100/resolution,res=600/resolution)
print(p2)
dev.off()


#### Figure 2

imgA = load.image(paste(path_require,"n1_n2_n3.png",sep=""))
imgB = load.image(paste(path_figure,"svr_busco_longevity.jpg",sep=""))
imgC = load.image(paste(path_figure,"C_Moreira_svr.jpg",sep=""))

{
  pdf(file=paste(path_pannel,"Figure2.pdf",sep=""), width=6.75, height=4)
  
  m=matrix(rep(1,10*11), nrow=11)
  
  for(i in 5:11){
    m[i,]=c(rep(2,5),rep(3,5))
  }
  m
  layout(m)
  
  par(mar=c(0, 0, 1, 0))
  plot(imgA, axes = F)
  mtext("A", side=2,at=100,adj=-12, line=1, font=2, cex=1.2,las=2)
  par(mar=c(0, 3, 0, 1))
  plot(imgB, axes = F)
  mtext("B",at=50,adj=1, side=2, line=1, font=2, cex=1.2,las=2)
  par(mar=c(0, 0, 0, 0))
  plot(imgC, axes = F)
  mtext("C", side=2,at=60,adj=1, line=1, font=2, cex=1.2,las=2)
  dev.off()
}
