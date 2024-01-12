source("figure/figure_refeLife_generator/library_path.R")

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <-
    lapply(sheets, function(X)
      readxl::read_excel(filename, sheet = X))
  if (!tibble)
    x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}


mysheets <- read_excel_allsheets(paste("data/metazoa_69species.xls",sep=""))
sp_studied = names(mysheets)

data_1 = data_1[data_1$clade == "Hymenoptera",]
data_1$socialite = sapply(data_1$species, function(x)  mysheets[[x]]$Socialite[1])
table(data_1$socialite)
table(data_1$species,data_1$socialite)

############## eLife referee Pannel 1 A
ylabel="dNdS_200k"
xlabel="socialite"

wilcox.test(data_1[data_1$socialite == "Solitary",ylabel],data_1[data_1$socialite == "Eusocial",ylabel])
t.test(data_1[data_1$socialite == "Solitary",ylabel],data_1[data_1$socialite == "Eusocial",ylabel])

stat.test <- data_1 %>%
  wilcox_test(dNdS_200k ~ socialite) %>%
  add_significance()
stat.test


bxp <- ggboxplot(data_1, x = "socialite", y = "dNdS_200k", fill = "#ba8e18")
stat.test <- stat.test %>% add_xy_position(x = "socialite")

p1A = bxp + stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01,
                               size = 10) + xlab("Degree of sociality")+theme_bw()+ theme(
  axis.title.x = element_text(color="black", size=31,margin = margin(t = 15, r = 0, b = 0, l = 0),family="serif"),
  axis.title.y = element_text(color="black", size=31,margin = margin(t = 0, r = 20, b = 0, l = 0), family="serif"),
  axis.text.y =  element_text(color="black", size=26, family="serif"),
  axis.text.x =  element_text(color="black", size=26, family="serif"),
  title =  element_text(color="black", size=26 ,family="serif"),
  text =  element_text(color="black", size=31, family="serif"),
  legend.text =  element_text(color="black", size=26, family="serif"),
  plot.caption = element_text(hjust = 0, face= "italic", size=23),
  plot.caption.position =  "plot"
) + geom_point(shape=21,size=5,alpha=0.7,fill="#ba8e18") + ylab("dN/dS") + ylim(0.09,0.12)

p1A


jpeg(paste(path_figure,"refeLife_p1A.jpg",sep=""), width = 4000/resolution, height = 4700/resolution,res=700/resolution)
print(p1A)
dev.off()



############## eLife referee Pannel 1 B
ylabel="mean_as_busco_high_as"
xlabel="socialite"
data_1$mean_as_busco_low_as
wilcox.test(data_1[data_1$socialite == "Solitary",ylabel],data_1[data_1$socialite == "Eusocial",ylabel])
t.test(data_1[data_1$socialite == "Solitary",ylabel],data_1[data_1$socialite == "Eusocial",ylabel])

stat.test <- data_1 %>%
  wilcox_test(mean_as_busco_high_as ~ socialite) %>%
  add_significance()
stat.test


bxp <- ggboxplot(data_1, x = "socialite", y = "mean_as_busco_high_as", fill = "#ba8e18")
stat.test <- stat.test %>% add_xy_position(x = "socialite")

p1B = bxp + stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01,
                               size = 10) + xlab("Degree of sociality")+theme_bw()+ggtitle("High-AS major-isoform introns\n(BUSCO genes)") + theme(
                                 axis.title.x = element_text(color="black", size=31,margin = margin(t = 15, r = 0, b = 0, l = 0),family="serif"),
                                 axis.title.y = element_text(color="black", size=31,margin = margin(t = 0, r = 20, b = 0, l = 0), family="serif"),
                                 axis.text.y =  element_text(color="black", size=26, family="serif"),
                                 axis.text.x =  element_text(color="black", size=26, family="serif"),
                                 title =  element_text(color="black", size=18 ,family="serif"),
                                 text =  element_text(color="black", size=31, family="serif"),
                                 legend.text =  element_text(color="black", size=26, family="serif"),
                                 plot.caption = element_text(hjust = 0, face= "italic", size=23),
                                 plot.caption.position =  "plot"
                               ) + geom_point(shape=21,size=5,alpha=0.7,fill="#ba8e18") +
  scale_y_continuous(breaks=seq(10,25,2), labels=paste(seq(10,25,2),"%"),limits = c(15,24)) + ylab("dN/dS")+labs(y=expression(paste("Average AS rate ",italic("per")," intron")))

p1B


jpeg(paste(path_figure,"refeLife_p1B.jpg",sep=""), width = 4000/resolution, height = 5500/resolution,res=700/resolution)
print(p1B)
dev.off()



############## eLife Figure 1
imgA = load.image(paste(path_figure,"refeLife_p1A.jpg",sep=""))
imgB = load.image(paste(path_figure,"refeLife_p1B.jpg",sep=""))

{
  pdf(file= paste(path_pannel,"Figure1.pdf",sep=""), width=4.5*1/1, height=2.75)
  
  m=matrix(rep(1,10*2), nrow=2)
  
  m[1,]=c(rep(1,5),rep(2,5))
  m[2,]=c(rep(1,5),rep(2,5))
  m
  layout(m)
  
  
  par(mar=c(0, 0.5, 2, 0.5))
  plot(imgA, axes=F)
  mtext("A",at=-20,adj=-.5, side=2, line=1, font=2, cex=1.4,las=2)
  
  par(mar=c(0, 0.5, 0.5, 0.5))
  plot(imgB, axes=F)
  mtext("B",at=80,adj=-.5, side=2, line=1, font=2, cex=1.4,las=2)
  dev.off()
}
