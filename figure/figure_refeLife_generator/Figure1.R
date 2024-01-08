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

############## eLife referee Pannel 1 A

ylabel = "mean_as_busco_low_as"
xlabel="socialite"

p1A =ggplot(  data_1,aes(data_1[,xlabel],data_1[,ylabel], fill=clade) ) + 
  geom_boxplot( aes(group=data_1[,xlabel] ))+ 
  geom_point(shape=21,size=7,alpha=0.7)+
  scale_fill_manual(values=vectorColor)+ ggtitle("Low-AS major introns (BUSCO genes)")+
  theme_bw() +labs(y=expression(paste("Average AS rate ",italic("per")," intron")))+
  xlab("Longevity (days, log scale)") +
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
  )
p1A

jpeg(paste(path_figure,"refeLife_p1A_.jpg",sep=""), width = 6100/resolution, height = 5500/resolution,res=700/resolution)
print(p1A)
dev.off()