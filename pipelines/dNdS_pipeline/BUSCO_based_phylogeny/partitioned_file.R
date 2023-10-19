options(stringsAsFactors = F, scipen = 999)
library(seqinr)

args = (commandArgs(TRUE))
infoPath = args[1]
aaspart_output = args[2]
cdspart_output = args[3]
aaspart_cons_output = args[4]






table = read.delim(infoPath)

df = data.frame(model="GTR+G",
              gene=paste(table[1,"gene.id"],"=",sep=""),
              tss=1,
              tes=table[1,"gene.size"])
for (i in 2:(nrow(table)-1)){print(i)
  df=rbind(df,data.frame(model="GTR+G",
                         gene=paste(table[i,"gene.id"],"=",sep=""),
                         tss=df[nrow(df),"tes"]+1,
                         tes=df[nrow(df),"tes"]+table[i,"gene.size"]
  ))
}

df$name=paste(df$model,", ",df$gene,df$tss,'-',df$tes,sep="")


write.table(df$name,cdspart_output, row.names=F, col.names=F, sep="\t", quote=F)


## AAS
table=read.delim(infoPath)

table$gene.size/3

df=data.frame(model="LG+G8+F",
              gene=paste(table[1,"gene.id"],"=",sep=""),
              tss=1,
              tes=table[1,"gene.size"]/3)
for (i in 2:(nrow(table)-1)){print(i)
  df=rbind(df,data.frame(model="LG+G8+F",
                         gene=paste(table[i,"gene.id"],"=",sep=""),
                         tss=df[nrow(df),"tes"]+1,
                         tes=df[nrow(df),"tes"]+table[i,"gene.size"]/3
  ))
}

df$name=paste(df$model,", ",df$gene,df$tss,'-',df$tes,sep="")


write.table(df$name,aaspart_output, row.names=F, col.names=F, sep="\t", quote=F)


## AAS
table = read.delim(infoPath)

table$gene_cons.size / 3

df=data.frame(model="LG+G8+F",
              gene=paste(table[1,"gene.id"],"=",sep=""),
              tss=1,
              tes=table[1,"gene_cons.size"]/3)
for (i in 2:(nrow(table)-1)){print(i)
  df=rbind(df,data.frame(model="LG+G8+F",
                         gene=paste(table[i,"gene.id"],"=",sep=""),
                         tss=df[nrow(df),"tes"]+1,
                         tes=df[nrow(df),"tes"]+table[i,"gene_cons.size"]/3
  ))
}

df$name=paste(df$model,", ",df$gene,df$tss,'-',df$tes,sep="")


write.table(df$name,aaspart_cons_output, row.names=F, col.names=F, sep="\t", quote=F)

