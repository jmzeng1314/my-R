suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(reshape2))
suppressMessages(library(annotate))

gpl_info=read.csv("GPL_info.csv",stringsAsFactors = F)
gpl_info=gpl_info[grepl("Mus|Rattus|Homo",gpl_info[,3]),]

### first download all of the annotation packages from bioconductor
for (i in 1:nrow(gpl_info)){
  print(i)
  platform=gpl_info[i,4]
  platform=gsub('^ ',"",platform) ##主要是因为我处理包的字符串前面有空格
  #platformDB='hgu95av2.db'
  platformDB=paste(platform,".db",sep="")
  if( platformDB  %in% rownames(installed.packages()) == FALSE) {
    BiocInstaller::biocLite(platformDB)
    #source("http://bioconductor.org/biocLite.R");
    #biocLite(platformDB )
  } 
}
#下载完了所有的包， 就可以进行批量导出芯片探针与gene的对应关系！
for (i in 1:nrow(gpl_info)){
  print(i)
  platform=gpl_info[i,4]
  platform=gsub('^ ',"",platform)
  #platformDB='hgu95av2.db'
  platformDB=paste(platform,".db",sep="")

  if( platformDB  %in% rownames(installed.packages()) != FALSE) {
    library(platformDB,character.only = T)
    #tmp=paste('head(mappedkeys(',platform,'ENTREZID))',sep='')
    #eval(parse(text = tmp))
    ###重点在这里，把字符串当做命令运行
    all_probe=eval(parse(text = paste('mappedkeys(',platform,'ENTREZID)',sep='')))
    EGID <- as.numeric(lookUp(all_probe, platformDB, "ENTREZID"))
    ##自己把内容写出来即可
  } 
}
