## GPL1352// U133_X3P//GSE43837
#exprSet_url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE43nnn/GSE43837/matrix/GSE43837_series_matrix.txt.gz"
#a=read.table(exprSet_url,sep = '\t',header = T,comment.char = "!")
library(GEOquery)
GSE43837 <- getGEO('GSE43837', destdir=".",getGPL = F)
## geneNames/sampleNames/pData/exprs
exprSet=exprs(GSE43837[[1]])

library("annotate") 
platformDB='u133x3p.db'
library(platformDB, character.only=TRUE)
probeset <- featureNames(GSE43837[[1]])
#EGID <- as.numeric(lookUp(probeset, platformDB, "ENTREZID"))
SYMBOL <-  lookUp(probeset, platformDB, "SYMBOL")
a=cbind(SYMBOL,exprSet)
## remove the duplicated probeset 
rmDupID <-function(a=matrix(c(1,1:5,2,2:6,2,3:7),ncol=6)){
  exprSet=a[,-1]
  rowMeans=apply(exprSet,1,function(x) mean(as.numeric(x),na.rm=T))
  a=a[order(rowMeans,decreasing=T),]
  exprSet=a[!duplicated(a[,1]),]
  #
  exprSet=exprSet[!is.na(exprSet[,1]),]
  rownames(exprSet)=exprSet[,1]
  exprSet=exprSet[,-1]
  return(exprSet)
}
exprSet=rmDupID(a)
rn=rownames(exprSet)
exprSet=apply(exprSet,2,as.numeric)
rownames(exprSet)=rn
exprSet=log(exprSet)

pdata=pData(GSE43837[[1]])
group_list=unlist(lapply(as.character(pdata$title),function(x) strsplit(x,"_")[[1]][2]))
library(limma) 
design=model.matrix(~factor(group_list))
fit=lmFit(exprSet,design)
fit=eBayes(fit)
#options(digits = 4)
DEG_limma=topTable(fit,coef=2,adjust='BH',n=Inf) 


