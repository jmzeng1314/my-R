## HuGene-1_0-st // GPL6244 // obesity
library(GEOquery)
library(limma) 
GSE70529 <- getGEO('GSE70529', destdir=".",getGPL = F)
exprSet=exprs(GSE70529[[1]])
library("annotate")
GSE70529[[1]]

platformDB='hugene10sttranscriptcluster.db'
library(platformDB, character.only=TRUE)
probeset <- featureNames(GSE70529[[1]])
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

#exprSet=log(exprSet)
#boxplot(exprSet,las=2)

pdata=pData(GSE70529[[1]])
## check the metadat and find the correct group information
# we must tell R that group should be interpreted as a factor !
individuals=factor(unlist(lapply(pdata$title,function(x) strsplit(as.character(x),"-")[[1]][1])))
treatment=unlist(lapply(pdata$title,function(x) strsplit(as.character(x),"-")[[1]][2]))
treatment=factor(sub('2','',treatment))

## if only take the treatment into accont 
design=model.matrix(~treatment)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
#vennDiagram(fit,include=c("up","down"))
write.csv(topTable(fit,coef='treatmentB',n=Inf,adjust='BH'),"treatmentBvsA.DEG.csv")
write.csv(topTable(fit,coef='treatmentC',n=Inf,adjust='BH'),"treatmentCvsA.DEG.csv")
write.csv(topTable(fit,coef='treatmentD',n=Inf,adjust='BH'),"treatmentDvsA.DEG.csv")
## and we don't care about the difference between individuals

