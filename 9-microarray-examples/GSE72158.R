## GPL10558 Illumina HumanHT-12 V4.0 expression beadchip 
##human subcutaneous adipose tissue of 42 subjects at 2 different time points (0, 1 year), after bariatric surgery
rm(list=ls())
library(GEOquery)
library(limma) 
 GSE72158 <- getGEO('GSE72158', destdir=".",getGPL = F)
exprSet=exprs(GSE72158[[1]])
library("annotate")
GSE72158[[1]]
 
probeset <- featureNames(GSE72158[[1]])
library(Biobase)
library(GEOquery)
#Download GPL file, put it in the current directory, and load it:
gpl <- getGEO('GPL10558', destdir=".")
colnames(Table(gpl)) ## [1] 41108    17
head(Table(gpl)[,c(1,10,13)])  ## you need to check this , which column do you need 
probe2symbol=Table(gpl)[,c(1,13)]

symbol=probe2symbol[match(probeset,probe2symbol[,1]),2]
symbol= as.character(symbol)
exprSet=exprs(GSE72158[[1]])

a=cbind(symbol,exprSet)
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
exprSet[1:4,1:4]
exprSet=log(exprSet)
boxplot(exprSet,las=2)

pdata=pData(GSE72158[[1]])
## check the metadat and find the correct group information
# we must tell R that group should be interpreted as a factor !
individuals=factor(unlist(lapply(pdata$characteristics_ch1.3,function(x) strsplit(as.character(x),":")[[1]][2])))
treatment=factor(unlist(lapply(pdata$characteristics_ch1.2,function(x) strsplit(as.character(x),":")[[1]][2])))

## paired analysis :
design=model.matrix(~individuals+treatment)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
#vennDiagram(decideTests(fit))
pair_DEG=topTable(fit, coef="treatment 1 year",n=Inf,adjust='BH')

## non paired analysis :
design=model.matrix(~ treatment)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
#vennDiagram(decideTests(fit))
non_pair_DEG=topTable(fit, coef="treatment 1 year",n=Inf,adjust='BH')


t1=rownames(pair_DEG[pair_DEG$P.Value<0.05,])
t2=rownames(non_pair_DEG[non_pair_DEG$P.Value<0.05,])
length(setdiff(t1,t2))   ## 1051
length(setdiff(t2,t1))   ## 159 
length(intersect(t1,t2))  ##7707

write.csv( pair_DEG,"GSE72158-pair_DEG.DEG.csv")
write.csv(non_pair_DEG ,"GSE72158-non_pair_DEG.DEG.csv") 
