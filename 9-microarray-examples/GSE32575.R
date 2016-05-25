##  obese GPL6102 Illumina human-6 v2.0 expression beadchip 
##  Transcriptome analysis of circulating monocytes in obese patients before and three months after bariatric surgery
##   
rm(list=ls())
library(GEOquery)
library(limma) 
GSE32575 <- getGEO('GSE32575', destdir=".",getGPL = F)
exprSet=exprs(GSE32575[[1]])
library("annotate")
GSE32575[[1]]

#platformDB='hgu133plus2.db'
#library(platformDB, character.only=TRUE)
probeset <- featureNames(GSE32575[[1]])
library(Biobase)
library(GEOquery)
#Download GPL file, put it in the current directory, and load it:
gpl <- getGEO('GPL6102', destdir=".")
colnames(Table(gpl)) ## [1] 41108    17
head(Table(gpl)[,c(1,10,13)])  ## you need to check this , which column do you need 
probe2symbol=Table(gpl)[,c(1,13)]

symbol=probe2symbol[match(probeset,probe2symbol[,1]),2]
symbol= as.character(symbol)
exprSet=exprs(GSE32575[[1]])
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
exprSet=log(exprSet) ## based on e 
boxplot(exprSet,las=2)


pdata=pData(GSE32575[[1]])
# we only want to compare the obese patients before and three months after bariatric surgery
exprSet=exprSet[,-c(1:12)]
pair <- factor(rep(1:18, each=2))
treatment=factor(rep(c('before','after'), times = 18))
treatment=relevel(treatment,'before')

## paired analysis :
design=model.matrix(~pair+treatment)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
#vennDiagram(decideTests(fit))
pair_DEG=topTable(fit, coef="treatmentafter",n=Inf,adjust='BH')

## non paired analysis :
design=model.matrix(~ treatment)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
#vennDiagram(decideTests(fit))
non_pair_DEG=topTable(fit, coef="treatmentafter",n=Inf,adjust='BH')


t1=rownames(pair_DEG[pair_DEG$P.Value<0.05,])
t2=rownames(non_pair_DEG[non_pair_DEG$P.Value<0.05,])
length(setdiff(t1,t2))   ## 960
length(setdiff(t2,t1))   ## 42 
length(intersect(t1,t2))  ##565

write.csv( pair_DEG,"GSE32575-pair_DEG.DEG.csv")
write.csv(non_pair_DEG ,"GSE32575-non_pair_DEG.DEG.csv") 
