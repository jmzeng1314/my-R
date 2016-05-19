library(GEOquery)
library(limma) 
GSE62832 <- getGEO('GSE62832', destdir=".",getGPL = F)
exprSet=exprs(GSE62832[[1]])
library("annotate")
GSE62832[[1]]

platformDB='hugene10sttranscriptcluster.db'
library(platformDB, character.only=TRUE)
probeset <- featureNames(GSE62832[[1]])
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

pdata=pData(GSE62832[[1]])
## check the metadat and find the correct group information
# we must tell R that group should be interpreted as a factor !
individuals=factor(unlist(lapply(pdata$title,function(x) strsplit(as.character(x),"-")[[1]][1])))
treatment=factor(unlist(lapply(pdata$title,function(x) strsplit(as.character(x),"-")[[1]][2])))
metabolic_status=factor(unlist(lapply(pdata$characteristics_ch1.1,function(x) strsplit(as.character(x)," ")[[1]][3])))
treatment <- relevel(treatment, "Pre")

## if only take the treatment into accont 
design=model.matrix(~treatment)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
topTable(fit,coef=2,n=20)

## if only take the metabolic_status into accont 
design=model.matrix(~metabolic_status)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
topTable(fit,coef=2,n=20)

## if take both treatment and metabolic_status into account :
library(limma) 
#design=model.matrix(~treatment+metabolic_status+individuals)
design=model.matrix(~treatment+metabolic_status)
fit=lmFit(exprSet,design)
fit=eBayes(fit)
options(digits = 4)
vennDiagram(fit,include=c("up","down"))
topTable(fit,coef=2,adjust='BH')  ## equal with below code 
topTable(fit,coef="treatmentPost",n=20)
topTable(fit,coef="metabolic_statusMNO",n=20)
## In fact, there's little different between them.

write.csv(topTable(fit,coef="treatmentPost",n=Inf,adjust='BH'),"treatmentPostVsPre.DEG.csv")
write.csv(topTable(fit,coef="metabolic_statusMNO",n=Inf,adjust='BH'),"metabolic_statusMNOVsMAO.DEG.csv")

## http://genomicsclass.github.io/book/pages/expressing_design_formula.html
## http://kasperdanielhansen.github.io/genbioconductor/html/limma.html 

## we will produce all possible pairwise comparisons for treatment and metabolic_status column
## also the design will  take into account the paired samples (individuals column) 
## please read the "8.7 Multi-level experiments" chapter of limma user guide if intesesting .
## technical replicates VS experimental replicates.

#https://support.bioconductor.org/p/26850/ 
