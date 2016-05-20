getwd()
library(ArrayExpress)
## getAE download all raw data for this study
## ArrayExpress just donwload the expressionSet 
rawset = ArrayExpress("E-MEXP-1422")    
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MEXP-1422/ 
## it's a AffyBatch object  
library(affy)
AEsetnorm = rma(rawset)
exprSet=exprs(AEsetnorm)
dim(exprSet)
#[1] 22277     6
length(featureNames(rawset))
#[1] 22277

#Investigation description  	E-MEXP-1422.idf.txt
#Sample and data relationship	E-MEXP-1422.sdrf.txt
#Raw data (1)               	E-MEXP-1422.raw.1.zip
#Processed data (1)         	E-MEXP-1422.processed.1.zip
#Array design	                A-AFFY-37.adf.txt
#R ExpressionSet	            E-MEXP-1422.eSet.r

library("annotate") 
platformDB='hgu133plus2.db'
library(platformDB, character.only=TRUE)
probeset <- featureNames(rawset)
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
dim(exprSet)
#[1] 12500     6
pdata=pData(rawset)
group_list=factor(c("PROX1 siRNA 1","PROX1 siRNA 1","GFP siRNA 1","GFP siRNA 1","GFP siRNA 2","GFP siRNA 2"))

#exprSet=log(exprSet)
#boxplot(exprSet,las=2)
library(limma) 
design=model.matrix(~factor(group_list))
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
vennDiagram(fit,include="up")
vennDiagram(fit,include="down")

topTable(fit,coef=2,n=20)
topTable(fit,coef=3,n=20)

write.csv(topTable(fit,coef=2,n=Inf,adjust='BH'),"GFP_siRNA_2_vS_1.DEG.csv")


