suppressPackageStartupMessages(library(airway))
suppressPackageStartupMessages(library(DESeq))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR)) 
suppressPackageStartupMessages(library(pasilla))
suppressPackageStartupMessages(library(pasilla))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(CLL))
data(sCLLex)
data(pasillaGenes) 
data(airway)

write.table(exprs(sCLLex),file='sCLLex.expression.txt',quote = F,sep = '\t',row.names = T)
exprSet=read.table('sCLLex.expression.txt',stringsAsFactors = F,header = T)
write.table(pData(sCLLex)[,c(2,1)],file='sCLLex.group.txt',quote = F,sep = '\t',row.names = T)
group=read.table('sCLLex.group.txt',stringsAsFactors = F,header = T)


write.table(counts(pasillaGenes),file='pasillaGenes.expression.txt',quote = F,sep = '\t',row.names = T)
exprSet=read.table('pasillaGenes.expression.txt',stringsAsFactors = F,header = T)
write.table(pData(pasillaGenes)[,c(2,3)],file='pasillaGenes.group.txt',quote = F,sep = '\t',row.names = T)
group=read.table('pasillaGenes.group.txt',stringsAsFactors = F,header = T)
 
write.table(assays(airway)$counts,file='airway.expression.txt',quote = F,sep = '\t',row.names = T)
exprSet=read.table('airway.expression.txt',stringsAsFactors = F,header = T)
write.table(colData(airway) [,c(3,1,2)],file='airway.group.txt',quote = F,sep = '\t',row.names = T)
group=read.table('airway.group.txt',stringsAsFactors = F,header = T)






