#suppressMessages(library(affy))
#suppressMessages(library(limma))
library(DESeq)
library(limma)
library(DESeq2)
library(edgeR)

suppressMessages(library(optparse))
txt_output=T ## whether need txt or not ?

option_list = list(
	make_option(c("-g", "--gct"),help="The full path of gct format file "),
	make_option(c("-c", "--cls"),help="The full path of cls format file "),
	make_option(c("-m", "--method"),default="limma",help="The full path of cls format file "),
)
opt = parse_args(OptionParser(option_list=option_list))
if(is.null(opt$gct) || is.null(opt$cls)){
   stop(sprintf("The gct and cls files is required !!!"))
}
cls_file=opt$gct
gct_file=opt$cls
DEG_method=tolower(opt$method)
TEST=T
if (TEST){
	cls_file="P53.cls"
	gct_file="P53_hgu95av2.gct"
	DEG_method='limma'
}
studyID=sub(".gct","",gct_file)
rmDupID <-function(a=matrix(c(1,1:5,2,2:6,2,3:7),ncol=6)){
  exprSet=a[,-1]
  rowMeans=apply(exprSet,1,mean)
  a=a[order(rowMeans,decreasing=T),]
  return(a[!duplicated(a[,1]),])
}

cat(paste("step1:",Sys.time(),"start to do DEG !!! \n"))
a=read.table(gct_file,stringsAsFactors=F,skip=2,header=T,sep="\t")
## probeID/geneID/value1/value2~~~~~~
exprSet=a[,-1]
exprSet=rmDupID(exprSet)
rownames(exprSet)=exprSet[,1]
exprSet=exprSet[,-1]
sample_list=colnames(exprSet)
b=read.table(cls_file,stringsAsFactors=F,skip=2)
group_list=as.character(b)
sample_group=data.frame(sample_list,group_list)
names(sample_group)=c('sample_id','compare')

##liver; males and femlaes human pubmed ID: 	20009012
phenodata=read.table("http://bowtie-bio.sourceforge.net/recount/phenotypeTables/gilad_phenodata.txt",header=T)
countdata=read.table("http://bowtie-bio.sourceforge.net/recount/countTables/gilad_count_table.txt",header = T,sep = "\t")
rmDupID <-function(a=matrix(c(1,1:5,2,2:6,2,3:7),ncol=6)){
  exprSet=a[,-1]
  rowMeans=apply(exprSet,1,mean)
  a=a[order(rowMeans,decreasing=T),]
  return(a[!duplicated(a[,1]),])
}
exprSet=rmDupID(countdata)
rownames(exprSet)=exprSet[,1]
exprSet=exprSet[,-1]
sample_list=colnames(exprSet)
group_list=as.character(phenodata[,3])

DEG_results=DEG_voom(exprSet,group_list)
DEG_results=DEG_DESeq2(exprSet,group_list)


if(DEG_method=='limma'){
	DEG_results=DEG_limma(exprSet,group_list)		
}else if (DEG_method=='voom'){
	DEG_results=DEG_voom(exprSet,group_list)
}else if (DEG_method=='deseq'){
	DEG_results=DEG_DESeq(exprSet,group_list)
}else if (DEG_method=='deseq2'){
	DEG_results=DEG_DESeq2(exprSet,group_list)
}else if (DEG_method=='edger'){
	DEG_results=DEG_edger(exprSet,group_list)
}else{
	stop(sprintf("we don't offer the method you want"))
}


DEG_edger <- function(exprSet=exprSet,group_list=group_list){
	suppressMessages(library(edgeR)	
	d <- DGEList(counts=exprSet,group=factor(group_list))
	d.full <- d # keep the old one in case we mess up
	apply(d$counts, 2, sum) # total gene counts per samplekeep <- rowSums(cpm(d)>100) >= 2
	d <- d[keep,]
	d$samples$lib.size <- colSums(d$counts)
	d <- calcNormFactors(d)
	png("MDS.png")
	plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
	legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)
	dev.off()
	d1 <- estimateCommonDisp(d, verbose=T)
	d1 <- estimateTagwiseDisp(d1)
	et12 <- exactTest(d1, pair=c(1,2)) 
	png("MA.png")
	de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
	de1tags12 <- rownames(d1)[as.logical(de1)] 
	plotSmear(et12, de.tags=de1tags12)
	abline(h = c(-2, 2), col = "blue")
	dev.off()
	topTags(et12, n=10)
}

DEG_limma <- function(exprSet=exprSet,group_list=group_list){
	suppressMessages(library(limma))
	design <- model.matrix(~0+factor(group_list))
	colnames(design)=levels(factor(group_list))
	rownames(design)=colnames(exprSet)
	fit <- lmFit(exprSet,design)
	png("MA.png")
	plotMA(fit)
	abline(0,0,col="blue")
	dev.off()
	fit2 <- eBayes(fit)
	tempOutput = topTable(fit2, coef=1, n=Inf)
	nrDEG2 = na.omit(tempOutput) 
}

DEG_voom <- function(exprSet=exprSet,group_list=group_list){
	suppressMessages(library(limma))
	design <- model.matrix(~0+factor(group_list))
	colnames(design)=levels(factor(group_list))
	rownames(design)=colnames(exprSet)
	v <- voom(exprSet,design,normalize="quantile")
	png("MDS.png")
	plotMDS(v, labels=1:ncol(exprSet),col=rainbow(ncol(exprSet)))
	#abline(0,0,col="blue")
	dev.off()
	fit <- lmFit(v,design)
	png("MA.png")
	plotMA(fit)
	#abline(0,0,col="blue")
	dev.off()
	fit2 <- eBayes(fit)
	tempOutput = topTable(fit2, coef=1, n=Inf)
	nrDEG2 = na.omit(tempOutput)	
}

DEG_DESeq2 <- function(exprSet=exprSet,group_list=group_list){
	suppressMessages(library(DESeq2))
	exprSet=ceiling(exprSet)
	(colData <- data.frame(row.names=colnames(exprSet), group_list=group_list))
	dds <- DESeqDataSetFromMatrix(countData = exprSet,
								colData = colData,
								design = ~ group_list)
	dds <- DESeq(dds)
	res <- results(dds)
	png("MA.png")
	#plotMA(res, main="DESeq2", ylim=c(-2,2))
	dev.off()
	resOrdered <- res[order(res$padj),]
	resOrdered=as.data.frame(resOrdered)
}

DEG_DESeq <- function(exprSet=exprSet,group_list=group_list){
	suppressMessages(library(DESeq))
	exprSet=ceiling(exprSet)
	de=newCountDataSet(exprSet,group_list)
	de=estimateSizeFactors(de)
	de=estimateDispersions(de)
	res=nbinomTest(de,"case","control")
	rownames(res)=res[,1]
	res=res[,-1]
}
