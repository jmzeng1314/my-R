#!/usr/bin/env Rscript

#############################################
#Author: Jianming Zeng
#email:  jmzeng1314@163.com
#Creat Time: 2017-10-12 11:32:45
#URL1:    http://www.bio-info-trainee.com/
#URL2:    https://github.com/jmzeng1314
### CAFS-->SUSTC-->LCRDC-->university of MACAU
#############################################

# Rscript run_DEG.R -e airway.expression.txt -g airway.group.txt -c 'trt-untrt' -s counts -m DESeq2
# Rscript run_DEG.R -e airway.expression.txt -g airway.group.txt -c 'trt-untrt' -s counts -m edgeR
# Rscript run_DEG.R -e sCLLex.expression.txt -g sCLLex.group.txt -c 'progres.-stable'
# Rscript run_DEG.R -e sCLLex.expression.txt -g sCLLex.group.txt -c 'progres.-stable' -m t.test

library("optparse")
 
option_list = list(
	make_option(c("-e", "--expression"), type="character", default=NULL, 
              help="expression matrix file name", metavar="character"),
	make_option(c("-g", "--group"), type="character", default=NULL, 
              help="group information file name", metavar="character"),
	make_option(c("-c", "--contrast"), type="character", default=NULL, 
              help="How to compare them,such as case-control", metavar="character"),
	make_option(c("-s", "--style"), type="character", default='signals', 
              help="What kind of data,read counts or signals?[default = signals]", metavar="character"),
	make_option(c("-m", "--method"), type="character", default="limma", 
              help="which method to use for the DEG[default = limma].For signals, you can choose limma or t.test.For read counts , you can choose DESeq2 or edgeR", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$expression) | is.null(opt$group) ){
  print_help(opt_parser)
  stop("expression matrix file and group information file must be supplied", call.=FALSE)
}
if (is.null(opt$contrast)  ){
  print_help(opt_parser)
  stop("The comparison information must be supplied", call.=FALSE)
}

## program...


for (pkg in c("limma","DESeq2","edgeR","pi0",'UpSetR')){
  if (! require(pkg,character.only=T) ) {
    #source("https://bioconductor.org/biocLite.R")
    BiocInstaller::biocLite(pkg,ask = F,suppressUpdates = T)
    require(pkg,character.only=T) 
    }
}
 

exprSet=read.table(opt$expression,stringsAsFactors = F,header = T)
group_info=read.table(opt$group,stringsAsFactors = F,header = T)

head(group_info)
group_list=group_info[,1]

## please make sure that the first column is the group information 
## Also please make sure that the order of exprSet and group are identical.


run_limma <- function(exprSet,group_list,contrast){
	suppressPackageStartupMessages(library(limma))
	design=model.matrix(~0+factor(group_list))
	colnames(design)=levels(factor(group_list))
	fit=lmFit(exprSet,design)
	cont.matrix=makeContrasts(contrasts=contrast ,levels = design)
	fit2=contrasts.fit(fit,cont.matrix)
	fit2=eBayes(fit2)
	results=topTable(fit2,adjust='BH',n=Inf)
	write.table(results,paste0("limma_DEG.",contrast,".txt"),quote = F,sep = '\t',row.names = T)

}

run_t.test <- function(exprSet,group_list,contrast){

	keep = group_list   %in%  strsplit(contrast,'-')[[1]]
	dat=exprSet[,keep]
	group_list=group_list[keep]

	group1 = which(group_list == strsplit(contrast,'-')[[1]][1] )
	group2 = which(group_list == strsplit(contrast,'-')[[1]][2] )

	dat1 = dat[, group1]
	dat2 = dat[, group2]
	dat = cbind(dat1, dat2)   
	library(pi0)
	pvals = matrix.t.test(dat, 1, 
		   length(group1), length(group2))
	p.adj = p.adjust(pvals, method = "BH")
	avg_1 = rowMeans(dat1)
	avg_2 = rowMeans(dat2)
	FC = avg_2/avg_1
	results = cbind(avg_1, avg_2, FC, pvals, p.adj)
	colnames(results) = c("avg_1", "avg_2", "logFC", "P.Value", "adj.P.Val") 
	DEG_t_test=results[order(results[,4]),]
	write.table(DEG_t_test,paste0("t_test_DEG.",contrast,".txt"),quote = F,sep = '\t',row.names = T)

}

run_DESeq2 <- function(exprSet,group_list,contrast){
	suppressPackageStartupMessages(library(DESeq2))

	geneLists=row.names(exprSet)
	keepGene=rowSums(edgeR::cpm(exprSet)>0) >=2
	table(keepGene);dim(exprSet)
	dim(exprSet[keepGene,])
	exprSet=exprSet[keepGene,]
	rownames(exprSet)=geneLists[keepGene]

	(colData <- data.frame(row.names=colnames(exprSet), group_list=group_list) )
	dds <- DESeqDataSetFromMatrix(countData = exprSet,
								  colData = colData,
								  design = ~ group_list)
	dds <- DESeq(dds) 
	plotDispEsts(dds, main="Dispersion plot") 

	rld <- rlogTransformation(dds)
	exprMatrix_rlog=assay(rld) 

	normalizedCounts1 <- t( t(counts(dds)) / sizeFactors(dds) )
	exprMatrix_rpm=as.data.frame(normalizedCounts1) 
	head(exprMatrix_rpm) 

	#res <- results(dds, contrast=c("group_list","trt","untrt"))
	res <- results(dds, contrast=c("group_list", strsplit(contrast,'-')[[1]] ))
	resOrdered <- res[order(res$padj),]
	head(resOrdered)
	DESeq2_DEG=as.data.frame(resOrdered)

	write.table(exprMatrix_rpm,file='DESeq2.exprMatrix_rpm.txt',quote = F,sep = '\t',row.names = T)
	write.table(exprMatrix_rlog,file='DESeq2.exprMatrix_rlog.txt',quote = F,sep = '\t',row.names = T)
	write.table(DESeq2_DEG,file=paste0("DESeq2_DEG.",contrast,".txt"),quote = F,sep = '\t',row.names = T)

}


run_edgeR <- function(exprSet,group_list,contrast){
		suppressPackageStartupMessages(library(edgeR)) 


		keep = group_list   %in%  strsplit(contrast,'-')[[1]]
		exprSet=exprSet[,keep]
		group_list=group_list[keep]

		d <- DGEList(counts=exprSet,group=factor(group_list))
		#d$samples$lib.size <- colSums(d$counts)
		d <- calcNormFactors(d)
		d$samples 
		dge=d

		design <- model.matrix(~factor(group_list))
		rownames(design)<-colnames(dge)
		colnames(design)<-levels(factor(group_list))
		dge <- estimateDisp(dge,design)

		#To perform quasi-likelihood F-tests:
		fit <- glmQLFit(dge,design)
		qlf <- glmQLFTest(fit,coef=2)
		# To perform likelihood ratio tests:
		fit <- glmFit(dge,design)
		lrt <- glmLRT(fit,coef=2)

		nrDEG=topTags(lrt, n=nrow(exprSet))
		nrDEG=as.data.frame(nrDEG)
		head(nrDEG)
		edgeR_lrt_DEG=nrDEG

		nrDEG=topTags(qlf, n=nrow(exprSet))
		nrDEG=as.data.frame(nrDEG)
		head(nrDEG)
		edgeR_qlf_DEG=nrDEG

		write.table(edgeR_lrt_DEG,file=paste0("edgeR_lrt.",contrast,".txt"),quote = F,sep = '\t',row.names = T)
		write.table(edgeR_qlf_DEG,file=paste0("edgeR_qlf_DEG.",contrast,".txt"),quote = F,sep = '\t',row.names = T)

		
}




if (opt$style == "counts") {

	if (opt$method == "DESeq2") {
	    
			suppressPackageStartupMessages(library(DESeq2))
			run_DESeq2(exprSet,group_list,opt$contrast)
	    }

	if (opt$method == "edgeR") {
	
		
		suppressPackageStartupMessages(library(edgeR)) 
		run_edgeR(exprSet,group_list,opt$contrast)

	}

}

if (opt$style == "signals") {

	if (opt$method == "t.test") {
		suppressPackageStartupMessages(library(pi0)) 
		run_t.test(exprSet,group_list,opt$contrast)
		}

	if (opt$method == "limma") {
		suppressPackageStartupMessages(library(limma)) 
		run_limma(exprSet,group_list,opt$contrast)
		
		
	}

}
