## the enrichment function "hyperGTest" comes from R package :GOstats
annotationPKG='org.Hs.eg.db'
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(GOstats))
#suppressMessages(library(XML))
suppressMessages(library(GO.db))
suppressMessages(library(KEGG.db))

if (F){
	f=list.files("./")
	tmp=sapply(f,function(diff_file){
		a=read.table(diff_file,header=T,stringsAsFactors=F)
		#return(nrow(a[abs(a$logFC)>2 & a$fdr<0.05,]))
		a$fdr=p.adjust(a$pvals,"fdr")
		diff_result=a
		write.table(diff_result[order(diff_result[,5]),],diff_file,row.names=F,sep="\t",quote=F)	
	})
	#[1] 1140    5
}

for (i in 2:nrow(a)){
	diff_probe=a[i,2]
	diff_genename=a[i,1]
	outdir=paste("enrichment_result/",diff_genename,"_",diff_probe,sep="")
	dir.create(outdir, recursive =T)
	b=read.table(file.path("diff_result",paste(diff_probe,".diff",sep="")),header=T,stringsAsFactors=F)
	b$fdr=p.adjust(b$pvals,"fdr")
	### b is a table , which produce by limma DEG analysis !
	b=b[order(b[,5]),]
	if (nrow(b[b$fdr<0.05,]) >1000){
		tmp=b[b$fdr<0.05,]
		tmp=tmp[order(abs(tmp[,3]),decreasing=T),]
		diff_gene_list=tmp[1:1000,1]
	}else if (nrow(b[b$pvals<0.05,]) >1000){
		tmp=b[b$pvals<0.05,]
		tmp=tmp[order(abs(tmp[,3]),decreasing=T),]
		diff_gene_list=tmp[1:1000,1]
	}else {
		diff_gene_list=b[1:1000,1]
	}

	diff_gene_list=intersect(diff_gene_list,mappedkeys(org.Hs.egSYMBOL2EG))
	diff_gene_list=as.numeric(unlist(as.list(org.Hs.egSYMBOL2EG[diff_gene_list])))
	GOES = c('BP','CC', 'MF');
	for (ontology in GOES) {
	 GO.hyperG.params = new("GOHyperGParams", geneIds=diff_gene_list, universeGeneIds=NULL, annotation=annotationPKG, 
							ontology=ontology, pvalueCutoff=1, conditional = FALSE, testDirection = "over")
	 GO.hyperG.results = hyperGTest(GO.hyperG.params);
	 outHTMLname=paste(outdir,"/GO_",ontology,".enrichment.html",sep="")
	 htmlReport(GO.hyperG.results, file=outHTMLname, summary.args=list("htmlLinks"=TRUE))
	 summary(GO.hyperG.results)
	}	
	options(digits=4);	
	hyperG.params = new("KEGGHyperGParams", geneIds=diff_gene_list, universeGeneIds=NULL, annotation=annotationPKG, 
						categoryName="KEGG", pvalueCutoff=1, testDirection = "over")
	KEGG.hyperG.results = hyperGTest(hyperG.params); 
	outHTMLname=paste(outdir,"/kegg.enrichment.html",sep="")
	htmlReport(KEGG.hyperG.results, file=outHTMLname, summary.args=list("htmlLinks"=TRUE))
}