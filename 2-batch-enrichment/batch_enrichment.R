##http://www.genome.jp/kegg-bin/get_htext?hsa00001+3101
##download files :hsa00001.keg ,about 2.7M (2016/02/16)
##perl -lane '{if(/^C/){$t=$F[1]};if(/^D/){print "$t\t$F[1]"}}' hsa00001.keg >kegg2geneID.txt
##perl -lane '{if(/^A<b>(.*?)<\/b>/){$a=$1};if(/^B.*?<b>(.*?)<\/b>/){$b=$1};if(/^C\s+\d+\s+(.*?)\s+\[/){print "$a\t$b\t$F[1]\t$1"};}' hsa00001.keg >kegg_hierarchical.txt

annotationPKG='org.Hs.eg.db'
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(GOstats))
#suppressMessages(library(XML))
tmp=read.table("kegg2geneID.txt",sep="\t")
#tmp=toTable(org.Hs.egPATH)
# first column is kegg ID, second column is entrez ID 
GeneID2kegg<<- tapply(tmp[,1],as.factor(tmp[,2]),function(x) x)
kegg2GeneID<<- tapply(tmp[,2],as.factor(tmp[,1]),function(x) x)
kegg2name=read.delim("kegg_hierarchical.txt",header=F,sep="\t",stringsAsFactors =F)
colnames(kegg2name)=c('parent1','parent2','pathway_id','pathway_name')
kegg2name$pathway_id=as.numeric(kegg2name$pathway_id)
rownames(kegg2name)=kegg2name$pathway_id

suppressMessages(library(GO.db))
suppressMessages(library(KEGG.db))
hyperGtest <- function(GeneID2Path,Path2GeneID,diff_gene){
  diff_gene_has_path=intersect(diff_gene,names(GeneID2Path))
  n=length(diff_gene_has_path) #306
  N=length(GeneID2Path) #5870
  options(digits = 4)
  results=c()
  for (i in names(Path2GeneID)){
    M=length(Path2GeneID[[i]]) 
    exp_count=n*M/N
    k=0
    for (j in diff_gene_has_path){
      if (i %in% GeneID2Path[[j]]) k=k+1
    }
    OddsRatio=k/exp_count
    if (k==0) next
    p=phyper(k-1,M, N-M, n, lower.tail=F)
    print(paste(i,p,OddsRatio,exp_count,k,M,sep="    "))
    results=rbind(results,c(i,p,OddsRatio,exp_count,k,M))
  }
  colnames(results)=c("PathwayID","Pvalue","OddsRatio","ExpCount","Count","Size")
  return(results)
}

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
	kegg_result=hyperGtest(GeneID2kegg,kegg2GeneID,diff_gene_list)
	kegg_result=as.data.frame(kegg_result)
	kegg_result$pathway_name=kegg2name[match(as.numeric(as.character(kegg_result[,1])),kegg2name[,'pathway_id']),'pathway_name']
	kegg_result=kegg_result[order(kegg_result[,2]),]
	write.csv(kegg_result,paste(outdir,"/update_kegg.enrichment.csv",sep=""),row.names = F)

}