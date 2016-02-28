##http://www.genome.jp/kegg-bin/get_htext?hsa00001+3101
##download files :hsa00001.keg ,about 2.7M (2016/02/16)
##perl -lane '{if(/^C/){$t=$F[1]};if(/^D/){print "$t\t$F[1]"}}' hsa00001.keg >kegg2geneID.txt
##perl -lane '{if(/^A<b>(.*?)<\/b>/){$a=$1};if(/^B.*?<b>(.*?)<\/b>/){$b=$1};if(/^C\s+\d+\s+(.*?)\s+\[/){print "$a\t$b\t$F[1]\t$1"};}' hsa00001.keg >kegg_hierarchical.txt
suppressMessages(library(optparse))
option_list = list(
	make_option(c("-p", "--path2gene"), action="store", help="The full path of path2gene file eg:kegg2geneID.txt"),
	make_option(c("-g", "--path2name"), action="store", help="The full path of path2name file eg:kegg_hierarchical.txt"),
	make_option(c("-d", "--diff_gene_list"), action="store",help="The full path of diff_gene_list file eg:diff_gene_list.txt")
)
opt = parse_args(OptionParser(option_list=option_list))

if(is.null(opt$path2gene)){
	path2gene_file='kegg2geneID.txt'
}else{
	path2gene_file=opt$path2gene
}
if(is.null(opt$path2name)){
	path2name_file='kegg_hierarchical.txt'
}else{
	path2name_file=opt$path2name
}
if(is.null(opt$diff_gene_list)){
	diff_gene_list_file='diff_gene_list.txt'
}else{
	diff_gene_list_file=opt$diff_gene_list
}

#setwd("D:\\test_analysis\\my_github\\my-R\\7-enrichment-with-newest-kegg")
hyperGtest_jimmy <- function(GeneID2Path,Path2GeneID,diff_gene){
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

if (file.exists(path2gene_file)){
	tmp=read.table(path2gene_file,sep="\t",colClasses=c('character'))
	#tmp=toTable(org.Hs.egPATH)
	# first column is kegg ID, second column is entrez ID 
	GeneID2kegg<<- tapply(tmp[,1],as.factor(tmp[,2]),function(x) x)
	kegg2GeneID<<- tapply(tmp[,2],as.factor(tmp[,1]),function(x) x)
}else{stop("we can not find the file:path2gene_file")}
if (file.exists(path2name_file)){
	kegg2name<<- read.delim(path2name_file,header=F,sep="\t",colClasses=c('character'),stringsAsFactors =F)
	colnames(kegg2name)=c('parent1','parent2','pathway_id','pathway_name')
	###kegg2name$pathway_id=as.numeric(kegg2name$pathway_id)
	rownames(kegg2name)=kegg2name$pathway_id
}else{stop("we can not find the file:path2name_file")}
if (file.exists(diff_gene_list_file)){
	a=read.table(diff_gene_list_file)
	diff_gene_list<<- a[,1]
}else{stop("we can not find the file:diff_gene_list_file")}

kegg_result=hyperGtest_jimmy(GeneID2kegg,kegg2GeneID,diff_gene_list)
kegg_result=as.data.frame(kegg_result)
kegg_result$pathway_name=kegg2name[match(kegg_result[,1],kegg2name[,'pathway_id']),'pathway_name']
kegg_result=kegg_result[order(kegg_result[,2]),]
write.csv(kegg_result,"update_kegg.enrichment.csv",row.names = F)