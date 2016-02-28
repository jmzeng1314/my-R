suppressMessages(library(GOstats))
suppressMessages(library(GSEABase))
suppressMessages(library(optparse))
option_list = list(
  make_option(c("-p", "--kegg2geneID"), action="store", help="The full path of kegg2geneID file eg:kegg2geneID.txt"),
  make_option(c("-g", "--GO2geneID"), action="store", help="The full path of GO2geneID file eg:GO2geneID.txt"),
  make_option(c("-u", "--universe_gene_list"), action="store", help="The full path of universe_gene_list file eg:universe_gene_list.txt"),
  make_option(c("-d", "--diff_gene_list"), action="store",help="The full path of diff_gene_list file eg:diff_gene_list.txt")
)
opt = parse_args(OptionParser(option_list=option_list))

## first check the params!!!
if(is.null(opt$kegg2geneID)){
  kegg2geneID_file='kegg2geneID.txt'
}else{
  kegg2geneID_file=opt$kegg2geneID
}
if(is.null(opt$GO2geneID)){
  GO2geneID_file='GO2geneID.txt'
}else{
  GO2geneID_file=opt$GO2geneID
}
if(is.null(opt$diff_gene_list)){
  diff_gene_list_file='diff_gene_list.txt'
}else{
  diff_gene_list_file=opt$diff_gene_list
}
if(is.null(opt$universe_gene_list)){
  universe_gene_list_file='universe_gene_list.txt'
}else{
  universe_gene_list_file=opt$universe_gene_list
}
## then read the files !!!
if (file.exists(kegg2geneID_file)){
	keggframeData=read.table(kegg2geneID_file,colClasses=c('character'))
	### two columns, first is keggID and second is geneID
	GeneID2kegg<<- tapply(tmp[,1],as.factor(tmp[,2]),function(x) x)
	kegg2GeneID<<- tapply(tmp[,2],as.factor(tmp[,1]),function(x) x)
}else{stop("we can not find the file:kegg2geneID_file")}
if (file.exists(GO2geneID_file)){
  goframeData=read.table(GO2geneID_file,header = F,stringsAsFactors = F)
  ### three columns, first is GO ID and second is evidence codes, last is gene ID
}else{stop("we can not find the file:GO2geneID_file")}
if (file.exists(diff_gene_list_file)){
  a=read.table(diff_gene_list_file)
  diff_gene_list<<- a[,1]
}else{stop("we can not find the file:diff_gene_list_file")}
if (file.exists(universe_gene_list_file)){
  a=read.table(universe_gene_list_file)
  universe=as.character(a[,1])
}else{stop("we can not find the file:universe_gene_list_file")}
library(KEGG.db)
library(GO.db)
all_go=mappedkeys(GOTERM)
go2name=data.frame(go_id=all_go,term=as.character(Term(all_go)))
keggID2name=toTable(KEGGPATHID2NAME)
kegg2name=keggID2name
##[1] "path_id"   "path_name"

#setwd("D:\\test_analysis\\my_github\\my-R\\6-enrichmentForUnsupportedOrganisms")
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

kegg_result=hyperGtest_jimmy(GeneID2kegg,kegg2GeneID,diff_gene_list)
kegg_result=as.data.frame(kegg_result)
kegg_result$pathway_name=kegg2name[match(kegg_result[,1],kegg2name[,'pathway_id']),'pathway_name']
kegg_result=kegg_result[order(kegg_result[,2]),]
write.csv(kegg_result,"update_kegg.enrichment.csv",row.names = F)


if (F){
  ## then do GO enrichment analysis
  #frame = toTable(org.Hs.egGO)
  #goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)
  goFrame=GOFrame(goframeData,organism="unknown")
  goAllFrame=GOAllFrame(goFrame)
  gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
  GOES = c('BP','CC', 'MF');
  for (ontology in GOES) {
    GO.hyperG.params = GSEAGOHyperGParams(name="My Custom GSEA based annot Params", 
                                          geneSetCollection=gsc,
                                          geneIds=diff_gene_list, 
                                          universeGeneIds=universe,
                                          ontology = ontology,
                                          pvalueCutoff=1, 
                                          conditional = FALSE,
                                          testDirection = "over"
    )
    GO.hyperG.results = hyperGTest(GO.hyperG.params);
    outHTMLname=paste("GO_",ontology,".enrichment.csv",sep="")
    htmlReport(GO.hyperG.results, file=outHTMLname, summary.args=list("htmlLinks"=TRUE))
    summary(GO.hyperG.results)
  }
}	
