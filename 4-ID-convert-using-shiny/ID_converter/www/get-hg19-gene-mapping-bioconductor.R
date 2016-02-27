library(RSQLite)
sqlite    <- dbDriver("SQLite")
con <- dbConnect(sqlite,"hg19_bioconductor.sqlite") # makes a new file
suppressMessages(library(org.Hs.eg.db))
kegg2ID=toTable(org.Hs.egPATH)
#[1] "gene_id" "path_id"
dbWriteTable(con,'keggID2geneID',kegg2ID,row.name=F,overwrite=T)
go2id=toTable(org.Hs.egGO2ALLEGS)
## gene_id      go_id Evidence Ontology
dbWriteTable(con,'goID2geneID',go2id,row.name=F,overwrite=T)
library(KEGG.db)
library(GO.db)
#ls("package:KEGG.db")
#ls("package:GO.db")
keggID2name=toTable(KEGGPATHID2NAME)
##[1] "path_id"   "path_name"
dbWriteTable(con,'keggID2name',keggID2name,row.name=F,overwrite=T)
all_go=mappedkeys(GOTERM)
go2name=data.frame(go_id=all_go,term=as.character(Term(all_go)))
dbWriteTable(con,'go2name',go2name,row.name=F,overwrite=T)

suppressMessages(library("org.Hs.eg.db"))
all_EG=mappedkeys(org.Hs.egSYMBOL)

tmp=unlist(as.list(org.Hs.egSYMBOL))
EG2Symbol=data.frame(EGID=names(tmp),symbol=as.character(tmp))

tmp=unlist(as.list(org.Hs.egENSEMBL))
EG2ENSEMBL=data.frame(EGID=names(tmp),ENSEMBL=as.character(tmp))

tmp=unlist(as.list(org.Hs.egGENENAME))
EG2name=data.frame(EGID=names(tmp),name=as.character(tmp))

tmp=unlist(as.list(org.Hs.egMAP))
EG2MAP=data.frame(EGID=names(tmp),MAP=as.character(tmp))


tmp=merge(EG2Symbol,EG2MAP,by='EGID',all=TRUE)
tmp=merge(tmp,EG2ENSEMBL,by='EGID',all=TRUE)
my_gene_mapping=merge(tmp,EG2name,by='EGID',all=TRUE)
##[1] "EGID"    "symbol"  "MAP"     "ENSEMBL" "name"
apply(my_gene_mapping,2,function(x) length(unique(x)))
dbWriteTable(con,'my_gene_mapping',my_gene_mapping,row.name=F,overwrite=T)
dbDisconnect(con)