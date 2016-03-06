# my-R
#This is a small collection for the R scripts 
##1-get-all-probeset-info/
    it will download all of the microarray plantform R bioconductor packages, and extract the probeset information and corresponding gene ID
    like hgu95av2,hgu133a,hgu133plus2 and so on
##2-batch-enrichment/
    it will get do batch enrichment for the interested gene lists stored as the files.
    I still need to modify it, don't just copy it.
##3-get-hg19-gene-mapping/
    it just get all of the interesting object in org.Hs.eg.db which relevant for all kinds of gene ID.
    like entrez ID, gene Symbol, gene name, ENSEMBL ID, and also locus
    also I will connect the gene with the KEGG pathway ID and the GO ontology term.
##4-ID-convert-using-shiny/
    This is a shiny APP which is able to convert the ID inputed by user to the other default ID
    We accept 2 types of ID, gene Symbol or entrez ID, the convert to  gene name, ENSEMBL ID,locus,KEGG pathway ID and the GO ontology term.
##5-wego-web-app/
    This is shiny APP to do enrichment analysis , which is a imitative of [WEGO](http://wego.genomics.org.cn/cgi-bin/wego/index.pl)
    it's far from perfect, I'm still working on it.
##6-enrichmentForUnsupportedOrganisms    
    We can use this scripts to do enrichment for the organisms which are not support by Bioconductor package
    one of them use the GOstats package to do enrichment ,the other use the function writen by me.
##7-enrichment-with-newest-kegg/
    most of the kegg database are to old , for example:org.Hs.egPATH has 5869 entrez genes and 229 pathways
    On Augest 2015, there are 6901 entrez genes and 295 pathways
    And now , there are 299 pathways and 6992 genes.
    Address the sustainable growth of kegg database, it's essential to update it.
    I have a blog to introduct how to update the kegg information by Chinese [link](www.bio-info-trainee.com/1188.html)
    
##8-DEG
    there will be a detailed instruction in the folder, please go ahead .
    just for differential expression analysis !!!
