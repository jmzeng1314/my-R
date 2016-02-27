library(shiny)
suppressMessages(library(RSQLite))
trimws <- function(x=' trimws '){
  x=sub("^\\s","",x)
  x=sub("\\s$","",x)
  x
}
sqlite    <- dbDriver("SQLite")
#con <- dbConnect(sqlite,"www/hg19_bioconductor.sqlite")
#dbDisconnect(con)
# Define server logic for random distribution application

shinyServer(function(input, output) {
  gene_info <- reactive({
    con <- dbConnect(sqlite,"www/hg19_bioconductor.sqlite")
    ID=input$id
    ID=trimws(ID)
    ##[1] "EGID"    "symbol"  "MAP"     "ENSEMBL"  "name" 
    query=paste0(
      "select * from my_gene_mapping where EGID=",shQuote(ID) ,
      " or symbol=upper(",shQuote(ID),")"
    )
    tmp=dbGetQuery(con,query)
    dbDisconnect(con)
    tmp
  })
  output$Entrez <- renderUI({
    tmp=gene_info()
    Entrez<<-tmp[1,'EGID']
    link=paste0("http://www.ncbi.nlm.nih.gov/gene/",Entrez)
    tags$a(href=link,Entrez)
  })
  #http://www.ncbi.nlm.nih.gov/gene?term=BRCA1
  output$Symbol <- renderUI({
    tmp=gene_info()
    Symbol<<-tmp[1,'symbol']
    link=paste0("http://www.ncbi.nlm.nih.gov/gene?term=",Symbol)
    tags$a(href=link,Symbol)
  })
  output$Locus      <- renderText({
    tmp=gene_info()
    Locus<<-tmp[1,'MAP']
  })
  #http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000141510
  output$ENSEMBL <- renderUI({
    tmp=gene_info()
    ENSEMBL<<-tmp[1,'ENSEMBL']
    link=paste0("http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",ENSEMBL)
    tags$a(href=link,ENSEMBL)  
  })
  output$Name      <- renderText({
    tmp=gene_info()
    Name<<-tmp[1,'name']   
  })
  output$kegg     <- renderTable({
    tmp=gene_info()
    Entrez<<-tmp[1,'EGID']
    con <- dbConnect(sqlite,"www/hg19_bioconductor.sqlite")
    query=paste0("select b.path_id,b.path_name from keggID2geneID a", 
                 "  left join keggID2name b on (a.path_id=b.path_id)", 
                 "  where a.gene_id=",shQuote(Entrez)
                 )
    tmp=dbGetQuery(con,query)
    dbDisconnect(con)
    tmp
  })
  output$go_bp     <- renderTable({
    tmp=gene_info()
    Entrez<<-tmp[1,'EGID']
    con <- dbConnect(sqlite,"www/hg19_bioconductor.sqlite")
    query=paste0("select a.go_id,a.Evidence,a.Ontology,b.term from goID2geneID a ",
            "  left join go2name b on (a.go_id=b.go_id)", 
            "  where a.Ontology='BP'  and a.gene_id=",shQuote(Entrez)
            )
    tmp=dbGetQuery(con,query)
    dbDisconnect(con)
    tmp
  })
  output$go_cc     <- renderTable({
    tmp=gene_info()
    Entrez<<-tmp[1,'EGID']
    con <- dbConnect(sqlite,"www/hg19_bioconductor.sqlite")   
    query=paste0("select a.go_id,a.Evidence,a.Ontology,b.term from goID2geneID a ",
            "  left join go2name b on (a.go_id=b.go_id)", 
            "  where a.Ontology='CC'  and a.gene_id=",shQuote(Entrez)
            )
    tmp=dbGetQuery(con,query)
    dbDisconnect(con)
    tmp
  })
  output$go_mf     <- renderTable({
    tmp=gene_info()
    Entrez<<-tmp[1,'EGID']
    con <- dbConnect(sqlite,"www/hg19_bioconductor.sqlite")  
    query=paste0("select a.go_id,a.Evidence,a.Ontology,b.term from goID2geneID a ",
            "  left join go2name b on (a.go_id=b.go_id)", 
            "  where a.Ontology='MF'  and a.gene_id=",shQuote(Entrez)
            )
    tmp=dbGetQuery(con,query)
    dbDisconnect(con)
    tmp
  })
  ## need to change
  if(T){
        output$batch_infor  <- renderTable({
      	inFile <- input$file1
          if (is.null(inFile))
            return(NULL)   
          gene_list=read.table(inFile$datapath, header=input$header)
          #gene_list=gene_list[,1]
          names(gene_list)[1]='ID';
          con <- dbConnect(sqlite,"www/hg19_bioconductor.sqlite")  
          tmp=dbGetQuery(con,"select * from my_gene_mapping")
          dbDisconnect(con)
          if ( !is.na(match(gene_list[1,1],tmp[,'symbol'])) ){
            t=merge(gene_list,tmp,by.x='ID',by.y='symbol',all.x = T)
          }
          if ( !is.na(match(gene_list[1,1],tmp[,'EGID'])) ){
            t=merge(gene_list,tmp,by.x='ID',by.y='EGID',all.x = T)
          }
          if ( !is.na(match(gene_list[1,1],tmp[,'ENSEMBL'])) ){
            t=merge(gene_list,tmp,by.x='ID',by.y='ENSEMBL',all.x = T)
          }
          t
        })
  
   }
  
  
})