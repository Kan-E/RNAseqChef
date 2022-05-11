library(shiny)
library(DT)
library(gdata)
library(rstatix)
library(multcomp)
library(tidyverse)
library(tools)
library(ggpubr)
library(tidyverse)
library(ggrepel)
library(ggdendro)
library(ggplotify)
library(gridExtra)
library(cowplot)
library(DESeq2)
library(EBSeq)
library(edgeR)
library(IHW)
library(qvalue)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(genefilter)
library(ComplexHeatmap)
library("shinyAce",verbose=FALSE) # for showing text files, code
library(shinyBS,verbose=FALSE) # for popup figures
library(plotly,verbose=FALSE)
library('shinyjs', verbose = FALSE)
library('reactable', verbose = FALSE)

ui<-navbarPage(
  footer=p(hr(),p("ShinyApp created by Kan Etoh",align="center",width=4),
           p(("Copyright (C) 2022, code licensed under MIT"),align="center",width=4),
           p(("Code available on Github:"),a("https://github.com/Kan-E/Automate_shiny",href="https://github.com/Kan-E/Automate_shiny"),align="center",width=4)),
  
  "RNAseqchef",
  id='navBar',
  tabPanel("Pair-wise",
           # titlePanel(h5("Upload Files")),
           sidebarLayout(
             # pair-wise -------------------------------------
             # sidebar---------------------------------
             sidebarPanel(
               radioButtons('data_file_type','Input:',
                            c('Row_count_matrix'="Row1",
                              'Option: Row_count_matrix + Metadata'="Row2"
                            ),selected = "Row1"),
               # Conditional panels appear based on input.data_file_type selection
               conditionalPanel(condition="input.data_file_type=='Row1'",
                                strong("Count matrix format: "),br(),
                                "The replication number is represented by the underbar.",br(),
                                "Do not use it for anything else.", br(),
                                img(src="input_format1.png", height = 125, width = 204), br(),
                                fileInput("file3",
                                          label = "Select a row count matrix file (txt, csv)",
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%")
               ),
               conditionalPanel(condition="input.data_file_type=='Row2'",
                                strong("Count matrix format: "),br(),
                                "You can use the matrix file whose column name is accession number, and extract the colums you want to analyze by using",
                                "the metadata.", br(),
                                "The replication number is represented by the underbar in the characteristics of metadata.",br(),br(),
                                img(src="input_format2.png", height = 302, width = 253), br(),
                                fileInput("file1",
                                          label = "Select a row count matrix file (txt, csv)",
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%"),
                                fileInput("file2",
                                          label = "Select a metadata file to define samples for the following analysis",
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%")
               ),
               radioButtons('DEG_method','DEG analysis method:',
                            c('DESeq2'="DESeq2",
                              'EBSeq'="EBSeq",
                              'edgeR'="edgeR"
                            ),selected = "DESeq2"),
               conditionalPanel(condition="input.DEG_method=='DESeq2'",
                                selectInput("FDR_method", "FDR method", c("BH", "Qvalue", "IHW"), selected = "BH")
               ),
               conditionalPanel(condition="input.DEG_method=='edgeR'",
                                selectInput("FDR_method", "FDR method", c("BH", "Qvalue", "IHW"), selected = "BH")
               ),
               conditionalPanel(condition="input.DEG_method=='edgeR'"),
               fluidRow(
                 column(6, selectInput("Species", "Species", c("not selected", "human", "mouse"), selected = "not selected"))),
               h4("Cut-off conditions:"),
               fluidRow(
                 column(4, numericInput("fc", "Fold Change", min   = 0, max   = NA, value = 2)),
                 column(4, numericInput("fdr", "FDR", min   = 0, max   = NA, value = 0.05)),
                 column(4, numericInput("basemean", "Basemean", min   = 0, max   = NA, value = 0))
               ),
               actionButton("goButton", "example data"),
               tags$head(tags$style("#goButton{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                         tags$style("
          body {
            padding: 0 !important;
          }"
                         ))
             ), #sidebarPanel
             
             # Main Panel -------------------------------------
             mainPanel(
               tabsetPanel(
                 type = "tabs",
                 tabPanel("Input Data",
                          bsCollapse(id="input_collapse_panel",open="Row_count_panel",multiple = FALSE,
                                     bsCollapsePanel(title="Row_count_matrix:",
                                                     value="Row_count_panel",
                                                     dataTableOutput('Row_count_matrix')
                                     ),
                                     bsCollapsePanel(title="Metadata:",
                                                     value="Metadata_panel",
                                                     dataTableOutput('Metadata')
                                     ),
                                     bsCollapsePanel(title="Defined_row_count_matrix:",
                                                     value="D_row_count_matrix_panel",
                                                     dataTableOutput('D_Row_count_matrix')
                                     )
                          )
                 ),
                 tabPanel("Result overview",
                          plotOutput("PCA"),
                          plotOutput("MA"),
                          dataTableOutput("DEG_result"),
                          dataTableOutput("Normalized_Count_matrix")),
                 tabPanel("Box plot",
                          plotOutput("box")),
                 tabPanel("Enrichment analysis",
                          plotOutput("enrichment1"),
                          plotOutput("enrichment2"))
               )
             ) # main panel
           ) #sidebarLayout
  ), #tabPanel
  
  tabPanel("3 conditions",
           # titlePanel(h5("Upload Files")),
           sidebarLayout(
             # 3conditions -------------------------------------
             # sidebar_3conditions---------------------------------
             sidebarPanel(
               radioButtons('data_file_type2','Input:',
                            c('Row_count_matrix'="Row3",
                              'Option: Row_count_matrix + Metadata'="Row4"
                            ),selected = "Row3"),
               # Conditional panels appear based on input.data_file_type selection
               conditionalPanel(condition="input.data_file_type2=='Row3'",
                                strong("Count matrix format: "),br(),
                                "The replication number is represented by the underbar.",br(),
                                "Do not use it for anything else.", br(),
                                img(src="input_format1.png", height = 125, width = 204), br(),
                                fileInput("file4",
                                          label = "Select a row count matrix file (txt, csv)",
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%")
               ),
               conditionalPanel(condition="input.data_file_type2=='Row4'",
                                strong("Count matrix format: "),br(),
                                "You can use the matrix file whose column name is accession number, and extract the colums you want to analyze by using",
                                "the metadata.", br(),
                                "The replication number is represented by the underbar in the characteristics of metadata.",br(),br(),
                                img(src="input_format2.png", height = 302, width = 253), br(),
                                fileInput("file5",
                                          label = "Select a row count matrix file (txt, csv)",
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%"),
                                fileInput("file6",
                                          label = "Select a metadata file to define samples for the following analysis",
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%")
               ),
               fluidRow(
                 column(6, selectInput("Species2", "Species", c("not selected", "human", "mouse"), selected = "not selected"))),
               h4("Cut-off conditions:"),
               fluidRow(
                 column(4, numericInput("fc2", "Fold Change", min   = 0, max   = NA, value = 2)),
                 column(4, numericInput("fdr2", "FDR", min   = 0, max   = NA, value = 0.05)),
                 column(4, numericInput("basemean2", "Basemean", min   = 0, max   = NA, value = 0))
               ),
               actionButton("goButton2", "example data"),
               tags$head(tags$style("#goButton{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                         tags$style("
          body {
            padding: 0 !important;
          }"
                         ))
             ), #sidebarPanel
             
             # Main Panel -------------------------------------
             mainPanel(
               tabsetPanel(
                 type = "tabs",
                 tabPanel("Input 3 conditions Data",
                          bsCollapse(id="input_collapse_panel2",open="Row_count_panel2",multiple = FALSE,
                                     bsCollapsePanel(title="Row_count_matrix:",
                                                     value="Row_count_panel2",
                                                     dataTableOutput('Row_count_matrix2')
                                     ),
                                     bsCollapsePanel(title="Metadata:",
                                                     value="Metadata_panel2",
                                                     dataTableOutput('Metadata2')
                                     ),
                                     bsCollapsePanel(title="Defined_row_count_matrix:",
                                                     value="D_row_count_matrix_panel2",
                                                     dataTableOutput('D_Row_count_matrix2')
                                     )
                          )
                 ),
                 tabPanel("Result overview",
                          plotOutput("PCA2"),
                          plotOutput("scatter_1"),
                          dataTableOutput("DEG_result2_1"),
                          plotOutput("scatter_2"),
                          dataTableOutput("DEG_result2_2"),
                          plotOutput("scatter_3"),
                          dataTableOutput("DEG_result2_3"),
                          dataTableOutput("Normalized_Count_matrix2")),
                 tabPanel("Box plot",
                          plotOutput("box2")),
                 tabPanel("Enrichment analysis",
                          plotOutput("keggenrichment2_1"),
                          plotOutput("keggenrichment2_2"),
                          plotOutput("keggenrichment2_3"))
               )
             ) # main panel
           ) #sidebarLayout
  ) #tabPanel
)
# server ---------------------------------
server <- function(input, output, session) {
  # pair-wise ------------------------------------------------------------------------------
  row_count_matrix <- reactive({
    if (input$data_file_type == "Row1"){
      tmp <- input$file3$datapath
      if(is.null(tmp)) {
        return(NULL)
      }else{
        if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
        if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
        if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
        return(df)
      }
    }else{
      tmp <- input$file1$datapath
      if(is.null(input$file1) && input$goButton == 0) return(NULL)
      if(is.null(input$file1) && input$goButton > 0 )  tmp = "data/example1.xlsx"
      if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
      if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
      if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
      return(df)
    }
  })
  metadata <- reactive({
    tmp <- input$file2$datapath
    if(is.null(input$file2) && input$goButton == 0) return(NULL)
    if(is.null(input$file2) && input$goButton > 0 )  tmp = "data/example2.xlsx"
    if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
    if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
    if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
    return(df)
  })
  norm_count_matrix <- reactive({
    data <- input$file4$datapath
    if(is.null(input$file4) && input$goButton == 0) return(NULL)
    if(is.null(input$file4) && input$goButton > 0 )  tmp = "data/example3.xlsx"
    if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
    if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
    if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
    return(df)
  })
  d_row_count_matrix <- reactive({
    row <- row_count_matrix()
    if (input$data_file_type == "Row1"){
      if(is.null(row)) {
        return(NULL)
      }else{
        return(row)
      }
    }else{
      meta <- metadata()
      if (is.null(row) || is.null(meta)){
        return(NULL)
      } else {
        row_t <- t(row)
        colname <- colnames(meta)
        data <- merge(meta, row_t, by=0, sort = F)
        rownames(data) <- data$characteristics
        data2 <- data[, - which(colnames(data) %in% c("Row.names", colname))]
        data2_t <- t(data2)
        data3 <- apply(data2_t, 2, as.numeric)
        rownames(data3) <- rownames(data2_t)
        return(data3)
      }
    }
  })
  
  
  # pair-wise DEG ------------------------------------------------------------------------------
  dds <- reactive({
    count <- d_row_count_matrix()
    file_name <- gsub("\\..+$", "", input$file1)
    collist <- gsub("\\_.+$", "", colnames(count))
    if (input$DEG_method == "DESeq2") {
      group <- data.frame(con = factor(collist))
      dds<- DESeqDataSetFromMatrix(countData = round(count),colData = group, design = ~ con)
      dds$con <- factor(dds$con, levels = unique(collist))
      dds <- DESeq(dds)
    }
    if (input$DEG_method == "edgeR") {
      group <- factor(collist)
      dds <- DGEList(counts = count, group = group)
      keep <- filterByExpr(dds)
      dds = dds[keep, , keep.lib.sizes=FALSE]
      dds <- calcNormFactors(dds)
      dds <- estimateCommonDisp(dds)
      dds <- estimateTagwiseDisp(dds)
    }
    return(dds)
  })
  
  deg_result <- reactive({
    count <- d_row_count_matrix()
    file_name <- gsub("\\..+$", "", input$file1)
    collist <- gsub("\\_.+$", "", colnames(count))
    if (input$DEG_method == "DESeq2") {
      dds <- dds()
      contrast <- c("con", unique(collist))
      res <- results(dds,  contrast = contrast)
      if(input$FDR_method == "IHW") {
        ihw_res <- ihw(pvalue ~ baseMean,  data=as.data.frame(res), alpha = 0.1)
        res$padj <- IHW::as.data.frame(ihw_res)$adj_pvalue
      }
      if(input$FDR_method == "Qvalue") {
        res <- results(dds,  contrast = contrast)
        qvalue <- qvalue::qvalue(res$pvalue)
        res$padj <- qvalue$qvalues
      }
    }
    if(input$DEG_method == "edgeR"){
      dds <- dds()
      group <- factor(collist)
      result <- exactTest(dds, pair = c(unique(group)[2],unique(group)[1]))
      res <- as.data.frame(topTags(result, n = nrow(count)))
      qvalue <- qvalue::qvalue(res$PValue)
      res$padj <- qvalue$qvalues
      ihw_res <- ihw(PValue ~ 2^logCPM,  data=res, alpha = 0.1)
      ihw_res_df <- IHW::as.data.frame(ihw_res)
      res$ihw_padj <- ihw_res_df$adj_pvalue
      if(input$FDR_method == "BH"){label <- c("log2FoldChange", "log2CPM", "PValue","padj", "Qvalue", "IHW_FDR")}
      if(input$FDR_method == "Qvalue"){label <- c("log2FoldChange", "log2CPM", "PValue","BH_FDR", "padj", "IHW_FDR")}
      if(input$FDR_method == "IHW"){label <- c("log2FoldChange", "log2CPM", "PValue","BH_FDR", "Qvalue", "padj")}
      colnames(res) <- label
    }
    if(input$DEG_method == "EBSeq"){
      count <- data.matrix(count)
      vec <- c()
      for (i in 1:length(unique(collist))) {
        num <- length(collist[collist == unique(collist)[i]])
        vec <- c(vec, num)
      }
      ngvector <- NULL
      conditions <- as.factor(rep(paste("C", 1:length(unique(collist)), sep=""), times = vec))
      Sizes <- MedianNorm(count)
      NormMat <- GetNormalizedMat(count, Sizes)
      EBOut <- NULL
      EBOut <- EBTest(Data = count, NgVector = ngvector, Conditions = conditions, sizeFactors = Sizes, maxround = 5)
      stopifnot(!is.null(EBOut))
      PP <- as.data.frame(GetPPMat(EBOut))
      fc_res <- PostFC(EBOut)
      results <- cbind(PP, fc_res$PostFC, fc_res$RealFC,unlist(EBOut$C1Mean)[rownames(PP)], unlist(EBOut$C2Mean)[rownames(PP)])
      colnames(results) <- c("PPEE", "PPDE", "PostFC", "RealFC","C1Mean","C2Mean")
      res <- results[order(results[,"PPDE"], decreasing = TRUE),]
    }
    res <- as.data.frame(res)
    if(input$Species != "not selected"){
      if(str_detect(rownames(count)[1], "ENS")){
        switch (input$Species,
                "mouse" = org <- org.Mm.eg.db,
                "human" = org <- org.Hs.eg.db)
        switch (input$Species,
                "mouse" = org_code <- "mmu",
                "human" = org_code <- "hsa")
        my.symbols <- rownames(res)
        gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                        keytype = "ENSEMBL",
                                        columns = c("ENSEMBL","SYMBOL"))
        colnames(gene_IDs) <- c("Row.names","SYMBOL")
        res$Row.names <- rownames(res)
        gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
        data2 <- merge(res, gene_IDs, by="Row.names")
        rownames(data2) <- data2$Row.names
        res <- data2[,-1]
      }
    }
    return(res)
  })
  
  deg_norm_count <- reactive({
    count <- d_row_count_matrix()
    file_name <- gsub("\\..+$", "", input$file1)
    collist <- gsub("\\_.+$", "", colnames(count))
    group <- data.frame(con = factor(collist))
    if (input$DEG_method == "DESeq2") {
      dds <- dds()
      contrast <- c("con", unique(collist))
      normalized_counts <- counts(dds, normalized=TRUE)
    }
    if (input$DEG_method == "edgeR") {
      dds <- dds()
      normalized_counts <- t(t(dds$pseudo.counts)*(dds$samples$norm.factors))
    }
    if(input$DEG_method == "EBSeq"){
      count <- data.matrix(count)
      vec <- c()
      for (i in 1:length(unique(collist))) {
        num <- length(collist[collist == unique(collist)[i]])
        vec <- c(vec, num)
      }
      ngvector <- NULL
      conditions <- as.factor(rep(paste("C", 1:length(unique(collist)), sep=""), times = vec))
      Sizes <- MedianNorm(count)
      normalized_counts <- GetNormalizedMat(count, Sizes)
    }
    if(input$Species != "not selected"){
      if(str_detect(rownames(count)[1], "ENS")){
        switch (input$Species,
                "mouse" = org <- org.Mm.eg.db,
                "human" = org <- org.Hs.eg.db)
        switch (input$Species,
                "mouse" = org_code <- "mmu",
                "human" = org_code <- "hsa")
        normalized_counts <- as.data.frame(normalized_counts)
        my.symbols <- rownames(normalized_counts)
        gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                        keytype = "ENSEMBL",
                                        columns = c("ENSEMBL","SYMBOL"))
        colnames(gene_IDs) <- c("Row.names","SYMBOL")
        normalized_counts$Row.names <- rownames(normalized_counts)
        gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
        data2 <- merge(normalized_counts, gene_IDs, by="Row.names")
        rownames(data2) <- data2$Row.names
        normalized_counts <- data2[,-1]
      }
    }
    return(normalized_counts)
  })
  
  observeEvent(input$file1, ({
    updateCollapse(session,id =  "input_collapse_panel", open="Row_count_matrix_panel")
  }))
  observeEvent(input$file2, ({
    updateCollapse(session,id =  "input_collapse_panel", open="Metadata_panel")
  }))
  output$Row_count_matrix <- DT::renderDataTable({
    row_count_matrix()
  }, extensions = c('Buttons'), options = list(dom = 'Blfrtip', buttons = c('csv', 'excel', 'pdf'))
  )
  output$Metadata <- DT::renderDataTable({
    metadata()
  }, extensions = c('Buttons'), options = list(dom = 'Blfrtip', buttons = c('csv', 'excel', 'pdf'))
  )
  output$D_Row_count_matrix <- DT::renderDataTable({
    d_row_count_matrix()
  }, extensions = c('Buttons'), options = list(dom = 'Blfrtip', buttons = c('csv', 'excel', 'pdf'))
  )
  output$DEG_result <- DT::renderDataTable({
    deg_result()
  }, extensions = c('Buttons'), options = list(dom = 'Blfrtip', buttons = c('csv', 'excel', 'pdf'))
  )
  output$DEG_result <- DT::renderDataTable({
    deg_result()
  }, extensions = c('Buttons'), options = list(dom = 'Blfrtip', buttons = c('csv', 'excel', 'pdf'))
  )
  #pair-wise DEG vis---------------------------
  data_degcount <- reactive({
    data <- deg_result()
    count <- deg_norm_count()
    if(str_detect(rownames(data)[1], "ENS")){
      if(length(grep("SYMBOL", colnames(data))) != 0){
        data <- data[, - which(colnames(data) == "SYMBOL")]
        count <- count[, - which(colnames(count) == "SYMBOL")]
      }
    }
    collist <- factor(gsub("\\_.+$", "", colnames(count)))
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    Cond_1 <- vec[1]
    Cond_2 <- vec[2]
    Row.names <- NULL
    log2FoldChange <- NULL
    value <- NULL
    data <- merge(data,count, by=0)
    Type <- input$DEG_method
    data <- dplyr::filter(data, apply(data[,8:(7 + Cond_1 + Cond_2)],1,mean) > input$basemean)
    
    if(input$Species != "not selected"){
      switch (input$Species,
              "mouse" = org <- org.Mm.eg.db,
              "human" = org <- org.Hs.eg.db)
      switch (input$Species,
              "mouse" = org_code <- "mmu",
              "human" = org_code <- "hsa")
    }
    if (Type == "EBSeq"){
      data$padj <- data$PPEE
      data$log2FoldChange <- -1 * log2(data$PostFC)
      baseMean <- (data$C1Mean + data$C2Mean)*(1/2)
      data <- cbind(data, baseMean)
    }
    
    if (Type == "DESeq2"){
      data$log2FoldChange <- -1 * data$log2FoldChange
    }
    if(Type == "edgeR") {
      colnames(data)[3] <- "baseMean"
      data$baseMean <- 2^data$baseMean
      data$log2FoldChange <- -1 * data$log2FoldChange
    }
    
    if(str_detect(rownames(count)[1], "ENS")){
      if(input$Species != "not selected"){
        my.symbols <- data$Row.names
        gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                        keytype = "ENSEMBL",
                                        columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
        colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
        gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
        data <- merge(data, gene_IDs, by="Row.names")
        genenames <- as.vector(data$SYMBOL)
      }else{
        genenames=NULL
      }
    }else{
      if(input$Species != "not selected"){
        my.symbols <- data$Row.names
        gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                        keytype = "SYMBOL",
                                        columns = c("SYMBOL", "ENTREZID"))
        colnames(gene_IDs) <- c("Row.names", "ENTREZID")
        gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
        data <- merge(data, gene_IDs, by="Row.names")
      }
      genenames <- as.vector(data$Row.names)
    }
    
    
    return(data)
  })
  
  data_degcount2 <- reactive({
    data <- data_degcount()
    count <- deg_norm_count()
    collist <- factor(gsub("\\_.+$", "", colnames(count)))
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    Cond_1 <- vec[1]
    Cond_2 <- vec[2]
    data2 <- dplyr::filter(data, abs(data$log2FoldChange) > log(input$fc, 2))
    if(nrow(data2) != 0){
      data2$group <- "upregulated"
      data2$group[data2$log2FoldChange < 0] <- "downregulated"
      data3 <- dplyr::filter(data2, abs(data2$padj) < input$fdr)
      return(data3)
    }else{return(NULL)}
  })
  # pair-wise MA ------------------------------------------------------------------------------
  output$MA <- renderPlot({
    data <- data_degcount()
    count <- deg_norm_count()
    if(str_detect(rownames(data)[1], "ENS")){
      if(length(grep("SYMBOL", colnames(data))) != 0){
        count <- count[, - which(colnames(count) == "SYMBOL")]
      }
    }
    collist <- factor(gsub("\\_.+$", "", colnames(count)))
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    Cond_1 <- vec[1]
    Cond_2 <- vec[2]
    if(str_detect(rownames(count)[1], "ENS")){
      if(input$Species != "not selected"){
        genenames <- as.vector(data$SYMBOL)
      }else{ genenames=NULL }
    }else{
      genenames <- as.vector(data$Row.names)
    }
    if(input$FDR_method == "edgeR") {
      xlab <- "LogCPM"
    }else{
      xlab <- "Log2 mean expression"
    }
    m1 <- as.grob(ggmaplot(data, fdr = input$fdr, fc = input$fc, size = 0.4,
                           palette = c("#B31B21", "#1465AC", "darkgray"),
                           genenames = genenames,
                           legend = "top", top = 20,
                           font.label = c("bold", 6),font.legend = "bold",
                           font.main = "bold",xlab = xlab,
                           ggtheme = ggplot2::theme_minimal(),
                           select.top.method = "fc"))
    data2 <- data_degcount2()
    if(is.null(data2)){
      ht <- NULL
    }else{
      data.z <- genescale(data2[,8:(7 + Cond_1 + Cond_2)], axis=1, method="Z")
      ht <- as.grob(Heatmap(data.z, name = "z-score",column_order = colnames(data.z),
                            clustering_method_columns = 'ward.D2',
                            show_row_names = F, show_row_dend = F))
    }
    suppressWarnings(print(plot_grid(m1, ht, rel_widths = c(2, 1))))
  })
  
  # pair-wise PCA ------------------------------------------------------------------------------
  output$PCA <- renderPlot({
    data <- deg_norm_count()
    if(length(grep("SYMBOL", colnames(data))) != 0){
      data <- data[, - which(colnames(data) == "SYMBOL")]
    }
    pca <- prcomp(data, scale. = T)
    label<- colnames(data)
    label<- gsub("\\_.+$", "", label)
    lab_x <- paste(summary(pca)$importance[2,1]*100,
                   "% of variance)", sep = "")
    lab_x <- paste("PC1 (", lab_x, sep = "")
    lab_y <- paste(summary(pca)$importance[2,2]*100,
                   "% of variance)", sep = "")
    lab_y <- paste("PC2 (", lab_y, sep = "")
    pca$rotation <- as.data.frame(pca$rotation)
    g1 <- ggplot(pca$rotation,aes(x=pca$rotation[,1],
                                  y=pca$rotation[,2],
                                  col=label, label = label)) +
      geom_point()+
      theme(panel.background =element_rect(fill=NA,color=NA),
            panel.border = element_rect(fill = NA)) +
      xlab(lab_x) + ylab(lab_y) + geom_text_repel()  +
      theme(legend.position="none", aspect.ratio=1)
    rho <- cor(data,method="spearman")
    d <- dist(1-rho)
    mds <- as.data.frame(cmdscale(d))
    label<-colnames(data)
    label<-gsub("\\_.+$", "", label)
    g2 <- ggplot(mds, aes(x = mds[,1], y = mds[,2],
                          col = label, label = label)) +
      geom_point()+
      theme(panel.background =element_rect(fill=NA,color=NA),
            panel.border = element_rect(fill = NA)) +
      xlab("dim 1") + ylab("dim 2") +
      geom_text_repel() + theme(legend.position="none", aspect.ratio=1)
    x <- NULL
    y <- NULL
    xend <- NULL
    yend <- NULL
    data.t <- t(data)
    hc <- hclust(dist(data.t), "ward.D2")
    dendr <- dendro_data(hc, type="rectangle")
    g3 <- ggplot() +
      geom_segment(data=segment(dendr),
                   aes(x=x, y=y, xend=xend, yend=yend)) +
      geom_text(data=label(dendr),
                aes(x, y, label=label, hjust=0), size=3) +
      coord_flip()+ scale_y_reverse(expand=c(0.2, 0)) +
      theme(axis.line.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank(), aspect.ratio=1)
    gridExtra::grid.arrange(g1, g2, g3, nrow = 1)
  })
  # pair-wise enrichment ------------------------------------------------------------------------------
  output$enrichment1 <- renderPlot({
    data <- data_degcount()
    data3 <- data_degcount2()
    count <- deg_norm_count()
    if(input$Species != "not selected"){
      switch (input$Species,
              "mouse" = org <- org.Mm.eg.db,
              "human" = org <- org.Hs.eg.db)
      switch (input$Species,
              "mouse" = org_code <- "mmu",
              "human" = org_code <- "hsa")
      formula_res <- try(compareCluster(ENTREZID~group, data=data3,
                                        fun="enrichKEGG", organism=org_code), silent = T)
      if (class(formula_res) == "try-error") {
        formula_res <- NA
        p1 <- NULL
      } else{
        if ((length(as.data.frame(formula_res)) == 0) ||
            is.na(unique(as.data.frame(formula_res)$qvalue))) {
          p1 <- NULL
        } else{
          p1 <- as.grob(dotplot(formula_res, color ="qvalue", font.size = 10))
        }
      }
      data <- na.omit(data)
      geneList <- data$log2FoldChange
      names(geneList) = as.character(data$ENTREZID)
      geneList <- sort(geneList, decreasing = TRUE)
      kk3 <- gseKEGG(geneList = geneList, pvalueCutoff = 0.05,
                     organism = org_code, keyType = "kegg",
                     exponent = 1, eps = 0, pAdjustMethod = "BH",
                     minGSSize = 50, maxGSSize = 500, by = "fgsea",
                     use_internal_data = FALSE, verbose = F)
      if (length(kk3$ID) == 0) {
        p4 <- NULL
      } else{
        kk3 <- setReadable(kk3, org, 'ENTREZID')
        if (length(kk3$ID) >= 5){
          p4 <- as.grob(gseaplot2(kk3, 1:5, pvalue_table = F))
        }else{
          p4 <- as.grob(gseaplot2(kk3, 1:length(kk3$ID), pvalue_table = F))
        }
      }
      print(plot_grid(p1, p4, nrow = 1))
    }
  })
  
  output$enrichment2 <- renderPlot({
    data <- data_degcount()
    data3 <- data_degcount2()
    count <- deg_norm_count()
    if(input$Species != "not selected"){
      switch (input$Species,
              "mouse" = org <- org.Mm.eg.db,
              "human" = org <- org.Hs.eg.db)
      switch (input$Species,
              "mouse" = org_code <- "mmu",
              "human" = org_code <- "hsa")
      formula_res <- try(compareCluster(ENTREZID~group, data=data3,
                                        fun="enrichKEGG", organism=org_code), silent = T)
      if (class(formula_res) == "try-error") {
        formula_res <- NA
        p1 <- NULL
      } else{
        if ((length(as.data.frame(formula_res)) == 0) ||
            is.na(unique(as.data.frame(formula_res)$qvalue))) {
          p1 <- NULL
        } else{
          p1 <- as.grob(dotplot(formula_res, color ="qvalue", font.size = 10))
        }
      }
      upgene <- data3[data3$log2FoldChange > log(input$fc, 2),]
      geneList_up <- upgene$log2FoldChange
      names(geneList_up) = as.character(upgene$ENTREZID)
      kk1 <- enrichKEGG(upgene$ENTREZID, organism =org_code,
                        pvalueCutoff = 0.05, minGSSize = 50, maxGSSize = 500)
      if(is.null(kk1)){
        cnet1 <- NULL
      } else cnet1 <- setReadable(kk1, org, 'ENTREZID')
      if (length(cnet1$ID) == 0) {
        p2 <- NULL
      } else{
        p2 <- as.grob(cnetplot(cnet1, foldChange=geneList_up,
                               cex_label_gene = 0.75, cex_label_category = 1,
                               cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none"))
        downgene <- data3[data3$log2FoldChange < log(1/input$fc, 2),]
        geneList_down <- downgene$log2FoldChange
        names(geneList_down) = as.character(downgene$ENTREZID)
        kk2 <- enrichKEGG(downgene$ENTREZID, organism =org_code,
                          pvalueCutoff = 0.05, minGSSize = 50, maxGSSize = 500)
        if(is.null(kk2)){
          cnet2 <- NULL
        } else cnet2 <- setReadable(kk2, org, 'ENTREZID')
        if (length(cnet2$ID) == 0) {
          p3 <- NULL
        } else{
          p3 <- as.grob(cnetplot(cnet2, foldChange=geneList_down,
                                 cex_label_gene = 0.75, cex_label_category = 1,
                                 cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none"))
        }
      }
      print(plot_grid(p2, p3, nrow = 1))
    }
  })
  
  
  
  output$Normalized_Count_matrix <- DT::renderDataTable({
    deg_norm_count()
  }, extensions = c('Buttons'), options = list(dom = 'Blfrtip', buttons = c('csv', 'excel', 'pdf'))
  )
  
  # 3 conditions ------------------------------------------------------------------------------
  row_count_matrix2 <- reactive({
    if (input$data_file_type2 == "Row3"){
      tmp <- input$file4$datapath
      if(is.null(tmp)) {
        return(NULL)
      }else{
        if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
        if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
        if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
        return(df)
      }
    }else{
      tmp <- input$file5$datapath
      if(is.null(input$file5) && input$goButton == 0) return(NULL)
      if(is.null(input$file5) && input$goButton > 0 )  tmp = "data/example3.xlsx"
      if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
      if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
      if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
      return(df)
    }
  })
  metadata2 <- reactive({
    tmp <- input$file6$datapath
    if(is.null(input$file6) && input$goButton == 0) return(NULL)
    if(is.null(input$file6) && input$goButton > 0 )  tmp = "data/example4.xlsx"
    if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
    if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
    if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
    return(df)
  })
  norm_count_matrix2 <- reactive({
    data <- input$file10$datapath
    if(is.null(input$file10) && input$goButton == 0) return(NULL)
    if(is.null(input$file10) && input$goButton > 0 )  tmp = "data/example3.xlsx"
    if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
    if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
    if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
    return(df)
  })
  d_row_count_matrix2 <- reactive({
    row <- row_count_matrix2()
    if (input$data_file_type2 == "Row3"){
      if(is.null(row)) {
        return(NULL)
      }else{
        return(row)
      }
    }else{
      meta <- metadata2()
      if (is.null(row) || is.null(meta)){
        return(NULL)
      } else {
        row_t <- t(row)
        colname <- colnames(meta)
        data <- merge(meta, row_t, by=0, sort = F)
        rownames(data) <- data$characteristics
        data2 <- data[, - which(colnames(data) %in% c("Row.names", colname))]
        data2_t <- t(data2)
        data3 <- apply(data2_t, 2, as.numeric)
        rownames(data3) <- rownames(data2_t)
        return(data3)
      }
    }
  })
  
  observeEvent(input$file5, ({
    updateCollapse(session,id =  "input_collapse_panel2", open="Row_count_matrix_panel2")
  }))
  observeEvent(input$file6, ({
    updateCollapse(session,id =  "input_collapse_panel2", open="Metadata_panel2")
  }))
  output$Row_count_matrix2 <- DT::renderDataTable({
    row_count_matrix2()
  }, extensions = c('Buttons'), options = list(dom = 'Blfrtip', buttons = c('csv', 'excel', 'pdf'))
  )
  output$Metadata2 <- DT::renderDataTable({
    metadata2()
  }, extensions = c('Buttons'), options = list(dom = 'Blfrtip', buttons = c('csv', 'excel', 'pdf'))
  )
  output$D_Row_count_matrix2 <- DT::renderDataTable({
    d_row_count_matrix2()
  }, extensions = c('Buttons'), options = list(dom = 'Blfrtip', buttons = c('csv', 'excel', 'pdf'))
  )
  
  gene_ID <- reactive({
    res <- d_row_count_matrix2()
    if(is.null(res)){
      return(NULL)
    }else{
      if(input$Species2 != "not selected"){
        if(str_detect(rownames(count)[1], "ENS")){
          switch (input$Species2,
                  "mouse" = org <- org.Mm.eg.db,
                  "human" = org <- org.Hs.eg.db)
          switch (input$Species2,
                  "mouse" = org_code <- "mmu",
                  "human" = org_code <- "hsa")
          my.symbols <- rownames(res)
          gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL")
          res$Row.names <- rownames(res)
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          return(gene_IDs)
        }
      }else{ return(NULL) }
    }
  })
  
  # 3 conditions DEG ------------------------------------------------------------------------------
  MultiOut <- reactive({
    count <- d_row_count_matrix2()
    collist <- gsub("\\_.+$", "", colnames(count))
    count <- data.matrix(count)
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    ngvector <- NULL
    conditions <- as.factor(rep(paste("C", 1:length(unique(collist)), sep=""), times = vec))
    Sizes <- MedianNorm(count)
    NormMat <- GetNormalizedMat(count, Sizes)
    patterns <- GetPatterns(conditions)
    eename <- rownames(patterns)[which(rowSums(patterns) == length(unique(collist)))]
    stopifnot(length(eename) == 1)
    MultiOut <- NULL
    MultiOut <- EBMultiTest(Data = count, NgVector = ngvector, Conditions = conditions, AllParti = patterns, sizeFactors = Sizes, maxround = 5)
    stopifnot(!is.null(MultiOut))
    return(MultiOut)
  })
  
  deg_result2 <- reactive({
    count <- d_row_count_matrix2()
    collist <- gsub("\\_.+$", "", colnames(count))
    count <- data.matrix(count)
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    ngvector <- NULL
    conditions <- as.factor(rep(paste("C", 1:length(unique(collist)), sep=""), times = vec))
    Sizes <- MedianNorm(count)
    patterns <- GetPatterns(conditions)
    eename <- rownames(patterns)[which(rowSums(patterns) == length(unique(collist)))]
    MultiOut <- MultiOut()
    MultiPP <- GetMultiPP(MultiOut)
    PP <- as.data.frame(MultiPP$PP)
    pos <- which(names(PP) == eename)
    probs <- rowSums(PP[,-pos])
    results <- cbind(PP, MultiPP$MAP[rownames(PP)], probs)
    colnames(results) <- c(colnames(PP), "MAP", "PPDE")
    ord <- order(results[,"PPDE"], decreasing = TRUE)
    results <- results[ord,]
    MultiFC <- GetMultiFC(MultiOut)
    res <- results
    
    if(input$Species2 != "not selected"){
      if(str_detect(rownames(count)[1], "ENS")){
        gene_IDs  <- gene_ID()
        data2 <- merge(res, gene_IDs, by="Row.names")
        rownames(data2) <- data2$Row.names
        res <- data2[,-1]
      }
    }
    return(as.data.frame(res))
  })
  deg_result2_pattern <- reactive({
    count <- d_row_count_matrix2()
    collist <- gsub("\\_.+$", "", colnames(count))
    count <- data.matrix(count)
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    ngvector <- NULL
    conditions <- as.factor(rep(paste("C", 1:length(unique(collist)), sep=""), times = vec))
    Sizes <- MedianNorm(count)
    patterns <- GetPatterns(conditions)
    eename <- rownames(patterns)[which(rowSums(patterns) == length(unique(collist)))]
    MultiOut <- MultiOut()
    MultiPP <- GetMultiPP(MultiOut)
    PP <- as.data.frame(MultiPP$PP)
    pos <- which(names(PP) == eename)
    probs <- rowSums(PP[,-pos])
    results <- cbind(PP, MultiPP$MAP[rownames(PP)], probs)
    colnames(results) <- c(colnames(PP), "MAP", "PPDE")
    ord <- order(results[,"PPDE"], decreasing = TRUE)
    results <- results[ord,]
    MultiFC <- GetMultiFC(MultiOut)
    res <- MultiPP$Pattern
    return(as.data.frame(res))
  })
  deg_result2_condmean <- reactive({
    count <- d_row_count_matrix2()
    collist <- gsub("\\_.+$", "", colnames(count))
    count <- data.matrix(count)
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    ngvector <- NULL
    conditions <- as.factor(rep(paste("C", 1:length(unique(collist)), sep=""), times = vec))
    Sizes <- MedianNorm(count)
    patterns <- GetPatterns(conditions)
    eename <- rownames(patterns)[which(rowSums(patterns) == length(unique(collist)))]
    MultiOut <- MultiOut()
    MultiPP <- GetMultiPP(MultiOut)
    PP <- as.data.frame(MultiPP$PP)
    pos <- which(names(PP) == eename)
    probs <- rowSums(PP[,-pos])
    results <- cbind(PP, MultiPP$MAP[rownames(PP)], probs)
    colnames(results) <- c(colnames(PP), "MAP", "PPDE")
    ord <- order(results[,"PPDE"], decreasing = TRUE)
    results <- results[ord,]
    MultiFC <- GetMultiFC(MultiOut)
    res <- MultiFC$CondMeans[ord,]
    return(as.data.frame(res))
  })
  
  deg_norm_count2 <- reactive({
    count <- d_row_count_matrix2()
    collist <- gsub("\\_.+$", "", colnames(count))
    group <- data.frame(con = factor(collist))
    count <- data.matrix(count)
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    ngvector <- NULL
    conditions <- as.factor(rep(paste("C", 1:length(unique(collist)), sep=""), times = vec))
    Sizes <- MedianNorm(count)
    normalized_counts <- GetNormalizedMat(count, Sizes)
    if(input$Species2 != "not selected"){
      if(str_detect(rownames(count)[1], "ENS")){
        gene_IDs  <- gene_ID()
        data2 <- merge(normalized_counts, gene_IDs, by="Row.names")
        rownames(data2) <- data2$Row.names
        normalized_counts <- data2[,-1]
      }
    }
    return(as.data.frame(normalized_counts))
  })
  
  output$DEG_result2_1 <- DT::renderDataTable({
    data_3degcount2_1()
  }, extensions = c('Buttons'), options = list(dom = 'Blfrtip', buttons = c('csv', 'excel', 'pdf'))
  )
  output$DEG_result2_2 <- DT::renderDataTable({
    data_3degcount2_2()
  }, extensions = c('Buttons'), options = list(dom = 'Blfrtip', buttons = c('csv', 'excel', 'pdf'))
  )
  output$DEG_result2_3 <- DT::renderDataTable({
    data_3degcount2_3()
  }, extensions = c('Buttons'), options = list(dom = 'Blfrtip', buttons = c('csv', 'excel', 'pdf'))
  )
  output$Normalized_Count_matrix2 <- DT::renderDataTable({
    deg_norm_count2()
  }, extensions = c('Buttons'), options = list(dom = 'Blfrtip', buttons = c('csv', 'excel', 'pdf'))
  )
  
  #3conditions DEG vis------------------------
  
  #3conditions DEG_1------------------------
  data_3degcount1_1 <- reactive({
    data <- deg_norm_count2()
    collist <- factor(gsub("\\_.+$", "", colnames(data)))
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    Cond_1 <- vec[1]
    Cond_2 <- vec[2]
    Cond_3 <- vec[3]
    collist <- unique(collist)
    result_Condm <- deg_result2_condmean()
    result_FDR <- deg_result2()
    if(str_detect(rownames(data)[1], "ENS")){
      if(length(grep("SYMBOL", colnames(data))) != 0){
        data <- data[, - which(colnames(data) == "SYMBOL")]
      }
    }
    if(input$Species2 != "not selected"){
      switch (input$Species2,
              "mouse" = org <- org.Mm.eg.db,
              "human" = org <- org.Hs.eg.db)
      switch (input$Species2,
              "mouse" = org_code <- "mmu",
              "human" = org_code <- "hsa")
    }
    specific = collist[1]
    FC_xlab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[2], ")"))
    FC_ylab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[3], ")"))
    result_Condm$FC_x <- log2((result_Condm$C1 + 0.01)/(result_Condm$C2 + 0.01))
    result_Condm$FC_y <- log2((result_Condm$C1 + 0.01)/(result_Condm$C3 + 0.01))
    Pattern1 <- "Pattern4"
    Pattern2 <- "Pattern5"
    result_FDR$FDR <- 1 - result_FDR$PPDE
    result <- merge(result_Condm, result_FDR, by=0)
    data$Row.names <- rownames(data)
    data2 <- merge(data, result, by="Row.names")
    result <- dplyr::filter(data2, apply(data2[,2:(Cond_1 + Cond_2 + Cond_3)],1,mean) > input$basemean2)
    sig <- rep(3, nrow(result))
    sig[which(result$FDR <= input$fdr2 & result$FC_x < log2(1/input$fc2) & result$FC_y < log2(1/input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 2
    sig[which(result$FDR <= input$fdr2 & result$FC_x > log2(input$fc2) & result$FC_y > log2(input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 1
    data3 <- data.frame(Row.names = result$Row.names, FC_x = result$FC_x,
                        FC_y = result$FC_y, padj = result$FDR, sig = sig, FC_xy = result$FC_x * result$FC_y)
    if((sum(sig == 1) >= 1) && (sum(sig == 2) >= 1)){
      new.levels <- c( paste0(paste0(specific,"_up: "), sum(sig == 1)), paste0(paste0(specific,"_down: "), sum(sig == 2)), "NS" )
      col = c("red","blue", "darkgray")}
    if((sum(sig == 1) >= 1) && (sum(sig == 2) == 0)){
      new.levels <- c(paste0(paste0(specific,"_up: "), sum(sig == 1)), "NS" )
      col = c("red", "darkgray")}
    if((sum(sig == 1) == 0) && (sum(sig == 2) >= 1)){
      new.levels <- c(paste0(paste0(specific,"_down: "), sum(sig == 2)), "NS" )
      col = c("blue", "darkgray")}
    if((sum(sig == 1) == 0) && (sum(sig == 2) == 0)){
      new.levels <- c("NS")
      col = "darkgray"}
    
    data3$sig <- factor(data3$sig, labels = new.levels)
    return(data3)
  })
  
  data_3degcount2_1 <- reactive({
    data3 <- data_3degcount1_1()
    if(input$Species2 != "not selected"){
      switch (input$Species2,
              "mouse" = org <- org.Mm.eg.db,
              "human" = org <- org.Hs.eg.db)
      switch (input$Species2,
              "mouse" = org_code <- "mmu",
              "human" = org_code <- "hsa")}
    if(length(unique(data3$sig)) == 1){
      data4 <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
      return(data4)
    } else {
      data4 <- dplyr::filter(data3, sig != "NS")
      if(str_detect(data4$Row.names[1], "ENS")){
        if(input$Species2 != "not selected"){
          my.symbols <- data4$Row.names
          gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data4 <- merge(data4, gene_IDs, by="Row.names")
        }
      }else{
        if(input$Species2 != "not selected"){
          my.symbols <- data4$Row.names
          gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                          keytype = "SYMBOL",
                                          columns = c("SYMBOL", "ENTREZID"))
          colnames(gene_IDs) <- c("Row.names", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data4 <- merge(data4, gene_IDs, by="Row.names")
        }
      }
      return(data4)
    }
  })
  
  #3conditions scatter + heatmap_1
  output$scatter_1 <- renderPlot({
    data3 <- data_3degcount1_1()
    data <- deg_norm_count2()
    if(str_detect(rownames(data)[1], "ENS")){
      if(length(grep("SYMBOL", colnames(data))) != 0){
        data <- data[, - which(colnames(data) == "SYMBOL")]
      }
    }
    collist <- factor(gsub("\\_.+$", "", colnames(data)))
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    Cond_1 <- vec[1]
    Cond_2 <- vec[2]
    Cond_3 <- vec[3]
    collist <- unique(collist)
    FC_xlab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[2], ")"))
    FC_ylab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[3], ")"))
    complete_data <- stats::na.omit(data3)
    labs_data <- subset(complete_data, padj <= input$fdr2 & Row.names !=
                          "" & (FC_x)*(FC_y) >= log2(input$fc2))
    labs_data<-  labs_data[sort(labs_data$FC_xy, decreasing = T, index=T)$ix,]
    labs_data <- dplyr::filter(labs_data, sig != "NS")
    labs_data <- utils::head(labs_data, 20)
    font.label <- data.frame(size=5, color="black", face = "plain")
    set.seed(42)
    FC_x <- FC_y <- sig <- Row.names <- padj <- NULL
    p <- ggplot(data3, aes(x = FC_x, y = FC_y)) + geom_point(aes(color = sig),size = 0.1)
    p <- p  + geom_hline(yintercept = c(-log2(input$fc2), log2(input$fc2)), linetype = c(2, 2), color = c("black", "black"))+
      geom_vline(xintercept = c(-log2(input$fc2), log2(input$fc2)),linetype = c(2, 2), color = c("black", "black"))
    p <- p + ggrepel::geom_text_repel(data = labs_data, mapping = aes(label = Row.names),
                                      box.padding = unit(0.35, "lines"), point.padding = unit(0.3,
                                                                                              "lines"), force = 1, fontface = font.label$face,
                                      size = font.label$size/2, color = font.label$color)
    p <- p +
      theme_bw()+ scale_color_manual(values = col)+
      theme(legend.position = "top" , legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 10),
            axis.text.y= ggplot2::element_text(size = 10),
            text = ggplot2::element_text(size = 10),
            title = ggplot2::element_text(size = 10)) +
      xlab(FC_xlab) + ylab(FC_ylab)
    
    data4 <- data_3degcount2_1()
    if(length(unique(data3$sig)) == 1){
      ht <- NULL
    }else{
      data$Row.names <- rownames(data)
      data5 <- merge(data, data4, by ="Row.names")
      rownames(data5) <- data5$Row.names
      data5 <- data5[,-1]
      data5 <- data5[,1: (Cond_1 + Cond_2 + Cond_3)]
      data.z <- genescale(data5, axis=1, method="Z")
      ht <- as.grob(Heatmap(data.z, name = "z-score", column_order = colnames(data.z),
                            clustering_method_columns = 'ward.D2',
                            show_row_names = F, show_row_dend = T))
    }
    suppressWarnings(print(plot_grid(p, ht, rel_widths = c(2, 1))))
  })
  
  
  #3conditions DEG_2------------------------
  data_3degcount1_2 <- reactive({
    data <- deg_norm_count2()
    collist <- factor(gsub("\\_.+$", "", colnames(data)))
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    Cond_1 <- vec[1]
    Cond_2 <- vec[2]
    Cond_3 <- vec[3]
    collist <- unique(collist)
    result_Condm <- deg_result2_condmean()
    result_FDR <- deg_result2()
    if(str_detect(rownames(data)[1], "ENS")){
      if(length(grep("SYMBOL", colnames(data))) != 0){
        data <- data[, - which(colnames(data) == "SYMBOL")]
      }
    }
    if(input$Species2 != "not selected"){
      switch (input$Species2,
              "mouse" = org <- org.Mm.eg.db,
              "human" = org <- org.Hs.eg.db)
      switch (input$Species2,
              "mouse" = org_code <- "mmu",
              "human" = org_code <- "hsa")
    }
    specific = collist[2]
    FC_xlab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[1], ")"))
    FC_ylab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[3], ")"))
    result_Condm$FC_x <- log2((result_Condm$C2 + 0.01)/(result_Condm$C1 + 0.01))
    result_Condm$FC_y <- log2((result_Condm$C2 + 0.01)/(result_Condm$C3 + 0.01))
    Pattern1 <- "Pattern4"
    Pattern2 <- "Pattern5"
    result_FDR$FDR <- 1 - result_FDR$PPDE
    result <- merge(result_Condm, result_FDR, by=0)
    data$Row.names <- rownames(data)
    data2 <- merge(data, result, by="Row.names")
    result <- dplyr::filter(data2, apply(data2[,2:(Cond_1 + Cond_2 + Cond_3)],1,mean) > input$basemean2)
    sig <- rep(3, nrow(result))
    sig[which(result$FDR <= input$fdr2 & result$FC_x < log2(1/input$fc2) & result$FC_y < log2(1/input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 2
    sig[which(result$FDR <= input$fdr2 & result$FC_x > log2(input$fc2) & result$FC_y > log2(input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 1
    data3 <- data.frame(Row.names = result$Row.names, FC_x = result$FC_x,
                        FC_y = result$FC_y, padj = result$FDR, sig = sig, FC_xy = result$FC_x * result$FC_y)
    if((sum(sig == 1) >= 1) && (sum(sig == 2) >= 1)){
      new.levels <- c( paste0(paste0(specific,"_up: "), sum(sig == 1)), paste0(paste0(specific,"_down: "), sum(sig == 2)), "NS" )
      col = c("red","blue", "darkgray")}
    if((sum(sig == 1) >= 1) && (sum(sig == 2) == 0)){
      new.levels <- c(paste0(paste0(specific,"_up: "), sum(sig == 1)), "NS" )
      col = c("red", "darkgray")}
    if((sum(sig == 1) == 0) && (sum(sig == 2) >= 1)){
      new.levels <- c(paste0(paste0(specific,"_down: "), sum(sig == 2)), "NS" )
      col = c("blue", "darkgray")}
    if((sum(sig == 1) == 0) && (sum(sig == 2) == 0)){
      new.levels <- c("NS")
      col = "darkgray"}
    
    data3$sig <- factor(data3$sig, labels = new.levels)
    return(data3)
  })
  
  data_3degcount2_2 <- reactive({
    data3 <- data_3degcount1_2()
    if(input$Species2 != "not selected"){
      switch (input$Species2,
              "mouse" = org <- org.Mm.eg.db,
              "human" = org <- org.Hs.eg.db)
      switch (input$Species2,
              "mouse" = org_code <- "mmu",
              "human" = org_code <- "hsa")}
    if(length(unique(data3$sig)) == 1){
      data4 <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
      return(data4)
    } else {
      data4 <- dplyr::filter(data3, sig != "NS")
      if(str_detect(data4$Row.names[1], "ENS")){
        if(input$Species2 != "not selected"){
          my.symbols <- data4$Row.names
          gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data4 <- merge(data4, gene_IDs, by="Row.names")
        }
      }else{
        if(input$Species2 != "not selected"){
          my.symbols <- data4$Row.names
          gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                          keytype = "SYMBOL",
                                          columns = c("SYMBOL", "ENTREZID"))
          colnames(gene_IDs) <- c("Row.names", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data4 <- merge(data4, gene_IDs, by="Row.names")
        }
      }
      return(data4)
    }
  })
  
  #3conditions scatter + heatmap_2
  output$scatter_2 <- renderPlot({
    data3 <- data_3degcount1_2()
    data <- deg_norm_count2()
    if(str_detect(rownames(data)[1], "ENS")){
      if(length(grep("SYMBOL", colnames(data))) != 0){
        data <- data[, - which(colnames(data) == "SYMBOL")]
      }
    }
    collist <- factor(gsub("\\_.+$", "", colnames(data)))
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    Cond_1 <- vec[1]
    Cond_2 <- vec[2]
    Cond_3 <- vec[3]
    collist <- unique(collist)
    FC_xlab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[1], ")"))
    FC_ylab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[3], ")"))
    complete_data <- stats::na.omit(data3)
    labs_data <- subset(complete_data, padj <= input$fdr2 & Row.names !=
                          "" & (FC_x)*(FC_y) >= log2(input$fc2))
    labs_data<-  labs_data[sort(labs_data$FC_xy, decreasing = T, index=T)$ix,]
    labs_data <- dplyr::filter(labs_data, sig != "NS")
    labs_data <- utils::head(labs_data, 20)
    font.label <- data.frame(size=5, color="black", face = "plain")
    set.seed(42)
    FC_x <- FC_y <- sig <- Row.names <- padj <- NULL
    p <- ggplot(data3, aes(x = FC_x, y = FC_y)) + geom_point(aes(color = sig),size = 0.1)
    p <- p  + geom_hline(yintercept = c(-log2(input$fc2), log2(input$fc2)), linetype = c(2, 2), color = c("black", "black"))+
      geom_vline(xintercept = c(-log2(input$fc2), log2(input$fc2)),linetype = c(2, 2), color = c("black", "black"))
    p <- p + ggrepel::geom_text_repel(data = labs_data, mapping = aes(label = Row.names),
                                      box.padding = unit(0.35, "lines"), point.padding = unit(0.3,
                                                                                              "lines"), force = 1, fontface = font.label$face,
                                      size = font.label$size/2, color = font.label$color)
    p <- p +
      theme_bw()+ scale_color_manual(values = col)+
      theme(legend.position = "top" , legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 10),
            axis.text.y= ggplot2::element_text(size = 10),
            text = ggplot2::element_text(size = 10),
            title = ggplot2::element_text(size = 10)) +
      xlab(FC_xlab) + ylab(FC_ylab)
    
    data4 <- data_3degcount2_2()
    if(length(unique(data3$sig)) == 1){
      ht <- NULL
    }else{
      data$Row.names <- rownames(data)
      data5 <- merge(data, data4, by ="Row.names")
      rownames(data5) <- data5$Row.names
      data5 <- data5[,-1]
      data5 <- data5[,1: (Cond_1 + Cond_2 + Cond_3)]
      data.z <- genescale(data5, axis=1, method="Z")
      ht <- as.grob(Heatmap(data.z, name = "z-score", column_order = colnames(data.z),
                            clustering_method_columns = 'ward.D2',
                            show_row_names = F, show_row_dend = T))
    }
    
    suppressWarnings(print(plot_grid(p, ht, rel_widths = c(2, 1))))
  })
  
  
  
  #3conditions DEG_3------------------------
  data_3degcount1_3 <- reactive({
    data <- deg_norm_count2()
    collist <- factor(gsub("\\_.+$", "", colnames(data)))
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    Cond_1 <- vec[1]
    Cond_2 <- vec[2]
    Cond_3 <- vec[3]
    collist <- unique(collist)
    result_Condm <- deg_result2_condmean()
    result_FDR <- deg_result2()
    if(str_detect(rownames(data)[1], "ENS")){
      if(length(grep("SYMBOL", colnames(data))) != 0){
        data <- data[, - which(colnames(data) == "SYMBOL")]
      }
    }
    if(input$Species2 != "not selected"){
      switch (input$Species2,
              "mouse" = org <- org.Mm.eg.db,
              "human" = org <- org.Hs.eg.db)
      switch (input$Species2,
              "mouse" = org_code <- "mmu",
              "human" = org_code <- "hsa")
    }
    specific = collist[3]
    FC_xlab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[1], ")"))
    FC_ylab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[2], ")"))
    result_Condm$FC_x <- log2((result_Condm$C3 + 0.01)/(result_Condm$C1 + 0.01))
    result_Condm$FC_y <- log2((result_Condm$C3 + 0.01)/(result_Condm$C2 + 0.01))
    Pattern1 <- "Pattern4"
    Pattern2 <- "Pattern5"
    result_FDR$FDR <- 1 - result_FDR$PPDE
    result <- merge(result_Condm, result_FDR, by=0)
    data$Row.names <- rownames(data)
    data2 <- merge(data, result, by="Row.names")
    result <- dplyr::filter(data2, apply(data2[,2:(Cond_1 + Cond_2 + Cond_3)],1,mean) > input$basemean2)
    sig <- rep(3, nrow(result))
    sig[which(result$FDR <= input$fdr2 & result$FC_x < log2(1/input$fc2) & result$FC_y < log2(1/input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 2
    sig[which(result$FDR <= input$fdr2 & result$FC_x > log2(input$fc2) & result$FC_y > log2(input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 1
    data3 <- data.frame(Row.names = result$Row.names, FC_x = result$FC_x,
                        FC_y = result$FC_y, padj = result$FDR, sig = sig, FC_xy = result$FC_x * result$FC_y)
    if((sum(sig == 1) >= 1) && (sum(sig == 2) >= 1)){
      new.levels <- c( paste0(paste0(specific,"_up: "), sum(sig == 1)), paste0(paste0(specific,"_down: "), sum(sig == 2)), "NS" )
      col = c("red","blue", "darkgray")}
    if((sum(sig == 1) >= 1) && (sum(sig == 2) == 0)){
      new.levels <- c(paste0(paste0(specific,"_up: "), sum(sig == 1)), "NS" )
      col = c("red", "darkgray")}
    if((sum(sig == 1) == 0) && (sum(sig == 2) >= 1)){
      new.levels <- c(paste0(paste0(specific,"_down: "), sum(sig == 2)), "NS" )
      col = c("blue", "darkgray")}
    if((sum(sig == 1) == 0) && (sum(sig == 2) == 0)){
      new.levels <- c("NS")
      col = "darkgray"}
    
    data3$sig <- factor(data3$sig, labels = new.levels)
    return(data3)
  })
  
  data_3degcount2_3 <- reactive({
    data3 <- data_3degcount1_3()
    if(input$Species2 != "not selected"){
      switch (input$Species2,
              "mouse" = org <- org.Mm.eg.db,
              "human" = org <- org.Hs.eg.db)
      switch (input$Species2,
              "mouse" = org_code <- "mmu",
              "human" = org_code <- "hsa")}
    if(length(unique(data3$sig)) == 1){
      data4 <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
      return(data4)
    } else {
      data4 <- dplyr::filter(data3, sig != "NS")
      if(str_detect(data4$Row.names[1], "ENS")){
        if(input$Species2 != "not selected"){
          my.symbols <- data4$Row.names
          gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data4 <- merge(data4, gene_IDs, by="Row.names")
        }
      }else{
        if(input$Species2 != "not selected"){
          my.symbols <- data4$Row.names
          gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                          keytype = "SYMBOL",
                                          columns = c("SYMBOL", "ENTREZID"))
          colnames(gene_IDs) <- c("Row.names", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data4 <- merge(data4, gene_IDs, by="Row.names")
        }
      }
      return(data4)
    }
  })
  
  #3conditions scatter + heatmap_3
  output$scatter_3 <- renderPlot({
    data3 <- data_3degcount1_3()
    complete_data <- stats::na.omit(data3)
    labs_data <- subset(complete_data, padj <= input$fdr2 & Row.names !=
                          "" & (FC_x)*(FC_y) >= log2(input$fc2))
    labs_data<-  labs_data[sort(labs_data$FC_xy, decreasing = T, index=T)$ix,]
    labs_data <- dplyr::filter(labs_data, sig != "NS")
    labs_data <- utils::head(labs_data, 20)
    font.label <- data.frame(size=5, color="black", face = "plain")
    data <- deg_norm_count2()
    if(str_detect(rownames(data)[1], "ENS")){
      if(length(grep("SYMBOL", colnames(data))) != 0){
        data <- data[, - which(colnames(data) == "SYMBOL")]
      }
    }
    collist <- factor(gsub("\\_.+$", "", colnames(data)))
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    Cond_1 <- vec[1]
    Cond_2 <- vec[2]
    Cond_3 <- vec[3]
    collist <- unique(collist)
    FC_xlab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[1], ")"))
    FC_ylab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[2], ")"))
    set.seed(42)
    FC_x <- FC_y <- sig <- Row.names <- padj <- NULL
    p <- ggplot(data3, aes(x = FC_x, y = FC_y)) + geom_point(aes(color = sig),size = 0.1)
    p <- p  + geom_hline(yintercept = c(-log2(input$fc2), log2(input$fc2)), linetype = c(2, 2), color = c("black", "black"))+
      geom_vline(xintercept = c(-log2(input$fc2), log2(input$fc2)),linetype = c(2, 2), color = c("black", "black"))
    p <- p + ggrepel::geom_text_repel(data = labs_data, mapping = aes(label = Row.names),
                                      box.padding = unit(0.35, "lines"), point.padding = unit(0.3,
                                                                                              "lines"), force = 1, fontface = font.label$face,
                                      size = font.label$size/2, color = font.label$color)
    p <- p +
      theme_bw()+ scale_color_manual(values = col)+
      theme(legend.position = "top" , legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 10),
            axis.text.y= ggplot2::element_text(size = 10),
            text = ggplot2::element_text(size = 10),
            title = ggplot2::element_text(size = 10)) +
      xlab(FC_xlab) + ylab(FC_ylab)
    
    data4 <- data_3degcount2_3()
    if(length(unique(data3$sig)) == 1){
      ht <- NULL
    }else{
      data$Row.names <- rownames(data)
      data5 <- merge(data, data4, by ="Row.names")
      rownames(data5) <- data5$Row.names
      data5 <- data5[,-1]
      data5 <- data5[,1: (Cond_1 + Cond_2 + Cond_3)]
      data.z <- genescale(data5, axis=1, method="Z")
      ht <- as.grob(Heatmap(data.z, name = "z-score", column_order = colnames(data.z),
                            clustering_method_columns = 'ward.D2',
                            show_row_names = F, show_row_dend = T))
    }
    
    suppressWarnings(print(plot_grid(p, ht, rel_widths = c(2, 1))))
  })
  
  
  
  
  #3conditions PCA--------------
  output$PCA2 <- renderPlot({
    data <- deg_norm_count2()
    if(length(grep("SYMBOL", colnames(data))) != 0){
      data <- data[, - which(colnames(data) == "SYMBOL")]
    }
    pca <- prcomp(data, scale. = T)
    label<- colnames(data)
    label<- gsub("\\_.+$", "", label)
    lab_x <- paste(summary(pca)$importance[2,1]*100,
                   "% of variance)", sep = "")
    lab_x <- paste("PC1 (", lab_x, sep = "")
    lab_y <- paste(summary(pca)$importance[2,2]*100,
                   "% of variance)", sep = "")
    lab_y <- paste("PC2 (", lab_y, sep = "")
    pca$rotation <- as.data.frame(pca$rotation)
    g1 <- ggplot(pca$rotation,aes(x=pca$rotation[,1],
                                  y=pca$rotation[,2],
                                  col=label, label = label)) +
      geom_point()+
      theme(panel.background =element_rect(fill=NA,color=NA),
            panel.border = element_rect(fill = NA)) +
      xlab(lab_x) + ylab(lab_y) + geom_text_repel()  +
      theme(legend.position="none", aspect.ratio=1)
    rho <- cor(data,method="spearman")
    d <- dist(1-rho)
    mds <- as.data.frame(cmdscale(d))
    label<-colnames(data)
    label<-gsub("\\_.+$", "", label)
    g2 <- ggplot(mds, aes(x = mds[,1], y = mds[,2],
                          col = label, label = label)) +
      geom_point()+
      theme(panel.background =element_rect(fill=NA,color=NA),
            panel.border = element_rect(fill = NA)) +
      xlab("dim 1") + ylab("dim 2") +
      geom_text_repel() + theme(legend.position="none", aspect.ratio=1)
    x <- NULL
    y <- NULL
    xend <- NULL
    yend <- NULL
    data.t <- t(data)
    hc <- hclust(dist(data.t), "ward.D2")
    dendr <- dendro_data(hc, type="rectangle")
    g3 <- ggplot() +
      geom_segment(data=segment(dendr),
                   aes(x=x, y=y, xend=xend, yend=yend)) +
      geom_text(data=label(dendr),
                aes(x, y, label=label, hjust=0), size=3) +
      coord_flip()+ scale_y_reverse(expand=c(0.2, 0)) +
      theme(axis.line.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank(), aspect.ratio=1)
    gridExtra::grid.arrange(g1, g2, g3, nrow = 1)
  })
  
  #3conditions enrichment ------------------------------------------------------------------------------
  output$keggenrichment2_1 <- renderPlot({
    
    data4 <- data_3degcount2_1()
    data3 <- data_3degcount1_1()
    cnetkegg_list <- list()
    if(input$Species2 != "not selected"){
      switch (input$Species2,
              "mouse" = org <- org.Mm.eg.db,
              "human" = org <- org.Hs.eg.db)
      switch (input$Species2,
              "mouse" = org_code <- "mmu",
              "human" = org_code <- "hsa")
      formula_res <- try(compareCluster(ENTREZID~sig, data=data4,fun="enrichKEGG", organism=org_code, universe = universe), silent = T)
      if (class(formula_res) == "try-error") {
        formula_res <- NA
        d <- NULL
      } else {
        if ((length(as.data.frame(formula_res)) == 0) ||
            is.na(unique(as.data.frame(formula_res)$qvalue))) {
          d <- NULL
        } else{
          d <- as.grob(dotplot(formula_res, showCategory=5, color ="qvalue" ,font.size=10))
        }
      }
      for (name in unique(data3$sig)) {
        if (name != "NS"){
          if (input$Species2 != "not selected"){
            kk1 <- enrichKEGG(data4$ENTREZID[data4$sig == name], organism =org_code)
            if (is.null(kk1)) {
              cnet1 <- NULL
            } else cnet1 <- setReadable(kk1, org, 'ENTREZID')
            if ((length(cnet1$ID) == 0) || is.na(unique(cnet1$qvalue))) {
              c <- NULL
            } else{
              c <- cnetplot(cnet1, cex_label_gene = 0.75, cex_label_category = 1,
                            cex_category = 0.75, colorEdge = TRUE)
              c <- as.grob(c + guides(edge_color = "none"))
              cnetkegg_list[[name]] = c
            }
          }
        }
      }
      if (input$Species2 != "not selected"){
        if (length(cnetkegg_list) == 2){
          cnetkegg1 <- cnetkegg_list[[1]]
          cnetkegg2 <- cnetkegg_list[[2]]}
        if (length(cnetkegg_list) == 1){
          cnetkegg1 <- cnetkegg_list[[1]]
          cnetkegg2 <- NULL}
        if (length(cnetkegg_list) == 0){
          cnetkegg1 <- NULL
          cnetkegg2 <- NULL}
      }
      print(plot_grid(d, cnetkegg1, cnetkegg2, ncol =3, nrow = 1))
      
    }
  })
  
  
  
  
}











# Run the app
shinyApp(ui, server)
