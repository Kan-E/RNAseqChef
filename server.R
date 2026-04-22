popoverTempate <- 
  '<div class="popover popover-lg" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content"></div></div>'
shinyServer(function(input, output, session) {
  options(shiny.maxRequestSize=3000*1024^2)
  safe_ordered_matrix <- function(count, order) {
    if (is.null(count) || is.null(dim(count)) || NCOL(count) == 0) {
      return(NULL)
    }
    available <- colnames(count)
    if (is.null(order) || !length(order)) {
      order <- available
    } else {
      order <- intersect(order, available)
      if (!length(order)) {
        order <- available
      }
    }
    count[, order, drop = FALSE]
  }
  safe_cluster_genes <- function(data, idx) {
    if (is.null(data) || is.null(dim(data)) || NROW(data) == 0) {
      return(character(0))
    }
    idx <- idx[!is.na(idx)]
    if (!length(idx)) {
      return(character(0))
    }
    rownames(data[idx, , drop = FALSE])
  }
  rename_limma_display_columns <- function(res) {
    if (is.null(res)) {
      return(NULL)
    }
    base_cols <- c("log2FoldChange", "baseMean", "t", "pval", "padj", "B")
    limit <- min(length(base_cols), ncol(res))
    colnames(res)[seq_len(limit)] <- base_cols[seq_len(limit)]
    res
  }

  observeEvent(input$Species,({
    if(sum(is.element(no_orgDb_plants,input$Species)) == 1){
      updateSelectInput(session,inputId = "Ortholog", 
                        "Ortholog",
                        "Arabidopsis thaliana", selected = "Arabidopsis thaliana")
      updateSelectInput(session,inputId = "Biomart_archive", "Biomart host", ensembl_archive_plants)
    }
    if(sum(is.element(no_orgDb_fungi,input$Species)) == 1){
      updateSelectInput(session,inputId = "Ortholog", 
                        "Ortholog",
                        "Saccharomyces cerevisiae", selected = "Saccharomyces cerevisiae")
      updateSelectInput(session,inputId = "Biomart_archive", "Biomart host", ensembl_archive_fungi)
    }
    if(sum(is.element(no_orgDb_metazoa,input$Species)) == 1){
      updateSelectInput(session,inputId = "Ortholog", 
                        "Ortholog",
                        c("Drosophila melanogaster","Caenorhabditis elegans"), selected = "Drosophila melanogaster")
      updateSelectInput(session,inputId = "Biomart_archive", "Biomart host", ensembl_archive_metazoa)
    }
    if(sum(is.element(no_orgDb_animals,input$Species)) == 1){
      updateSelectInput(session,inputId = "Ortholog", 
                        "Ortholog", orgDb_list, selected = "Mus musculus")
      updateSelectInput(session,inputId = "Biomart_archive", "Biomart host", ensembl_archive)
    }
  }))
  observeEvent(input$Species2,({
    if(sum(is.element(no_orgDb_plants,input$Species2)) == 1){
      updateSelectInput(session,inputId = "Ortholog2", 
                        "Ortholog",
                        "Arabidopsis thaliana", selected = "Arabidopsis thaliana")
      updateSelectInput(session,inputId = "Biomart_archive2", "Biomart host", ensembl_archive_plants)
    }
    if(sum(is.element(no_orgDb_fungi,input$Species2)) == 1){
      updateSelectInput(session,inputId = "Ortholog2", 
                        "Ortholog",
                        "Saccharomyces cerevisiae", selected = "Saccharomyces cerevisiae")
      updateSelectInput(session,inputId = "Biomart_archive2", "Biomart host", ensembl_archive_fungi)
    }
    if(sum(is.element(no_orgDb_metazoa,input$Species2)) == 1){
      updateSelectInput(session,inputId = "Ortholog2", 
                        "Ortholog",
                        c("Drosophila melanogaster","Caenorhabditis elegans"), selected = "Drosophila melanogaster")
      updateSelectInput(session,inputId = "Biomart_archive2", "Biomart host", ensembl_archive_metazoa)
    }
    if(sum(is.element(no_orgDb_animals,input$Species2)) == 1){
      updateSelectInput(session,inputId = "Ortholog2", 
                        "Ortholog", orgDb_list, selected = "Mus musculus")
      updateSelectInput(session,inputId = "Biomart_archive2", "Biomart host", ensembl_archive)
    }
  }))
  observeEvent(input$Species3,({
    if(sum(is.element(no_orgDb_plants,input$Species3)) == 1){
      updateSelectInput(session,inputId = "Ortholog3", 
                        "Ortholog",
                        "Arabidopsis thaliana", selected = "Arabidopsis thaliana")
      updateSelectInput(session,inputId = "Biomart_archive3", "Biomart host", ensembl_archive_plants)
    }
    if(sum(is.element(no_orgDb_fungi,input$Species3)) == 1){
      updateSelectInput(session,inputId = "Ortholog3", 
                        "Ortholog",
                        "Saccharomyces cerevisiae", selected = "Saccharomyces cerevisiae")
      updateSelectInput(session,inputId = "Biomart_archive3", "Biomart host", ensembl_archive_fungi)
    }
    if(sum(is.element(no_orgDb_metazoa,input$Species3)) == 1){
      updateSelectInput(session,inputId = "Ortholog3", 
                        "Ortholog",
                        c("Drosophila melanogaster","Caenorhabditis elegans"), selected = "Drosophila melanogaster")
      updateSelectInput(session,inputId = "Biomart_archive3", "Biomart host", ensembl_archive_metazoa)
    }
    if(sum(is.element(no_orgDb_animals,input$Species3)) == 1){
      updateSelectInput(session,inputId = "Ortholog3", 
                        "Ortholog", orgDb_list, selected = "Mus musculus")
      updateSelectInput(session,inputId = "Biomart_archive3", "Biomart host", ensembl_archive)
    }
  }))
  observeEvent(input$Species4,({
    if(sum(is.element(no_orgDb_plants,input$Species4)) == 1){
      updateSelectInput(session,inputId = "Ortholog4", 
                        "Ortholog",
                        "Arabidopsis thaliana", selected = "Arabidopsis thaliana")
      updateSelectInput(session,inputId = "Biomart_archive4", "Biomart host", ensembl_archive_plants)
    }
    if(sum(is.element(no_orgDb_fungi,input$Species4)) == 1){
      updateSelectInput(session,inputId = "Ortholog4", 
                        "Ortholog",
                        "Saccharomyces cerevisiae", selected = "Saccharomyces cerevisiae")
      updateSelectInput(session,inputId = "Biomart_archive4", "Biomart host", ensembl_archive_fungi)
    }
    if(sum(is.element(no_orgDb_metazoa,input$Species4)) == 1){
      updateSelectInput(session,inputId = "Ortholog4", 
                        "Ortholog",
                        c("Drosophila melanogaster","Caenorhabditis elegans"), selected = "Drosophila melanogaster")
      updateSelectInput(session,inputId = "Biomart_archive4", "Biomart host", ensembl_archive_metazoa)
    }
    if(sum(is.element(no_orgDb_animals,input$Species4)) == 1){
      updateSelectInput(session,inputId = "Ortholog4", 
                        "Ortholog", orgDb_list, selected = "Mus musculus")
      updateSelectInput(session,inputId = "Biomart_archive4", "Biomart host", ensembl_archive)
    }
  }))
  observeEvent(input$Species5,({
    if(sum(is.element(no_orgDb_plants,input$Species5)) == 1){
      updateSelectInput(session,inputId = "Ortholog5", 
                        "Ortholog",
                        "Arabidopsis thaliana", selected = "Arabidopsis thaliana")
      updateSelectInput(session,inputId = "Biomart_archive5", "Biomart host", ensembl_archive_plants)
    }
    if(sum(is.element(no_orgDb_fungi,input$Species5)) == 1){
      updateSelectInput(session,inputId = "Ortholog5", 
                        "Ortholog",
                        "Saccharomyces cerevisiae", selected = "Saccharomyces cerevisiae")
      updateSelectInput(session,inputId = "Biomart_archive5", "Biomart host", ensembl_archive_fungi)
    }
    if(sum(is.element(no_orgDb_metazoa,input$Species5)) == 1){
      updateSelectInput(session,inputId = "Ortholog5", 
                        "Ortholog",
                        c("Drosophila melanogaster","Caenorhabditis elegans"), selected = "Drosophila melanogaster")
      updateSelectInput(session,inputId = "Biomart_archive5", "Biomart host", ensembl_archive_metazoa)
    }
    if(sum(is.element(no_orgDb_animals,input$Species5)) == 1){
      updateSelectInput(session,inputId = "Ortholog5", 
                        "Ortholog", orgDb_list, selected = "Mus musculus")
      updateSelectInput(session,inputId = "Biomart_archive5", "Biomart host", ensembl_archive)
    }
  }))
  observeEvent(input$Species6,({
    if(sum(is.element(no_orgDb_plants,input$Species6)) == 1){
      updateSelectInput(session,inputId = "Ortholog6", 
                        "Ortholog",
                        "Arabidopsis thaliana", selected = "Arabidopsis thaliana")
      updateSelectInput(session,inputId = "Biomart_archive6", "Biomart host", ensembl_archive_plants)
    }
    if(sum(is.element(no_orgDb_fungi,input$Species6)) == 1){
      updateSelectInput(session,inputId = "Ortholog6", 
                        "Ortholog6",
                        "Saccharomyces cerevisiae", selected = "Saccharomyces cerevisiae")
      updateSelectInput(session,inputId = "Biomart_archive6", "Biomart host", ensembl_archive_fungi)
    }
    if(sum(is.element(no_orgDb_metazoa,input$Species6)) == 1){
      updateSelectInput(session,inputId = "Ortholog6", 
                        "Ortholog",
                        c("Drosophila melanogaster","Caenorhabditis elegans"), selected = "Drosophila melanogaster")
      updateSelectInput(session,inputId = "Biomart_archive6", "Biomart host", ensembl_archive_metazoa)
    }
    if(sum(is.element(no_orgDb_animals,input$Species6)) == 1){
      updateSelectInput(session,inputId = "Ortholog6", 
                        "Ortholog", orgDb_list, selected = "Mus musculus")
      updateSelectInput(session,inputId = "Biomart_archive6", "Biomart host", ensembl_archive)
    }
  }))
  observeEvent(input$Species7,({
    if(sum(is.element(no_orgDb_plants,input$Species7)) == 1){
      updateSelectInput(session,inputId = "Ortholog7", 
                        "Ortholog",
                        "Arabidopsis thaliana", selected = "Arabidopsis thaliana")
      updateSelectInput(session,inputId = "Biomart_archive7", "Biomart host", ensembl_archive_plants)
    }
    if(sum(is.element(no_orgDb_fungi,input$Species7)) == 1){
      updateSelectInput(session,inputId = "Ortholog7", 
                        "Ortholog",
                        "Saccharomyces cerevisiae", selected = "Saccharomyces cerevisiae")
      updateSelectInput(session,inputId = "Biomart_archive7", "Biomart host", ensembl_archive_fungi)
    }
    if(sum(is.element(no_orgDb_metazoa,input$Species7)) == 1){
      updateSelectInput(session,inputId = "Ortholog7", 
                        "Ortholog",
                        c("Drosophila melanogaster","Caenorhabditis elegans"), selected = "Drosophila melanogaster")
      updateSelectInput(session,inputId = "Biomart_archive7", "Biomart host", ensembl_archive_metazoa)
    }
    if(sum(is.element(no_orgDb_animals,input$Species7)) == 1){
      updateSelectInput(session,inputId = "Ortholog7", 
                        "Ortholog", orgDb_list, selected = "Mus musculus")
      updateSelectInput(session,inputId = "Biomart_archive7", "Biomart host", ensembl_archive)
    }
  }))
  observeEvent(input$Species_ens,({
    if(sum(is.element(no_orgDb_plants,input$Species_ens)) == 1){
      updateSelectInput(session,inputId = "Ortholog_ens", 
                        "Ortholog",
                        "Arabidopsis thaliana", selected = "Arabidopsis thaliana")
      updateSelectInput(session,inputId = "Biomart_archive_ens", "Biomart host", ensembl_archive_plants)
    }
    if(sum(is.element(no_orgDb_fungi,input$Species_ens)) == 1){
      updateSelectInput(session,inputId = "Ortholog_ens", 
                        "Ortholog",
                        "Saccharomyces cerevisiae", selected = "Saccharomyces cerevisiae")
      updateSelectInput(session,inputId = "Biomart_archive_ens", "Biomart host", ensembl_archive_fungi)
    }
    if(sum(is.element(no_orgDb_metazoa,input$Species_ens)) == 1){
      updateSelectInput(session,inputId = "Ortholog_ens", 
                        "Ortholog",
                        c("Drosophila melanogaster","Caenorhabditis elegans"), selected = "Drosophila melanogaster")
      updateSelectInput(session,inputId = "Biomart_archive_ens", "Biomart host", ensembl_archive_metazoa)
    }
    if(sum(is.element(no_orgDb_animals,input$Species_ens)) == 1){
      updateSelectInput(session,inputId = "Ortholog_ens", 
                        "Ortholog", orgDb_list, selected = "Mus musculus")
      updateSelectInput(session,inputId = "Biomart_archive_ens", "Biomart host", ensembl_archive)
    }
  }))
  # pair-wise ------------------------------------------------------------------------------
  org1 <- reactive({
    return(org(Species = input$Species,Ortholog = input$Ortholog))
  })
  ortholog1 <- reactive({
    return(no_org_ID(count = row_count_matrix(),Species = input$Species,Ortholog = input$Ortholog,Biomart_archive=input$Biomart_archive))
  })
  isoform1 <- reactive({
    return(isoform_ID(count = row_count_matrix(),Species = input$Species,Ortholog = input$Ortholog,Biomart_archive=input$Biomart_archive,RNA_type=input$Level_pair))
  })
  org_code1 <- reactive({
    return(org_code(Species = input$Species, Ortholog= input$Ortholog))
  })
  gene_type1 <- reactive({
    return(gene_type(my.symbols=rownames(row_count_matrix()),org=org1(),Species=input$Species,RNA_type=input$Level_pair))
  })
  gene_type1_batch <- reactive({
    return(gene_type(my.symbols=rownames(batch_files()[[1]]),org=org1(),Species=input$Species,RNA_type=input$Level_pair))
  })
  ortholog1_batch <- reactive({
    return(no_org_ID(count = batch_files()[[1]],Species = input$Species,Ortholog = input$Ortholog,Biomart_archive=input$Biomart_archive))
  })
  isoform1_batch <- reactive({
    return(isoform_ID(count = batch_files()[[1]],Species = input$Species,Ortholog = input$Ortholog,Biomart_archive=input$Biomart_archive,RNA_type=input$Level_pair))
  })
  observeEvent(pre_d_row_count_matrix(),({
    if(!is.null(pre_d_row_count_matrix())){
      updateSelectizeInput(session,inputId = "sample_order","Select samples:",
                           choices = colnames(pre_d_row_count_matrix()),selected = colnames(pre_d_row_count_matrix()))
    }
  }))
  output$paired_sample_file <- renderUI({
    if(input$paired_sample == "Yes"){
      fileInput("paired_sample_file","Select a paired-sample file.",accept = c("txt", "csv","xlsx"),
                multiple = FALSE, width = "100%")
    }
  })
  d_paired_sample_file <- reactive({
    if(!is.null(input$sample_order)){
      tmp <- input$paired_sample_file$datapath
      if(is.null(input$paired_sample_file) && input$goButton > 0 )  tmp = example_data_path("data/241011_test_paired.csv", "https://raw.githubusercontent.com/Kan-E/RNAseqChef_imeg/main/data/241011_test_paired.csv")
      print(tmp)
      df <- read_df(tmp = tmp)
      if(!is.null(df)) {
        rownames(df) <- gsub("-",".",rownames(df))
        df2 <- data.frame(sample=rownames(df),pair=df[,1])
        df <- df2 %>% dplyr::filter(sample %in% input$sample_order) %>% dplyr::arrange(factor(sample, levels=input$sample_order))
      }
      return(df)
    }
  })
  output$paired_table <- renderDataTable({
    if(input$paired_sample == "Yes" && !is.null(d_row_count_matrix())){
      if(input$DEG_method == "EBSeq") validate("The EBSeq method is not available for paired-sample analysis. Please use a different method.")
      data.frame(pair=d_paired_sample_file()$pair,row.names = d_paired_sample_file()$sample)
    }
  })
  observeEvent(d_row_count_matrix(), ({
    updateCollapse(session,id =  "input_collapse_panel", open="D_row_count_matrix_panel")
  }))
  d_row_count_matrix <- reactive({
    count <- pre_d_row_count_matrix()
    order <- input$sample_order
    data <- try(safe_ordered_matrix(count, order), silent = TRUE)
    if(inherits(data, "try-error")) validate("")
    return(data)
  })
  row_count_matrix <- reactive({
    withProgress(message = "Importing row count matrix, please wait",{
      if (input$data_file_type == "Row1"){
        tmp <- input$file3$datapath
        if(is.null(input$file3) && input$goButton > 0 )  tmp = example_data_path("data/example1.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example1.txt")
        return(read_df(tmp = tmp))
      }
      if (input$data_file_type == "Row2"){
        tmp <- input$file1$datapath
        if(is.null(input$file1) && input$goButton > 0 )  tmp = example_data_path("data/example2.csv", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example2.csv")
        return(read_df(tmp = tmp))
      }
    })
  })
  metadata <- reactive({
    if (input$data_file_type != "Row2"){
      return(NULL)
    }else{
      tmp <- input$file2$datapath
      if(is.null(input$file2) && input$goButton > 0 )  tmp = example_data_path("data/example3.csv", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example3.csv")
      df <- read_df(tmp = tmp)
      if(!is.null(df)) rownames(df) <- gsub("-",".",rownames(df))
      return(df)
    }
  })
  norm_count_matrix <- reactive({
    if(is.null(input$norm_file1)){
      return(NULL)
    }else{
      if(length(input$norm_file1[, 1]) == 1){
        if(!is.null(input$sample_order)){
        upload <- read_df(input$norm_file1[[1, 'datapath']])
        if (input$data_file_type == "Row2"){
          meta <- anno_rep_meta(metadata())
          if (!is.null(upload) && !is.null(meta)){
            row_t <- t(upload)
            meta <- data.frame(characteristics = meta[,1], row.names = rownames(meta))
            colname <- colnames(meta)
            data <- merge(meta, row_t, by=0, sort = F)
            if(dim(data)[1] == 0) {
              rownames(meta) <- gsub("\\.","-",rownames(meta))
              data <- merge(meta, row_t, by=0, sort = F)
              if(dim(data)[1] == 0) {
                rownames(row_t) <- gsub("\\.","-",rownames(row_t))
                data <- merge(meta, row_t, by=0, sort = F)
                validate("Error: failed to merge count data with metadata. Please check row names of matadata.")
              }
            }
            rownames(data) <- data$characteristics
            data2 <- data[, - which(colnames(data) %in% c("Row.names", colname))]
            data2_t <- t(data2)
            upload <- apply(data2_t, 2, as.numeric)
            rownames(upload) <- rownames(data2_t)
            upload <- as.data.frame(upload)
          }else return(NULL)
        }else upload <- anno_rep(upload)
        order <- input$sample_order
        upload <- try(upload[,order])
        if(length(upload) == 1){
          if(inherits(upload, "try-error")) validate("")
        }
        }else return(NULL)
      }else{
        upload = list()
        for(nr in 1:length(input$norm_file1[, 1])){
          df <- anno_rep(read_df(input$norm_file1[[nr, 'datapath']]))
          upload[gsub("\\..+$", "", input$norm_file1[nr,]$name)] <- list(df)
        }
      } 
      return(upload)
    }
  })
  pre_d_row_count_matrix <- reactive({
    withProgress(message = "Creating defined count matrix, please wait",{
      row <- row_count_matrix()
      if (input$data_file_type == "Row1"){
        if(is.null(row)) {
          return(NULL)
        }else{
          return(anno_rep(row = row))
        }
      }
      if (input$data_file_type == "Row2"){
        meta <- anno_rep_meta(metadata())
        if (is.null(row) || is.null(meta)){
          return(NULL)
        } else {
          row_t <- t(row)
          meta <- data.frame(characteristics = meta[,1], row.names = rownames(meta))
          colname <- colnames(meta)
          data <- merge(meta, row_t, by=0, sort = F)
          if(dim(data)[1] == 0) {
            rownames(meta) <- gsub("\\.","-",rownames(meta))
            data <- merge(meta, row_t, by=0, sort = F)
            if(dim(data)[1] == 0) {
              rownames(row_t) <- gsub("\\.","-",rownames(row_t))
              data <- merge(meta, row_t, by=0, sort = F)
              validate("Error: failed to merge count data with metadata. Please check row names of matadata.")
            }
          }
          rownames(data) <- data$characteristics
          data2 <- data[, - which(colnames(data) %in% c("Row.names", colname))]
          data2_t <- t(data2)
          data3 <- apply(data2_t, 2, as.numeric)
          rownames(data3) <- rownames(data2_t)
          return(data3)
        }
      }
    })
  })
  
  #-----------------
  output$not_cond2 <- renderText({
    count <- d_row_count_matrix()
    if(!is.null(count)){
      collist <- gsub("\\_.+$", "", colnames(count))
      if(length(unique(collist)) != 2) print("Uploaded count data is in an inappropriate format. Please refer to the RNAseqChef manual for guidance and make the necessary corrections.")
    }
  })
  output$not_cond2_pair <- renderText({
    count <- d_row_count_matrix()
    if(!is.null(count)){
      collist <- gsub("\\_.+$", "", colnames(count))
      if(length(unique(collist)) != 2) print(paste0("Your count file contains data for ", length(unique(collist))," conditions. Please select samples to narrow it down to 2 conditions."))
    }
  })
  # pair-wise DEG ------------------------------------------------------------------------------
  observeEvent(input$DEG_method,({
    if(input$DEG_method == "edgeR") {
      updateNumericInput(session,"pair_prefilter", "Minimum count required for at least some samples", value = 10)
      updateNumericInput(session,"pair_prefilterTotal", "Minimum total count required", value = 15)
    }
    if(input$DEG_method == "limma") {
      updateNumericInput(session,"pair_prefilter", "Minimum count required for at least some samples", value = 0)
      updateNumericInput(session,"pair_prefilterTotal", "Minimum total count required", value = 0)
    }
  }))
  pair_analysis_requested <- reactiveVal(FALSE)

  observeEvent(input$pair_tabs, {
    if(!is.null(input$pair_tabs) && input$pair_tabs != "pair_input_data_tab"){
      pair_analysis_requested(TRUE)
    }else{
      pair_analysis_requested(FALSE)
    }
  }, ignoreInit = FALSE)

    dds <- reactive({
    if(!isTRUE(pair_analysis_requested())){
      return(NULL)
    }
    count <- d_row_count_matrix()
    file_name <- gsub("\\..+$", "", input$file1)
    collist <- gsub("\\_.+$", "", colnames(count))
    if(length(unique(collist)) == 2){
      if (input$DEG_method == "DESeq2") {
        withProgress(message = "DESeq2",{
          if(input$paired_sample == "Yes"){
            if(!is.null(d_paired_sample_file())) {
              group <- data.frame(con = factor(collist),pair=factor(d_paired_sample_file()$pair))
              dds<- DESeqDataSetFromMatrix(countData = round(count),colData = group, design = ~ pair + con)
            }else validate("Please upload a paired-sample file")
          }else{
            group <- data.frame(con = factor(collist))
            dds<- DESeqDataSetFromMatrix(countData = round(count),colData = group, design = ~ con)
          }
          dds$con <- factor(dds$con, levels = unique(collist))
          dds <- DESeq(dds)
          incProgress(1)
        })
      }
      if (input$DEG_method == "edgeR") {
        withProgress(message = "edgeR",{
          group <- factor(collist)
          dds <- DGEList(counts = count, group = group)
          if(input$pair_prefilterON == "ON"){
          keep <- filterByExpr(dds,min.count=input$pair_prefilter,min.total.count=input$pair_prefilterTotal)
          dds = dds[keep, , keep.lib.sizes=FALSE]
          }
          dds <- calcNormFactors(dds)
          if(input$paired_sample == "Yes"){
            if(!is.null(d_paired_sample_file())) {
              con = factor(collist)
              pair=factor(d_paired_sample_file()$pair)
              design <- model.matrix(~ pair + con)
              dds <- estimateGLMCommonDisp(dds, design)
              dds <- estimateGLMTrendedDisp(dds, design)
              dds <- estimateGLMTagwiseDisp(dds, design)
            }else validate("Please upload a paired-sample file")
          }else{
          dds <- estimateCommonDisp(dds)
          dds <- estimateTagwiseDisp(dds)
          }
          incProgress(1)
        })
      }
      return(dds)
    }
  })
  
  deg_result_raw <- reactive({
    if(!isTRUE(pair_analysis_requested())){
      return(NULL)
    }
    if(is.null(d_row_count_matrix())){
      return(NULL)
    }else{
      count <- d_row_count_matrix()
      collist <- gsub("\\_.+$", "", colnames(count))
      if(length(unique(collist)) == 2){
        if (input$DEG_method == "DESeq2") {
          dds <- dds()
          contrast <- c("con", unique(collist))
          res <- DESeq2::results(dds, contrast = contrast)
          if(input$FDR_method == "IHW") {
            ihw_res <- IHW::ihw(pvalue ~ baseMean, data = as.data.frame(res), alpha = 0.1)
            res$padj <- IHW::as.data.frame(ihw_res)$adj_pvalue
          }
          if(input$FDR_method == "Qvalue") {
            qvalue <- qvalue::qvalue(res$pvalue)
            res$padj <- qvalue$qvalues
          }
        }
        if(input$DEG_method == "edgeR"){
          dds <- dds()
          group <- factor(collist)
          if(input$paired_sample == "Yes") {
            con = factor(collist)
            pair=factor(d_paired_sample_file()$pair)
            design <- model.matrix(~ pair + con)
            fit <- glmFit(dds, design)
            result <- glmLRT(fit, coef = dim(design)[2])
          }else result <- exactTest(dds, pair = c(unique(group)[2],unique(group)[1]))
          res <- as.data.frame(topTags(result, n = nrow(count)))
          if(input$paired_sample == "Yes") res <- res[,-3]
          qvalue <- qvalue::qvalue(res$PValue)
          res$padj <- qvalue$qvalues
          ihw_res <- try(IHW::ihw(PValue ~ 2^logCPM, data = res, alpha = 0.1))
          if(inherits(ihw_res, "try-error")){
            res$ihw_padj <- NA
          }else{
            ihw_res_df <- IHW::as.data.frame(ihw_res)
            res$ihw_padj <- ihw_res_df$adj_pvalue
          }
          if(input$FDR_method == "BH"){label <- c("log2FoldChange", "log2CPM", "PValue","padj", "Qvalue", "IHW_FDR")}
          if(input$FDR_method == "Qvalue"){label <- c("log2FoldChange", "log2CPM", "PValue","BH_FDR", "padj", "IHW_FDR")}
          if(input$FDR_method == "IHW"){label <- c("log2FoldChange", "log2CPM", "PValue","BH_FDR", "Qvalue", "padj")}
          colnames(res) <- label
        }
        if(input$DEG_method == "limma"){
          collist <- gsub(" ", ".", collist)
          collist <- gsub("\\+", ".", collist)
          group <- factor(collist)

          if(input$pair_prefilterON == "ON"){
            dds <- DGEList(counts = count, group = group)
            keep <- filterByExpr(dds,min.count=input$pair_prefilter,min.total.count=input$pair_prefilterTotal)
            dds = dds[keep, , keep.lib.sizes=FALSE]
            count <- dds$counts
          }

          count <- log(count + 1,2)
          eset = new("ExpressionSet", exprs=as.matrix(count))

          if(input$paired_sample == "Yes"){
            if(!is.null(d_paired_sample_file())) {
              con = factor(collist)
              pair=factor(d_paired_sample_file()$pair,levels = unique(d_paired_sample_file()$pair))
              design <- model.matrix(~ pair + con)
              fit <- lmFit(eset, design)
              if(input$limma_trend == TRUE) fit2 <- eBayes(fit,trend = TRUE,robust = input$regression_mode)
              if(input$limma_trend == FALSE) fit2 <- eBayes(fit,trend = FALSE,robust = input$regression_mode)
              res <- topTable(fit2,coef = dim(design)[2], number = 1e12)
            }
          }else{
            design <- model.matrix(~0+collist)
            colnames(design) <- gsub("^collist","",colnames(design))
            fit <- lmFit(eset, design)
            comparisons <-  paste(unique(collist)[1],"-",unique(collist)[2],sep="")
            cont.matrix <- makeContrasts(contrasts=comparisons, levels=design)
            fit <- contrasts.fit(fit, cont.matrix)
            if(input$limma_trend == TRUE) fit2 <- eBayes(fit,trend = TRUE,robust = input$regression_mode)
            if(input$limma_trend == FALSE) fit2 <- eBayes(fit,trend = FALSE,robust = input$regression_mode)
            res =topTable(fit2, number = 1e12)
          }
          if(input$cutoff_limma == "fdr"){
            colnames(res) <- c("log2FoldChange","baseMean","t","pval","padj","B")
          }else{
            colnames(res) <- c("log2FoldChange","baseMean","t","padj","padj2","B")
          }
          res$baseMean <- 2^res$baseMean
        }
        if(input$DEG_method == "EBSeq"){
          withProgress(message = "EBSeq",{
            count <- data.matrix(count)
            vec <- c()
            for (i in 1:length(unique(collist))) {
              num <- length(collist[collist == unique(collist)[i]])
              vec <- c(vec, num)
            }
            ngvector <- NULL
            conditions <- as.factor(rep(paste("C", 1:length(unique(collist)), sep=""), times = vec))
            Sizes <- MedianNorm(count)
            EBOut <- NULL
            EBOut <- EBTest(Data = count, NgVector = ngvector, Conditions = conditions, sizeFactors = Sizes, maxround = 5)
            stopifnot(!is.null(EBOut))
            PP <- as.data.frame(GetPPMat(EBOut))
            fc_res <- PostFC(EBOut)
            results <- cbind(PP, fc_res$PostFC, fc_res$RealFC,EBOut$Mean[rownames(PP),1], EBOut$Mean[rownames(PP),2])
            colnames(results) <- c("PPEE", "PPDE", "PostFC", "RealFC","C1Mean","C2Mean")
            res <- results[order(results[,"PPDE"], decreasing = TRUE),]
            incProgress(1)
          })
        }
        as.data.frame(res)
      }
    }
  })

  deg_result <- reactive({
    res <- deg_result_raw()
    if(is.null(res)){
      return(NULL)
    }
    if(input$Species != "not selected"){
      res <- ensembl2symbol(gene_type=gene_type1(),data = res, input$Species,Ortholog=ortholog1(),Isoform=isoform1(),org = org1())
    }
    return(res)
  })
  
  deg_norm_count_raw <- reactive({
    if(!isTRUE(pair_analysis_requested())){
      return(NULL)
    }
    if(is.null(d_row_count_matrix())){
      return(NULL)
    }else{
      count <- d_row_count_matrix()
      collist <- gsub("\\_.+$", "", colnames(count))
      if(length(unique(collist)) == 2){
        if(!is.null(norm_count_matrix())){
          return(norm_count_matrix())
        }else {
          if (input$DEG_method == "DESeq2") {
            dds <- dds()
            normalized_counts <- counts(dds, normalized=TRUE)
          }
          if (input$DEG_method == "edgeR") {
            dds <- dds()
            if(input$paired_sample == "Yes"){
              if(!is.null(d_paired_sample_file())) {
                dds <- estimateCommonDisp(dds)
                dds <- estimateTagwiseDisp(dds)
              }else validate("Please upload a paired-sample file")
            }
            normalized_counts <- t(t(dds$pseudo.counts)*(dds$samples$norm.factors))
          }
          if (input$DEG_method == "limma") {
            normalized_counts <- count
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
          return(normalized_counts)
        }
      }
    }
  })

  deg_norm_count <- reactive({
    normalized_counts <- deg_norm_count_raw()
    if(is.null(normalized_counts)){
      return(NULL)
    }
    if(input$Species != "not selected"){
      normalized_counts <- ensembl2symbol(gene_type=gene_type1(),data = normalized_counts,
                                          Species=input$Species,Ortholog=ortholog1(),Isoform=isoform1(),org = org1())
    }
    return(normalized_counts)
  })
  
  observeEvent(input$goButton,({
    updateSelectInput(session,inputId = "Species","Species",species_list, selected = "Homo sapiens")
  }))
  
  observeEvent(input$file1, ({
    updateCollapse(session,id =  "input_collapse_panel", open="Row_count_matrix_panel")
  }))
  observeEvent(input$file2, ({
    updateCollapse(session,id =  "input_collapse_panel", open="Metadata_panel")
  }))
  output$Row_count_matrix <- DT::renderDataTable({
    if(input$data_file_type == "Row11"){
      uploaded_files = names(batch_files())
      as.data.frame(uploaded_files)
    }else{
      row_count_matrix()
    }
  })
  output$Metadata <- DT::renderDataTable({
    metadata()
  })
  output$D_Row_count_matrix <- DT::renderDataTable({
    d_row_count_matrix()
  })
  deg_result_limma <- reactive({
    if (!is.null(deg_result()) && input$DEG_method == "limma") {
      res <- deg_result()
      return(rename_limma_display_columns(res))
    }
  })
  output$DEG_result <- DT::renderDataTable({
    if (input$DEG_method == "limma") {
      deg_result_limma()
    }else{
      deg_result()
    }
  })
  output$download_pair_d_row_count = downloadHandler(
    filename = function() {
      if (input$data_file_type == "Row1"){
        paste(gsub("\\..+$", "", input$file3), paste0(input$FDR_method,".txt"), sep ="-")
      }else{
        paste(gsub("\\..+$", "", input$file1), paste0(gsub("\\..+$", "", input$file2),".txt"), sep ="-")
      }},
    content = function(file){write.table(d_row_count_matrix(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  output$download_pair_DEG_result = downloadHandler(
    filename = function() {
      if (input$data_file_type == "Row1"){
        paste(paste(gsub("\\..+$", "", input$file3), input$DEG_method, sep = "-"), paste0(input$FDR_method,".txt"), sep ="-")
      }else{
        paste(paste(gsub("\\..+$", "", input$file1), input$DEG_method, sep = "-"), paste0(input$FDR_method,".txt"), sep ="-")
      }},
    content = function(file){write.table(deg_result(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  output$download_pair_norm_count = downloadHandler(
    filename = function() {
      if (input$data_file_type == "Row1"){
        paste(paste(gsub("\\..+$", "", input$file3), input$DEG_method, sep = "-"), "normalized.txt", sep ="-")
      }else{
        paste(paste(gsub("\\..+$", "", input$file1), input$DEG_method, sep = "-"), "normalized.txt", sep ="-")
      }},
    content = function(file){write.table(deg_norm_count(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  download_pair_overview_dir <-reactive({
    if (input$data_file_type == "Row1"){
      dir_name <- paste(paste(gsub("\\..+$", "", input$file3), input$DEG_method, sep = "-"), input$FDR_method , sep ="-")
    }else{
      dir_name <- paste(paste(gsub("\\..+$", "", input$file1), input$DEG_method, sep = "-"), input$FDR_method , sep ="-")
    }
    dir_name <- paste0(dir_name, paste0("_fc", input$fc))
    dir_name <- paste0(dir_name, paste0("_fdr", input$fdr))
    dir_name <- paste0(dir_name, paste0("_basemean", input$basemean))
    return(dir_name)
  })
  #pair-wise DEG vis---------------------------
  pair_gene_id_map <- reactive({
    data <- deg_result()
    if(is.null(data) || input$Species == "not selected"){
      return(NULL)
    }

    my.symbols <- rownames(data)
    if(length(my.symbols) == 0){
      return(NULL)
    }

    if(gene_type1() != "SYMBOL"){
      if(sum(is.element(no_orgDb, input$Species)) == 1){
        gene_IDs <- ortholog1()
      }else if(gene_type1() == "isoform"){
        gene_IDs <- isoform1()
      }else{
        key <- gene_primary_key(my.symbols)
        if(org1()$packageName == "org.Sc.sgd.db") SYMBOL <- "GENENAME" else SYMBOL <- "SYMBOL"
        gene_IDs <- AnnotationDbi::select(org1(),keys = my.symbols,
                                          keytype = key,
                                          columns = c(key,SYMBOL, "ENTREZID"))
      }
      colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
    }else{
      if(sum(is.element(no_orgDb, input$Species)) == 1){
        gene_IDs <- ortholog1()[,-1]
      }else if(gene_type1() == "isoform"){
        gene_IDs <- isoform1()[,-1]
      }else{
        if(org1()$packageName == "org.Sc.sgd.db") SYMBOL <- "GENENAME" else SYMBOL <- "SYMBOL"
        gene_IDs <- AnnotationDbi::select(org1(),keys = my.symbols,
                                          keytype = SYMBOL,
                                          columns = c(SYMBOL, "ENTREZID"))
      }
      colnames(gene_IDs) <- c("Row.names", "ENTREZID")
    }

    gene_IDs %>% distinct(Row.names, .keep_all = T)
  })

  gene_ID_pair <- reactive({
    gene_IDs <- pair_gene_id_map()
    if(is.null(gene_IDs) || gene_type1() == "SYMBOL" || !"SYMBOL" %in% colnames(gene_IDs)){
      return(NULL)
    }
    dplyr::select(gene_IDs, Row.names, SYMBOL)
  })
  
  
  data_degcount <- reactive({
    data <- deg_result()
    count <- deg_norm_count()
    if(is.null(count) || is.null(data)){
      return(NULL)
    }else{
      if(gene_type1() != "SYMBOL"){
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
      if (Type == "limma"){
        data$log2FoldChange <- -1 * data$log2FoldChange
      }
      if(input$Species != "not selected"){
        gene_IDs <- pair_gene_id_map()
        if(!is.null(gene_IDs)){
          data <- merge(data, gene_IDs, by="Row.names")
          if(gene_type1() != "SYMBOL" && "SYMBOL" %in% colnames(gene_IDs)){
            data$Unique_ID <- paste(data$SYMBOL,data$Row.names, sep = "\n- ")
          }
        }
      }
      return(data)
    }
  })
  
  data_degcount2 <- reactive({
    data <- data_degcount()
    count <- deg_norm_count()
    if(is.null(count) || is.null(data)){
      return(NULL)
    }else{
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
        data2$group <- "Up"
        data2$group[data2$log2FoldChange < 0] <- "Down"
        data3 <- dplyr::filter(data2, abs(data2$padj) < input$fdr)
        return(data3)
      }else{return(NULL)}
    }
  })
  
  data_degcount_up <- reactive({
    data <- data_degcount()
    count <- deg_norm_count()
    if(is.null(count) || is.null(data)){
      return(NULL)
    }else{
      collist <- factor(gsub("\\_.+$", "", colnames(count)))
      vec <- c()
      for (i in 1:length(unique(collist))) {
        num <- length(collist[collist == unique(collist)[i]])
        vec <- c(vec, num)
      }
      Cond_1 <- vec[1]
      Cond_2 <- vec[2]
      data2 <- data_degcount2()
      up_all <- dplyr::filter(data2, log2FoldChange > 0)
      rownames(up_all) <- up_all$Row.names
      up_all <- up_all[,8:(7 + Cond_1 + Cond_2)]
      if(input$Species != "not selected"){
        if(gene_type1() != "SYMBOL"){
          up_all <- merge(up_all, gene_ID_pair(), by=0)
          rownames(up_all) <- up_all$Row.names
          up_all <- up_all[,-1]
          up_all <- up_all[, - which(colnames(up_all) == "Row.names.y")]
        }
      }
      return(up_all)
    }
  })
  
  data_degcount_down <- reactive({
    data <- data_degcount()
    count <- deg_norm_count()
    if(is.null(count) || is.null(data)){
      return(NULL)
    }else{
      collist <- factor(gsub("\\_.+$", "", colnames(count)))
      vec <- c()
      for (i in 1:length(unique(collist))) {
        num <- length(collist[collist == unique(collist)[i]])
        vec <- c(vec, num)
      }
      Cond_1 <- vec[1]
      Cond_2 <- vec[2]
      data2 <- data_degcount2()
      down_all <- dplyr::filter(data2, log2FoldChange < 0)
      rownames(down_all) <- down_all$Row.names
      down_all <- down_all[,8:(7 + Cond_1 + Cond_2)]
      if(input$Species != "not selected"){
        if(gene_type1() != "SYMBOL"){
          down_all <- merge(down_all, gene_ID_pair(), by=0)
          rownames(down_all) <- down_all$Row.names
          down_all <- down_all[,-1]
          down_all <- down_all[, - which(colnames(down_all) == "Row.names.y")]
        }
      }
      return(down_all)
    }
  })
  
  output$pair_deg_up <- DT::renderDataTable({
    data_degcount_up()
  })
  output$pair_deg_down <- DT::renderDataTable({
    data_degcount_down()
  })
  output$download_pair_deg_count_up = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), "_DEG_count_up.txt")
    },
    content = function(file) {write.table(data_degcount_up(), file, quote = F, row.names = T, col.names=NA, sep = "\t")})
  
  output$download_pair_deg_count_down = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), "_DEG_count_down.txt")
    },
    content = function(file) {write.table(data_degcount_down(), file, quote = F, row.names = T, col.names=NA, sep = "\t")})
  
  # pair-wise MA ------------------------------------------------------------------------------
  output$MA <- renderPlot({
    withProgress(message = "MA plot and heatmap",{
      if(is.null(d_row_count_matrix()) || is.null(ma_heatmap_plot())){
        return(NULL)
      }else{
        plot(ma_heatmap_plot())
      }
      incProgress(1)
    })
  })
  ma_heatmap_plot <- reactive({
    data <- data_degcount()
    count <- deg_norm_count()
    if(!is.null(count)){
      if(gene_type1() != "SYMBOL"){
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
      if(gene_type1() != "SYMBOL"){
        if(input$Species != "not selected"){
          genenames <- as.vector(data$SYMBOL)
        }else{ genenames=as.vector(data$Row.names) }
      }else{
        genenames <- as.vector(data$Row.names)
      }
      if(input$FDR_method == "edgeR") {
        xlab <- "LogCPM"
      }else{
        xlab <- "Log2 mean expression"
      }
      m1 <- as.grob(ggmaplot(data, fdr = input$fdr, fc = input$fc, size = 0.1,
                             palette = c("#B31B21", "#1465AC", "darkgray"),
                             genenames = genenames,
                             legend = "top", top = 20,
                             font.label = c("bold.italic", 10),font.legend = "bold",
                             font.main = c("bold", 15),xlab = xlab,
                             ggtheme = ggplot2::theme_minimal(base_size = 15),
                             select.top.method = "fc"))
      data2 <- data_degcount2()
      if(is.null(data2)){
        ht <- NULL
      }else{
        data.z <- genefilter::genescale(data2[,8:(7 + Cond_1 + Cond_2)], axis=1, method="Z")
        ht <- as.grob(GOIheatmap(data.z,show_row_names = FALSE))
      }
      p <- plot_grid(m1, ht, rel_widths = c(2, 1))
      return(p)
    }
  })
  
  output$download_pair_MA = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), "_MAplot-heatmap.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 4
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(ma_heatmap_plot())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  
  # pair-wise volcano--------------
  normalize_goi_choices <- function(x) {
    unique(stats::na.omit(as.character(x)))
  }

  goi_selectize_options <- list(
    delimiter = " ",
    create = TRUE,
    plugins = list("remove_button"),
    persist = FALSE
  )

  GOI_list <- reactive({
    data <- data_degcount()
    if(is.null(data) || is.null(d_row_count_matrix())){
      return(NULL)
    }else{
      if(input$Level_pair != "gene_level" && "padj" %in% colnames(data)) {
        data <- data %>% dplyr::filter(!is.na(padj))
      }
      if(is.null(input$GOI_color_type)) validate("")
      if(input$GOI_color_type != "default" && !is.null(pair_pathway_color_gene())){
        data <- data %>% dplyr::filter(Row.names %in% rownames(pair_pathway_color_gene()))
      }
      if(gene_type1() != "SYMBOL"){
        if(input$Species != "not selected"){
          GOI <- data$Unique_ID
        }else GOI <- data$Row.names
      }else{
        if(input$Species != "not selected"){
          GOI <- data$Row.names
        }else GOI <- data$Row.names
      }
      return(normalize_goi_choices(GOI))
    }
  })
  
  observeEvent(GOI_list(), {
    choices <- GOI_list()
    if(is.null(choices)){
      choices <- character(0)
    }
    selected <- intersect(isolate(input$GOI), choices)
    updateSelectizeInput(session, "GOI", choices = choices, selected = selected,
                         options = goi_selectize_options, server = TRUE)
  }, ignoreNULL = FALSE)

  observeEvent(input$GOIreset_pair, {
    choices <- GOI_list()
    if(is.null(choices)){
      choices <- character(0)
    }
    updateSelectizeInput(session, "GOI", choices = choices,
                         selected = character(0),
                         options = goi_selectize_options, server = TRUE)
  })
  output$GOI_color_type <- renderUI({
    if(input$Species == "not selected"){
      radioButtons("GOI_color_type","Filter", c("All genes" = "default"),selected="default")
    }else{
      radioButtons("GOI_color_type","Filter", c("All genes" = "default","Pathway of interest"="pathway"),selected="default")
    }
  })
  observeEvent(input$Species, ({
    if(input$Species == "not selected"){
      updateSelectInput(session, "GOI_color_pathway1","Select a gene set for gene extraction","")
    }else  if(input$Species != "Xenopus laevis" && input$Ortholog != "Arabidopsis thaliana" && input$Species != "Arabidopsis thaliana"){
      updateSelectInput(session, "GOI_color_pathway1","Select a gene set for gene extraction",gene_set_list) 
    }else {
      updateSelectInput(session, "GOI_color_pathway1","Select a gene set for gene extraction",c("KEGG", "GO biological process", 
                                                                                                "GO cellular component","GO molecular function")) 
    }
  }))
  GOI_color_pathway_list <- reactive({
    list <-pair_pathway_color_gene_list()
    if(!is.null(input$GOI_color_pathway1)){
    if(input$Species == "Xenopus laevis" || input$Ortholog == "Arabidopsis thaliana" || input$Species == "Arabidopsis thaliana"){
      print(head(list))
      if(input$GOI_color_pathway1 == "GO biological process") list <- list %>% dplyr::filter(Ontology == "BP")
      if(input$GOI_color_pathway1 == "GO cellular component") list <- list %>% dplyr::filter(Ontology == "CC")
      if(input$GOI_color_pathway1 == "GO molecular function") list <- list %>% dplyr::filter(Ontology == "MF") 
    }}
    list <- unique(list$gs_name)
    return(list)
  })
  observeEvent(input$GOI_color_pathway1, ({
    if(is.null(input$GOI_color_pathway1)){
      updateSelectInput(session, "GOI_color_pathway2","","")
    }else{
      updateSelectInput(session, "GOI_color_pathway2","",GOI_color_pathway_list())
    }
  }))
  output$GOIreset_pair <- renderUI({
    actionButton("GOIreset_pair", "GOI reset")
  })
  
  
  
  output$volcano_x <- renderUI({
    if(!is.null(data_degcount())){
      data <- as.data.frame(data_degcount())
      if(input$GOI_plot_select == "Volcano plot"){
        min <- floor(min(data$log2FoldChange))
        max <- ceiling(max(data$log2FoldChange))
        sliderInput("xrange","X_axis range:",min = min-1,
                    max=max+1, step = 0.5,
                    value = c(min, max))
      }else{
        max <- ceiling(max(log(data$baseMean,2)))
        sliderInput("xrange","X_axis range:",min = 0,
                    max=max+1, step = 0.5, value = max)
      }
    }
  })
  output$volcano_y <- renderUI({
    if(!is.null(data_degcount())){
      data <- as.data.frame(data_degcount())
      if(input$GOI_plot_select == "Volcano plot"){
        data$padj[data$padj == 0] <- 10^(-300)
        data <- na.omit(data)
        max <- ceiling(max(-log10(data$padj)))
        sliderInput("yrange","Y_axis range:",min = 0, max= max+1, step = 1,
                    value = max)
      }else{
        min <- floor(min(data$log2FoldChange))
        max <- ceiling(max(data$log2FoldChange))
        sliderInput("yrange","Y_axis range:",min = min-1,
                    max=max+1, step = 0.5,
                    value = c(min, max))
      }
    }
  })
  pair_pathway_color_gene_list <- reactive({
    genes <- GeneList_for_enrichment(Species = input$Species, Ortholog = input$Ortholog,
                                     Gene_set=input$GOI_color_pathway1, org = org1(), 
                                     Biomart_archive=input$Biomart_archive,gene_type=gene_type1())
    return(genes)
  })
  pair_pathway_color_gene <- reactive({
    ##extract pathway genes
    if(is.null(input$GOI_color_pathway1) || is.null(input$GOI_color_pathway2) || 
       is.null(gene_type1()) || is.null(org1()) || input$GOI_color_pathway1 == "" || 
       input$GOI_color_pathway2 == "" || !input$GOI_color_pathway2 %in% GOI_color_pathway_list()) validate("")
    genes <- pair_pathway_color_gene_list()
    genes <- try(dplyr::filter(genes, gs_name == input$GOI_color_pathway2))
    if(length(genes) == 1) if(class(genes)=="try-error") validate("")
    
    my.symbols <- as.character(genes$entrez_gene)
    if(gene_type1() == "non-model organism"){
      gene_IDs <-  try(dplyr::filter(ortholog1(), ENTREZID %in% my.symbols))
      if(length(gene_IDs) == 1) if(class(gene_IDs)=="try-error") validate("")
      df <- data.frame(gene = gene_IDs$ENSEMBL, row.names = gene_IDs$ENSEMBL)
    }else if(gene_type1() == "isoform"){
      gene_IDs <-  try(dplyr::filter(isoform1(), ENTREZID %in% my.symbols))
      if(length(gene_IDs) == 1) if(class(gene_IDs)=="try-error") validate("")
      df <- data.frame(gene = gene_IDs$Transcript_ID, row.names = gene_IDs$Transcript_ID)
    }else{
      if(gene_type1() == "SYMBOL") columns <- c("ENTREZID", sgd_symbol_column(org1())) else columns <- c("ENTREZID","ENSEMBL")
      if(org1()$packageName == "org.At.tair.db") {
        if(gene_type1() == "SYMBOL"){
          gene_IDs <- AnnotationDbi::select(org1(), keys = my.symbols,
                                            keytype = "TAIR",
                                            columns = c("TAIR","SYMBOL"))
          colnames(gene_IDs) <- c("TAIR","GeneID")
        }else{
          gene_IDs <- data.frame(TAIR = my.symbols, GeneID = my.symbols)
        }
      } else {
        gene_IDs <- AnnotationDbi::select(org1(), keys = my.symbols,
                                          keytype = "ENTREZID",
                                          columns = columns)
        colnames(gene_IDs) <- c("entrezid","GeneID")
        gene_IDs <- na.omit(gene_IDs) 
      }
      gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
      df <- data.frame(gene = gene_IDs$GeneID, row.names = gene_IDs$GeneID)
      if(dim(df)[1] == 0) validate("No filtered genes.")
    }
    return(df)
  })
  output$uniqueID_cut <- renderUI({
    if(gene_type1() != "SYMBOL" && input$Species != "not selected") 
      radioButtons("uniqueID_cut","Short unique ID",c("ON"=TRUE,"OFF"=FALSE),selected = TRUE)
  })
  pair_volcano <- reactive({
    if(!is.null(input$xrange) && !is.null(input$yrange)){
      data <- as.data.frame(data_degcount())
      if(dim(brush_info())[1]!=0){
        if(gene_type1() != "SYMBOL"){
          if(input$Species != "not selected"){
            label_data <- brush_info()$Unique_ID
          }else{
            label_data <- brush_info()$Row.names
          }
        }else{
          label_data <- brush_info()$Row.names
        }
      }else{
        if(!is.null(input$GOI)){
          label_data <- input$GOI
        }else label_data <- NULL
      }
      data$padj[data$padj == 0] <- 10^(-300)
      if(input$GOI_color_type == "default"){
        data$color <- "NS"
        data$color[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr] <- "down"
        data$color[data$log2FoldChange > log2(input$fc) & data$padj < input$fdr] <- "up"
        if(!is.null(label_data)) {
          Color <- c(down = "blue", GOI = "green", NS = "darkgray", up = "red")
          for(name in label_data){
            if(gene_type1() != "SYMBOL"){
              if(input$Species != "not selected"){
                data$color[data$Unique_ID == name] <- "GOI"
              }else{
                data$color[data$Row.names == name] <- "GOI"
              }
            }else{
              data$color[data$Row.names == name] <- "GOI"
            }
          }
          legend_levels <- c("down", "GOI", "NS", "up")
          data$color <- factor(data$color, levels = legend_levels)
        }else{
          Color <- c(down = "blue", NS = "darkgray", up = "red")
          legend_levels <- c("down", "NS", "up")
          data$color <- factor(data$color, levels = legend_levels)
        }
      }else{
        df <- pair_pathway_color_gene()
        data$color <- "others"
        Color <- c(others = "lightgray", down = "blue", NS = "gray1", up = "red")
        for(name in rownames(df)){
          data$color[data$Row.names == name] <- "NS"
          data$color[data$Row.names == name & data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr] <- "down"
          data$color[data$Row.names == name & data$log2FoldChange > log2(input$fc) & data$padj < input$fdr] <- "up"
        }
        data$color <- factor(data$color, levels = c("others", "NS","down","up"))
        
        if(!is.null(label_data)) {
          Color <- c(others = "lightgray", down = "blue", GOI = "green", NS = "gray1", up = "red")
          for(name2 in label_data){
            if(gsub(".+\\s", "", name2) %in% rownames(df)){
              if(gene_type1() != "SYMBOL"){
                if(input$Species != "not selected"){
                  data <- data %>% dplyr::mutate(color=if_else(Unique_ID==name2, "GOI", color))
                }else{
                  data <- data %>% dplyr::mutate(color=if_else(Row.names==name2, "GOI", color))
                }
              }else{
                data <- data %>% dplyr::mutate(color=if_else(Row.names==name2, "GOI", color))
              }
            }
          }
          legend_levels <- c("others", "down", "GOI", "NS", "up")
          data$color <- factor(data$color, levels = legend_levels)
        }else{
          Color <- c(others = "lightgray", down = "blue", NS = "gray1", up = "red")
          legend_levels <- c("others", "down", "NS", "up")
          data$color <- factor(data$color, levels = legend_levels)
        }
        
        
        
      }
      legend_breaks <- legend_levels[legend_levels %in% unique(as.character(stats::na.omit(data$color)))]
      data$minusLog10padj<--log10(data$padj)
      lab_y <- "-log10(padj)"
      if(input$DEG_method == "limma"){
        if(input$cutoff_limma == "pval") lab_y <- "-log10(pval)"
      }
      if(input$GOI_plot_select == "Volcano plot"){
        v <- try(ggplot(data, aes(x = log2FoldChange, y = minusLog10padj)) + 
                   geom_vline(xintercept = c(-log2(input$fc), log2(input$fc)), linetype = c(2, 2), color = c("black", "black")) +
                   geom_hline(yintercept = c(-log10(input$fdr)), linetype = 2, color = c("black")) +
                   xlab("log2 fold change") + ylab(lab_y) +
                   xlim(input$xrange)+
                   ylim(c(0, input$yrange)))
      }else{
        data$log2baseMean <- log(data$baseMean,2)
        v <- try(ggplot(data, aes(x = log2baseMean, y = log2FoldChange)) +
                   geom_hline(yintercept = c(-log2(input$fc), 0 ,log2(input$fc)), linetype = c(2,1,2), color = c("black", "black","black")) +
                   xlab("log2 mean expression") + ylab("log2 fold change") +
                   xlim(c(0, input$xrange))+
                   ylim(c(input$yrange)))
      }
      if(inherits(v, "try-error")) validate("")
      v <- v + ggrastr::geom_point_rast(aes(color = color),size = 0.4)  +
        theme_bw()+ scale_color_manual(values = Color, breaks = legend_breaks)+
        theme(legend.position = "top" , legend.title = element_blank(),
              axis.text.x= ggplot2::element_text(size = 12),
              axis.text.y= ggplot2::element_text(size = 12),
              text = ggplot2::element_text(size = 12),
              title = ggplot2::element_text(size = 12))+
        guides(color = guide_legend(override.aes = list(size=2))) +guides(shape = guide_legend(override.aes = list(size=4))) 
      if(input$GOI_color_type == "pathway") {
        v <- v + geom_point(data=dplyr::filter(data, color == "NS"),color="gray1", size=1)
        v <- v + geom_point(data=dplyr::filter(data, color == "up"),color="red", size=1)
        v <- v + geom_point(data=dplyr::filter(data, color == "down"),color="blue", size=1)
      }
      if(!is.null(label_data)) {
        if(gene_type1() != "SYMBOL"){
          if(input$Species != "not selected"){
            if(input$uniqueID_cut) {
              id_list <- gsub("\\\n.+$", "", data$Unique_ID)
              dup_list <- unique(id_list[duplicated(id_list)])
              for(i in 1:length(data$Unique_ID)){
                if(gsub("\\\n.+$", "", data$Unique_ID[i]) == "NA") {
                  data$Unique_ID[i] <- gsub(".+\\s", "", data$Unique_ID[i])
                }else if(! gsub("\\\n.+$", "", data$Unique_ID[i]) %in% dup_list) {
                  data$Unique_ID[i] <- gsub("\\\n.+$", "", data$Unique_ID[i])
                }
              }
            }
            v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
            v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Unique_ID),
                                              box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"),
                                              fontface = "bold.italic",
                                              bg.color = grDevices::adjustcolor("white", alpha.f = 0.6), bg.r = 0.15)
          }else{
            v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
            v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Row.names),
                                              box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"),
                                              fontface = "bold.italic",
                                              bg.color = grDevices::adjustcolor("white", alpha.f = 0.6), bg.r = 0.15)
          }
        }else{
          v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
          v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Row.names),
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"),
                                            fontface = "bold.italic",
                                            bg.color = grDevices::adjustcolor("white", alpha.f = 0.6), bg.r = 0.15)
        }
      }
      return(v)
    }else return(NULL)
  })
  
  
  brush_info <- reactive({
    if(!is.null(input$xrange) && !is.null(data_degcount())){
      data <- as.data.frame(data_degcount())
      if(input$GOI_plot_select == "Volcano plot"){
        data$padj[data$padj == 0] <- 10^(-300)
        data$minusLog10padj<--log10(data$padj)
        brush <- brushedPoints(data, input$plot1_brush,xvar = "log2FoldChange",yvar="minusLog10padj")
      }else{
        data$log2baseMean <- log(data$baseMean,2)
        brush <- brushedPoints(data, input$plot1_brush,xvar = "log2baseMean",yvar="log2FoldChange")
      }
      return(brush)
    }
  })
  
  output$volcano1 <- renderPlot({
    if(!is.null(input$xrange)){
      if(is.null(d_row_count_matrix())){
        return(NULL)
      }else{
        pair_volcano()
      }
    }
  })
  
  output$download_pair_volcano = downloadHandler(
    filename = function(){
      if(input$GOI_color_type == "default") paste0(download_pair_overview_dir(), "_volcano.pdf") else paste0(download_pair_overview_dir(), "_volcano_",input$GOI_color_pathway2,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 5
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(pair_volcano())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  pair_GOIheatmap <- reactive({
    count <- deg_norm_count()
    if(gene_type1() != "SYMBOL"){
      if(length(grep("SYMBOL", colnames(count))) != 0){
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
    data2 <- pair_GOI_count()
    if(is.null(data2)){
      ht <- NULL
    }else{
      data.z <- genefilter::genescale(data2[,8:(7 + Cond_1 + Cond_2)], axis=1, method="Z")
      ht <- GOIheatmap(data.z)
    }
    return(ht)
  })
  
  output$GOIheatmap <- renderPlot({
    if(is.null(d_row_count_matrix())){
      return(NULL)
    }else{
      if(!is.null(brush_info())){
        if(!is.null(input$GOI) || dim(brush_info())[1] != 0){
          withProgress(message = "Heatmap",{
            suppressWarnings(print(pair_GOIheatmap()))
            incProgress(1)
          })
        }}
    }
  })
  
  output$download_pair_GOIheatmap = downloadHandler(
    filename = function(){
      if(input$GOI_color_type == "default") paste0(download_pair_overview_dir(), "_GOIheatmap.pdf") else paste0(download_pair_overview_dir(), "_GOIheatmap_",input$GOI_color_pathway2,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 10
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(pair_GOIheatmap())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  pair_GOI_count <- reactive({
    data <- data_degcount()
    count <- deg_norm_count()
    if(is.null(data) || is.null(count)){
      return(NULL)
    }
    if(gene_type1() != "SYMBOL"){
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
    if(dim(brush_info())[1] == 0){
      if(gene_type1() != "SYMBOL"){
        if(input$Species != "not selected"){
          Unique_ID <- input$GOI
          label_data <- as.data.frame(Unique_ID, row.names = Unique_ID)
          data2 <- merge(data, label_data, by="Unique_ID")
          if(input$uniqueID_cut) {
            id_list <- gsub("\\\n.+$", "", data2$Unique_ID)
            dup_list <- unique(id_list[duplicated(id_list)])
            for(i in 1:length(data2$Unique_ID)){
              if(gsub("\\\n.+$", "", data2$Unique_ID[i]) == "NA") {
                data2$Unique_ID[i] <- gsub(".+\\s", "", data2$Unique_ID[i])
              }else if(! gsub("\\\n.+$", "", data2$Unique_ID[i]) %in% dup_list) {
                data2$Unique_ID[i] <- gsub("\\\n.+$", "", data2$Unique_ID[i])
              }
            }
          }
          rownames(data2) <- data2$Unique_ID
          data2 <- data2[, - which(colnames(data2) == "Row.names")]
        }else{
          Row.names <- input$GOI
          label_data <- as.data.frame(Row.names, row.names = Row.names)
          data2 <- merge(data, label_data, by="Row.names")
          rownames(data2) <- data2$Row.names}
      }else{
        Row.names <- input$GOI
        label_data <- as.data.frame(Row.names, row.names = Row.names)
        data2 <- merge(data, label_data, by="Row.names")
        rownames(data2) <- data2$Row.names
      }
    }else{
      data2<-brush_info()
      if(input$GOI_plot_select == "Volcano plot"){
        data2 <- data2[, - which(colnames(data2) == "minusLog10padj")]
      }else{
        data2 <- data2[, - which(colnames(data2) == "log2baseMean")]
      }
      if(gene_type1() != "SYMBOL"){
        if(input$Species != "not selected"){
          if(input$uniqueID_cut) {
            id_list <- gsub("\\\n.+$", "", data2$Unique_ID)
            dup_list <- unique(id_list[duplicated(id_list)])
            for(i in 1:length(data2$Unique_ID)){
              if(gsub("\\\n.+$", "", data2$Unique_ID[i]) == "NA") {
                data2$Unique_ID[i] <- gsub(".+\\s", "", data2$Unique_ID[i])
              }else if(! gsub("\\\n.+$", "", data2$Unique_ID[i]) %in% dup_list) {
                data2$Unique_ID[i] <- gsub("\\\n.+$", "", data2$Unique_ID[i])
              }
            }
          }
          rownames(data2) <- data2$Unique_ID
        }else{
          rownames(data2) <- data2$Row.names
        }
      }else{
        rownames(data2) <- data2$Row.names
      }
      if(input$GOI_color_type != "default" && !is.null(pair_pathway_color_gene())){
        data2 <- data2 %>% dplyr::filter(Row.names %in% rownames(pair_pathway_color_gene()))
      }
    }
    return(data2)
  })
  output$paired_for_GOItype <- renderUI({
    if(input$paired_sample == "Yes"){
      selectInput("paired_for_GOItype","Plot type",c("Boxplot","without boxplot"))
    }
  })
  pair_GOIbox <- reactive({
    count <- deg_norm_count()
    if(gene_type1() != "SYMBOL"){
      if(length(grep("SYMBOL", colnames(count))) != 0){
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
    data2 <- pair_GOI_count()
    if(dim(data2)[1] == 0) validate("")
    if(is.null(data2)){
      p <- NULL
    }else{
      data3 <- data2[,8:(7 + Cond_1 + Cond_2)]
      if(input$paired_sample != "Yes") {
        pair= NULL
        plottype="Boxplot"
      }else {
        pair = d_paired_sample_file()
        plottype = input$paired_for_GOItype
      }
      p <- GOIboxplot(data = data3, pair=pair,plottype = plottype) + scale_fill_manual(values=c("gray", "#ff8082"))
    }
    return(p)
  })
  
  output$GOIbox <- renderPlot({
    if(is.null(d_row_count_matrix())){
      return(NULL)
    }else{
      if(!is.null(brush_info())){
        if(!is.null(input$GOI) || dim(brush_info())[1] != 0){
          withProgress(message = "Boxplot",{
            suppressWarnings(print(pair_GOIbox()))
            incProgress(1)
          })
        }
      }}
  })
  
  output$download_pair_GOIbox = downloadHandler(
    filename = function(){
      if(input$GOI_color_type == "default") paste0(download_pair_overview_dir(), "_GOIboxplot.pdf") else paste0(download_pair_overview_dir(), "_GOIboxplot_",input$GOI_color_pathway2,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- data_degcount()
        count <- deg_norm_count()
        if(gene_type1() != "SYMBOL"){
          if(length(grep("SYMBOL", colnames(data))) != 0){
            count <- count[, - which(colnames(count) == "SYMBOL")]
          }
        }
        if(dim(brush_info())[1] == 0){
          if(gene_type1() != "SYMBOL"){
            if(input$Species != "not selected"){
              Unique_ID <- input$GOI
              label_data <- as.data.frame(Unique_ID, row.names = Unique_ID)
              data2 <- merge(data, label_data, by="Unique_ID")
              rownames(data2) <- data2$Unique_ID
              data2 <- data2[, - which(colnames(data2) == "Row.names")]
            }else{
              Row.names <- input$GOI
              label_data <- as.data.frame(Row.names, row.names = Row.names)
              data2 <- merge(data, label_data, by="Row.names")
              rownames(data2) <- data2$Row.names}
          }else{
            Row.names <- input$GOI
            label_data <- as.data.frame(Row.names, row.names = Row.names)
            data2 <- merge(data, label_data, by="Row.names")
            rownames(data2) <- data2$Row.names
          }
        }else{
          data2<-brush_info()
          if(input$GOI_plot_select == "Volcano plot"){
            data2 <- data2[, - which(colnames(data2) == "minusLog10padj")]
          }else{
            data2 <- data2[, - which(colnames(data2) == "log2baseMean")]
          }
          if(gene_type1() != "SYMBOL"){
            if(input$Species != "not selected"){
              rownames(data2) <- data2$Unique_ID
            }else{
              rownames(data2) <- data2$Row.names
            }
          }else{
            rownames(data2) <- data2$Row.names
          }
          if(input$GOI_color_type != "default" && !is.null(pair_pathway_color_gene())){
            data2 <- data2 %>% dplyr::filter(Row.names %in% rownames(pair_pathway_color_gene()))
          }
        }
        rowlist <- rownames(data2)
        if(input$pair_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(pair_GOIbox())
        dev.off()
        incProgress(1)
      })
    }
  )
  #pair-wise transcript-------------
  output$Transcript_top <- renderUI({
    select <-c('manual'="manual")
    selected <- "manual"
    if(input$Transcript_group != "DET_manual"){
      select <- c(select,'Top10'="Top10",'Top20'="Top20",'Top40'="Top40")
      selected <- "Top10"
    }
    radioButtons('Transcript_top','Mode:',select,selected = selected)
  })
  
  GOI_list_transcript <- reactive({
    data <- data_degcount()
    if(!is.null(d_row_count_matrix()) && !is.null(deg_norm_count()) &&
       input$Level_pair != "gene_level" && input$Species != "not selected"){
      GOI <- unique(deg_norm_count()$SYMBOL)
      return(GOI)
    }
  })
  
  output$transcript_manual  <- renderUI({
    if(is.null(GOI_list_transcript())){
      return(NULL)
    }else{
      if(!is.null(input$Transcript_top)){
      if(input$Transcript_top == "manual" && 
         input$Level_pair != "gene_level" && input$Species != "not selected" && 
         !is.null(d_row_count_matrix())){
      withProgress(message = "Preparing Gene list (about 10 sec)",{
        selectizeInput("transcript_manual", "GOI for transcripts", c(GOI_list_transcript()),multiple = TRUE, 
                       options = list(delimiter = " ", create = T))
      })
      }
      }else validate("")
    }
  })
  output$GOIreset_transcript_manual <- renderUI({
    actionButton("GOIreset_transcript_manual", "GOI reset")
  })
  observeEvent(input$GOIreset_transcript_manual, {
    withProgress(message = "Preparing Gene list (about 10 sec)",{
      updateSelectizeInput(session, "transcript_manual", choices = c(GOI_list_transcript()), 
                           selected = character(0),
                           options = list(delimiter = " ", create=TRUE, 'plugins' = list('remove_button'), persist = FALSE))
    })
  })
  
  DETs <- reactive({
    if(input$Level_pair != "gene_level" && input$Species != "not selected" && 
       !is.null(d_row_count_matrix())){
    up <- data_degcount_up()
    down <- data_degcount_down()
    if(input$Transcript_group == "DET"){
    data <- rbind(up,down)
    }else if(input$Transcript_group == "DET_both"){
      up <- up %>% tibble::rownames_to_column() %>% dplyr::filter(!is.na(SYMBOL), SYMBOL != "")
      down <- down %>% tibble::rownames_to_column() %>% dplyr::filter(!is.na(SYMBOL), SYMBOL != "")
      merge<-merge(up,down,by="SYMBOL")
      sample_n <- dim(up)[2]
      data <- merge  %>% dplyr::filter(SYMBOL!="") %>%
        dplyr::select(c(1:sample_n)) %>% dplyr::select(-starts_with("rowname"))
    }else if(input$Transcript_group == "DET_manual"){
      gene <- data.frame(SYMBOL=input$transcript_manual)
      norm_count <- deg_norm_count()
      data <- merge(gene,norm_count,by="SYMBOL")
    }
    data <- isoform_balance(data,norm_count = deg_norm_count())
    return(data)
    }
  })
  DETs_plot <- reactive({
    if(input$Transcript_top == "Top10") gene <- DETs()[["res"]]$SYMBOL[1:10]
    if(input$Transcript_top == "Top20") gene <- DETs()[["res"]]$SYMBOL[1:20]
    if(input$Transcript_top == "Top40") gene <- DETs()[["res"]]$SYMBOL[1:40]
    if(input$Transcript_top == "manual") if(!is.null(input$transcript_manual)) gene <- input$transcript_manual else validate("")
    p <- barplot_forTranscript(table=deg_norm_count(),name=gene)
    return(p)
  })
  output$transcript_barplot <- renderPlot({
    if(!is.null(input$Transcript_top)&& 
       input$Level_pair != "gene_level" && !is.null(d_row_count_matrix())){
      DETs_plot() 
    }
  })
  output$transcript_barplot_table <- DT::renderDT({
    if(!is.null(input$Transcript_top)&& 
       input$Level_pair != "gene_level" && input$Species != "not selected" && 
       !is.null(d_row_count_matrix())){
    DETs()[["res"]] %>%
        datatable(
          selection = "single",
          filter = "top")
    }
  })
  transcript_GOIbox <- reactive({
    if(!is.null(input$transcript_barplot_table_rows_selected)){
      gene <- DETs()[["res"]][input$transcript_barplot_table_rows_selected,]$SYMBOL
      data <- deg_norm_count() %>% dplyr::filter(SYMBOL == gene) %>%
        dplyr::select(-SYMBOL)
      return(data)
    }
  })
  
  output$transcript_GOIboxplot <- renderPlot({
    if(!is.null(input$transcript_barplot_table_rows_selected)){
      GOIboxplot(data = transcript_GOIbox())+ scale_fill_manual(values=c("gray", "#ff8082"))
    }
  })
  output$download_transcript_barplot = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), "_transcript_barplot.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$Transcript_top == "Top10") rowlist <- 10
        if(input$Transcript_top == "Top20") rowlist <- 20
        if(input$Transcript_top == "Top40") rowlist <- 40
        if(input$Transcript_top == "manual") rowlist <- length(input$transcript_manual)
        if(input$pair_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)*2
        }else pdf_height <- input$pair_pdf_height*2
        if(input$pair_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)*3
        }else pdf_width <- input$pair_pdf_width*3
        pdf(file, height = pdf_height, width = pdf_width)
        print(DETs_plot())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_transcript_barplot_table = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), "_transcript_analysis.txt")
    },
    content = function(file) {
      write.table(DETs()[["res"]],file,row.names = F,col.names = T,quote = F,sep = "\t")
    }
  )
  output$download_transcript_GOIboxplot = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), 
             "_transcript_",DETs()[["res"]][input$transcript_barplot_table_rows_selected,]$SYMBOL,
             ".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        gene <- DETs()[["res"]][input$transcript_barplot_table_rows_selected,]$SYMBOL
        rowlist <- deg_norm_count() %>% dplyr::filter(SYMBOL == gene)
        rowlist <- dim(rowlist)[1]
        if(input$pair_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)*2
        }else pdf_height <- input$pair_pdf_height*2
        if(input$pair_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)*3
        }else pdf_width <- input$pair_pdf_width*3
        pdf(file, height = pdf_height, width = pdf_width)
        print(GOIboxplot(data = transcript_GOIbox())+ scale_fill_manual(values=c("gray", "#ff8082")))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  #pair-wise DTU-------------
  output$pair_DTU <- renderText({
    if(input$Level_pair == "gene_level") print("Input data must be a 'Transcript-level' count matrix")
  })
  output$DTU_top <- renderUI({
    select <-c('manual'="manual")
    select <- c(select,'Top10'="Top10",'Top20'="Top20",'Top40'="Top40")
    selected <- "Top10"
    radioButtons('DTU_top','Mode:',select,selected = selected)
  })
  
  GOI_list_DTU <- reactive({
    data <- data_degcount()
    if(!is.null(d_row_count_matrix()) && !is.null(deg_norm_count()) && !is.null(DRIMSeq()) &&
       input$Level_pair != "gene_level" && input$Species != "not selected"){
      GOI <- unique(DRIMSeq()$gene_id)
      ##more than two transcripts per geneに変更
      return(GOI)
    }
  })

  observeEvent(GOI_list_DTU(), {
    choices <- normalize_goi_choices(GOI_list_DTU())
    selected <- intersect(isolate(input$DTU_manual), choices)
    freezeReactiveValue(input, "DTU_manual")
    updateSelectizeInput(
      session,
      "DTU_manual",
      choices = choices,
      selected = selected,
      options = goi_selectize_options,
      server = TRUE
    )
  }, ignoreNULL = FALSE)

  observeEvent(input$GOIreset_DTU_manual, {
    choices <- normalize_goi_choices(GOI_list_DTU())
    freezeReactiveValue(input, "DTU_manual")
    updateSelectizeInput(
      session,
      "DTU_manual",
      choices = choices,
      selected = character(0),
      options = goi_selectize_options,
      server = TRUE
    )
  })
  
  DRIMSeq <- reactive({
    if(!is.null(d_row_count_matrix()) && !is.null(deg_norm_count()) && !is.null(isoform1()) &&
       input$Level_pair != "gene_level" && input$Species != "not selected"){
    isoform <- isoform1()
    isoform[is.na(isoform$SYMBOL),]$SYMBOL <- isoform[is.na(isoform$SYMBOL),]$Transcript_ID
    count <- d_row_count_matrix()
    collist <- gsub("\\_.+$", "", colnames(count))
    group <- data.frame(sample_id = colnames(count),group = factor(collist))
    count$Transcript_ID <- rownames(count)
    count <- merge(isoform[,1:2],count,by="Transcript_ID")
    colnames(count)[1:2]<-c("feature_id","gene_id")
    print(head(count))
    d = DRIMSeq::dmDSdata(counts = count, samples = group)
    n = nrow(group)
    n.small = min(table(group$group))
    d = DRIMSeq::dmFilter(d,
                 min_samps_feature_expr = n.small, min_feature_expr = 10,
                 min_samps_feature_prop = n.small, min_feature_prop = 0.1,
                 min_samps_gene_expr = n, min_gene_expr = 10)
    design = model.matrix(~group, data = group)
    withProgress(message = "DRIMSeq (about 15 minutes)",{
    set.seed(12345)
      d <- DRIMSeq::dmPrecision(d, design = design)    # Estimate the precision (Higher dispersion is associated with lower precision)
      d <- DRIMSeq::dmFit(d, design = design)          # Fit regression coefficients
      d <- DRIMSeq::dmTest(d, coef = colnames(design)[2])     # Perform null hypothesis testing on the coefficient of interest
    })
      res <- DRIMSeq::results(d)
      res.txp <- DRIMSeq::results(d, level="feature")
      no.na <- function(x) ifelse(is.na(x), 1, x)
      res$pvalue <- no.na(res$pvalue)
      res.txp$pvalue <- no.na(res.txp$pvalue)
      ##Step-wise analysis
      pScreen <- res$pvalue
      strp <- function(x) substr(x,1,15)
      names(pScreen) <- strp(res$gene_id)
      pConfirmation <- matrix(res.txp$pvalue, ncol=1)
      rownames(pConfirmation) <- strp(res.txp$feature_id)
      tx2gene <- res.txp[,c("feature_id", "gene_id")]
      for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
      stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                            pScreenAdjusted=FALSE, tx2gene=tx2gene)
      stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05) 
      suppressWarnings({
        drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                        onlySignificantGenes=TRUE)
      })
        #Post-hoc filtering on the standard deviation in proportions
        res.txp.filt <- DRIMSeq::results(d, level="feature")
        getSampleProportions <- function(d) {
          cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
          gene.cts <-rowsum(cts,  counts(d)$gene_id)
          total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
          cts/total.cts
        }
        prop.d <- getSampleProportions(d)
        res.txp.filt$prop.sd <- sqrt(rowVars(prop.d))
        res.txp.filt$pvalue[res.txp.filt$prop.sd < .1] <- 1
        res.txp.filt$adj_pvalue <- p.adjust(res.txp.filt$pvalue, method="BH")
        res.txp.filt <- res.txp.filt %>% dplyr::arrange(adj_pvalue)
        return(res.txp.filt)
    }
  })

  DTU_plot <- reactive({
    if(input$DTU_top == "Top10") gene <- unique(DRIMSeq()$gene_id)[1:10]
    if(input$DTU_top == "Top20") gene <- unique(DRIMSeq()$gene_id)[1:20]
    if(input$DTU_top == "Top40") gene <- unique(DRIMSeq()$gene_id)[1:40]
    if(input$DTU_top == "manual") if(!is.null(input$DTU_manual)) gene <- input$DTU_manual else validate("")
    p <- barplot_forTranscript(table=deg_norm_count(),name=gene)
    return(p)
  })
  output$DRIMSeq_barplot <- renderPlot({
    if(!is.null(input$DTU_top) && !is.null(DTU_plot( )) &&
       input$Level_pair != "gene_level" && !is.null(d_row_count_matrix())){
      DTU_plot() 
    }
  })
  output$DRIMSeq_table <- DT::renderDT({
    if(!is.null(input$DTU_top)&& !is.null(DRIMSeq()) &&
       input$Level_pair != "gene_level" && input$Species != "not selected" && 
       !is.null(d_row_count_matrix())){
      DRIMSeq() %>%
        datatable(
          selection = "single",
          filter = "top")
    }
  })
  output$download_DRIMSeq_table = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), "_DRIMSeq.txt")
    },
    content = function(file) {
      write.table(DRIMSeq(),file,row.names = F,col.names = T,quote = F,sep = "\t")
    }
  )
  
  DTU_GOIbox <- reactive({
    if(!is.null(input$DRIMSeq_table_rows_selected)){
      gene <- DRIMSeq()[input$DRIMSeq_table_rows_selected,]$gene_id
      data <- deg_norm_count() %>% dplyr::filter(SYMBOL == gene) %>%
        dplyr::select(-SYMBOL)
      return(data)
    }
  })
  
  output$DRIMSeq_GOIboxplot <- renderPlot({
    if(!is.null(input$DRIMSeq_table_rows_selected)){
      GOIboxplot(data = DTU_GOIbox())+ scale_fill_manual(values=c("gray", "#ff8082"))
    }
  })
  output$download_DRIMSeq_barplot = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), "_DTU_barplot.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$DTU_top == "Top10") rowlist <- 10
        if(input$DTU_top == "Top20") rowlist <- 20
        if(input$DTU_top == "Top40") rowlist <- 40
        if(input$DTU_top == "manual") rowlist <- length(input$DTU_manual)
        if(input$pair_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)*2
        }else pdf_height <- input$pair_pdf_height*2
        if(input$pair_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)*3
        }else pdf_width <- input$pair_pdf_width*3
        pdf(file, height = pdf_height, width = pdf_width)
        print(DTU_plot())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_DRIMSeq_GOIboxplot = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), 
             "_transcript_",DRIMSeq()[input$DRIMSeq_table_rows_selected,]$gene_id,
             ".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        gene <- DRIMSeq()[input$DRIMSeq_table_rows_selected,]$gene_id
        rowlist <- deg_norm_count() %>% dplyr::filter(SYMBOL == gene)
        rowlist <- dim(rowlist)[1]
        if(input$pair_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)*2
        }else pdf_height <- input$pair_pdf_height*2
        if(input$pair_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)*3
        }else pdf_width <- input$pair_pdf_width*3
        pdf(file, height = pdf_height, width = pdf_width)
        print(GOIboxplot(data = DTU_GOIbox())+ scale_fill_manual(values=c("gray", "#ff8082")))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  gtf_input <- reactive({
    if(is.null(input$GTF_file)){
      if(input$goButton > 0 ){
        gtf <- rtracklayer::import("/data/hg19.ncbiRefSeq.gtf")
        return(gtf)
      }
      return(NULL)
    }else{
      file <- input$GTF_file$datapath
      gtf <- rtracklayer::import(file)
      return(gtf)
    }
  })
  transcript_plot <- reactive({
    gtf <- gtf_input() %>% dplyr::as_tibble()
    gene_of_interest <- DRIMSeq()[input$DRIMSeq_table_rows_selected,]$gene_id
    data <- deg_norm_count() %>% dplyr::filter(SYMBOL == gene_of_interest)
    id <- rownames(data)
    annotation_from_gtf <- gtf %>%
      dplyr::filter(
        !is.na(gene_name),
        gene_name == gene_of_interest
      )
    annotation_from_gtf$transcript_id <- gsub("\\..+$", "", annotation_from_gtf$transcript_id )
    exons <- annotation_from_gtf %>% 
      dplyr::filter(type == "exon") %>%
      dplyr::filter(transcript_id %in% id)
    p <- exons %>%
      ggplot(aes(
        xstart = start,
        xend = end,
        y = transcript_id
      )) +
      geom_range(
        aes(fill = transcript_id)
      ) +
      geom_intron(
        data = to_intron(exons, "transcript_id"),
        aes(strand = strand)
      )+theme_bw()+
      theme(legend.position = "none" , legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 12),
            axis.text.y= ggplot2::element_text(size = 12),
            text = ggplot2::element_text(size = 15),
            title = ggplot2::element_text(size = 15))
    return(p)
  })
  output$transcript_structure <- renderPlot({
    if(!is.null(input$DRIMSeq_table_rows_selected)){
    transcript_plot()
    }
  })
  # pair-wise PCA ------------------------------------------------------------------------------
  pair_pca_plot <- reactive({
    data <- deg_norm_count()
    if(is.null(data)){
      return(NULL)
    }
    PCAplot(data = data,legend = input$PCA_legend)
  })

  pair_pca_data <- reactive({
    row_count <- d_row_count_matrix()
    deg_norm_count_value <- deg_norm_count()
    if(is.null(row_count) || is.null(deg_norm_count_value)){
      return(NULL)
    }
    PCAdata(row_count = row_count, deg_norm_count = deg_norm_count_value)
  })

  output$download_pair_PCA = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), "_PCA-MDS-dendrogram.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$pair_pdf_height == 0){
          pdf_height <- 3.5
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 9
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(pair_pca_plot())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$PCA <- renderPlot({
    plot_obj <- pair_pca_plot()
    if(is.null(plot_obj)){
      return(NULL)
    }else{
      print(plot_obj)
    }
  })
  
  output$pair_PCA_data <- DT::renderDataTable({
    pair_pca_data()
  })
  
  output$download_pair_PCA_table = downloadHandler(
    filename = function() {
      paste0(download_pair_overview_dir(), "_PCA_table.txt")
    },
    content = function(file){write.table(pair_pca_data(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  
  output$download_pair_report = downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%m%d_"),download_pair_overview_dir(), "_Pair-wiseDEG",".zip")
    },
    content = function(fname){
      withProgress(message = "Preparing download, please wait",{
        pair_analysis_requested(TRUE)
        on.exit({
          if(is.null(input$pair_tabs) || input$pair_tabs == "pair_input_data_tab"){
            pair_analysis_requested(FALSE)
          }
        }, add = TRUE)
        fs <- c()
        setwd(tempdir())
        print(fs)
        dir.create("DEG_result/",showWarnings = FALSE)
        dir.create("Clustering/",showWarnings = FALSE)
        DEG <- "DEG_result/DEG_result.txt"
        up <- "DEG_result/up.txt"
        down <- "DEG_result/down.txt"
        count <- "DEG_result/normalized_count.txt"
        PCA <- "Clustering/clustering.pdf"
        PCA_table <- "Clustering/pca.txt"
        MAplot <- "DEG_result/MAplot.pdf"
        fs <- c(DEG, up,down,count,PCA,PCA_table,MAplot)
        write.table(deg_result(), DEG, row.names = T, col.names=NA, sep = "\t", quote = F)
        write.table(deg_norm_count(), count, row.names = T, col.names=NA, sep = "\t", quote = F)
        write.table(data_degcount_down(), down, quote = F, row.names = T, col.names=NA, sep = "\t")
        write.table(data_degcount_up(), up, quote = F, row.names = T, col.names=NA, sep = "\t")
        write.table(PCAdata(row_count = d_row_count_matrix(), deg_norm_count = deg_norm_count()), PCA_table, row.names = T, col.names=NA, sep = "\t", quote = F)
        pdf(PCA, height = 3.5, width = 9)
        print(PCAplot(data = deg_norm_count(),legend = input$PCA_legend))
        dev.off()
        pdf(MAplot, height = 4, width = 7)
        print(ma_heatmap_plot())
        dev.off()
        if(!is.null(input$xrange)){
          dir.create("GOI_profiling/",showWarnings = FALSE)
          if(input$GOI_color_type == "default") volcano <- "GOI_profiling/volcano_plot.pdf" else volcano <- paste0("GOI_profiling/volcano_plot_",input$GOI_color_pathway2,".pdf")
          fs <- c(fs,volcano)
          pdf(volcano, height = 5, width = 5)
          print(pair_volcano())
          dev.off()
          if(!is.null(brush_info())){
            if(!is.null(input$GOI) || dim(brush_info())[1] != 0){
                if(input$GOI_color_type == "default") {
                  boxplot <- "GOI_profiling/boxplot.pdf"
                  heat <- "GOI_profiling/heatmap.pdf"
                }else{
                  boxplot <- paste0("GOI_profiling/boxplot_",input$GOI_color_pathway2,".pdf")
                  heat <- paste0("GOI_profiling/heatmap_",input$GOI_color_pathway2,".pdf")
                } 
              fs <- c(fs,boxplot,heat)
              data <- data_degcount()
              count <- deg_norm_count()
              if(gene_type1() != "SYMBOL"){
                if(length(grep("SYMBOL", colnames(data))) != 0){
                  count <- count[, - which(colnames(count) == "SYMBOL")]
                }
              }
              if(dim(brush_info())[1] == 0){
                if(gene_type1() != "SYMBOL"){
                  if(input$Species != "not selected"){
                    Unique_ID <- input$GOI
                    label_data <- as.data.frame(Unique_ID, row.names = Unique_ID)
                    data2 <- merge(data, label_data, by="Unique_ID")
                    rownames(data2) <- data2$Unique_ID
                    data2 <- data2[, - which(colnames(data2) == "Row.names")]
                  }else{
                    Row.names <- input$GOI
                    label_data <- as.data.frame(Row.names, row.names = Row.names)
                    data2 <- merge(data, label_data, by="Row.names")
                    rownames(data2) <- data2$Row.names}
                }else{
                  Row.names <- input$GOI
                  label_data <- as.data.frame(Row.names, row.names = Row.names)
                  data2 <- merge(data, label_data, by="Row.names")
                  rownames(data2) <- data2$Row.names
                }
              }else{
                data2<-brush_info()
                if(input$GOI_plot_select == "Volcano plot"){
                  data2 <- data2[, - which(colnames(data2) == "minusLog10padj")]
                }else{
                  data2 <- data2[, - which(colnames(data2) == "log2baseMean")]
                }
                if(gene_type1() != "SYMBOL"){
                  if(input$Species != "not selected"){
                    rownames(data2) <- data2$Unique_ID
                  }else{
                    rownames(data2) <- data2$Row.names
                  }
                }else{
                  rownames(data2) <- data2$Row.names
                }
                if(input$GOI_color_type != "default" && !is.null(pair_pathway_color_gene())){
                  data2 <- data2 %>% dplyr::filter(Row.names %in% rownames(pair_pathway_color_gene()))
                }
              }
              rowlist <- rownames(data2)
              pdf_height <- pdf_h(rowlist)
              pdf_width <- pdf_w(rowlist)
              pdf(boxplot, height = pdf_height, width = pdf_width)
              print(pair_GOIbox())
              dev.off()
              pdf(heat, height = 10, width = 7)
              print(pair_GOIheatmap())
              dev.off()
            }
          }
        }
        if(input$Species != "not selected" && !is.null(input$Gene_set)){
          dir.create("Enrichment_analysis/",showWarnings = FALSE)
          enrichplot <- paste0("Enrichment_analysis/Enrichment_analysis_",input$Gene_set,".pdf")
          enrich_table <- paste0("Enrichment_analysis/Enrichment_analysis_",input$Gene_set,".txt")
          enrich_gseatable <- paste0("Enrichment_analysis/GSEA_",input$Gene_set,".txt")
          fs <- c(fs,enrichplot,enrich_table,enrich_gseatable)
          p1 <- pair_enrich1_H()
          p2 <- pair_enrich2()
          pdf(enrichplot, height = 10, width = 12)
          print(plot_grid(p1, p2, nrow =2))
          dev.off()
          write.table(pair_enrich_table(), enrich_table, row.names = F, sep = "\t", quote = F)
          write.table(pair_gsea_table(), enrich_gseatable, row.names = F, sep = "\t", quote = F)
        }
        report <- paste0(format(Sys.time(), "%Y%m%d_"),"pairwise_report",".docx")
        fs <- c(fs,report)
        rmarkdown::render("pair_report.Rmd", output_format = "word_document", output_file = report,
                          params = list(raw_count = d_row_count_matrix(),
                                        norm_count = norm_count_matrix(),
                                        input = input,
                                        deg_norm_count = deg_norm_count(),
                                        data_degcount = data_degcount(),
                                        data_degcount2 = data_degcount2(),
                                        enrichment_enricher = enrichment_enricher(),
                                        enrichment_1_gsea = enrichment_1_gsea(),
                                        pair_enrich_table = pair_enrich_table(),
                                        pair_gsea_table = pair_gsea_table(),
                                        metadata = metadata()), 
                          envir = new.env(parent = globalenv()),intermediates_dir = tempdir(),encoding="utf-8"
        )
        zip(zipfile=fname, files=fs)
      })
    },
    contentType = "application/zip"
  )
  # pair-wise enrichment ------------------------------------------------------------------------------
  output$Spe1 <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  Hallmark_set <- reactive({
    return(GeneList_for_enrichment(Species = input$Species, Ortholog=input$Ortholog,Custom_gene_list = Custom_input_pair(),
                                   Biomart_archive=input$Biomart_archive,Gene_set = input$Gene_set, org = org1(),gene_type = gene_type1()))
  })
  
  enrichment_1_1 <- reactive({
    df <- enrichment_enricher()
    if(!is.null(input$Gene_set) && input$Species != "not selected" && !is.null(df)){
      data <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
      colnames(data) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "Group")
      for(name in names(df)){
        if(!is.null(df[[name]])) {
          group1 <- as.data.frame(df[[name]])
        }else group1 <- NULL
        group1$Group <- name
        data <- rbind(data, group1)
      }
      if(length(data$Description) != 0) {
        data["Description"] <- lapply(data["Description"], gsub, pattern="HALLMARK_", replacement = "")
        data$GeneRatio <- DOSE::parse_ratio(data$GeneRatio)
        return(data)
      }else return(NULL)
    }else{return(NULL)}
  })
  
  enrichment_enricher <- reactive({
    data3 <- data_degcount2()
    if(!is.null(input$Gene_set) && input$Species != "not selected" && !is.null(data3)){
      withProgress(message = "enrichment analysis",{
        if(input$Species != "Xenopus laevis" && input$Ortholog != "Arabidopsis thaliana" && input$Species != "Arabidopsis thaliana"){
          H_t2g <- Hallmark_set()
          H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
          em_up <- try(clusterProfiler::enricher(dplyr::filter(data3, group == "Up")$ENTREZID, TERM2GENE=H_t2g2, pvalueCutoff = 0.05))
          em_down <- try(clusterProfiler::enricher(dplyr::filter(data3, group == "Down")$ENTREZID, TERM2GENE=H_t2g2, pvalueCutoff = 0.05))
        }else{
          if(input$Gene_set == "KEGG"){
            em_up <- try(clusterProfiler::enrichKEGG(dplyr::filter(data3, group == "Up")$ENTREZID, organism = org_code(input$Species, Ortholog= input$Ortholog), pvalueCutoff = 0.05,keyType = "ncbi-geneid")) 
            em_down <- try(clusterProfiler::enrichKEGG(dplyr::filter(data3, group == "Down")$ENTREZID, organism = org_code(input$Species, Ortholog= input$Ortholog), pvalueCutoff = 0.05,keyType = "ncbi-geneid"))
          }
          if(input$Gene_set == "GO biological process"){
            em_up <- try(clusterProfiler::enrichGO(dplyr::filter(data3, group == "Up")$ENTREZID, OrgDb = org(input$Species, Ortholog= input$Ortholog), ont = "BP",pvalueCutoff = 0.05)) 
            em_down <- try(clusterProfiler::enrichGO(dplyr::filter(data3, group == "Down")$ENTREZID, OrgDb = org(input$Species, Ortholog= input$Ortholog), ont = "BP",pvalueCutoff = 0.05))
          }
          if(input$Gene_set == "GO cellular component"){
            em_up <- try(clusterProfiler::enrichGO(dplyr::filter(data3, group == "Up")$ENTREZID, OrgDb= org(input$Species, Ortholog= input$Ortholog), ont = "CC",pvalueCutoff = 0.05)) 
            em_down <- try(clusterProfiler::enrichGO(dplyr::filter(data3, group == "Down")$ENTREZID,OrgDb= org(input$Species, Ortholog= input$Ortholog), ont = "CC",pvalueCutoff = 0.05))
          }
          if(input$Gene_set == "GO molecular function"){
            em_up <- try(clusterProfiler::enrichGO(dplyr::filter(data3, group == "Up")$ENTREZID, OrgDb = org(input$Species, Ortholog= input$Ortholog), ont = "MF",pvalueCutoff = 0.05)) 
            em_down <- try(clusterProfiler::enrichGO(dplyr::filter(data3, group == "Down")$ENTREZID, OrgDb = org(input$Species, Ortholog= input$Ortholog), ont = "MF",pvalueCutoff = 0.05))
          }
        }
        df <- list()
        df[["Up"]] <- em_up
        df[["Down"]] <- em_down
        for(name in names(df)){
          if (length(as.data.frame(df[[name]])$ID) == 0) {
            df[[name]] <- NULL
          } else{
            df[[name]] <- try(clusterProfiler::setReadable(df[[name]], org1(), 'ENTREZID'))
            if(inherits(df[[name]], "try-error")) df[[name]] <- NULL
          }
        }
        incProgress(1)
        return(df)
      })
    }else{return(NULL)}
  })
  
  enrichment_1_gsea <- reactive({
    data <- data_degcount()
    if(!is.null(input$Gene_set) && input$Species != "not selected" &&
       !is.null(data) && !is.null(data_degcount2())){
      data <- na.omit(data)
      geneList <- data$log2FoldChange
      names(geneList) = as.character(data$ENTREZID)
      geneList <- sort(geneList, decreasing = TRUE)
      withProgress(message = "GSEA",{
        if(input$Species != "Xenopus laevis" && input$Ortholog != "Arabidopsis thaliana" && input$Species != "Arabidopsis thaliana"){
          H_t2g <- Hallmark_set()
          H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
          set.seed(123)
          em3 <- try(clusterProfiler::GSEA(geneList, TERM2GENE = H_t2g2,pvalueCutoff = 0.05,exponent = 1, eps = 0, pAdjustMethod = "BH",
                          minGSSize = 50, maxGSSize = 5000,by = "fgsea",verbose = F))
        }else{
          set.seed(123)
          if(input$Gene_set == "KEGG"){
            em3 <- try(clusterProfiler::gseKEGG(geneList, organism = org_code(input$Species, Ortholog= input$Ortholog),pvalueCutoff = 0.05,exponent = 1, eps = 0, pAdjustMethod = "BH",
                               minGSSize = 50, maxGSSize = 500,by = "fgsea",verbose = F,keyType = "ncbi-geneid"))
          }
          if(input$Gene_set == "GO biological process"){
            em3 <- try(clusterProfiler::gseGO(geneList, OrgDb = org(input$Species, Ortholog= input$Ortholog),ont = "BP",pvalueCutoff = 0.05,exponent = 1, eps = 0, pAdjustMethod = "BH",
                             minGSSize = 50, maxGSSize = 500,by = "fgsea",verbose = F))
          }
          if(input$Gene_set == "GO cellular component"){
            em3 <- try(clusterProfiler::gseGO(geneList, OrgDb = org(input$Species, Ortholog= input$Ortholog),ont = "CC",pvalueCutoff = 0.05,exponent = 1, eps = 0, pAdjustMethod = "BH",
                             minGSSize = 50, maxGSSize = 500,by = "fgsea",verbose = F))
          }
          if(input$Gene_set == "GO molecular function"){
            em3 <- try(clusterProfiler::gseGO(geneList, OrgDb = org(input$Species, Ortholog= input$Ortholog),ont = "MF",pvalueCutoff = 0.05,exponent = 1, eps = 0, pAdjustMethod = "BH",
                             minGSSize = 50, maxGSSize = 500,by = "fgsea",verbose = F))
          }
        }
        if(inherits(em3, "try-error")) validate(em3)
        if (length(as.data.frame(em3)$ID) == 0) {
          em4 <- NA
        } else{
          em4 <- clusterProfiler::setReadable(em3, org1(), 'ENTREZID')
        }
        return(em4)
        incProgress(1)
      })
    }else return(NULL)
  })
  
  
  # pair-wise enrichment plot ------------------------------------------------------------------------------
  Custom_input_pair <- reactive({
    tmp <- input$custom_input_pair$datapath
    data <- read_gene_list(tmp)
    df <- gene_list_convert_for_enrichment(gene_type=gene_type1(),data= data, org=org1(),Species = input$Species,Ortholog=ortholog1(),Isoform=isoform1())
    return(df)
  })
  output$Custom_input_pair <- renderUI({
    if(is.null(input$Gene_set)){
      return(NULL)
    }else{
      if(input$Gene_set == "Custom gene set"){
        fileInput("custom_input_pair",
                  "Select a file (txt, csv, xlsx)",
                  accept = c("txt", "csv","xlsx"),
                  multiple = FALSE,
                  width = "80%")
      }else return(NULL)
    }
  })
  
  pair_enrich1_H <- reactive({
    if(!is.null(input$Gene_set) && input$Species != "not selected"){
      df <- enrichment_enricher()
      data3 <- data_degcount2()
      if (is.null(df[["Up"]]) && is.null(df[["Down"]]))  {
        p1 <- NULL
      } else{
        data <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
        colnames(data) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "Group")
        for(name in names(df)){
          if(!is.null(df[[name]])) {
            group1 <- as.data.frame(df[[name]])
            group1$Group <- paste(name, "\n(",length(dplyr::filter(data3, group == name)$ENTREZID),")",sep = "")
            if (length(group1$pvalue) > 5){
              group1 <- group1[sort(group1$pvalue, decreasing = F, index=T)$ix,]
              group1 <- group1[1:5,]
            }
          }else group1 <- NULL
          data <- rbind(data, group1)
        }
        if(length(data$Description) != 0){
          data["Description"] <- lapply(data["Description"], gsub, pattern="HALLMARK_", replacement = "")
          data$GeneRatio <- DOSE::parse_ratio(data$GeneRatio)
          if ((length(data$Description) == 0) || length(which(!is.na(unique(data$qvalue))))==0) {
            p1 <- NULL
          } else{
            data$Description <- gsub("_", " ", data$Description)
            data <- dplyr::mutate(data, x = paste0(Group, 1/(-log10(eval(parse(text = "qvalue"))))))
            data$x <- gsub(":","", data$x)
            data <- dplyr::arrange(data, x)
            idx <- order(data[["x"]], decreasing = FALSE)
            data$Description <- factor(data$Description,
                                       levels=rev(unique(data$Description[idx])))
            p1 <- as.grob(ggplot(data, aes_string(x = "Group",y= "Description",color="qvalue",size="GeneRatio"))+
                            geom_point() +
                            scale_color_continuous(low="red", high="blue",
                                                   guide=guide_colorbar(reverse=TRUE)) +
                            scale_size(range=c(1, 6))+ DOSE::theme_dose(font.size=12)+ylab(NULL)+xlab(NULL) + 
                            scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
          }}else p1 <- NULL
      }
      em3 <- enrichment_1_gsea()
      if (length(as.data.frame(em3)$ID) == 0) {
        p4 <- NULL
      } else{
        if (length(as.data.frame(em3)$ID) >= 5){
          p4 <- enrichplot::gseaplot2(em3, 1:5, pvalue_table = F,base_size = 14)
        }else{
          p4 <- enrichplot::gseaplot2(em3, 1:length(em3$ID), pvalue_table = F,base_size = 14)
        }
      }
      p <- plot_grid(p1, print(p4), nrow = 1)
      return(p)
    }else return(NULL)
  })
  
  output$enrichment1 <- renderPlot({
    dotplot_for_output(data = d_row_count_matrix(),
                       plot_genelist = pair_enrich1_H(), Gene_set = input$Gene_set, 
                       Species = input$Species)
  })
  
  pair_enrich2 <- reactive({
    if(!is.null(input$Gene_set) && input$Species != "not selected"){
      data3 <- data_degcount2()
      df <- enrichment_enricher()
      upgene <- data3[data3$log2FoldChange > log(input$fc, 2),]
      downgene <- data3[data3$log2FoldChange < log(1/input$fc, 2),]
      p <- list()
      for(name in names(df)){
        if(length(as.data.frame(df[[name]])$ID) == 0){
          cnet1 <- NULL
        } else {
          cnet1 <- df[[name]]
        }
        if (length(as.data.frame(cnet1)$ID) == 0) {
          p2 <- NULL
        } else{
          if(name == "Up") {
            genes <- upgene
            color <- c("white", "red")
            limits <- c(0,NA)
          }
          if(name == "Down") {
            genes <- downgene
            color <- c("blue","white")
            limits <- c(NA,0)
          }
          geneList <- genes$log2FoldChange
          names(geneList) = as.character(genes$ENTREZID)
          p2_plot <- build_cnetplot(cnet1, foldChange = geneList,
                                    cex_label_gene = 0.7,
                                    cex_label_category = 0.75,
                                    cex_category = 0.75)
          if(!is.null(p2_plot)) {
            p2_plot <- try(p2_plot + guides(edge_colour = "none") +
                             scale_color_gradientn(name = "log2FC",
                                                   colours = color,
                                                   limits = limits),
                           silent = TRUE)
            if(inherits(p2_plot, "try-error")) {
              p2_plot <- try(build_cnetplot(cnet1, foldChange = geneList,
                                            cex_label_gene = 0.7,
                                            cex_label_category = 0.75,
                                            cex_category = 0.75) +
                               guides(edge_color = "none") +
                               scale_color_gradientn(name = "log2FC",
                                                     colours = color,
                                                     limits = limits),
                             silent = TRUE)
            }
            if(inherits(p2_plot, "try-error")) {
              p2_plot <- build_cnetplot(cnet1, foldChange = geneList,
                                        cex_label_gene = 0.7,
                                        cex_label_category = 0.75,
                                        cex_category = 0.75)
            }
          }
          p2 <- as_cnet_grob(p2_plot, remove_all_legend = FALSE)
        }
        p[[name]] <- p2
      }
      p <- plot_grid(p[["Up"]], p[["Down"]], nrow = 1)
      return(p)
    }else return(NULL)
  })
  
  output$enrichment2 <- renderPlot({
    cnet_for_output(data = d_row_count_matrix(), plot_data = pair_enrich2(), 
                    Gene_set = input$Gene_set, Species = input$Species)
  })
  
  output$download_pair_enrichment = downloadHandler(
    filename = function(){
      paste(download_pair_overview_dir(), paste0(input$Gene_set,"-enrichment.pdf"), sep="_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- pair_enrich1_H()
        p2 <- pair_enrich2()
        if(input$pair_pdf_height == 0){
          pdf_height <- 10
        }else pdf_height <- input$pair_pdf_height
        if(input$pair_pdf_width == 0){
          pdf_width <- 12
        }else pdf_width <- input$pair_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(p1, p2, nrow =2))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$Gene_set <- renderUI({
    if(input$Species != "Xenopus laevis" && input$Ortholog != "Arabidopsis thaliana" && input$Species != "Arabidopsis thaliana"){
      selectInput('Gene_set', 'Gene Set', gene_set_list)
    }else selectInput('Gene_set', 'Gene Set', c("KEGG", "GO biological process", 
                                                "GO cellular component","GO molecular function"))
  })
  
  
  pair_enrich_table <- reactive({
    if(!is.null(input$Gene_set) && input$Species != "not selected"){
      if(input$Species != "Xenopus laevis" && input$Ortholog != "Arabidopsis thaliana" && input$Species != "Arabidopsis thaliana"){
        table <- enrich_for_table(data = as.data.frame(enrichment_1_1()), H_t2g = Hallmark_set(), Gene_set = input$Gene_set)
      }else table <- as.data.frame(enrichment_1_1())
      return(table)
    }
  })
  
  output$pair_enrichment_result <- DT::renderDataTable({
    pair_enrich_table()
  })
  
  pair_gsea_table <- reactive({
    if(!is.null(input$Gene_set) && input$Species != "not selected"){
      data <- as.data.frame(enrichment_1_gsea())
      if(input$Species != "Xenopus laevis" && input$Ortholog != "Arabidopsis thaliana" && input$Species != "Arabidopsis thaliana"){
        H_t2g <- Hallmark_set()
        if(length(as.data.frame(data)$Description) == 0 || is.null(H_t2g)){
          return(NULL)
        }else{
          colnames(data)[1] <- "gs_name"
          H_t2g <- H_t2g %>% distinct(gs_name, .keep_all = T)
          data2 <- left_join(data, H_t2g, by="gs_name")  %>% as.data.frame()
          if(input$Gene_set == "DoRothEA regulon (activator)" || input$Gene_set == "DoRothEA regulon (repressor)"){
            data3 <- data.frame(Gene_set_name = data2$gs_name, Confidence = data2$confidence,
                                setSize = data2$setSize, enrichmentScore = data2$enrichmentScore, NES = data2$NES, 
                                pvalue = data2$pvalue, p.adjust = data2$p.adjust, qvalue = data2$qvalue, 
                                rank = data2$rank, leading_edge = data2$leading_edge, core_enrichment = data2$core_enrichment)
          }else{
            if(input$Gene_set == "Custom gene set"){
              data3 <- data.frame(Gene_set_name = data2$gs_name,
                                  setSize = data2$setSize, enrichmentScore = data2$enrichmentScore, NES = data2$NES, 
                                  pvalue = data2$pvalue, p.adjust = data2$p.adjust, qvalue = data2$qvalue, 
                                  rank = data2$rank, leading_edge = data2$leading_edge, core_enrichment = data2$core_enrichment)
            }else{
              data3 <- data.frame(Gene_set_name = data2$gs_name, Description = data2$gs_description,
                                  setSize = data2$setSize, enrichmentScore = data2$enrichmentScore, NES = data2$NES, 
                                  pvalue = data2$pvalue, p.adjust = data2$p.adjust, qvalue = data2$qvalue, 
                                  rank = data2$rank, leading_edge = data2$leading_edge, core_enrichment = data2$core_enrichment)
            }
            return(data3) 
          }
        }
      }else return(data)
    }
  })
  
  
  output$pair_GSEA_result <- DT::renderDataTable({
    pair_gsea_table()
  })
  
  output$download_pair_enrichment_table = downloadHandler(
    filename = function() {
      paste(download_pair_overview_dir(), paste0(input$Gene_set,"-enrichment.txt"), sep="_")
    },
    content = function(file){write.table(pair_enrich_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  output$download_pair_GSEA_table = downloadHandler(
    filename = function() {
      paste(download_pair_overview_dir(), paste0(input$Gene_set,"-GSEA.txt"), sep="_")
    },
    content = function(file){write.table(pair_gsea_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  
  output$Normalized_Count_matrix <- DT::renderDataTable({
    deg_norm_count()
  })
  
  batch_files <- reactive({
    upload = list()
    name = c()
    if(is.null(input$file11)){
      if(input$goButton > 0 ){
        df <- list()
        df["day0"] <- list(read.table(example_data_path("data/iwat-raw/day0.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/iwat-raw/day0.txt"), header = T, row.names = 1))
        df["day1"] <- list(read.table(example_data_path("data/iwat-raw/day1.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/iwat-raw/day1.txt"), header = T, row.names = 1))
        df["day5"] <- list(read.table(example_data_path("data/iwat-raw/day5.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/iwat-raw/day5.txt"), header = T, row.names = 1))
        return(df)
      }
      return(NULL)
    }else{
      for(nr in 1:length(input$file11[, 1])){
        df <- read_df(input$file11[[nr, 'datapath']])
        upload[gsub("\\..+$", "", input$file11[nr,]$name)] <- list(df)
      }
      name_list <- names(upload) %>% sort()
      upload <- upload[name_list]
      return(upload)
    }
  })
  
  dds_batch <- reactive({
    count_files <- batch_files()
    count_list <- list()
    withProgress(message = "DEG analysis",{
      file_num <- length(names(count_files))
      for (name in names(count_files)) {
        count <- count_files[[name]]
        collist <- gsub("\\_.+$", "", colnames(count))
        if (input$DEG_method == "DESeq2") {
          group <- data.frame(con = factor(collist))
          dds<- DESeqDataSetFromMatrix(countData = round(count),colData = group, design = ~ con)
          dds$con <- factor(dds$con, levels = unique(collist))
          dds <- DESeq(dds)
          incProgress(1/file_num, message = paste("DEG analysis of", name))
        }
        if (input$DEG_method == "edgeR") {
          group <- factor(collist)
          dds <- DGEList(counts = count, group = group)
          if(input$pair_prefilterON == "ON"){
          keep <- filterByExpr(dds,min.count=input$pair_prefilter,min.total.count=input$pair_prefilterTotal)
          dds = dds[keep, , keep.lib.sizes=FALSE]
          }
          dds <- calcNormFactors(dds)
          dds <- estimateCommonDisp(dds)
          dds <- estimateTagwiseDisp(dds)
          incProgress(1/file_num, message = paste("DEG analysis of", name))
        }
        count_list[name] <- list(dds)
      }
    })
    return(count_list)
  })
  
  
  deg_result_batch_raw <- reactive({
    count_files <- batch_files()
    if(is.null(count_files)){
      return(NULL)
    }
    count_list <- list()
    file_num <- length(names(count_files))
    for (name in names(count_files)) {
      if(name != "combined"){
        count <- count_files[[name]]
        collist <- gsub("\\_.+$", "", colnames(count))
        if (input$DEG_method == "DESeq2") {
          dds <- dds_batch()[[name]]
          contrast <- c("con", unique(collist))
          res <- DESeq2::results(dds,  contrast = contrast)
          if(input$FDR_method == "IHW") {
            ihw_res <- IHW::ihw(pvalue ~ baseMean,  data=as.data.frame(res), alpha = 0.1)
            res$padj <- IHW::as.data.frame(ihw_res)$adj_pvalue
          }
          if(input$FDR_method == "Qvalue") {
            res <- DESeq2::results(dds,  contrast = contrast)
            qvalue <- qvalue::qvalue(res$pvalue)
            res$padj <- qvalue$qvalues
          }
        }
        if (input$DEG_method == "edgeR") {
          group <- factor(collist)
          dds <- dds_batch()[[name]]
          result <- exactTest(dds, pair = c(unique(group)[2],unique(group)[1]))
          res <- as.data.frame(topTags(result, n = nrow(count)))
          qvalue <- qvalue::qvalue(res$PValue)
          res$padj <- qvalue$qvalues
          ihw_res <- try(IHW::ihw(PValue ~ 2^logCPM,  data=res, alpha = 0.1))
          if(inherits(ihw_res, "try-error")){
            res$ihw_padj <- NA
          }else{
            ihw_res_df <- IHW::as.data.frame(ihw_res)
            res$ihw_padj <- ihw_res_df$adj_pvalue
          }
          if(input$FDR_method == "BH"){label <- c("log2FoldChange", "log2CPM", "PValue","padj", "Qvalue", "IHW_FDR")}
          if(input$FDR_method == "Qvalue"){label <- c("log2FoldChange", "log2CPM", "PValue","BH_FDR", "padj", "IHW_FDR")}
          if(input$FDR_method == "IHW"){label <- c("log2FoldChange", "log2CPM", "PValue","BH_FDR", "Qvalue", "padj")}
          colnames(res) <- label
        }
        if(input$DEG_method == "limma"){
          withProgress(message = "limma",{
            collist <- gsub(" ", ".", collist)
            collist <- gsub("\\+", ".", collist)
            group <- factor(collist)
            if(input$pair_prefilterON == "ON"){
              dds <- DGEList(counts = count, group = group)
              keep <- filterByExpr(dds,min.count=input$pair_prefilter,min.total.count=input$pair_prefilterTotal)
              dds = dds[keep, , keep.lib.sizes=FALSE]
              count <- dds$counts
            }
            count <- log(count + 1,2)
            eset = new("ExpressionSet", exprs=as.matrix(count))
            design <- model.matrix(~0+collist)
            colnames(design) <- gsub("^collist","",colnames(design))
            fit <- lmFit(eset, design)
            comparisons <-  paste(unique(collist)[1],"-",unique(collist)[2],sep="")
            cont.matrix <- makeContrasts(contrasts=comparisons, levels=design)
            fit <- contrasts.fit(fit, cont.matrix)
            if(input$limma_trend == TRUE) fit2 <- eBayes(fit,trend = TRUE,robust = input$regression_mode)
            if(input$limma_trend == FALSE) fit2 <- eBayes(fit,trend = FALSE,robust = input$regression_mode)
            res =topTable(fit2, number = 1e12)
            colnames(res) <- c("log2FoldChange","baseMean","t","pval","padj","B")
            res$baseMean <- 2^res$baseMean
          })
        }
        if(input$DEG_method == "EBSeq"){
          withProgress(message = paste("DEG analysis of", name),{
            count <- data.matrix(count)
            vec <- c()
            for (i in 1:length(unique(collist))) {
              num <- length(collist[collist == unique(collist)[i]])
              vec <- c(vec, num)
            }
            ngvector <- NULL
            conditions <- as.factor(rep(paste("C", 1:length(unique(collist)), sep=""), times = vec))
            Sizes <- MedianNorm(count)
            EBOut <- NULL
            EBOut <- EBTest(Data = count, NgVector = ngvector, Conditions = conditions, sizeFactors = Sizes, maxround = 5)
            stopifnot(!is.null(EBOut))
            PP <- as.data.frame(GetPPMat(EBOut))
            fc_res <- PostFC(EBOut)
            results <- cbind(PP, fc_res$PostFC, fc_res$RealFC,EBOut$Mean[rownames(PP),1], EBOut$Mean[rownames(PP),2])
            colnames(results) <- c("PPEE", "PPDE", "PostFC", "RealFC","C1Mean","C2Mean")
            res <- results[order(results[,"PPDE"], decreasing = TRUE),]
            incProgress(1)
          })
        }
        count_list[name] <- list(as.data.frame(res))
      }
    }
    return(count_list)
  })

  deg_result_batch <- reactive({
    count_list <- deg_result_batch_raw()
    if(is.null(count_list) || input$Species == "not selected"){
      return(count_list)
    }
    for (name in names(count_list)) {
      count_list[[name]] <- list(
        ensembl2symbol(gene_type=gene_type1_batch(),Species=input$Species,
                       Ortholog=ortholog1_batch(),Isoform=isoform1_batch(),
                       data = count_list[[name]], org = org1())
      )
    }
    return(count_list)
  })
  
  deg_result_batch_limma <- reactive({
    if (!is.null(deg_result()) && input$DEG_method == "limma") {
      res <- deg_result()
      return(rename_limma_display_columns(res))
    }
  })
  
  deg_norm_count_batch_raw <- reactive({
    if(!is.null(norm_count_matrix())){
      files <- norm_count_matrix()
      name_list <- names(files) %>% sort()
      files <- files[name_list]
      names(files) <- names(batch_files())
      files["combined"] <- list(batch_count_combined())
      return(files)
    }
    count_files <- batch_files()
    if(is.null(count_files)){
      return(NULL)
    }
    count_files["combined"] <- list(batch_count_combined())
    count_list <- list()
    for (name in names(count_files)) {
      count <- count_files[[name]]
      collist <- gsub("\\_.+$", "", colnames(count))
      if (input$DEG_method == "DESeq2") {
        if(name != "combined"){
          dds <- dds_batch()[[name]]
        }else{
          meta <- data.frame(condition = factor(collist))
          dds<- DESeqDataSetFromMatrix(countData = round(count),colData = meta, design = ~ condition)
          dds$meta <- factor(paste0(dds[[colnames(meta)[1]]], dds[[colnames(meta)[2]]]))
          design(dds) <- ~ meta
          dds <- DESeq(dds,test = "LRT", full = ~ meta, reduced = ~ 1)
        }
        normalized_counts <- counts(dds, normalized=TRUE)
      }
      if (input$DEG_method == "edgeR") {
        if(name != "combined"){
          dds <- dds_batch()[[name]]
          normalized_counts <- t(t(dds$pseudo.counts)*(dds$samples$norm.factors))
        }else{
          group <- factor(collist)
          dds <- DGEList(counts = count, group = group)
          if(input$pair_prefilterON == "ON"){
            keep <- filterByExpr(dds,min.count=input$pair_prefilter,min.total.count=input$pair_prefilterTotal)
            dds = dds[keep, , keep.lib.sizes=FALSE]
          }
          dds <- calcNormFactors(dds, method = "TMM")
          normalized_counts <- cpm(dds)
        }
      }
      if(input$DEG_method == "limma"){
        normalized_counts <- count
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
      count_list[name] <- list(normalized_counts)
    }
    return(count_list)
  })

  deg_norm_count_batch <- reactive({
    count_list <- deg_norm_count_batch_raw()
    if(is.null(count_list) || input$Species == "not selected"){
      return(count_list)
    }
    norm_list <- list()
    for(name in names(count_list)){
      norm_list[name] <- list(
        ensembl2symbol(gene_type=gene_type1_batch(),Species=input$Species,
                       Ortholog=ortholog1_batch(),Isoform=isoform1_batch(),
                       data = count_list[[name]], org = org1())
      )
    }
    return(norm_list)
  })
  
  batch_count_combined <- reactive({
    if(!is.null(batch_files())){
      if(!is.null(norm_count_matrix())){
        files <- norm_count_matrix()
      }else{
        files <- batch_files()
      }
      if(is.null(files)){
        return(NULL)
      }else{
        matrix_list <- list()
        for (name in names(files)) {
          matrix <- as.data.frame(files[name])
          matrix_2 <- matrix
          matrix_3 <- merge(matrix, matrix_2, by = 0)[,-2:-(1 + length(colnames(matrix)))]
          matrix_list[name] <- list(matrix_3)
        }
        base <- Reduce(function(x, y) merge(x, y, by = "Row.names"), matrix_list)
        rownames(base) <- base$Row.names
        base <- data.matrix(base[,-1])
        colnames(base) <- gsub("\\.y$", "", colnames(base))
        if(length(names(files)) == 1) {
          colnames(base) <- gsub(names(files),"",colnames(base))
          colnames(base) <- gsub("^\\.", "", colnames(base))
        }
        return(base)
      }
    }else return(NULL)
  })
  
  gene_ID_pair_batch <- reactive({
    res1 <- batch_files()
    reslist<-list()
    if(is.null(res1)){
      return(NULL)
    }else{
      if(input$Species != "not selected"){
        for (name in names(res1)) {
          res <- res1[[name]]
          if(gene_type1_batch() != "SYMBOL"){
            gene_IDs <- ensembl2symbol(gene_type=gene_type1_batch(),Species=input$Species,
                                       Ortholog=ortholog1_batch(),Isoform=isoform1_batch(),data = res, org = org1(),merge=FALSE)
            reslist[name] <- list(gene_IDs)
          }else{ reslist[name] <- list(NULL) }
        }
        return(reslist) 
      }
    }
  })
  
  
  data_degcount_batch <- reactive({
    data1 <- deg_result_batch()
    count1 <- deg_norm_count_batch()
    deglist <- list()
    if(is.null(count1) || is.null(data1)){
      return(NULL)
    }else{
      for (name in names(data1)) {
        if(name != "combined"){
          data <- data1[[name]]
          count <- count1[[name]]
          if(gene_type1_batch() != "SYMBOL"){
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
          if (Type == "limma"){
            data$log2FoldChange <- -1 * data$log2FoldChange
          }
          if(gene_type1_batch() != "SYMBOL"){
            if(input$Species != "not selected"){
              if(sum(is.element(no_orgDb, input$Species))==1){
                gene_IDs <- ortholog1_batch()
              }else if(gene_type1_batch() == "isoform"){
                gene_IDs <- isoform1_batch()
              }else{
                my.symbols <- data$Row.names
                if(str_detect(my.symbols[1], "^AT.G")) key = "TAIR" else key = "ENSEMBL"
                if(org1()$packageName == "org.Sc.sgd.db") SYMBOL <- "GENENAME" else SYMBOL <- "SYMBOL"
                gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                                keytype = key,
                                                columns = c(key,SYMBOL, "ENTREZID"))
              }
              colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
              gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
              data <- merge(data, gene_IDs, by="Row.names")
              data$Unique_ID <- paste(data$SYMBOL,data$Row.names, sep = "\n- ")
            }
          }else{
            if(input$Species != "not selected"){
              if(sum(is.element(no_orgDb, input$Species))==1){
                gene_IDs <- ortholog1_batch()[,-1]
              }else if(gene_type1_batch() == "isoform"){
                gene_IDs <- isoform1_batch()[,-1]
              }else{
                my.symbols <- data$Row.names
                if(org1()$packageName == "org.Sc.sgd.db") SYMBOL <- "GENENAME" else SYMBOL <- "SYMBOL"
                gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                                keytype = SYMBOL,
                                                columns = c(SYMBOL, "ENTREZID"))
              }
              colnames(gene_IDs) <- c("Row.names", "ENTREZID")
              gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
              data <- merge(data, gene_IDs, by="Row.names")
            }
          }
          deglist[name] <- list(data)
        }
      }
      return(deglist)
    }
  })
  
  data_degcount2_batch <- reactive({
    data1 <- data_degcount_batch()
    count1 <- deg_norm_count_batch()
    deglist <- list()
    if(is.null(count1) || is.null(data1)){
      return(NULL)
    }else{
      for (name in names(data1)) {
        if(name != "combined"){
          data <- data1[[name]]
          count <- count1[[name]]
          collist <- factor(gsub("\\_.+$", "", colnames(count)))
          vec <- c()
          for (i in 1:length(unique(collist))) {
            num <- length(collist[collist == unique(collist)[i]])
            vec <- c(vec, num)
          }
          Cond_1 <- vec[1]
          Cond_2 <- vec[2]
          data2 <- dplyr::filter(data, abs(data$log2FoldChange) > log(input$fc, 2))
          if(nrow(data2) == 0){
            deglist[name] <- list(NULL)
          }else{
            data2$group <- "Up"
            data2$group[data2$log2FoldChange < 0] <- "Down"
            data3 <- dplyr::filter(data2, abs(data2$padj) < input$fdr)
            deglist[name]<- list(data3)
          }
        }
      }
      return(deglist)
    }
  })
  
  data_degcount_up_batch <- reactive({
    data1 <- data_degcount_batch()
    count1 <- deg_norm_count_batch()
    deglist<-list()
    if(is.null(count1) || is.null(data1)){
      return(NULL)
    }else{
      for (name in names(data1)) {
        if(name != "combined"){
          data <- data1[[name]]
          count <- count1[[name]]
          collist <- factor(gsub("\\_.+$", "", colnames(count)))
          vec <- c()
          for (i in 1:length(unique(collist))) {
            num <- length(collist[collist == unique(collist)[i]])
            vec <- c(vec, num)
          }
          Cond_1 <- vec[1]
          Cond_2 <- vec[2]
          data2 <- as.data.frame(data_degcount2_batch()[[name]])
          colnames(data2) <- sub(paste0(name,"."), "", colnames(data2))
          up_all <- dplyr::filter(data2, log2FoldChange > 0)
          rownames(up_all) <- up_all$Row.names
          up_all <- up_all[,8:(7 + Cond_1 + Cond_2)]
          if(input$Species != "not selected"){
            if(gene_type1_batch() != "SYMBOL"){
              up_all <- merge(up_all, gene_ID_pair_batch(), by=0)
              rownames(up_all) <- up_all$Row.names
              up_all <- up_all[,2:(1 + Cond_1 + Cond_2)]
            }
          }
          deglist[name] <- list(up_all)
        }
      }
      return(deglist)
    }
  })
  
  data_degcount_down_batch <- reactive({
    data1 <- data_degcount_batch()
    count1 <- deg_norm_count_batch()
    deglist<-list()
    if(is.null(count1) || is.null(data1)){
      return(NULL)
    }else{
      for (name in names(data1)) {
        if(name != "combined"){
          data <- data1[[name]]
          count <- count1[[name]]
          collist <- factor(gsub("\\_.+$", "", colnames(count)))
          vec <- c()
          for (i in 1:length(unique(collist))) {
            num <- length(collist[collist == unique(collist)[i]])
            vec <- c(vec, num)
          }
          Cond_1 <- vec[1]
          Cond_2 <- vec[2]
          data2 <- as.data.frame(data_degcount2_batch()[[name]])
          colnames(data2) <- sub(paste0(name,"."), "", colnames(data2))
          down_all <- dplyr::filter(data2, log2FoldChange < 0)
          rownames(down_all) <- down_all$Row.names
          down_all <- down_all[,8:(7 + Cond_1 + Cond_2)]
          if(input$Species != "not selected"){
            if(gene_type1_batch() != "SYMBOL"){
              down_all <- merge(down_all, gene_ID_pair_batch(), by=0)
              rownames(down_all) <- down_all$Row.names
              down_all <- down_all[,2:(1 + Cond_1 + Cond_2)]
            }
          }
          deglist[name] <- list(down_all)
        }
      }
      return(deglist)
    }
  })
  
  ma_heatmap_plot_batch <- reactive({
    data1 <- data_degcount_batch()
    count1 <- deg_norm_count_batch()
    malist<-list()
    if(is.null(count1) || is.null(data1)){
      return(NULL)
    }else{
      for (name in names(data1)) {
        if(name != "combined"){
          data <- data1[[name]]
          count <- count1[[name]]
          if(gene_type1_batch() != "SYMBOL"){
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
          if(gene_type1_batch() != "SYMBOL"){
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
          m1 <- as.grob(ggmaplot(data, fdr = input$fdr, fc = input$fc, size = 0.1,
                                 palette = c("#B31B21", "#1465AC", "darkgray"),
                                 genenames = genenames,
                                 legend = "top", top = 20,
                                 font.label = c("bold", 6),font.legend = "bold",
                                 font.main = "bold",xlab = xlab,
                                 ggtheme = ggplot2::theme_minimal(),
                                 select.top.method = "fc"))
          data2 <- as.data.frame(data_degcount2_batch()[[name]])
          if(is.null(data2)){
            ht <- NULL
          }else{
            data.z <- genefilter::genescale(data2[,8:(7 + Cond_1 + Cond_2)], axis=1, method="Z")
            ht <- as.grob(Heatmap(data.z, name = "z-score",column_order = colnames(data.z),
                                  clustering_method_columns = 'ward.D2',use_raster = TRUE,
                                  show_row_names = F, show_row_dend = F,column_names_side = "top"))
          }
          p <- plot_grid(m1, ht, rel_widths = c(2, 1))
          malist[[name]] <- p
        }
      }
      return(malist)
    }
  })
  
  pair_pca_plot_batch <- reactive({
    data1 <- deg_norm_count_batch()
    pca_list <- list()
    for (name in names(data1)) {
      data <- data1[[name]]
      p2 <- PCAplot(data = data,legend = input$PCA_legend)
      pca_list[[name]] <- p2
    }
    return(pca_list)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%m%d_"),"pair-wise_batch_",
             input$DEG_method,"_fc",input$fc,"_fdr",input$fdr,"_basemean",input$basemean,".zip")
    },
    content = function(fname) {
      withProgress(message = "Preparing download, please wait",{
        fs <- c()
        setwd(tempdir())
        dir.create("DEG_result",showWarnings = FALSE)
        dir.create("normalized_count",showWarnings = FALSE)
        dir.create("up",showWarnings = FALSE)
        dir.create("down",showWarnings = FALSE)
        dir.create("MAplot",showWarnings = FALSE)
        dir.create("clustering",showWarnings = FALSE)
        files <- batch_files()
        if (input$DEG_method == "limma") {
          deg_r <- deg_result_batch_limma()
        }else{
          deg_r <- deg_result_batch()
        }
        norm_c <- deg_norm_count_batch()
        up_c <- data_degcount_up_batch()
        down_c <- data_degcount_down_batch()
        ma_r <- ma_heatmap_plot_batch()
        pca_r <- pair_pca_plot_batch()
        for(name in names(norm_c)){
          DEG <- paste0("DEG_result/" ,paste0(name, ".txt"))
          norm_count <- paste0("normalized_count/" ,paste0(name, ".txt"))
          UP <- paste0("up/" ,paste0(name, ".txt"))
          DOWN <- paste0("down/" ,paste0(name, ".txt"))
          MA <- paste0("MAplot/" ,paste0(name, ".pdf"))
          PCA <- paste0("clustering/" ,paste0(name, ".pdf"))
          fs <- c(fs, DEG, norm_count, UP, DOWN, MA, PCA)
          deg <- deg_r[[name]]
          norm <- norm_c[[name]]
          up <- up_c[[name]]
          down <- down_c[[name]]
          if(name != "combined"){
            write.table(deg, DEG, quote = F, row.names = T, col.names=NA, sep = "\t")
            write.table(up, UP, quote = F, row.names = T, col.names=NA, sep = "\t")
            write.table(down, DOWN, quote = F, row.names = T, col.names=NA, sep = "\t")
            pdf(MA, height = 4, width = 7)
            print(ma_r[[name]])
            dev.off() 
          }
          write.table(norm, norm_count, quote = F, row.names = T, col.names=NA, sep = "\t")
          pdf(PCA, height = 3.5, width = 9)
          print(pca_r[[name]])
          dev.off()
        }
        report <- paste0(format(Sys.time(), "%Y%m%d_"),"pairwise_report",".docx")
        fs <- c(fs,report)
        rmarkdown::render("pair_batch_report.Rmd", output_format = "word_document", output_file = report,
                          params = list(batch_files = batch_files(),
                                        deg_norm_count_batch = deg_norm_count_batch(),
                                        input = input,
                                        norm_count_matrix = norm_count_matrix()), 
                          envir = new.env(parent = globalenv()),intermediates_dir = tempdir(),encoding="utf-8"
        )
        
        zip(zipfile=fname, files=fs)
      })
    },
    contentType = "application/zip"
  )
  
  #multi DEG--------------------------
  org6 <- reactive({
    return(org(Species = input$Species6,Ortholog = input$Ortholog6))
  })
  ortholog6 <- reactive({
    return(no_org_ID(count = multi_row_count_matrix(),Species = input$Species6,Ortholog = input$Ortholog6,Biomart_archive=input$Biomart_archive6))
  })
  isoform6 <- reactive({
    return(isoform_ID(count = multi_row_count_matrix(),Species = input$Species6,Ortholog = input$Ortholog6,Biomart_archive=input$Biomart_archive6,RNA_type=input$Level_multi))
  })
  gene_type6 <- reactive({
    return(gene_type(my.symbols=rownames(multi_row_count_matrix()),org=org6(),Species=input$Species6,RNA_type=input$Level_multi))
  })
  org_code6 <- reactive({
    return(org_code(Species = input$Species6, Ortholog= input$Ortholog6))
  })
  
  multi_row_count_matrix <- reactive({
    withProgress(message = "Importing row count matrix, please wait",{
      if (input$multi_data_file_type == "Row1"){
        tmp <- input$multi_file1$datapath
        if(is.null(input$multi_file1) && input$goButton6 > 0 )  tmp = example_data_path("data/example4.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example4.txt")
        return(read_df(tmp = tmp))
      }
      if (input$multi_data_file_type == "Row2"){
        tmp <- input$multi_file2$datapath
        if(is.null(input$multi_file2) && input$goButton6 > 0 )  tmp = example_data_path("data/example8.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example8.txt")
        return(read_df(tmp = tmp))
      }
    })
  })
  multi_metadata <- reactive({
    if (input$multi_data_file_type != "Row2"){
      return(NULL)
    }else{
      tmp <- input$multi_file3$datapath
      if(is.null(input$multi_file3) && input$goButton6 > 0 )  tmp = example_data_path("data/example5.csv", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example5.csv")
      df <- read_df(tmp = tmp)
      if(!is.null(df)){
        rownames(df) <- gsub("-",".",rownames(df))
        df[,1] <- gsub("\\_.+$", "", df[,1])
        df[,2] <- gsub("\\_.+$", "", df[,2])
      }
      return(df)
    }
  })
  
  multi_d_row_count_matrix <- reactive({
    row <- multi_row_count_matrix()
    if (input$multi_data_file_type == "Row1"){
      if(is.null(row)) {
        return(NULL)
      }else{
        return(anno_rep(row))
      }
    }else{
      meta <- multi_metadata()
      if (is.null(row) || is.null(meta)){
        return(NULL)
      } else {
        row_t <- t(row)
        colname <- colnames(meta)
        data <- merge(meta, row_t, by=0, sort = F)
        if(dim(data)[1] == 0) {
          rownames(meta) <- gsub("\\.","-",rownames(meta))
          data <- merge(meta, row_t, by=0, sort = F)
          if(dim(data)[1] == 0) {
            rownames(row_t) <- gsub("\\.","-",rownames(row_t))
            data <- merge(meta, row_t, by=0, sort = F)
            validate("Error: failed to merge count data with metadata. Please check row names of matadata.")
          }
        }
        rownames(data) <- data[,1]
        data2 <- data[, - which(colnames(data) %in% c("Row.names", colname))]
        data2 <- data2[,1:length(rownames(row))]
        data2_t <- t(data2)
        data3 <- apply(data2_t, 2, as.numeric)
        rownames(data3) <- rownames(data2_t)
        if(input$Species6 != "not selected"){
          if(gene_type6() != "SYMBOL"){
            rownames(data3) < gsub("\\..*","", rownames(data3))
          }
        }
        return(data3)
      }
    }
  })
  
  multi_norm_count_matrix <- reactive({
    if(is.null(input$multi_norm_file1)){
      return(NULL)
    }else{
      if(length(input$multi_norm_file1[, 1]) == 1){
        upload <- anno_rep(read_df(input$multi_norm_file1[[1, 'datapath']]))
      }else{
        upload = list()
        for(nr in 1:length(input$multi_norm_file1[, 1])){
          df <- anno_rep(read_df(input$multi_norm_file1[[nr, 'datapath']]))
          upload[gsub("\\..+$", "", input$multi_norm_file1[nr,]$name)] <- list(df)
        }
      } 
      return(upload)
    }
  })
  
  
  # Multi DEG ------------------------------------------------------------------------------
  observeEvent(input$DEG_method,({
    if(input$DEG_method == "limma") updateNumericInput(session,"multi_prefilter", "Minimum count required for at least some samples", value = 0)
  }))
  multi_analysis_requested <- reactiveVal(FALSE)

  observeEvent(input$multi_main_tab, {
    if(!is.null(input$multi_main_tab) && input$multi_main_tab != "multi_input_data_tab"){
      multi_analysis_requested(TRUE)
    }
  }, ignoreInit = TRUE)

  observeEvent(
    list(
      input$multi_file1,
      input$multi_file2,
      input$multi_file3,
      input$multi_data_file_type,
      input$DEG_method_multi,
      input$Level_multi,
      input$FDR_method6,
      input$multi_prefilterON,
      input$multi_prefilter,
      input$multi_prefilterTotal,
      input$limma_trend_multi,
      input$regression_mode_multi,
      input$multi_norm_file1,
      input$goButton6
    ),
    {
      if(is.null(input$multi_main_tab) || input$multi_main_tab == "multi_input_data_tab"){
        multi_analysis_requested(FALSE)
      }
    },
    ignoreInit = TRUE
  )

  multi_dds <- reactive({
    if(!isTRUE(multi_analysis_requested())){
      return(NULL)
    }
    count <- multi_d_row_count_matrix()
    meta <- multi_metadata()
    if(is.null(count)){
      return(NULL)
    }else{
      withProgress(message = "DESeq2",{
        if (input$multi_data_file_type == "Row1"){
          collist <- gsub("\\_.+$", "", colnames(count))
          meta <- data.frame(condition = factor(collist))
        }else meta <- data.frame(condition=factor(meta[,1]), type=factor(meta[,2]))
        dds<- DESeqDataSetFromMatrix(countData = round(count),colData = meta, design = ~ condition)
        dds$meta <- factor(paste0(dds[[colnames(meta)[1]]], dds[[colnames(meta)[2]]]))
        design(dds) <- ~ meta
        dds <- DESeq(dds,test = "LRT", full = ~ meta, reduced = ~ 1)
        incProgress(1)
      })
      return(dds)
    }
  })
  
  multi_deg_result_raw <- reactive({
    if(!isTRUE(multi_analysis_requested())){
      return(NULL)
    }
    if(is.null(multi_d_row_count_matrix())){
      return(NULL)
    }else{
        count <- multi_d_row_count_matrix()
        meta <- multi_metadata()
        if(input$DEG_method_multi == "DESeq2"){
          withProgress(message = "Prepare a DEG result",{
        dds <- multi_dds()
        df <- list()
        for(i in 1:choose(n=length(unique(dds$meta)),k=2)){
          res <- DESeq2::results(dds, contrast = c("meta", as.character(unique(dds$meta)[combn(x=length(unique(dds$meta)),m=2)[1,i]]),as.character(unique(dds$meta)[combn(x=length(unique(dds$meta)),m=2)[2,i]])))
          foldchange <- data.frame(gene = rownames(res), log2FC = res$log2FoldChange)
          df[res@elementMetadata@listData[["description"]][2]] <- list(foldchange)
        }
        fc_matrix <- as.data.frame(df[names(df)[1]])
        colnames(fc_matrix) <- sub("log2.fold.change..MLE...meta.","",colnames(fc_matrix))
        colnames(fc_matrix)[1] <- "gene"
        for(i in 1:(length(names(df))-1)){
          matrix <- as.data.frame(df[names(df)[1+i]])
          colnames(matrix) <- sub("log2.fold.change..MLE...meta.","",colnames(matrix))
          colnames(matrix)[1] <- "gene"
          fc_matrix <- merge(fc_matrix, matrix, by="gene")
        }
        if(input$FDR_method6 == "IHW") {
          ihw_res <- IHW::ihw(pvalue ~ baseMean,  data=as.data.frame(res), alpha = 0.1)
          res$padj <- IHW::as.data.frame(ihw_res)$adj_pvalue
        }
        if(input$FDR_method6 == "Qvalue") {
          res <- DESeq2::results(dds)
          qvalue <- qvalue::qvalue(res$pvalue)
          res$padj <- qvalue$qvalues
        }
        res <- data.frame(gene=rownames(res), padj=res$padj)
        fc_matrix <- merge(fc_matrix, res, by="gene")
        rownames(fc_matrix)<- fc_matrix$gene
        fc_matrix <- fc_matrix[,-1]
        
        return(fc_matrix)
        })
        }else{
          collist <- gsub("\\_.+$", "", colnames(count))
          collist <- gsub(" ", ".", collist)
          group <- factor(collist)
          
          if(input$multi_prefilterON == "ON"){
          ##pre-filtering
          dds <- DGEList(counts = count, group = group)
          keep <- filterByExpr(dds,min.count=input$multi_prefilter,min.total.count=input$multi_prefilterTotal)
          dds = dds[keep, , keep.lib.sizes=FALSE]
          count <- dds$counts
          }
          
          count <- log(count + 1,2)
          if(is.element(TRUE, duplicated(collist)) == TRUE){
              print("limma")
            if (input$multi_data_file_type == "Row1"){
              meta <- data.frame(condition = factor(collist))
              pair <- collist
            }else {
              meta <- data.frame(condition=factor(meta[,1]), type=factor(meta[,2]))
              pair <- paste0(meta$condition,meta$type)
            }
            pair <- gsub("\\+", ".", pair)
            design <- model.matrix(~0+pair)
            colnames(design) <- factor(gsub("^pair","",colnames(design)),levels = unique(pair))
            cont <- c()
              for(i in 1:choose(n=length(unique(pair)),k=2)){
                contrast = paste0(as.character(unique(pair)[combn(x=length(unique(pair)),m=2)[1,i]]),"-",as.character(unique(pair)[combn(x=length(unique(pair)),m=2)[2,i]]))
                cont <- c(cont,contrast)
              }
              cont.matrix <- makeContrasts(contrasts=cont, levels=design)
              eset = new("ExpressionSet", exprs=as.matrix(count))
              fit1 <- lmFit(eset,design)
              fit2 <- contrasts.fit(fit1, cont.matrix)
              if(input$limma_trend_multi == TRUE) fit3 <- eBayes(fit2,trend = TRUE,robust = input$regression_mode_multi)
              if(input$limma_trend_multi == FALSE) fit3 <- eBayes(fit2,trend = FALSE,robust = input$regression_mode_multi)
              result <- topTable(fit3,coef=1:length(cont), number = 1e12)
              lab <- paste0("log2(",cont,")")
              lab <- gsub("-","/",lab)
              if(length(cont) != 1) label <- c(lab,"AveExpr","F","p_value","padj") else label <- c(lab,"AveExpr","F","p_value","padj","B")
              colnames(result) <- label
              print(label)
            return(result)
          }
          
        }
    }
  })
  
  multi_deg_result <- reactive({
    res <- multi_deg_result_raw()
    if(is.null(res)){
      return(NULL)
    }
    if(input$Species6 != "not selected") {
      res <- ensembl2symbol(gene_type=gene_type6(),Species=input$Species6,
                            Ortholog=ortholog6(),Isoform=isoform6(),data = res, org = org6())
    }
    return(res)
  })

  multi_deg_norm_count_raw <- reactive({
    if(!isTRUE(multi_analysis_requested())){
      return(NULL)
    }
    if(is.null(multi_d_row_count_matrix())){
      return(NULL)
    }else{
      if(input$DEG_method_multi == "DESeq2"){
        if(!is.null(multi_norm_count_matrix())){
          return(multi_norm_count_matrix())
        }else {
          dds <- multi_dds()
          normalized_counts <- counts(dds, normalized=TRUE)
          return(normalized_counts)
        }
      }else{
        return(multi_d_row_count_matrix())
      }
    }
  })

  multi_deg_norm_count <- reactive({
    normalized_counts <- multi_deg_norm_count_raw()
    if(is.null(normalized_counts)){
      return(NULL)
    }
    if(input$Species6 != "not selected"){
      normalized_counts <- ensembl2symbol(gene_type=gene_type6(),data = normalized_counts,
                                          Species=input$Species6,Ortholog=ortholog6(),
                                          Isoform=isoform6(),org = org6())
    }
    return(normalized_counts)
  })
  
  observeEvent(input$goButton6,({
    updateSelectInput(session,inputId = "Species6","Species",species_list, selected = "Homo sapiens")
  }))
  
  observeEvent(input$multi_file1, ({
    updateCollapse(session,id =  "multi_input_collapse_panel", open="multi_Row_count_panel")
  }))
  observeEvent(input$multi_file2, ({
    updateCollapse(session,id =  "multi_input_collapse_panel", open="multi_Metadata_panel")
  }))
  output$multi_Row_count_matrix <- DT::renderDataTable({
    multi_row_count_matrix()
  })
  output$multi_d_Row_count_matrix <- DT::renderDataTable({
    multi_d_row_count_matrix()
  })
  
  output$multi_Metadata <- DT::renderDataTable({
    multi_metadata()
  })
  output$multi_DEG_result <- DT::renderDataTable({
    multi_deg_result()
  })
  output$multi_Normalized_Count_matrix <- DT::renderDataTable({
    multi_deg_norm_count()
  })
  
  output$download_multi_DEG_result = downloadHandler(
    filename = function() {
      if (input$multi_data_file_type == "Row1"){
        paste(gsub("\\..+$", "", input$multi_file1), paste0(input$FDR_method6,".txt"), sep ="-")
      }else{
        paste(gsub("\\..+$", "", input$multi_file2), paste0(input$FDR_method6,".txt"), sep ="-")
      }},
    content = function(file){write.table(multi_deg_result(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  output$download_multi_norm_count = downloadHandler(
    filename = function() {
      if (input$multi_data_file_type == "Row1"){
        paste(gsub("\\..+$", "", input$multi_file1), "normalized.txt", sep ="-")
      }else{
        paste(gsub("\\..+$", "", input$multi_file2), "normalized.txt", sep ="-")
      }},
    content = function(file){write.table(multi_deg_norm_count(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  download_multi_overview_dir <-reactive({
    if (input$multi_data_file_type == "Row1"){
      dir_name <- paste(gsub("\\..+$", "", input$multi_file1), input$FDR_method6, sep ="-")
    }else{
      dir_name <- paste(gsub("\\..+$", "", input$multi_file2), input$FDR_method6, sep ="-")
    }
    dir_name <- paste0(dir_name, paste0("_fdr", input$fdr6))
    dir_name <- paste0(dir_name, paste0("_basemean", input$basemean6))
    dir_name <- paste0(dir_name, paste0("_fc", input$fc6))
    return(dir_name)
  })
  
  #multi DEG vis---------------------------
  updateCounter_multi <- reactiveValues(i = 0)
  
  observe({
    input$start_multi
    isolate({
      updateCounter_multi$i <- updateCounter_multi$i + 1
    })
  })
  
  
  #Restart
  observeEvent(multi_pattern1(), {
    isolate(updateCounter_multi$i == 0)
    updateCounter_multi$i <- 0
  }) 
  
  output$multi_DEG_total1 <- renderText({
    if(is.null(multi_pattern1())){
      return(NULL)
    }else{ 
      if(length(input$selectFC) == 2 && length(input$selectFC_2) == 2){
        print(paste0("The number of genes after the filtration (fdr < ", input$fdr6,", 
                     basemean > ",input$basemean6,", |log2(", 
                     input$selectFC[1],"/", input$selectFC[2],")| > ", log2(input$fc6),", |log2(",
                     input$selectFC_2[1],"/", input$selectFC_2[2],")| > ", log2(input$fc6),"): ", 
                     length(multi_pattern1()$gene)))
      }else if(length(input$selectFC) == 2 && length(input$selectFC_2) != 2){
        print(paste0("The number of genes after the filtration (fdr < ", input$fdr6,", 
                     basemean > ",input$basemean6,", |log2(", 
                     input$selectFC[1],"/", input$selectFC[2],")| > ", log2(input$fc6),"): ", 
                     length(multi_pattern1()$gene)))
      }else if(length(input$selectFC) != 2 && length(input$selectFC_2) != 2){
        print(paste0("The number of genes after the filtration (fdr < ", input$fdr6,", 
                     basemean > ",input$basemean6,"): ", 
                     length(multi_pattern1()$gene)))
      }
    }
  })
  
  multi_gene_id_map <- reactive({
    res <- multi_deg_result_raw()
    if(is.null(res) || input$Species6 == "not selected"){
      return(NULL)
    }
    my.symbols <- rownames(res)
    if(length(my.symbols) == 0){
      return(NULL)
    }
    if(gene_type6() != "SYMBOL"){
      if(sum(is.element(no_orgDb, input$Species6)) == 1){
        gene_IDs <- ortholog6()
      }else if(gene_type6() == "isoform"){
        gene_IDs <- isoform6()
      }else{
        if(str_detect(my.symbols[1], "^AT.G")) key = "TAIR" else key = "ENSEMBL"
        if(org6()$packageName == "org.Sc.sgd.db") SYMBOL <- "GENENAME" else SYMBOL <- "SYMBOL"
        gene_IDs<-AnnotationDbi::select(org6(),keys = my.symbols,
                                        keytype = key,
                                        columns = c(key,SYMBOL, "ENTREZID"))
      }
      colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
    }else{
      if(sum(is.element(no_orgDb, input$Species6)) == 1){
        gene_IDs <- ortholog6()[,-1]
      }else if(gene_type6() == "isoform"){
        gene_IDs <- isoform6()[,-1]
      }else{
        if(org6()$packageName == "org.Sc.sgd.db") SYMBOL <- "GENENAME" else SYMBOL <- "SYMBOL"
        gene_IDs<-AnnotationDbi::select(org6(),keys = my.symbols,
                                        keytype = SYMBOL,
                                        columns = c(SYMBOL, "ENTREZID"))
      }
      colnames(gene_IDs) <- c("Row.names", "ENTREZID")
    }
    gene_IDs %>% distinct(Row.names, .keep_all = T)
  })

  multi_gene_ID_pair <- reactive({
    gene_IDs <- multi_gene_id_map()
    if(is.null(gene_IDs) || gene_type6() == "SYMBOL" || !"SYMBOL" %in% colnames(gene_IDs)){
      return(NULL)
    }
    dplyr::select(gene_IDs, Row.names, SYMBOL)
  })
  
  multi_select <- reactive({
    if(input$DEG_method_multi == "DESeq2"){
      dds <- multi_dds()
      if(is.null(dds)){
        return(NULL)
      }
      pair <- dds$meta
    } else{
      count <- multi_d_row_count_matrix()
      meta <- multi_metadata()
      collist <- gsub("\\_.+$", "", colnames(count))
      collist <- gsub(" ", ".", collist)
        if (input$multi_data_file_type == "Row1"){
          meta <- data.frame(condition = factor(collist))
          pair <- collist
        }else {
          meta <- data.frame(condition=factor(meta[,1]), type=factor(meta[,2]))
          pair <- paste0(meta$condition,meta$type)
        }
      pair <- gsub("\\+",".",pair)
    }
    return(unique(pair))
  })

  multi_norm_count_filter_base <- reactive({
    data <- multi_deg_norm_count_raw()
    if(is.null(data)){
      return(NULL)
    }
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    data <- data[, !colnames(data) %in% c("SYMBOL", "Unique_ID", "ENTREZID"), drop = FALSE]
    keep <- rowMeans(as.matrix(data), na.rm = TRUE) > input$basemean6
    data[keep, , drop = FALSE]
  })

  multi_pair_column_info <- function(data, pair) {
    if(is.null(data) || length(pair) != 2){
      return(NULL)
    }
    candidates <- list(
      list(column = paste0(pair[1], ".vs.", pair[2], ".log2FC"), sign = 1),
      list(column = paste0(pair[2], ".vs.", pair[1], ".log2FC"), sign = -1),
      list(column = paste0("log2(", pair[1], "/", pair[2], ")"), sign = 1),
      list(column = paste0("log2(", pair[2], "/", pair[1], ")"), sign = -1)
    )
    for(candidate in candidates){
      if(candidate$column %in% colnames(data)){
        return(candidate)
      }
    }
    NULL
  }

  multi_pair_fc_table_raw <- function(data, pair) {
    info <- multi_pair_column_info(data, pair)
    if(is.null(info)){
      return(NULL)
    }
    data.frame(
      Row.names = rownames(data),
      log2FoldChange = info$sign * data[[info$column]],
      stringsAsFactors = FALSE
    )
  }

  multi_filtered_pair_result <- function(primary_pair = NULL, secondary_pair = NULL) {
    res <- multi_deg_result_raw()
    filtered_count <- multi_norm_count_filter_base()
    if(is.null(res) || is.null(filtered_count)){
      return(NULL)
    }
    res <- as.data.frame(res, stringsAsFactors = FALSE)
    out <- data.frame(
      gene = rownames(res),
      padj = res$padj,
      stringsAsFactors = FALSE
    )
    out <- dplyr::filter(out, !is.na(padj), padj < input$fdr6, gene %in% rownames(filtered_count))

    if(length(primary_pair) == 2){
      primary_tbl <- multi_pair_fc_table_raw(res, primary_pair)
      if(is.null(primary_tbl)){
        return(NULL)
      }
      out <- merge(out, primary_tbl, by.x = "gene", by.y = "Row.names")
      out <- dplyr::filter(out, !is.na(log2FoldChange), abs(log2FoldChange) > log2(input$fc6))
    } else {
      out$log2FoldChange <- NA_real_
    }

    if(length(secondary_pair) == 2){
      secondary_tbl <- multi_pair_fc_table_raw(res, secondary_pair)
      if(is.null(secondary_tbl)){
        return(NULL)
      }
      secondary_tbl <- dplyr::filter(
        secondary_tbl,
        !is.na(log2FoldChange),
        abs(log2FoldChange) > log2(input$fc6)
      )
      out <- merge(out, dplyr::select(secondary_tbl, Row.names), by.x = "gene", by.y = "Row.names")
    }

    tibble::as_tibble(out)
  }

  multi_cluster_transform <- reactive({
    if(input$DEG_method_multi == "DESeq2"){
      dds <- multi_dds()
      if(is.null(dds)){
        return(NULL)
      }
      withProgress(message = "Transform counts", {
        assay(rlogTransformation(dds))
      })
    }else{
      count <- multi_d_row_count_matrix()
      if(is.null(count)){
        return(NULL)
      }
      log(count + 1, 2)
    }
  })
  
  output$selectFC <- renderUI({
    if(is.null(multi_deg_result())){
      return(NULL)
    }else{
      selectizeInput("selectFC", "Select a pair for fold change cut-off", c(multi_select()),
                     selected = "", multiple = TRUE, 
                     options = list(maxItems = 2))
    }
  })
  output$selectFC_2 <- renderUI({
    if(is.null(multi_deg_result())){
      return(NULL)
    }else{
      selectizeInput("selectFC_2", "Option: select a pair for fold change cut-off", c(multi_select()),
                     selected = "", multiple = TRUE, 
                     options = list(maxItems = 2))
    }
  })
  output$topP <- renderUI({
    if(is.null(multi_deg_result())){
      return(NULL)
    }else{
      sliderInput("topP", "Most significant genes", min = 0,
                  max=5000, step = 100,
                  value = 2000)
    }
  })
  
  multi_pattern1 <- reactive({
    if(is.null(multi_deg_result_raw())){
      return(NULL)
    }
    withProgress(message = "Fold Change, FDR, and base mean cut-off", {
      multi_filtered_pair_result(input$selectFC, input$selectFC_2)
    })
  })
  
  multi_pattern1_2 <- reactive({
    sig_res_LRT <- multi_pattern1()
    if(length(sig_res_LRT$gene) == 0 || is.null(input$topP) || input$start_multi == 0 || updateCounter_multi$i == 0){
      return(NULL)
    }else{
      withProgress(message = "Select most significant genes",{
        clustering_sig_genes <- sig_res_LRT %>%
          arrange(padj) %>%
          head(n=input$topP)
        rld_mat <- multi_cluster_transform()
        if(is.null(rld_mat)){
          return(NULL)
        }
        cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
        return(cluster_rlog)
      })
    }
  })
  
  
  multi_pattern2 <- reactive({
    count <- multi_d_row_count_matrix()
    meta <- multi_metadata()
    cluster_rlog <- multi_pattern1_2()
    if(is.null(cluster_rlog)){
      return(NULL)
    }else{
      withProgress(message = "Performing a clustering analysis",{
        if (input$multi_data_file_type == "Row1"){
          collist <- gsub("\\_.+$", "", colnames(count))
          meta <- data.frame(condition = factor(collist,levels = unique(collist),ordered = TRUE), row.names = colnames(count))
          clusters <- degPatterns(cluster_rlog, metadata = meta, time = "condition", plot = FALSE,minc = 0)
        }else{
          meta <- data.frame(condition=factor(meta[,1]), type=factor(meta[,2]))
          meta[,2] <- factor(meta[,2], levels = unique(meta[,2]),ordered = TRUE)
          rownames(meta) <- colnames(count)
          clusters <- degPatterns(cluster_rlog, metadata = meta, time = "condition", col=colnames(meta)[2], plot = FALSE,minc = 0)
        }
        rownames(clusters$df) <- rownames(cluster_rlog)
        clusters$df$genes <- rownames(cluster_rlog)
        return(clusters)
      })
    }
  })
  
  multi_boxplot_reactive <- reactive({
    meta <- multi_metadata()
    count <- multi_d_row_count_matrix()
    cluster_rlog <- multi_pattern1_2()
    clusters <- multi_pattern2()
    if(is.null(clusters)){
      return(NULL)
    }else{
      withProgress(message = "Boxplot",{
        if (input$multi_data_file_type == "Row1"){
          collist <- gsub("\\_.+$", "", colnames(count))
          meta <- data.frame(condition = factor(collist,levels = unique(collist),ordered = TRUE), row.names = colnames(count))
        }else {
          meta <- data.frame(condition=factor(meta[,1]), type=factor(meta[,2]))
          rownames(meta) <- colnames(count)
          meta[,2] <- factor(meta[,2], levels = unique(meta[,2]),ordered = TRUE)
        }
        meta$condition <- factor(meta$condition, levels = unique(meta$condition),ordered = TRUE)
        table <- rownames_to_column(as.data.frame(cluster_rlog), "genes") %>%
          gather("sample", "normalized", -genes) %>%
          right_join(distinct(clusters$df[,c("genes", "cluster")]),
                     by = "genes") %>%
          left_join(rownames_to_column(as.data.frame(meta), "sample"),
                    by = "sample") %>%
          as.data.frame()
        colnames(table)[3] <- "expression"
        table <- na.omit(table)
        if (input$multi_data_file_type == "Row1"){
          p <- degPlotCluster(table, time = "condition", process = TRUE, min_genes = 0)+ 
            scale_color_brewer(palette = "Set1", direction=-1)+
            theme_bw(base_size = 15)+ theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
        }
        else{
          p <- degPlotCluster(table, time = "condition", color=colnames(meta)[2], process = TRUE, min_genes = 0)+ 
            scale_color_brewer(palette = "Set1", direction=-1)+
            theme_bw(base_size = 15)+ theme(legend.position = "top")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
        }
        return(p)
      })
    }
  })
  
  output$multi_boxplot <- renderPlot({
    plot_obj <- multi_boxplot_reactive()
    if(is.null(plot_obj)){
      return(NULL)
    }else{
      print(plot_obj)
    }
  })
  
  observeEvent(multi_pattern2(), {
    clusters_obj <- multi_pattern2()
    choices <- character(0)
    if (!is.null(clusters_obj) && !is.null(clusters_obj$df)) {
      choices <- unique(paste0("Group", clusters_obj$df$cluster))
    }
    selected <- isolate(input$multi_selectfile1)
    if (length(choices) == 0) {
      selected <- character(0)
    } else if (is.null(selected) || !selected %in% choices) {
      selected <- choices[[1]]
    }
    freezeReactiveValue(input, "multi_selectfile1")
    updateSelectInput(session, "multi_selectfile1",
                      choices = choices, selected = selected)
  }, ignoreInit = TRUE)
  
  multi_pattern_extract <- reactive({
    data <- as.data.frame(multi_deg_norm_count())
    clusters <- multi_pattern2()$df
    if(is.null(data) || is.null(clusters)){
      return(NULL)
    }else{
      if(is.null(input$multi_selectfile1) || !nzchar(input$multi_selectfile1)){
        return(NULL)
      }else{
        clusters$cluster <- paste0("Group",clusters$cluster)
        cluster_name <- input$multi_selectfile1
        clusterCount <- dplyr::filter(clusters, cluster == cluster_name)
        clusterCount <- merge(clusterCount, data, by=0)
        clusterCount <- clusterCount[,-2:-3]
        rownames(clusterCount) <- clusterCount$Row.names
        clusterCount <- clusterCount[,-1]
        return(clusterCount)
      }
    }
  })
  
  output$multi_pattern1_list <- DT::renderDataTable({
    if(is.null(multi_pattern2())){
      return(NULL)
    }else{
      clusters <- multi_pattern2()$df
      clusters$cluster <- paste0("Group",clusters$cluster)
      data <- as.data.frame(multi_deg_norm_count())
      if(input$Species6 != "not selected"){
        if(gene_type6() != "SYMBOL"){
          collist <- gsub("\\_.+$", "", colnames(data))
          data2 <-merge(clusters,data, by=0)
          data2 <- data2[,-1]
          clusters <- data2[,-3:-(1+length(collist))]
        }
      }
      clusters
    }
  })
  output$multi_pattern1_count <- DT::renderDT({
    res <- multi_pattern_extract()
    return(res)
  })
  output$download_deg_pattern_list = downloadHandler(
    filename = function(){
      paste0(download_multi_overview_dir(), "_DEG_pattern.txt")
    },
    content = function(file) {
      clusters <- multi_pattern2()$df
      clusters$cluster <- paste0("Group",clusters$cluster)
      write.table(clusters, file, quote = F, row.names = F, sep = "\t")})
  
  
  output$download_deg_pattern_count = downloadHandler(
    filename = function(){
      paste(paste(download_multi_overview_dir(), input$multi_selectfile1, sep = "_"), 
            "_pattern_extracted.txt", sep = "_")
    },
    content = function(file) {write.table(multi_pattern_extract(), file, quote = F, row.names = T, col.names=NA, sep = "\t")})
  
  multi_GOIbox <- reactive({
    if(!is.null(input$multi_pattern1_count_rows_selected)){
      data <- multi_pattern_extract()[input$multi_pattern1_count_rows_selected,]
      if(input$Species6 != "not selected"){
        if(gene_type6() != "SYMBOL"){
          rownames(data) <- paste(data$SYMBOL,rownames(data), sep = "\n- ")
          data <- data[, - which(colnames(data) == "SYMBOL")]
        }
      }
      return(data)
    }
  })
  
  output$multi_pattern_boxplot <- renderPlot({
    if(!is.null(input$multi_pattern1_count_rows_selected)){
      GOIboxplot(data = multi_GOIbox())
    }
  })
  
  output$download_deg_pattern_boxplot = downloadHandler(
    filename = function(){
      paste(paste(download_multi_overview_dir(), input$multi_selectfile1, sep = "_"), 
            "_boxplot.pdf", sep = "_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- multi_GOIbox()
        rowlist <- rownames(data)
        if(input$multi_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)
        }else pdf_height <- input$multi_pdf_height
        if(input$multi_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)
        }else pdf_width <- input$multi_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(GOIboxplot(data = data))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  observeEvent(input$multi_pattern1_count_rows_selected, ({
    updateCollapse(session,id =  "multi_collapse_panel1", open="multi_deg_pattern_boxplot_panel")
  }))
  
  output$download_multi_boxplot = downloadHandler(
    filename = function(){
      paste0(download_multi_overview_dir(), "_pattern.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(is.null(multi_boxplot_reactive)){
          return(NULL)
        }else{
          clusters <- multi_pattern2()$df
          clusterNumber <- length(unique(clusters$cluster))
          print(clusterNumber)
          if(input$multi_pdf_height == 0){
            pdf_height <- pdf_h(clusterNumber)+2
          }else pdf_height <- input$multi_pdf_height
          if(input$multi_pdf_width == 0){
            pdf_width <- pdf_w(clusterNumber)+2
          }else pdf_width <- input$multi_pdf_width
          pdf(file, height = pdf_height, width = pdf_width)
          print(multi_boxplot_reactive()+
                  theme(axis.text.x= element_text(size = 8),
                        axis.text.y= element_text(size = 8),
                        title = element_text(size = 8),text = element_text(size = 8)))
          dev.off()
        }
        incProgress(1)
      })
    }
  )
  
  #Multi kmeans-----------
  updateCounter_kmeans_multi <- reactiveValues(i = 0)
  multi_kmeans_running <- reactiveVal(FALSE)
  multi_kmeans_status <- reactiveVal("")
  set_multi_kmeans_status <- function(message = NULL) {
    if (is.null(message) || !nzchar(message)) {
      multi_kmeans_status("")
      removeNotification(id = "multi_kmeans_progress")
      return(invisible(NULL))
    }
    message <- paste(as.character(message), collapse = "\n")
    current <- isolate(multi_kmeans_status())
    if (!identical(current, message)) {
      multi_kmeans_status(message)
      showNotification(
        message,
        id = "multi_kmeans_progress",
        duration = NULL,
        closeButton = FALSE,
        type = "message"
      )
    }
    invisible(message)
  }
  reset_multi_kmeans_state <- function() {
    updateCounter_kmeans_multi$i <- 0
    multi_kmeans_running(FALSE)
    set_multi_kmeans_status(NULL)
    freezeReactiveValue(input, "kmeans_order_multi")
    updateSelectInput(session, "kmeans_order_multi",
                      choices = character(0), selected = character(0))
    freezeReactiveValue(input, "multi_selectfile2")
    updateSelectInput(session, "multi_selectfile2",
                      choices = character(0), selected = character(0))
  }
  build_multi_kmeans_order <- function(data.z, cl) {
    if (is.null(cl) || !length(cl)) {
      return(character(0))
    }
    cluster_ids <- sort(unique(cl))
    if (length(cluster_ids) <= 1) {
      return(as.character(cluster_ids))
    }
    cl_df <- data.frame(cluster = cl, stringsAsFactors = FALSE)
    data2 <- merge(cl_df, data.z, by = 0)
    rownames(data2) <- data2[, 1]
    data2 <- data2[, -1, drop = FALSE]
    df <- data.frame(matrix(nrow = 0, ncol = ncol(data2) - 1))
    if (ncol(data2) > 1) {
      colnames(df) <- colnames(data2[, -1, drop = FALSE])
    }
    for (i in cluster_ids) {
      data3 <- data2 %>% dplyr::filter(cluster == i)
      if (nrow(data3) == 0) {
        next
      }
      data4 <- colSums(data3[, -1, drop = FALSE], na.rm = TRUE)
      df <- rbind(df, data4)
    }
    if (nrow(df) <= 1) {
      return(as.character(cluster_ids))
    }
    as.character(cluster_ids[hclust(dist(df), "average")$order])
  }
  
  observeEvent(input$kmeans_start_multi, {
    multi_kmeans_running(TRUE)
    set_multi_kmeans_status("k-means clustering is running.\nPreparing clusters and heatmap...")
    updateCounter_kmeans_multi$i <- updateCounter_kmeans_multi$i + 1
  }, ignoreInit = TRUE)
  
  output$multi_kmeans_progress_text <- renderText({
    multi_kmeans_status()
  })
  
  
  #Restart
  observeEvent(list(
    input$multi_data_file_type,
    input$multi_file1,
    input$multi_file2,
    input$multi_file3,
    input$multi_norm_file1,
    input$goButton6,
    input$DEG_method_multi,
    input$Level_multi,
    input$FDR_method6,
    input$multi_prefilterON,
    input$multi_prefilter,
    input$multi_prefilterTotal,
    input$limma_trend_multi,
    input$regression_mode_multi,
    input$selectFC2,
    input$selectFC2_2,
    input$topP2,
    input$fdr6,
    input$basemean6,
    input$fc6,
    input$multi_kmeans_number,
    input$Species6,
    input$Ortholog6,
    input$Biomart_archive6
  ), {
    reset_multi_kmeans_state()
  }, ignoreInit = TRUE) 
  
  output$selectFC2 <- renderUI({
    if(is.null(multi_deg_result())){
      return(NULL)
    }else{
      selectizeInput("selectFC2", "Select a pair for fold change cut-off", c(multi_select()),
                     selected = "", multiple = TRUE, 
                     options = list(maxItems = 2))
    }
  })
  output$selectFC2_2 <- renderUI({
    if(is.null(multi_deg_result())){
      return(NULL)
    }else{
      if(length(input$selectFC2) == 2){
      selectizeInput("selectFC2_2", "Option: select a pair for fold change cut-off", c(multi_select()),
                     selected = "", multiple = TRUE, 
                     options = list(maxItems = 2))
      }
    }
  })
  output$topP2 <- renderUI({
    if(is.null(multi_deg_result())){
      return(NULL)
    }else{
      sliderInput("topP2", "Most significant genes", min = 0,
                  max=8000, step = 100,
                  value = 2000)
    }
  })
  
  output$multi_kmeans_num <- renderUI({
    if(is.null(multi_deg_norm_count())){
      return(NULL)
    }else{
      withProgress(message = "Preparing kmeans clustering",{
        sliderInput("multi_kmeans_number", "k-means number", min = 1,
                    max=20, step = 1,
                    value = 2)
      })
    }
  })
  
  output$multi_DEG_total2 <- renderText({
    if(is.null(multi_deg_count1())){
      return(NULL)
    }else{ 
      if(length(input$selectFC2) == 2 && length(input$selectFC2_2) == 2){
        print(paste0("The number of genes after the filtration (fdr < ", input$fdr6,", 
                     basemean > ",input$basemean6,", |log2(", 
                     input$selectFC2[1],"/", input$selectFC2[2],")| > ", log2(input$fc6),", |log2(",
                     input$selectFC2_2[1],"/", input$selectFC2_2[2],")| > ", log2(input$fc6),"): ", 
                     length(rownames(multi_deg_count1()))))
      }else if(length(input$selectFC2) == 2 && length(input$selectFC2_2) != 2){
        print(paste0("The number of genes after the filtration (fdr < ", input$fdr6,", 
                     basemean > ",input$basemean6,", |log2(", 
                     input$selectFC2[1],"/", input$selectFC2[2],")| > ", log2(input$fc6),"): ", 
                     length(rownames(multi_deg_count1()))))
      }else if(length(input$selectFC2) != 2 && length(input$selectFC2_2) != 2){
        print(paste0("The number of genes after the filtration (fdr < ", input$fdr6,", 
                     basemean > ",input$basemean6,"): ", 
                     length(rownames(multi_deg_count1()))))
      }
    }
  })
  
  multi_deg_count1 <- reactive({
    if(is.null(multi_deg_result_raw())){
      return(NULL)
    }
    multi_filtered_pair_result(input$selectFC2, input$selectFC2_2)
  })
  
  multi_deg_count <- reactive({
    data <- multi_norm_count_filter_base()
    sig_res_LRT <- multi_deg_count1()
    if(length(sig_res_LRT$gene) == 0 || is.null(input$topP2)){
      return(NULL)
    }else{
      sig_res_LRT <- sig_res_LRT %>%
        arrange(padj) %>%
        head(n=input$topP2)
      sig_res_LRT <- as.data.frame(sig_res_LRT)
      rownames(sig_res_LRT) <- sig_res_LRT$gene
      sig_res_LRT <- sig_res_LRT[,-1]

      collist <- gsub("\\_.+$", "", colnames(data))
      data2 <- merge(data, sig_res_LRT,by=0)
      rownames(data2) <- data2$Row.names
      data2<- data2[,-1]
      data2 <- data2[,1:length(collist)]
      return(data2)
    }
  })
  
  multi_data_z <- reactive({
    if(!is.null(input$topP2)){
      data <- multi_deg_count()
      if(is.null(data)){
        return(NULL)
      }else{
        data.z <- genefilter::genescale(data, axis = 1, method = "Z")
        data.z <- na.omit(data.z)
        if(input$Species6 != "not selected") {
          if(gene_type6() != "SYMBOL"){
            data.z <- ensembl2symbol(gene_type=gene_type6(),data = data.z, 
                                     Species=input$Species6,Ortholog=ortholog6(),Isoform=isoform6(),org = org6())
            rownames(data.z) <- paste0(data.z$SYMBOL,"\n- ",rownames(data.z))
            data.z <- data.z[, - which(colnames(data.z) == "SYMBOL")]
          }
        }
        return(data.z)
      }
    }
  })
  multi_kmeans_run <- eventReactive(input$kmeans_start_multi, {
    data <- multi_deg_count()
    data.z <- multi_data_z()
    if(is.null(data) || is.null(data.z) || is.null(input$multi_kmeans_number) || input$multi_kmeans_number < 1){
      multi_kmeans_running(FALSE)
      set_multi_kmeans_status(NULL)
      return(NULL)
    }else{
      tryCatch({
        withProgress(message = "k-means clustering", value = 0,{
          cluster_n <- min(input$multi_kmeans_number, nrow(data.z))
          if (cluster_n < 1) {
            multi_kmeans_running(FALSE)
            set_multi_kmeans_status(NULL)
            return(NULL)
          }
          set_multi_kmeans_status("k-means clustering is running.\nComputing consensus clusters...")
          incProgress(0.15, detail = "Computing consensus clusters")
          set.seed(123)
          cl = consensus_kmeans(data.z, cluster_n, 100)
          names(cl) <- rownames(data.z)
          if(length(unique(cl)) != cluster_n && cluster_n > 1){
            set.seed(123)
            cl = consensus_kmeans(data.z, cluster_n-1, 100)
            names(cl) <- rownames(data.z)
            if(length(unique(cl)) != cluster_n-1 && cluster_n > 2){
              set.seed(123)
              cl = consensus_kmeans(data.z, cluster_n-2, 100)
              names(cl) <- rownames(data.z)
            }
          }
          set_multi_kmeans_status("k-means clustering is running.\nPreparing cluster order...")
          incProgress(0.55, detail = "Preparing cluster order")
          order <- build_multi_kmeans_order(data.z, cl)
          incProgress(0.30, detail = "Finalizing k-means inputs")
        })
        return(list(
          data = data,
          data.z = data.z,
          cl = cl,
          order = order
        ))
      }, error = function(e) {
        multi_kmeans_running(FALSE)
        set_multi_kmeans_status(NULL)
        stop(e)
      })
    }
  }, ignoreInit = TRUE)
  pre_multi_kmeans <- reactive({
    run <- multi_kmeans_run()
    if (is.null(run)) {
      return(NULL)
    }
    run$cl
  })
  pre_multi_kmeans_order <- reactive({
    run <- multi_kmeans_run()
    if (is.null(run)) {
      return(character(0))
    }
    run$order
  })
  
  output$kmeans_order_multi <- renderUI({
    if(!is.null(multi_deg_count())){
      order <- as.character(pre_multi_kmeans_order())
      withProgress(message = "Draw heatmap",{
        selectInput("kmeans_order_multi","Order of clusters on heatmap", order,
                    selected = order, multiple = TRUE)
      })
    }
  })
  
  multi_kmeans <- reactive({
    run <- multi_kmeans_run()
    if(is.null(run) || updateCounter_kmeans_multi$i == 0){
      return(NULL)
    }else{
      cluster_order <- as.character(input$kmeans_order_multi)
      if (!length(cluster_order)) {
        cluster_order <- as.character(pre_multi_kmeans_order())
      }
      if(length(cluster_order) == length(unique(pre_multi_kmeans()))){
        set.seed(123)
        ht <- Heatmap(run$data.z, name = "z-score",
                      column_order = colnames(run$data.z),
                      clustering_method_columns = 'ward.D2',
                      cluster_row_slices = F, split = factor(pre_multi_kmeans(),levels = cluster_order),
                      show_row_names = F,column_names_side = "top",use_raster = TRUE)
      }else validate("Select all clusters from 'Order of clusters on heatmap'")
      return(ht)
    }
  })
  multi_kmeans_GOI <- reactive({
    ht <- multi_kmeans()
    run <- multi_kmeans_run()
    if(is.null(run) || is.null(ht) || is.null(pre_multi_kmeans()) || updateCounter_kmeans_multi$i == 0){
      return(NULL)
    }else{
      data.z <- run$data.z
      if(!is.null(input$multi_kmeans_count_table_rows_selected)){
        clusters <- multi_kmeans_cluster()
        if (is.null(clusters)) {
          return(ht)
        }
        data <- clusters[input$multi_kmeans_count_table_rows_selected, , drop = FALSE]
        if (is.null(data) || nrow(data) == 0) {
          return(ht)
        }
        lab <- rownames(data)
        if(input$Species6 != "not selected"){
          if(gene_type6() != "SYMBOL"){
            lab <- paste0(data$SYMBOL,"\n- ", rownames(data))
          }}
        indexes <- which(rownames(data.z) %in% lab)
        labels <- rownames(data.z)[indexes]
        set.seed(123)
        ht <- ht + rowAnnotation(
          link = anno_mark(at = indexes, labels = labels,which="row",link_width = unit(1, "cm"),
                           labels_gp = gpar(fontface = "italic")),
          width = unit(1, "cm") + max_text_width(labels))
      }
      return(ht)
    }
  })
  
  multi_kmeans_cluster <- reactive({
    ht <- multi_kmeans()
    run <- multi_kmeans_run()
    if (is.null(run)) {
      return(NULL)
    }
    data.z <- run$data.z
    data <- run$data
    if(is.null(ht) || is.null(data.z) || is.null(data)){
      return(NULL)
    }else{
      r.dend <- suppressWarnings(row_dend(ht))
      rcl.list <- suppressWarnings(row_order(ht))
      lapply(rcl.list, function(x) length(x))
      Cluster <- NULL
      expected_clusters <- length(unique(pre_multi_kmeans()))
      if(expected_clusters > 0){
        if(length(lapply(rcl.list, function(x) length(x))) != expected_clusters){
          return(NULL)
        }else{
          out <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
          cluster_order <- as.character(input$kmeans_order_multi)
          if (!length(cluster_order)) {
            cluster_order <- as.character(pre_multi_kmeans_order())
          }
          for (i in cluster_order){
            genes <- safe_cluster_genes(data.z, suppressWarnings(row_order(ht)[[i]]))
            if (!length(genes)) next
            clu <- data.frame(GeneID = genes, Cluster = paste("cluster", i, sep=""), stringsAsFactors = FALSE)
            out <- rbind(out, clu)
          }
          if (!nrow(out)) return(NULL)
          colnames(out) <- c("GeneID", "Cluster")
          out <- as.data.frame(out)
          rownames(out) <- out$GeneID
          if(input$Species6 != "not selected"){
            if(gene_type6() != "SYMBOL"){
              rownames(out) <- gsub(".+\\ ", "", out$GeneID)
              out$GeneID <- rownames(out)
              out <- ensembl2symbol(gene_type=gene_type6(),data = out, 
                                    Species=input$Species6,Ortholog=ortholog6(),Isoform=isoform6(),org = org6())
            }}
          clusterCount <- merge(out, data, by=0)
          rownames(clusterCount) <- clusterCount$GeneID
          clusterCount <- clusterCount[,-1:-2]
          return(clusterCount)
        }
      }else return(NULL)
    }
  })
  
  output$multi_kmeans_heatmap <- renderPlot({
    if(!is.null(multi_deg_count())){
      if(is.null(input$multi_kmeans_count_table_rows_selected)) ht <- multi_kmeans() else ht <- multi_kmeans_GOI()
      if(is.null(ht)){
        if (isTRUE(multi_kmeans_running())) {
          multi_kmeans_running(FALSE)
          set_multi_kmeans_status(NULL)
        }
        return(NULL)
      }else{
        if (isTRUE(multi_kmeans_running())) {
          set_multi_kmeans_status("k-means clustering is running.\nRendering heatmap...")
        }
        withProgress(message = "Draw heatmap", value = 0,{
          set.seed(123)
          incProgress(0.2, detail = "Rendering heatmap")
          draw(ht)
          incProgress(0.8)
        })
        if (isTRUE(multi_kmeans_running())) {
          multi_kmeans_running(FALSE)
          set_multi_kmeans_status(NULL)
        }
      }
    }
  })
  
  multi_kmeans_box <- reactive({
    res <- multi_kmeans_cluster()
    
    if (input$Species6 != "not selected") {
      if (gene_type6() != "SYMBOL") {
        res <- res[, -which(colnames(res) == "SYMBOL"), drop = FALSE]
      }
    }
    
    ma <- as.data.frame(multi_deg_count(), stringsAsFactors = FALSE)
    meta <- multi_metadata()
    
    if (is.null(ma) || is.null(res)) {
      return(NULL)
    } else {
      withProgress(message = "Boxplot", {
        
        if (input$multi_data_file_type == "Row1") {
          collist <- gsub("\\_.+$", "", colnames(ma))
          meta <- data.frame(
            condition = factor(collist),
            row.names = colnames(ma),
            stringsAsFactors = FALSE
          )
        } else {
          meta <- as.data.frame(meta, stringsAsFactors = FALSE)
          
          if (ncol(meta) == 1) {
            meta <- data.frame(
              condition = factor(meta[[1]]),
              row.names = colnames(ma),
              stringsAsFactors = FALSE
            )
          } else {
            meta <- data.frame(
              condition = factor(meta[[1]]),
              type = factor(meta[[2]]),
              row.names = colnames(ma),
              stringsAsFactors = FALSE
            )
            meta$type <- factor(meta$type, levels = unique(meta$type), ordered = TRUE)
          }
        }
        
        meta$condition <- factor(meta$condition, levels = unique(meta$condition), ordered = TRUE)
        
        res$Cluster <- gsub("cluster", "", res$Cluster)
        
        gene_ids <- rownames(res)
        if (is.null(gene_ids) || all(gene_ids %in% as.character(seq_len(nrow(res))))) {
          if ("genes" %in% colnames(res)) {
            gene_ids <- res$genes
          } else if ("SYMBOL" %in% colnames(res)) {
            gene_ids <- res$SYMBOL
          } else {
            validate(need(FALSE, "kmeans結果に遺伝子名列がありません。rownames(res) を確認してください。"))
          }
        }
        
        res <- data.frame(
          genes = as.character(gene_ids),
          cluster = as.character(res$Cluster),
          stringsAsFactors = FALSE
        )
        
        table <- rownames_to_column(as.data.frame(ma), "genes") %>%
          tidyr::pivot_longer(
            cols = -genes,
            names_to = "sample",
            values_to = "expression"
          ) %>%
          dplyr::inner_join(
            dplyr::distinct(res[, c("genes", "cluster")]),
            by = "genes"
          ) %>%
          dplyr::left_join(
            tibble::rownames_to_column(as.data.frame(meta), "sample"),
            by = "sample"
          ) %>%
          as.data.frame()
        
        validate(
          need(nrow(table) > 0, "boxplot用テーブルが空です。genes または sample の対応を確認してください。")
        )
        
        table$cluster <- gsub("[^0-9.-]", "", as.character(table$cluster))
        table$cluster[table$cluster == ""] <- NA
        table$cluster <- as.integer(table$cluster)
        
        table <- na.omit(table)
        
        validate(
          need(nrow(table) > 0, "cluster変換後に有効なデータが残っていません。Cluster列を確認してください。")
        )
        
        if (input$multi_data_file_type == "Row1" || ncol(meta) == 1) {
          p <- try(
            DEGreport::degPlotCluster(
              table,
              time = "condition",
              process = TRUE,
              min_genes = 0
            ) +
              scale_color_brewer(palette = "Set1", direction = -1) +
              theme_bw(base_size = 15) +
              theme(
                legend.position = "none",
                axis.text.x = element_text(angle = 90, hjust = 1)
              ),
            silent = TRUE
          )
        } else {
          p <- try(
            DEGreport::degPlotCluster(
              table,
              time = "condition",
              color = "type",
              process = TRUE,
              min_genes = 0
            ) +
              scale_color_brewer(palette = "Set1", direction = -1) +
              theme_bw(base_size = 15) +
              theme(
                legend.position = "top",
                axis.text.x = element_text(angle = 90, hjust = 1)
              ),
            silent = TRUE
          )
        }
        
        if (inherits(p, "try-error")) {
          validate(need(FALSE, as.character(p)))
        }
        
        return(p)
      })
    }
  })
  output$multi_kmeans_boxplot<- renderPlot({
    if(is.null(multi_kmeans_box())|| length(input$selectFC2) != 2){
      return(NULL)
    }else{
      print(multi_kmeans_box())
    }
  })
  
  output$download_multi_kmeans_boxplot = downloadHandler(
    filename = function(){
      paste0(download_multi_overview_dir(), paste(input$multi_kmeans_number,"kmeans_boxplot.pdf",sep = "_"))
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(is.null(multi_kmeans_box())|| length(input$selectFC2) != 2){
          return(NULL)
        }else{
          clusters <- multi_kmeans_cluster()
          clusterNumber <- length(unique(clusters$Cluster))
          if(input$multi_pdf_height == 0){
            pdf_height <- pdf_h(clusterNumber)+2
          }else pdf_height <- input$multi_pdf_height
          if(input$multi_pdf_width == 0){
            pdf_width <- pdf_w(clusterNumber)+2
          }else pdf_width <- input$multi_pdf_width
          pdf(file, height = pdf_height, width = pdf_width)
          print(multi_kmeans_box()+
                  theme(axis.text.x= element_text(size = 8),
                        axis.text.y= element_text(size = 8),
                        title = element_text(size = 8),text = element_text(size = 8)))
          dev.off()
        }
        incProgress(1)
      })
    }
  )
  
  output$multi_kmeans_count_table <- DT::renderDataTable({
    clusters <- multi_kmeans_cluster()
    clusters
  })
  
  output$download_multi_kmeans_cluster = downloadHandler(
    filename = function() {
      paste0(download_multi_overview_dir(), "kmeans_count_table.txt")
    },
    content = function(file){write.table(multi_kmeans_cluster(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  output$download_multi_kmeans_heatmap = downloadHandler(
    filename = function() {
      paste0(download_multi_overview_dir(), paste(input$multi_kmeans_number,"kmeans_heatmap.pdf",sep = "_"))
    },
    content = function(file){
      withProgress(message = "Preparing download",{
        if(input$multi_pdf_height == 0){
          pdf_height <- 10
        }else pdf_height <- input$multi_pdf_height
        if(input$multi_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$multi_pdf_width
        if(is.null(input$multi_kmeans_count_table_rows_selected)) ht <- multi_kmeans() else ht <- multi_kmeans_GOI()
        pdf(file, height = pdf_height, width = pdf_width)
        set.seed(123)
        draw(ht)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$multi_select_file2 <- renderUI({
    clusters <- multi_kmeans_cluster()
    if(is.null(clusters)){
      return(NULL)
    }else{
      selectInput("multi_selectfile2", "cluster_list", choices = c(unique(clusters$Cluster)), multiple = FALSE)
    }
  })
  
  multi_kmeans_pattern_extract <- reactive({
    count <- multi_d_row_count_matrix()
    clusters <- multi_kmeans_cluster()
    if(is.null(clusters)){
      return(NULL)
    }else{
      if(is.null(input$multi_selectfile2) || !nzchar(input$multi_selectfile2)){
        return(NULL)
      }else{
        cluster_name <- input$multi_selectfile2
        clusterCount <- dplyr::filter(clusters, Cluster == cluster_name)
        clusterCount <- clusterCount[,-1]
        return(clusterCount)
      }
    }
  })
  
  output$multi_pattern2_count <- DT::renderDT({
    if(is.null(input$multi_selectfile2) || !nzchar(input$multi_selectfile2)){
      return(NULL)
    }else{
      clusters <- multi_kmeans_pattern_extract()
      clusters
    }
  })
  
  
  multi_kmeans_GOIbox <- reactive({
    if(!is.null(input$multi_kmeans_count_table_rows_selected)){
      data <- multi_kmeans_cluster()[input$multi_kmeans_count_table_rows_selected,]
      data <- data[, - which(colnames(data) == "Cluster")]
      if(input$Species6 != "not selected"){
        if(gene_type6() != "SYMBOL"){
          rownames(data) <- paste(data$SYMBOL,rownames(data), sep = "\n- ")
          data <- data[, - which(colnames(data) == "SYMBOL")]
        }
      }
      return(data)
    }
  })
  
  output$multi_kmeans_GOIboxplot <- renderPlot({
    if(!is.null(input$multi_kmeans_count_table_rows_selected)){
      GOIboxplot(data = multi_kmeans_GOIbox())
    }
  })
  
  output$download_deg_kmeans_boxplot = downloadHandler(
    filename = function(){
      paste(paste(download_multi_overview_dir(), input$multi_selectfile2, sep = "_"), 
            "_boxplot.pdf", sep = "_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- multi_kmeans_GOIbox()
        rowlist <- rownames(data)
        if(input$multi_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)
        }else pdf_height <- input$multi_pdf_height
        if(input$multi_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)
        }else pdf_width <- input$multi_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(GOIboxplot(data = data))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  observeEvent(input$multi_kmeans_count_table_rows_selected, ({
    updateCollapse(session,id =  "multi_collapse_panel2", open="multi_deg_kmeans_boxplot_panel")
  }))
  
  
  output$download_deg_kmeans_pattern_count = downloadHandler(
    filename = function() {
      paste(paste(download_multi_overview_dir(),input$multi_selectfile2, sep = "_"), 
            "kmeans_count_table.txt", sep = "_")
    },
    content = function(file){write.table(multi_kmeans_pattern_extract(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  
  #Multi DEG enrichment------------
  output$multi_Spe1 <- renderText({
    if(input$Species6 == "not selected") print("Please select 'Species'")
  })
  multi_gsea_requested <- reactive({
    isTRUE(multi_analysis_requested()) && identical(input$multi_main_tab, "multi_gsea_tab")
  })
  multi_Hallmark_set <- reactive({
    if(!isTRUE(multi_gsea_requested()) || is.null(input$Gene_set6) || input$Gene_set6 == "" || input$Species6 == "not selected"){
      return(NULL)
    }
    return(GeneList_for_enrichment(Species = input$Species6, Ortholog=input$Ortholog6,
                                   Biomart_archive=input$Biomart_archive6, Gene_set = input$Gene_set6, org = org6()))
  })
  
  observe({
    choices <- multi_select()
    if(is.null(choices)){
      choices <- character(0)
    }
    selected <- intersect(isolate(input$selectEnrich_pair), choices)
    updateSelectizeInput(session, "selectEnrich_pair", choices = choices,
                         selected = selected, options = list(maxItems = 2),
                         server = TRUE)
  })

  multi_enrich_pairFC_raw <- reactive({
    res <- multi_deg_result_raw()
    if(is.null(res) || length(input$selectEnrich_pair) != 2){
      return(NULL)
    }
    pair_tbl <- multi_pair_fc_table_raw(as.data.frame(res, stringsAsFactors = FALSE), input$selectEnrich_pair)
    if(is.null(pair_tbl)){
      return(NULL)
    }
    data <- data.frame(
      Row.names = pair_tbl$Row.names,
      log2FoldChange = -1 * pair_tbl$log2FoldChange,
      stringsAsFactors = FALSE
    )
    if("padj" %in% colnames(res)){
      data$padj <- res[pair_tbl$Row.names, "padj"]
    }
    data
  })
  
  multi_enrich_pairFC <- reactive({
    data <- multi_enrich_pairFC_raw()
    if(is.null(data)){
      return(NULL)
    }
    if(input$Species6 != "not selected"){
      gene_IDs <- multi_gene_id_map()
      if(!is.null(gene_IDs)){
        data <- merge(data, gene_IDs, by = "Row.names")
        if("SYMBOL" %in% colnames(data)){
          data$Unique_ID <- paste(data$SYMBOL, data$Row.names, sep = "\n- ")
        }
      }
    }
    data
  })

  multi_gsea_input <- reactive({
    if(!isTRUE(multi_gsea_requested()) || is.null(input$Gene_set6) || input$Gene_set6 == "" || input$Species6 == "not selected"){
      return(NULL)
    }
    data <- multi_enrich_pairFC()
    if(is.null(data)){
      return(NULL)
    }
    data <- na.omit(data)
    data <- data %>% dplyr::filter(log2FoldChange != 0)
    if(!nrow(data) || !"ENTREZID" %in% colnames(data)){
      return(NULL)
    }
    geneList <- data$log2FoldChange
    names(geneList) <- as.character(data$ENTREZID)
    geneList <- sort(geneList, decreasing = TRUE)
    list(data = data, geneList = geneList)
  })

  multi_gsea_cache_key <- reactive({
    gsea_input <- multi_gsea_input()
    if(is.null(gsea_input)){
      return(NULL)
    }
    digest::digest(list(
      Species = input$Species6,
      Ortholog = input$Ortholog6,
      Biomart_archive = input$Biomart_archive6,
      Gene_set = input$Gene_set6,
      selectEnrich_pair = as.character(input$selectEnrich_pair),
      gsea_data = gsea_input$data[, intersect(c("Row.names", "log2FoldChange", "ENTREZID", "SYMBOL", "Unique_ID", "padj"), colnames(gsea_input$data)), drop = FALSE]
    ))
  })

  multi_gsea_cache <- reactiveValues(
    key = NULL,
    raw = NULL,
    readable = NULL,
    table = NULL
  )
  
  multi_enrichment_1_gsea_raw <- reactive({
    gsea_input <- multi_gsea_input()
    cache_key <- multi_gsea_cache_key()
    if(is.null(gsea_input) || is.null(cache_key)){
      return(NULL)
    }
    if(identical(isolate(multi_gsea_cache$key), cache_key) && !is.null(isolate(multi_gsea_cache$raw))){
      return(isolate(multi_gsea_cache$raw))
    }
    geneList <- gsea_input$geneList
    withProgress(message = "GSEA",{
      if(input$Species6 != "Xenopus laevis" && input$Ortholog6 != "Arabidopsis thaliana" && input$Species6 != "Arabidopsis thaliana"){
        H_t2g <- multi_Hallmark_set()
        H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
        set.seed(123)
        em3 <- try(clusterProfiler::GSEA(geneList, TERM2GENE = H_t2g2,pvalueCutoff = 0.05,exponent = 1, eps = 0, pAdjustMethod = "BH",
                        minGSSize = 50, maxGSSize = 500,by = "fgsea",verbose = F))
      }else{
        set.seed(123)
        if(input$Gene_set6 == "KEGG"){
          em3 <- try(clusterProfiler::gseKEGG(geneList, organism = org_code(input$Species6, Ortholog= input$Ortholog6),pvalueCutoff = 0.05,exponent = 1, eps = 0, pAdjustMethod = "BH",
                             minGSSize = 50, maxGSSize = 500,by = "fgsea",verbose = F,keyType = "ncbi-geneid"))
        }
        if(input$Gene_set6 == "GO biological process"){
          em3 <- try(clusterProfiler::gseGO(geneList, OrgDb = org6(),ont = "BP",pvalueCutoff = 0.05,exponent = 1, eps = 0, pAdjustMethod = "BH",
                           minGSSize = 50, maxGSSize = 500,by = "fgsea",verbose = F))
        }
        if(input$Gene_set6 == "GO cellular component"){
          em3 <- try(clusterProfiler::gseGO(geneList, OrgDb = org6(),ont = "CC",pvalueCutoff = 0.05,exponent = 1, eps = 0, pAdjustMethod = "BH",
                           minGSSize = 50, maxGSSize = 500,by = "fgsea",verbose = F))
        }
        if(input$Gene_set6 == "GO molecular function"){
          em3 <- try(clusterProfiler::gseGO(geneList, OrgDb = org6(),ont = "MF",pvalueCutoff = 0.05,exponent = 1, eps = 0, pAdjustMethod = "BH",
                           minGSSize = 50, maxGSSize = 500,by = "fgsea",verbose = F))
        }
      }
      if(inherits(em3, "try-error")) validate(em3)
      multi_gsea_cache$key <- cache_key
      multi_gsea_cache$raw <- em3
      multi_gsea_cache$readable <- NULL
      multi_gsea_cache$table <- NULL
      return(em3)
      incProgress(1)
    })
  })

  multi_enrichment_1_gsea <- reactive({
    em3 <- multi_enrichment_1_gsea_raw()
    cache_key <- multi_gsea_cache_key()
    if(is.null(em3) || is.null(cache_key)){
      return(NULL)
    }
    if(identical(isolate(multi_gsea_cache$key), cache_key) && !is.null(isolate(multi_gsea_cache$readable))){
      return(isolate(multi_gsea_cache$readable))
    }
    if (length(as.data.frame(em3)$ID) == 0) {
      em4 <- NA
    } else{
      em4 <- clusterProfiler::setReadable(em3, org6(), 'ENTREZID')
    }
    multi_gsea_cache$readable <- em4
    return(em4)
  })
  
  
  # Multi DEG enrichment plot ------------------------------------------------------------------------------
  multi_enrich1_H <- reactive({
    if(!is.null(input$Gene_set6) && input$Species6 != "not selected"){
      em3 <- multi_enrichment_1_gsea()
      if (length(as.data.frame(em3)$ID) == 0) {
        p4 <- NULL
      } else{
        if (length(as.data.frame(em3)$ID) >= 5){
          p4 <- enrichplot::gseaplot2(em3, 1:5, pvalue_table = F,base_size = 14)
        }else{
          p4 <- enrichplot::gseaplot2(em3, 1:length(em3$ID), pvalue_table = F,base_size = 14)
        }
      }
      return(p4)
    }else return(NULL)
  })
  
  output$multi_enrichment1 <- renderPlot({
    dotplot_for_output(data = multi_enrich_pairFC(),
                       plot_genelist = multi_enrich1_H(), Gene_set = input$Gene_set6, 
                       Species = input$Species6)
  })
  
  observe({
    if(input$Species6 != "Xenopus laevis" && input$Ortholog6 != "Arabidopsis thaliana" && input$Species6 != "Arabidopsis thaliana"){
      choices <- gene_set_list
    }else{
      choices <- c("KEGG", "GO biological process", "GO cellular component","GO molecular function")
    }
    selected <- isolate(input$Gene_set6)
    if(is.null(selected) || !selected %in% choices){
      selected <- if(length(choices)) choices[[1]] else ""
    }
    updateSelectInput(session, "Gene_set6", choices = choices, selected = selected)
  })
  
  multi_GSEA_table <- reactive({
    gsea_result <- multi_enrichment_1_gsea()
    cache_key <- multi_gsea_cache_key()
    if(is.null(gsea_result) || length(input$selectEnrich_pair) != 2){
      return(NULL)
    }else{
      if(identical(isolate(multi_gsea_cache$key), cache_key) && !is.null(isolate(multi_gsea_cache$table))){
        return(isolate(multi_gsea_cache$table))
      }
      data <- as.data.frame(gsea_result)
      if(input$Species6 != "Xenopus laevis" && input$Ortholog6 != "Arabidopsis thaliana" && input$Species6 != "Arabidopsis thaliana"){
        H_t2g <- multi_Hallmark_set()
        if(length(as.data.frame(data)$Description) == 0 || is.null(H_t2g)){
          return(NULL)
        }else{
          colnames(data)[1] <- "gs_name"
          H_t2g <- H_t2g %>% distinct(gs_name, .keep_all = T)
          data2 <- left_join(data, H_t2g, by="gs_name")  %>% as.data.frame()
          if(input$Gene_set6 == "DoRothEA regulon (activator)" || input$Gene_set6 == "DoRothEA regulon (repressor)"){
            data3 <- data.frame(Gene_set_name = data2$gs_name, Confidence = data2$confidence,
                                setSize = data2$setSize, enrichmentScore = data2$enrichmentScore, NES = data2$NES, 
                                pvalue = data2$pvalue, p.adjust = data2$p.adjust, qvalue = data2$qvalue, 
                                rank = data2$rank, leading_edge = data2$leading_edge, core_enrichment = data2$core_enrichment)
          }else{
            if(input$Gene_set6 == "Custom gene set"){
              data3 <- data.frame(Gene_set_name = data2$gs_name,
                                  setSize = data2$setSize, enrichmentScore = data2$enrichmentScore, NES = data2$NES, 
                                  pvalue = data2$pvalue, p.adjust = data2$p.adjust, qvalue = data2$qvalue, 
                                  rank = data2$rank, leading_edge = data2$leading_edge, core_enrichment = data2$core_enrichment)
            }else{
              data3 <- data.frame(Gene_set_name = data2$gs_name, Description = data2$gs_description,
                                  setSize = data2$setSize, enrichmentScore = data2$enrichmentScore, NES = data2$NES, 
                                pvalue = data2$pvalue, p.adjust = data2$p.adjust, qvalue = data2$qvalue, 
                                rank = data2$rank, leading_edge = data2$leading_edge, core_enrichment = data2$core_enrichment)
            }
          }
          multi_gsea_cache$table <- data3
          return(data3) 
        }
      }else{
        multi_gsea_cache$table <- data
        return(data)
      }
    }
  })
  
  output$multi_GSEA_result <- DT::renderDataTable({
    if(is.null(multi_enrichment_1_gsea()) || length(input$selectEnrich_pair) != 2){
      return(NULL)
    }else{
      multi_GSEA_table()
    }
  })
  
  output$download_multi_enrichment = downloadHandler(
    filename = function(){
      paste(download_multi_overview_dir(), paste0(input$Gene_set6,"-GSEA.pdf"), sep="_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- multi_enrich1_H()
        if(input$multi_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$multi_pdf_height
        if(input$multi_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <-input$multi_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(p1)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_multi_GSEA_table = downloadHandler(
    filename = function() {
      paste(download_multi_overview_dir(), paste0(input$Gene_set6,"-GSEA.txt"), sep="_")
    },
    content = function(file){write.table(multi_GSEA_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  # Multi ssGSEA-----
  Custom_input_ssGSEA <- reactive({
    tmp <- input$custom_input_ssGSEA$datapath
    data <- read_gene_list(tmp)
    df <- gene_list_convert_for_enrichment(gene_type=gene_type6(),data= data, org=org6(),Species = input$Species6,Ortholog=ortholog6(),Isoform=isoform6())
    return(df)
  })
  output$Custom_input_ssGSEA <- renderUI({
    if(is.null(input$Gene_set_ssGSEA)){
      return(NULL)
    }else{
      if(input$Gene_set_ssGSEA == "Custom gene set"){
        fileInput("custom_input_ssGSEA",
                  "Select a file (txt, csv, xlsx)",
                  accept = c("txt", "csv","xlsx"),
                  multiple = FALSE,
                  width = "80%")
      }else return(NULL)
    }
  })
  output$multi_Spe1_ssGSEA <- renderText({
    if(input$Species6 == "not selected") print("Please select 'Species'")
  })
  multi_ssgsea_custom_input_info <- reactive({
    if(is.null(input$custom_input_ssGSEA$datapath)){
      return(NULL)
    }
    info <- file.info(input$custom_input_ssGSEA$datapath)
    list(
      name = input$custom_input_ssGSEA$name,
      size = info$size[[1]],
      mtime = as.character(info$mtime[[1]])
    )
  })
  multi_ssgsea_gene_set_cache_key <- reactive({
    if(is.null(input$Gene_set_ssGSEA) || input$Gene_set_ssGSEA == "" || input$Species6 == "not selected"){
      return(NULL)
    }
    if(input$Gene_set_ssGSEA == "Custom gene set" && is.null(input$custom_input_ssGSEA$datapath)){
      return(NULL)
    }
    digest::digest(list(
      Species = input$Species6,
      Ortholog = input$Ortholog6,
      Biomart_archive = input$Biomart_archive6,
      Gene_set = input$Gene_set_ssGSEA,
      gene_type = gene_type6(),
      custom_input = multi_ssgsea_custom_input_info()
    ))
  })
  multi_ssgsea_gene_set_cache <- reactiveValues(
    key = NULL,
    value = NULL
  )
  multi_Hallmark_set_ssGSEA <- reactive({
    cache_key <- multi_ssgsea_gene_set_cache_key()
    if(is.null(cache_key)){
      return(NULL)
    }
    if(identical(isolate(multi_ssgsea_gene_set_cache$key), cache_key) &&
       !is.null(isolate(multi_ssgsea_gene_set_cache$value))){
      return(isolate(multi_ssgsea_gene_set_cache$value))
    }
    if(input$Gene_set_ssGSEA == "Custom gene set" && is.null(Custom_input_ssGSEA())) validate("")
    gene_set <- GeneList_for_enrichment(Species = input$Species6, Ortholog=input$Ortholog6,
                                        Biomart_archive=input$Biomart_archive6, gene_type = gene_type6(),
                                        Gene_set = input$Gene_set_ssGSEA, Custom_gene_list = Custom_input_ssGSEA(),org = org6())
    my.symbols <- as.character(gene_set$entrez_gene)
    if(gene_type6() != "SYMBOL"){
      if(sum(is.element(no_orgDb, input$Species6)) == 1){
        gene_IDs <- ortholog6()
      }else if(gene_type6() == "isoform"){
        gene_IDs <- isoform6()
      }else{
        if(str_detect(my.symbols[1], "^AT.G")) key = "TAIR" else key = "ENSEMBL"
        gene_IDs<-AnnotationDbi::select(org6(),keys = my.symbols,
                                        keytype = "ENTREZID",
                                        columns = c("ENTREZID",key,sgd_symbol_column(org6())))
      }
      colnames(gene_IDs) <- c("entrez_gene","GeneID","SYMBOL")
    }else{
      if(sum(is.element(no_orgDb, input$Species6)) == 1){
        gene_IDs <- ortholog6()[,-1]
      }else if(gene_type6() == "isoform"){
        gene_IDs <- isoform6()[,-1]
      }else{
        gene_IDs <- AnnotationDbi::select(org6(), keys = my.symbols,
                                          keytype = "ENTREZID",
                                          columns = c("ENTREZID",sgd_symbol_column(org6())))
      }
      colnames(gene_IDs) <- c("entrez_gene","GeneID")
    }
    gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
    gene_set <- merge(gene_set,gene_IDs, by="entrez_gene")
    gene_set <- gene_set %>% distinct(GeneID, .keep_all = T)
    multi_ssgsea_gene_set_cache$key <- cache_key
    multi_ssgsea_gene_set_cache$value <- gene_set
    return(gene_set)
  })
  observe({
    if(input$Species6 != "Xenopus laevis" && input$Ortholog6 != "Arabidopsis thaliana" && input$Species6 != "Arabidopsis thaliana"){
      choices <- c("", gene_set_list)
    }else{
      choices <- c("", "KEGG", "GO biological process", "GO cellular component","GO molecular function")
    }
    selected <- isolate(input$Gene_set_ssGSEA)
    if(is.null(selected) || !selected %in% choices){
      selected <- ""
    }
    updateSelectInput(session, "Gene_set_ssGSEA", choices = choices, selected = selected)
  })
  
  multi_enrichment_1_ssGSEA <- reactive({
    if(!is.null(input$Gene_set_ssGSEA) && input$Gene_set_ssGSEA != "" && input$Species6 != "not selected"){
      count <- multi_deg_norm_count()
      withProgress(message = "ssGSEA",{
        if(input$Gene_set_ssGSEA == "Custom gene set" && is.null(Custom_input_ssGSEA())) validate("")
        result <- ssGSEA(norm_count = count,gene_set = multi_Hallmark_set_ssGSEA(), org = org6(),
                       gene_type=gene_type6(),Species=input$Species6,Ortholog=input$Ortholog6)
        return(result)
        incProgress(1)
      })
    }
  })
  output$multi_ssGSEA_score <- DT::renderDataTable({
    if(is.null(multi_enrichment_1_ssGSEA())){
      return(NULL)
    }else{
      multi_enrichment_1_ssGSEA()
    }
  })
  output$download_multi_ssGSEA_score = downloadHandler(
    filename = function() {
      paste(download_multi_overview_dir(), paste0(input$Gene_set_ssGSEA,"-ssGSEA_score.txt"), sep="_")
    },
    content = function(file){write.table(multi_enrichment_1_ssGSEA(), file, row.names = T,col.names = NA, sep = "\t", quote = F)}
  )
  # Multi ssGSEA limma-----
  multi_ssGSEA_limma <- reactive({
      if(is.null(multi_enrichment_1_ssGSEA())){
        return(NULL)
      }else{
        count <- multi_enrichment_1_ssGSEA()
          collist <- gsub("\\_.+$", "", colnames(count))
          collist <- gsub(" ", ".", collist)
          collist <- gsub("\\+", ".", collist)
          if(is.element(TRUE, duplicated(collist)) == TRUE){
            if (input$multi_data_file_type == "Row1"){
              meta <- data.frame(condition = factor(collist))
              pair <- collist
            }else {
              meta <- multi_metadata()
              meta <- data.frame(condition=factor(meta[,1]), type=factor(meta[,2]))
              pair <- paste0(meta$condition,meta$type)
            }
            design <- model.matrix(~0+pair)
            colnames(design) <- factor(gsub("^pair","",colnames(design)),levels = unique(pair))
            cont <- c()
            for(i in 1:choose(n=length(unique(pair)),k=2)){
              contrast = paste0(as.character(unique(pair)[combn(x=length(unique(pair)),m=2)[1,i]]),"-",as.character(unique(pair)[combn(x=length(unique(pair)),m=2)[2,i]]))
              cont <- c(cont,contrast)
            }
            cont.matrix <- makeContrasts(contrasts=cont, levels=design)
            eset = new("ExpressionSet", exprs=as.matrix(count))
            fit1 <- lmFit(eset,design)
            fit2 <- contrasts.fit(fit1, cont.matrix)
            fit3 <- eBayes(fit2)
            result <- topTable(fit3,coef=1:length(cont), number = 1e12)
            lab <- cont
            if(length(cont) != 1) label <- c(lab,"AveExpr","F","p_value","padj") else label <- c(lab,"AveExpr","F","p_value","padj","B")
            colnames(result) <- label
            return(result)
          }
          
        }
  })
  output$multi_ssGSEA_limma <- DT::renderDataTable({
    if(is.null(multi_ssGSEA_limma())){
      return(NULL)
    }else{
      multi_ssGSEA_limma()
    }
  })
  output$download_multi_ssGSEA_limma = downloadHandler(
    filename = function() {
      paste(download_multi_overview_dir(), paste0(input$Gene_set_ssGSEA,"-ssGSEA_differential_analysis_all_result.txt"), sep="_")
    },
    content = function(file){write.table(multi_ssGSEA_limma(), file, row.names = T,col.names = NA, sep = "\t", quote = F)}
  )
  multi_ssGSEA_limma_dp <- reactive({
    result <- multi_ssGSEA_limma()
    if(!is.null(result)){
      if(dim(result)[1] != 0){
      res2 <- result %>% dplyr::filter(padj < input$ssGSEA_fdr)
      return(res2)
      }
    }
  })
  output$multi_ssGSEA_limma_dp <- DT::renderDataTable({
    if(is.null(multi_ssGSEA_limma_dp())){
      return(NULL)
    }else{
      multi_ssGSEA_limma_dp()
    }
  })
  output$download_multi_ssGSEA_limma_dp = downloadHandler(
    filename = function() {
      paste(download_multi_overview_dir(), paste0(input$Gene_set_ssGSEA,"-ssGSEA_differential_pathways.txt"), sep="_")
    },
    content = function(file){write.table(multi_ssGSEA_limma_dp(), file, row.names = T,col.names = NA, sep = "\t", quote = F)}
  )
  #ssGSEA GOI------------------------------------------------------
  observeEvent(multi_ssGSEA_limma_dp(), {
    choices <- rownames(multi_ssGSEA_limma_dp())
    if(is.null(choices)){
      choices <- character(0)
    }
    selected <- intersect(isolate(input$GOI_multi_ssGSEA), choices)
    updateSelectizeInput(session, "GOI_multi_ssGSEA", choices = choices,
                         selected = selected,
                         options = goi_selectize_options, server = TRUE)
  }, ignoreNULL = FALSE)
  observeEvent(input$GOIreset_multi_ssGSEA, {
    choices <- rownames(multi_ssGSEA_limma_dp())
    if(is.null(choices)){
      choices <- character(0)
    }
    updateSelectizeInput(session, "GOI_multi_ssGSEA", choices = choices,
                         selected = character(0),
                         options = goi_selectize_options, server = TRUE)
  })
  observeEvent(list(input$GOI_type_multi_ssGSEA, input$GOI_type_multi_ssGSEA_custom), {
    if(identical(input$GOI_type_multi_ssGSEA, "custom") &&
       !identical(input$GOI_type_multi_ssGSEA_custom, "Manual")){
      choices <- rownames(multi_ssGSEA_limma_dp())
      if(is.null(choices)){
        choices <- character(0)
      }
      freezeReactiveValue(input, "GOI_multi_ssGSEA")
      updateSelectizeInput(session, "GOI_multi_ssGSEA", choices = choices,
                           selected = character(0),
                           options = goi_selectize_options, server = TRUE)
    }
  }, ignoreNULL = FALSE)
  output$GOI_type_multi_ssGSEA_all <- renderUI({
    if(input$GOI_type_multi_ssGSEA == "ALL"){
      radioButtons('GOI_type_multi_ssGSEA_all','Label of heatmap:',
                   c('ALL'=TRUE,
                     'Select pathways from the below list'=FALSE
                   ),selected = FALSE)
    }
  })
  output$GOI_type_multi_ssGSEA_custom <- renderUI({
    if(input$GOI_type_multi_ssGSEA == "custom"){
      radioButtons('GOI_type_multi_ssGSEA_custom','Mode:',
                   c('top 10 pathways'="top10",
                     'top 20 pathways'="top20",
                     'top 40 pathways'="top40",
                     'Manual'="Manual"
                   ),selected = "top10")
    }
  })
  GOI_multi_ssGSEA_INPUT <- reactive({
    if(is.null(multi_ssGSEA_limma_dp())){
      return(NULL)
    }else{
      if(dim(multi_ssGSEA_limma_dp())[1] == 0) validate("There are no differential pathways. Please change FDR cut-off value.")
      if(input$GOI_type_multi_ssGSEA == "ALL") return(rownames(multi_ssGSEA_limma_dp()))
      if(input$GOI_type_multi_ssGSEA == "custom"){
        if(input$GOI_type_multi_ssGSEA_custom == "top10"){
          data <- multi_ssGSEA_limma_dp() %>% dplyr::arrange(padj) %>% head(n=10)
          return(rownames(data))
        }
        if(input$GOI_type_multi_ssGSEA_custom == "top20"){
          data <- multi_ssGSEA_limma_dp() %>% dplyr::arrange(padj) %>% head(n=20)
          return(rownames(data))
        }
        if(input$GOI_type_multi_ssGSEA_custom == "top40"){
          data <- multi_ssGSEA_limma_dp() %>% dplyr::arrange(padj) %>% head(n=40)
          return(rownames(data))
        }
        if(length(input$GOI_multi_ssGSEA) == 0) validate("")
        if(input$GOI_type_multi_ssGSEA_custom == "Manual"){
          return(input$GOI_multi_ssGSEA)
        }else {
          data <- multi_ssGSEA_limma_dp()
          df2 <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          for(i in rownames(data)){
            corr<-suppressWarnings(cor.test(x=as.numeric(data[input$GOI_multi_ssGSEA,]),y=as.numeric(data[i,]),method="spearman"))
            df <- data.frame(y_axis = i, x_axis = input$GOI_multi_ssGSEA, statistics=corr$statistic,corr_score = corr$estimate,
                             pvalue = corr$p.value)
            df2 <- rbind(df2,df)
          }
          padj <- p.adjust(df2$pvalue,method="BH")
          df2$padj <- padj
          colnames(df2) <- c("prey","bait","statistics","corr_score","pvalue","padj")
          df2 <- na.omit(df2)
          df2 <- df2%>% dplyr::filter(pvalue < 0.05)
          return(df2$prey)
        }
      } 
    }
  })
  multi_ssGSEA_GOIcount <- reactive({
    count <- multi_enrichment_1_ssGSEA()
    Row.names <- GOI_multi_ssGSEA_INPUT()
    label_data <- as.data.frame(Row.names, row.names = Row.names)
    data <- merge(label_data,count, by=0)
    rownames(data) <- data[,1]
    data <- data[, -1:-2]
    return(data)
  })
  
  multi_ssGSEA_GOIheat <- reactive({
    data <- multi_ssGSEA_GOIcount()
    data <- genefilter::genescale(data,axis=1, method="Z")
    if(is.null(data)){
      ht <- NULL
    }else{
      rownames(data) <- gsub("_"," ", rownames(data))
      for(i in 1:length(rownames(data))){
        rownames(data)[i] <- paste(strwrap(rownames(data)[i], width = 30),collapse = "\n")
      }
      ht <- GOIheatmap(data, type = input$GOI_type_multi_ssGSEA, GOI = input$GOI_multi_ssGSEA,all=input$GOI_type_multi_ssGSEA_all)
    }
    return(ht)
  })
  
  output$multi_ssGSEA_GOIheatmap <- renderPlot({
    if(is.null(multi_ssGSEA_limma_dp())){
      return(NULL)
    }else{
      if(!is.null(GOI_multi_ssGSEA_INPUT())){
        withProgress(message = "heatmap",{
          suppressWarnings(print(multi_ssGSEA_GOIheat()))
          incProgress(1)
        })
      }
    }
  })
  output$statistics_multi_ssGSEA <- renderUI({
    if(!is.null(multi_ssGSEA_limma_dp()) && !is.null(GOI_multi_ssGSEA_INPUT())){
      data <- multi_ssGSEA_GOIcount()
      if(!is.null(data)){
        collist <- gsub("\\_.+$", "", colnames(data))
        collist <- unique(collist)
        if(length(collist) == 2){
          selectInput("statistics_multi_ssGSEA","statistics",choices = c("not_selected","Welch's t-test","Wilcoxon test"),selected="not_selected",multiple = F)
        }else{
          selectInput("statistics_multi_ssGSEA","statistics",choices = c("not_selected","TukeyHSD","Dunnet's test","Dunn's test"),selected="not_selected", multiple = F)
        }
      }}
  })
  statistical_analysis_goi_multi_ssGSEA <- reactive({
    data <- multi_ssGSEA_GOIcount()
    if(is.null(data) || is.null(input$statistics_multi_ssGSEA)){
      p <- NULL
    }else{
      rownames(data) <- gsub("_"," ", rownames(data))
      for(i in 1:length(rownames(data))){
        rownames(data)[i] <- paste(strwrap(rownames(data)[i], width = 10),collapse = "\n")
      }
      p <- GOIboxplot(data = data,statistical_test =input$statistics_multi_ssGSEA,plottype=input$PlotType_multi_ssGSEA,ssGSEA=TRUE)
    }
    return(p)
  })
  
  multi_ssGSEA_GOIbox <- reactive({
    print("multi_gsca_GOIbox start")
    if(!is.null(statistical_analysis_goi_multi_ssGSEA())){
      if(input$statistics_multi_ssGSEA == "not_selected"){
        return(statistical_analysis_goi_multi_ssGSEA())
      }else return(statistical_analysis_goi_multi_ssGSEA()[["plot"]])
    }
  })
  
  
  output$multi_ssGSEA_GOIboxplot <- renderPlot({
    if(is.null(multi_ssGSEA_limma_dp())){
      return(NULL)
    }else{
      if(!is.null(GOI_multi_ssGSEA_INPUT())){
        if(length(rownames(multi_ssGSEA_GOIcount())) >200){
          validate("Unable to display more than 200 genes. Please adjust the threshold to narrow down the number of genes to less than 200, or utilize the 'Custom' mode.")
        }
        withProgress(message = "Boxplot",{
          suppressWarnings(print(multi_ssGSEA_GOIbox()))
          incProgress(1)
        })
      }
    }
  })
  
  multi_ssGSEA_GOIbox_statistic <- reactive({
    if(!is.null(statistical_analysis_goi_multi_ssGSEA())){
      if(input$statistics_multi_ssGSEA != "not_selected"){
        data <- as.data.frame(statistical_analysis_goi_multi_ssGSEA()[["statistical_test"]])
        colnames(data)[1] <- "gene"
        if(input$statistics_multi_ssGSEA == "TukeyHSD" || input$statistics_multi_ssGSEA == "Dunnet's test"){
          data <- data[, - which(colnames(data) == "term")]
        }else{
          data <- data[, - which(colnames(data) == ".y.")]
          data <- data[, - which(colnames(data) == "n1")]
          data <- data[, - which(colnames(data) == "n2")]
        }
        data <- data[, - which(colnames(data) == "y.position")]
        data <- data[, - which(colnames(data) == "groups")]
        data <- data[, - which(colnames(data) == "xmin")]
        data <- data[, - which(colnames(data) == "xmax")]
        return(data)
      }else return(NULL)}
  })
  
  output$statistical_table_multi_ssGSEA <- DT::renderDataTable({
    if(is.null(multi_ssGSEA_limma_dp())){
      return(NULL)
    }else{
      if(!is.null(GOI_multi_ssGSEA_INPUT())){
        if(length(rownames(multi_ssGSEA_GOIcount())) >200){
          validate("Cannot display more than 200 pathways.")
        }
        multi_ssGSEA_GOIbox_statistic()
      }}
  })
  
  output$download_multi_ssGSEA_statisics = downloadHandler(
    filename = function(){
      paste0(download_multi_overview_dir(), "GOIboxplot_",input$statistics_multi_ssGSEA,".txt")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        write.table(multi_ssGSEA_GOIbox_statistic(),file, row.names = F, col.names=TRUE, sep = "\t", quote = F)
        incProgress(1)
      })
    }
  )
  
  output$download_multi_ssGSEA_GOIbox = downloadHandler(
    filename = function(){
      paste0(download_multi_overview_dir(), "ssGSEA_boxplot.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- multi_ssGSEA_GOIcount()
        rowlist <- rownames(data)
        if(input$multi_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)
        }else pdf_height <- input$multi_pdf_height
        if(input$multi_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)
        }else pdf_width <- input$multi_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(multi_ssGSEA_GOIbox())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_multi_ssGSEA_GOIheat = downloadHandler(
    filename = function(){
      paste0(download_multi_overview_dir(), "GOIheatmap.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- multi_ssGSEA_GOIcount()
        rowlist <- rownames(data)
        if(input$multi_pdf_height == 0){
          pdf_height <- 10
        }else pdf_height <- input$multi_pdf_height
        if(input$multi_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$multi_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(multi_ssGSEA_GOIheat())
        dev.off()
        incProgress(1)
      })
    }
  )
  #ssGSEA contribute---------
  observeEvent(multi_ssGSEA_limma_dp(), {
    choices <- rownames(multi_ssGSEA_limma_dp())
    if(is.null(choices)){
      choices <- character(0)
    }
    selected <- isolate(input$selectssGSEA_contribute_pathway)
    if(is.null(selected) || !selected %in% choices){
      selected <- if(length(choices)) choices[[1]] else ""
    }
    updateSelectizeInput(session, "selectssGSEA_contribute_pathway",
                         choices = choices, selected = selected, server = TRUE)
  }, ignoreNULL = FALSE)
  multi_ssGSEA_df_selected_id <- reactive({
    data <- multi_ssGSEA_limma_dp()
    if(!is.null(data) && !is.null(input$selectssGSEA_contribute_pathway)){
      score <- multi_enrichment_1_ssGSEA() %>% 
        as.data.frame() %>%
        dplyr::filter(row.names(.) == input$selectssGSEA_contribute_pathway)
      return(score)
    }
  })
  multi_norm_count_for_ssGSEA_selected_id <- reactive({
    count <- multi_deg_norm_count()
    if(length(grep("SYMBOL", colnames(count))) != 0) count <- count[, - which(colnames(count) == "SYMBOL")]
    dp <- multi_ssGSEA_df_selected_id()
    if(!is.null(data) && !is.null(dp)){
      geneset <- multi_Hallmark_set_ssGSEA() %>% dplyr::filter(gs_name == input$selectssGSEA_contribute_pathway)
      rownames(geneset) <- geneset$GeneID
      count2 <- merge(geneset, count, by = 0)
      rownames(count2) <- count2$GeneID
      if(input$Gene_set_ssGSEA == "DoRothEA regulon (activator)" || 
        input$Gene_set_ssGSEA == "DoRothEA regulon (repressor)") {
        num <- -1
        }else if(input$Gene_set_ssGSEA =="Custom gene set"){
          num <- dim(geneset)[2]-5
          }else num <- 0
      count2 <- count2[,-1:-(6+num)]
      print(head(geneset))
      return(count2)
    }
  })
  
  multi_ssGSEA_contribute_cor <- reactive({
    count <- multi_norm_count_for_ssGSEA_selected_id()
    if(length(grep("SYMBOL", colnames(count))) != 0) norm_count <- count[, - which(colnames(count) == "SYMBOL")] else norm_count <- count
    score <- multi_ssGSEA_df_selected_id()
    print("score")
    print(head(as.numeric(score)))
    print(head(norm_count))
    print(dim(score))
    print(dim(norm_count))
    if(!is.null(norm_count)){
      df2 <- data.frame(matrix(rep(NA, 4), nrow=1))[numeric(0), ]
      for(i in rownames(norm_count)){
        corr<-suppressWarnings(try(cor.test(x=as.numeric(score[1,]),y=as.numeric(norm_count[i,]),method="spearman")))
        if(!inherits(corr, "try-error")){
          if(length(grep("SYMBOL", colnames(count))) != 0){
            symbol <- count %>% dplyr::filter(row.names(.) == i)
            df <- data.frame(Gene = i, statistics=corr$statistic,corr_score = corr$estimate,
                             pvalue = corr$p.value, UniqueID = paste0(i,"-",symbol$SYMBOL))
          }else{
        df <- data.frame(Gene = i, statistics=corr$statistic,corr_score = corr$estimate,
                         pvalue = corr$p.value)
          }
        df2 <- rbind(df2,df)
        }
      }
      print(head(df2))
      label <- c("Gene","statistics","corr_score","pvalue")
      if(length(grep("SYMBOL", colnames(count))) != 0) label <- c(label,"Unique_ID")
      colnames(df2) <- label
      df2 <- na.omit(df2)
      df2 <- df2%>% dplyr::arrange(-corr_score, pvalue) %>%
        dplyr::mutate(rank = row_number())
      rownames(df2) <- df2$rank
      
      df2$ssGSEAscore_vs_Expression <- "NS"
      df2$ssGSEAscore_vs_Expression[df2$pvalue < 0.05 & df2$corr_score > 0] <- "positive_correlation"
      df2$ssGSEAscore_vs_Expression[df2$pvalue < 0.05 & df2$corr_score < 0] <- "negative_correlation"
      df2 <- df2 %>% dplyr::select(Gene,ssGSEAscore_vs_Expression,everything())
      df2 <- df2%>% dplyr::arrange(-corr_score, pvalue)
      return(df2)
    }
  })
  multi_ssGSEA_contribute_cor_count <- reactive({
    if(!is.null(multi_ssGSEA_contribute_cor())){
      norm_count <- multi_norm_count_for_ssGSEA_selected_id()
      print("norm_count")
      print(head(norm_count))
      data <- multi_ssGSEA_contribute_cor() %>% 
        dplyr::filter(ssGSEAscore_vs_Expression == "positive_correlation")
      rownames(data) <- data$Gene
      print("data")
      print(head(data))
      count <- merge(data, norm_count, by=0)
      rownames(count) <- count$Gene
      print("count")
      print(head(count))
      if(length(grep("SYMBOL", colnames(count))) != 0) rownames(count) <- count$Unique_ID
      count <- count[,-1:-7]
      print("count")
      print(head(count))
      if(length(grep("SYMBOL", colnames(count))) != 0){
        count <- count[,-1]
      count <- count[, - which(colnames(count) == "SYMBOL")]
      }
      print(head(count))
      return(count)
    }
  })
  multi_ssGSEA_contribute_corplot <- reactive({
    if(!is.null(multi_ssGSEA_contribute_cor_count())){
      data <- multi_ssGSEA_contribute_cor_count()
      rownames(data) <- gsub("-","- ", rownames(data))
      for(i in 1:length(rownames(data))){
        rownames(data)[i] <- paste(strwrap(rownames(data)[i], width = 10),collapse = "\n")
      }
      print(head(data))
      p <- GOIboxplot(data = data)
      return(p)
    }
  })
  output$multi_ssGSEA_contribute <- renderPlot({
    if(!is.null(multi_ssGSEA_contribute_corplot()))
      multi_ssGSEA_contribute_corplot()
  })
  output$multi_ssGSEA_contribute_table <-DT::renderDataTable({
    if(!is.null(multi_ssGSEA_contribute_corplot()))
      multi_ssGSEA_contribute_cor()
  })
  output$download_multi_ssGSEA_contribute = downloadHandler(
    filename = function(){
      paste0(download_multi_overview_dir(), "ssGSEAcontributed_",input$selectssGSEA_contribute_pathway,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- multi_ssGSEA_contribute_cor_count()
        rowlist <- rownames(data)
        if(input$multi_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)
        }else pdf_height <- input$multi_pdf_height
        if(input$multi_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)
        }else pdf_width <- input$multi_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(multi_ssGSEA_contribute_corplot())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_multi_ssGSEA_contribute_table = downloadHandler(
    filename = function(){
      paste0(download_multi_overview_dir(), "ssGSEAcontributed_",input$selectssGSEA_contribute_pathway,".txt")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        write.table(multi_ssGSEA_contribute_cor(),file, row.names = F, col.names=TRUE, sep = "\t", quote = F)
        incProgress(1)
      })
    }
  )
  #ssGSEA dorothea---------
  multi_ssGSEA_limma_dp_entrez <- reactive({
    data <- multi_ssGSEA_limma_dp()
    if(input$Gene_set_ssGSEA != "DoRothEA regulon (activator)" && 
       input$Gene_set_ssGSEA != "DoRothEA regulon (repressor)") validate("You have to select the 'DoRothEA regulon' from the 'Gene sets' list.")
    if(!is.null(data)){
      score <- multi_enrichment_1_ssGSEA()
      Row.names <- rownames(data)
      label_data <- as.data.frame(Row.names, row.names = Row.names)
      score_dp <- merge(label_data,score, by=0)
      rownames(score_dp) <- score_dp[,1]
      score_dp <- score_dp[, -1:-2]
      gene_IDs <- geneid_convert_ssGSEA(norm_count=score_dp,org = org6(),
                                        gene_type=gene_type6(),Species=input$Species6,Ortholog=ortholog6(),Isoform=isoform6())
      if(gene_type6() != "SYMBOL"){
        score_dp$SYMBOL <- rownames(score_dp)
        data2 <- merge(gene_IDs,score_dp,by="SYMBOL")
      }else{
        data2 <- merge(gene_IDs,score_dp,by=0)
      }
      rownames(data2) <- data2$GeneID
      data2 <- data2[,-1:-3]
      if(length(grep("SYMBOL", colnames(data2))) != 0) data2 <- data2[, - which(colnames(data2) == "SYMBOL")]
      return(data2)
    }
  })
  multi_norm_count_for_ssGSEA_dorothea_entrez <- reactive({
    count <- multi_deg_norm_count()
    dp <- multi_ssGSEA_limma_dp_entrez()
    if(!is.null(data) && !is.null(dp)){
      Row.names <- rownames(dp)
      label_data <- as.data.frame(Row.names, row.names = Row.names)
      count_dp <- merge(label_data,count, by=0)
      rownames(count_dp) <- count_dp[,1]
      count_dp <- count_dp[, -1:-2]
      return(count_dp)
    }
  })
  
  multi_ssGSEA_TF_cor <- reactive({
    count <- multi_norm_count_for_ssGSEA_dorothea_entrez()
    if(length(grep("SYMBOL", colnames(count))) != 0) norm_count <- count[, - which(colnames(count) == "SYMBOL")] else norm_count <- count
    score <- multi_ssGSEA_limma_dp_entrez()
    if(!is.null(norm_count)){
      df2 <- data.frame(matrix(rep(NA, 4), nrow=1))[numeric(0), ]
      for(i in rownames(score)){
        corr<-suppressWarnings(try(cor.test(x=as.numeric(score[i,]),y=as.numeric(norm_count[i,]),method="spearman")))
        if(!inherits(corr, "try-error")){
          if(length(grep("SYMBOL", colnames(count))) != 0){
            symbol <- count %>% dplyr::filter(row.names(.) == i)
            df <- data.frame(TF = i, statistics=corr$statistic,corr_score = corr$estimate,
                             pvalue = corr$p.value, UniqueID = paste0(i,"-",symbol$SYMBOL))
          }else{
          df <- data.frame(TF = i, statistics=corr$statistic,corr_score = corr$estimate,
                           pvalue = corr$p.value)
          }
          df2 <- rbind(df2,df)
        }
      }
      label <- c("TF","statistics","corr_score","pvalue")
      if(length(grep("SYMBOL", colnames(count))) != 0) label <- c(label,"Unique_ID")
      colnames(df2) <- label
      df2 <- na.omit(df2)
      df2 <- df2%>% dplyr::arrange(-corr_score, pvalue) %>%
        dplyr::mutate(rank = row_number())
      rownames(df2) <- df2$rank
      
      
      df2$ssGSEAscore_vs_Expression <- "NS"
      df2$ssGSEAscore_vs_Expression[df2$pvalue < 0.05 & df2$corr_score > 0] <- "positive_correlation"
      df2$ssGSEAscore_vs_Expression[df2$pvalue < 0.05 & df2$corr_score < 0] <- "negative_correlation"
      df2 <- df2 %>% dplyr::select(TF,ssGSEAscore_vs_Expression,everything())
      df2 <- df2%>% dplyr::arrange(-corr_score, pvalue)
      return(df2)
    }
  })
  multi_ssGSEA_TF_cor_count <- reactive({
    if(!is.null(multi_ssGSEA_TF_cor())){
      norm_count <- multi_deg_norm_count()
      data <- multi_ssGSEA_TF_cor() %>% 
        dplyr::filter(ssGSEAscore_vs_Expression == "positive_correlation")
      rownames(data) <- data$TF
      count <- merge(data, norm_count, by=0)
      rownames(count) <- count$TF
      if(length(grep("SYMBOL", colnames(count))) != 0) rownames(count) <- count$Unique_ID
      count <- count[,-1:-7]
      if(length(grep("SYMBOL", colnames(count))) != 0){
        count <- count[,-1]
        count <- count[, - which(colnames(count) == "SYMBOL")]
      }
      return(count)
    }
  })
  multi_ssGSEA_TF_corplot <- reactive({
    if(!is.null(multi_ssGSEA_TF_cor())){
      data <- multi_ssGSEA_TF_cor_count()
      rownames(data) <- gsub("-","- ", rownames(data))
      for(i in 1:length(rownames(data))){
        rownames(data)[i] <- paste(strwrap(rownames(data)[i], width = 10),collapse = "\n")
      }
      p <- GOIboxplot(data = data)
      return(p)
    }
  })
  output$multi_ssGSEA_dorothea <- renderPlot({
    if(!is.null(multi_ssGSEA_TF_corplot()))
      multi_ssGSEA_TF_corplot()
  })
  output$multi_ssGSEA_dorothea_table <-DT::renderDataTable({
    if(!is.null(multi_ssGSEA_TF_corplot()))
      multi_ssGSEA_TF_cor()
  })
  output$download_multi_ssGSEA_dorothea = downloadHandler(
    filename = function(){
      paste0(download_multi_overview_dir(), "ssGSEAscore_vsExpression_level.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- multi_ssGSEA_TF_cor_count()
        rowlist <- rownames(data)
        if(input$multi_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)
        }else pdf_height <- input$multi_pdf_height
        if(input$multi_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)
        }else pdf_width <- input$multi_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(multi_ssGSEA_TF_corplot())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_multi_ssGSEA_dorothea_table = downloadHandler(
    filename = function(){
      paste0(download_multi_overview_dir(), "ssGSEAscore_vsExpression_level.txt")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        write.table(multi_ssGSEA_TF_cor(),file, row.names = F, col.names=TRUE, sep = "\t", quote = F)
        incProgress(1)
      })
    }
  )
  #multi DEG enrichment 2--------
  multi_Hallmark_set2 <- reactive({
    return(GeneList_for_enrichment(Species = input$Species6, Ortholog=input$Ortholog6,Biomart_archive=input$Biomart_archive6, Gene_set = input$Gene_set7, org = org6()))
  })
  
  multi_Hallmark_set3 <- reactive({
    return(GeneList_for_enrichment(Species = input$Species6, Ortholog=input$Ortholog6,Biomart_archive=input$Biomart_archive6, Gene_set = input$Gene_set8, org = org6()))
  })
  
  
  output$Gene_set7 <- renderUI({
    if(input$Species6 != "Xenopus laevis" && input$Ortholog6 != "Arabidopsis thaliana" && input$Species6 != "Arabidopsis thaliana"){
      selectInput('Gene_set7', 'Gene Set', gene_set_list)
    }else selectInput('Gene_set7', 'Gene Set', c("KEGG", "GO biological process", 
                                                 "GO cellular component","GO molecular function"))
  })
  
  output$Gene_set8 <- renderUI({
    if(input$Species6 != "Xenopus laevis" && input$Ortholog6 != "Arabidopsis thaliana" && input$Species6 != "Arabidopsis thaliana"){
      selectInput('Gene_set8', 'Gene Set', gene_set_list)
    }else selectInput('Gene_set8', 'Gene Set', c("KEGG", "GO biological process", 
                                                 "GO cellular component","GO molecular function"))
  })
  
  output$multi_Spe <- renderText({
    if(input$Species6 == "not selected") print("Please select 'Species'")
  })
  output$multi_Spe2 <- renderText({
    if(input$Species6 == "not selected") print("Please select 'Species'")
  })
  
  output$multi_whichGroup1_1 <- renderUI({
    clusters <- multi_pattern2()$df
    if(is.null(clusters)){
      return(NULL)
    }else{
      clusters$cluster <- paste0("Group",clusters$cluster)
      selectInput("multi_whichGroup1_1", "cluster_list", choices = c(sort(unique(clusters$cluster))), multiple = TRUE)
    }
  })
  
  output$multi_whichGroup2_1 <- renderUI({
    clusters <- multi_kmeans_cluster()
    if(is.null(clusters)){
      return(NULL)
    }else{
      selectInput("multi_whichGroup2_1", "cluster_list", choices = c(sort(unique(clusters$Cluster))),multiple = TRUE)
    }
  })
  
  output$multi_whichGroup1_2 <- renderUI({
    clusters <- multi_pattern2()$df
    if(is.null(clusters)){
      return(NULL)
    }else{
      clusters$cluster <- paste0("Group",clusters$cluster)
      selectInput("multi_whichGroup1_2", "cluster_list", choices = c("not selected", unique(clusters$cluster)), selected = "not selected" ,multiple = F)
    }
  })
  
  output$multi_whichGroup2_2 <- renderUI({
    clusters <- multi_kmeans_cluster()
    if(is.null(clusters)){
      return(NULL)
    }else{
      selectInput("multi_whichGroup2_2", "cluster_list", choices = c("not selected", unique(clusters$Cluster)), selected = "not selected", multiple = F)
    }
  })
  
  multi_enrich_input1 <- reactive({
    if(is.null(input$multi_whichGroup1_1)){
      return(NULL)
    }else{
      clusters <- multi_pattern2()$df
      clusters$cluster <- paste0("Group",clusters$cluster)
      cluster_name <- input$multi_whichGroup1_1
      clusterCount <- data.frame(GeneID = NA, Group=NA)
      for(name in cluster_name){
        clusterCount2 <- dplyr::filter(clusters, cluster == name)
        clusterCount2 <- data.frame(GeneID = rownames(clusterCount2), Group = clusterCount2$cluster)
        clusterCount <- rbind(clusterCount, clusterCount2) 
      }
      clusterCount <- na.omit(clusterCount)
      return(clusterCount)
    }
  })
  
  multi_enrich_input2 <- reactive({
    if(is.null(input$multi_whichGroup2_1)){
      return(NULL)
    }else{
      clusters <- multi_kmeans_cluster()
      if(input$Species6 != "not selected"){
        if(gene_type6() != "SYMBOL"){
          clusters <- clusters[, - which(colnames(clusters) == "SYMBOL")]
        }}
      cluster_name <- input$multi_whichGroup2_1
      clusterCount <- data.frame(GeneID = NA, Group=NA)
      for(name in cluster_name){
        clusterCount2 <- dplyr::filter(clusters, Cluster == name)
        clusterCount2 <- data.frame(GeneID = rownames(clusterCount2), Group = clusterCount2$Cluster)
        clusterCount <- rbind(clusterCount, clusterCount2) 
      }
      clusterCount <- na.omit(clusterCount)
      return(clusterCount)
    }
  })

  multi_enrich_div_source <- reactive({
    if(input$Species6 != "Xenopus laevis" && input$Ortholog6 != "Arabidopsis thaliana" && input$Species6 != "Arabidopsis thaliana"){
      enrich_viewer_forMulti1(
        gene_type = gene_type6(),
        df = multi_enrich_input1(),
        Species = input$Species6,
        Ortholog = ortholog6(),
        Isoform = isoform6(),
        org = org6()
      )
    }else{
      enrich_viewer_forMulti1(
        gene_type = gene_type6(),
        df = multi_enrich_input1(),
        Species = input$Species6,
        Ortholog = ortholog6(),
        Isoform = isoform6(),
        org = org6()
      )
    }
  })

  multi_enrich_k_source <- reactive({
    enrich_viewer_forMulti1(
      gene_type = gene_type6(),
      df = multi_enrich_input2(),
      Species = input$Species6,
      Ortholog = ortholog6(),
      Isoform = isoform6(),
      org = org6()
    )
  })

  multi_enrich_div_cnet_source <- reactive({
    enrich_viewer_forMulti1(
      gene_type = gene_type6(),
      df = multi_enrich_input3(),
      Species = input$Species6,
      Ortholog = ortholog6(),
      Isoform = isoform6(),
      org = org6()
    )
  })

  multi_enrich_k_cnet_source <- reactive({
    enrich_viewer_forMulti1(
      gene_type = gene_type6(),
      df = multi_enrich_input4(),
      Species = input$Species6,
      Ortholog = ortholog6(),
      Isoform = isoform6(),
      org = org6()
    )
  })
  
  multi_enrich_viewer2 <- reactive({
    if(input$Species6 != "Xenopus laevis" && input$Ortholog6 != "Arabidopsis thaliana" && input$Species6 != "Arabidopsis thaliana"){
      return(enrich_viewer_forMulti2(gene_type=gene_type6(),df = multi_enrich_input1(), Species = input$Species6,Ortholog = ortholog6(),Isoform=isoform6(), org = org6(),
                                     org_code = org_code6(),H_t2g = multi_Hallmark_set2(),Gene_set = input$Gene_set7))
    }else return(enrich_viewer_forMulti2_xenopus(gene_type=gene_type6(),df = multi_enrich_input1(), Species = input$Species6, org = org6(),
                                                 org_code = org_code6(),Ortholog = ortholog6(),Isoform=isoform6(),Gene_set = input$Gene_set7))
  })
  
  multi_enrich_viewer12 <- reactive({
    if(input$Species6 != "Xenopus laevis" && input$Ortholog6 != "Arabidopsis thaliana" && input$Species6 != "Arabidopsis thaliana"){
      return(enrich_viewer_forMulti2(gene_type=gene_type6(),df = multi_enrich_input2(), Species = input$Species6,Ortholog = ortholog6(),Isoform=isoform6(), org = org6(),
                                     org_code = org_code6(),H_t2g = multi_Hallmark_set3(),Gene_set = input$Gene_set8))
    }else return(enrich_viewer_forMulti2_xenopus(gene_type=gene_type6(),df = multi_enrich_input2(), Species = input$Species6, org = org6(),
                                                 org_code = org_code6(),Ortholog = ortholog6(),Isoform=isoform6(),Gene_set = input$Gene_set8))
  })
  multi_enrich_h <- reactive({
    source_data <- multi_enrich_div_source()
    if(input$Species6 != "Xenopus laevis" && input$Ortholog6 != "Arabidopsis thaliana" && input$Species6 != "Arabidopsis thaliana"){
      return(enrich_gene_list(data = source_data,
                              Gene_set = input$Gene_set7, org = org6(), H_t2g = multi_Hallmark_set2()))
    }else return(enrich_gene_list_xenopus(data = source_data,
                                          Gene_set = input$Gene_set7, org = org6(), org_code = org_code6()))
  })
  multi_enrich_H <- reactive({
    return(enrich_genelist(data = multi_enrich_div_source(),
                           enrich_gene_list = multi_enrich_h(),group_order = input$multi_whichGroup1_1))
  })
  multi_enrich_h2 <- reactive({
    source_data <- multi_enrich_k_source()
    if(input$Species6 != "Xenopus laevis" && input$Ortholog6 != "Arabidopsis thaliana" && input$Species6 != "Arabidopsis thaliana"){
      return(enrich_gene_list(data = source_data,
                              Gene_set = input$Gene_set8, org = org6(), H_t2g = multi_Hallmark_set3()))
    }else return(enrich_gene_list_xenopus(data = source_data,
                                          Gene_set = input$Gene_set8, org = org6(),org_code = org_code6()))
  })
  multi_enrich_H2 <- reactive({
    return(enrich_genelist(data = multi_enrich_k_source(),
                           enrich_gene_list = multi_enrich_h2(),group_order = input$multi_whichGroup2_1))
  })
  
  output$multi_enrichment3 <- renderPlot({
    dotplot_for_output(data = multi_enrich_viewer2(),
                       plot_genelist = multi_enrich_H(), Gene_set = input$Gene_set7, 
                       Species = input$Species6)
  })
  
  output$multi_enrichment5 <- renderPlot({
    dotplot_for_output(data = multi_enrich_viewer12(), 
                       plot_genelist = multi_enrich_H2(), Gene_set = input$Gene_set8, 
                       Species = input$Species6)
  })
  
  multi_enrich_input3 <- reactive({
    if(is.null(input$multi_whichGroup1_2) || input$multi_whichGroup1_2 == "not selected"){
      return(NULL)
    }else{
      clusters <- multi_pattern2()$df
      clusters$cluster <- paste0("Group",clusters$cluster)
      cluster_name <- input$multi_whichGroup1_2
      clusterCount <- data.frame(GeneID = NA, Group=NA)
      for(name in cluster_name){
        clusterCount2 <- dplyr::filter(clusters, cluster == name)
        clusterCount2 <- data.frame(GeneID = rownames(clusterCount2), Group = clusterCount2$cluster)
        clusterCount <- rbind(clusterCount, clusterCount2) 
      }
      clusterCount <- na.omit(clusterCount)
      return(clusterCount)
    }
  })
  
  multi_enrich_input4 <- reactive({
    if(is.null(input$multi_whichGroup2_2) || input$multi_whichGroup2_2 == "not selected"){
      return(NULL)
    }else{
      clusters <- multi_kmeans_cluster()
      cluster_name <- input$multi_whichGroup2_2
      clusterCount <- data.frame(GeneID = NA, Group=NA)
      for(name in cluster_name){
        clusterCount2 <- dplyr::filter(clusters, Cluster == name)
        clusterCount2 <- data.frame(GeneID = rownames(clusterCount2), Group = clusterCount2$Cluster)
        clusterCount <- rbind(clusterCount, clusterCount2) 
      }
      clusterCount <- na.omit(clusterCount)
      return(clusterCount)
    }
  })
  
  multi_enrich2 <- reactive({
    cnet_global(data = multi_enrich_div_cnet_source(), 
                group = input$multi_whichGroup1_2, enrich_gene_list = multi_enrich_h())
  })
  
  multi_enrich12 <- reactive({
    cnet_global(data = multi_enrich_k_cnet_source(), 
                group = input$multi_whichGroup2_2, enrich_gene_list = multi_enrich_h2())
  })
  
  output$multi_enrichment4 <- renderPlot({
    cnet_for_output(data = multi_enrich_input3(), plot_data = multi_enrich2(), 
                    Gene_set = input$Gene_set7, Species = input$Species6)
  })
  
  output$multi_enrichment6 <- renderPlot({
    cnet_for_output(data = multi_enrich_input4(), plot_data = multi_enrich12(), 
                    Gene_set = input$Gene_set8, Species = input$Species6)
  })
  
  output$download_multi_cluster_enrichment = downloadHandler(
    filename = function(){
      paste(download_multi_overview_dir(), paste0(input$Gene_set7,"_divisive_enrichment.pdf"), sep="_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- multi_enrich_H()
        if(input$multi_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$multi_pdf_height
        if(input$multi_pdf_width == 0){
          pdf_width <- 8
        }else pdf_width <- input$multi_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(p1)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_multi_cluster_enrichment2 = downloadHandler(
    filename = function(){
      paste(download_multi_overview_dir(), paste0(input$Gene_set8,"_kmeans_enrichment.pdf"), sep="_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- multi_enrich_H2()
        if(input$multi_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$multi_pdf_height
        if(input$multi_pdf_width == 0){
          pdf_width <- 8
        }else pdf_width <- input$multi_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(p1)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_multi_enrichment_cnet = downloadHandler(
    filename = function(){
      paste(download_multi_overview_dir(), paste0(input$Gene_set7,"_divisive_cnet.pdf"), sep="_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p <- multi_enrich2()
        if(input$multi_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$multi_pdf_height
        if(input$multi_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$multi_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(p)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_Multi_report = downloadHandler(
    filename = function() {
      if (input$multi_data_file_type == "Row1"){
        paste(format(Sys.time(), "%Y%m%d"),download_multi_overview_dir(),"MultiDEG.zip", sep ="_")
      }else{
        paste(format(Sys.time(), "%Y%m%d"),download_multi_overview_dir(), "MultiDEG.zip",sep ="_")
      }
    },
    content = function(fname){
      withProgress(message = "Preparing download, please wait",{
        multi_analysis_requested(TRUE)
        fs <- c()
        setwd(tempdir())
        print(fs)
        dir.create("DEG_result/",showWarnings = FALSE)
        dir.create("Clustering/",showWarnings = FALSE)
        DEG <- "DEG_result/DEG_result.txt"
        count <- "DEG_result/normalized_count.txt"
        PCA <- "Clustering/clustering.pdf"
        PCA_table <- "Clustering/pca.txt"
        fs <- c(DEG,count,PCA,PCA_table)
        write.table(multi_deg_result(), DEG, row.names = T, col.names=NA, sep = "\t", quote = F)
        write.table(multi_deg_norm_count(), count, row.names = T, col.names=NA, sep = "\t", quote = F)
        write.table(PCAdata(row_count = multi_d_row_count_matrix(), deg_norm_count = multi_deg_norm_count()), PCA_table, row.names = T, col.names=NA, sep = "\t", quote = F)
        pdf(PCA, height = 3.5, width = 9)
        print(multi_pca_plot())
        dev.off()
        if(length(multi_umap_plot()) != 1){
          umap <- "Clustering/umap.pdf"
          fs <- c(fs,umap)
          pdf(umap, height = 3.5, width = 4.7)
          print(multi_umap_plot())
          dev.off()
        }
        print("finish pca")
        if(length(input$selectFC) == 2){
          if(!is.null(multi_boxplot_reactive())){
            dir.create(paste0("Divisive clustering_",as.character(input$selectFC[1]),"_vs_",as.character(input$selectFC[2]),"/"),showWarnings = FALSE)
            dir.create(paste0("Divisive clustering_",as.character(input$selectFC[1]),"_vs_",as.character(input$selectFC[2]),"/group_list"),showWarnings = FALSE)
            DEG_pattern <- paste0("Divisive clustering_",as.character(input$selectFC[1]),"_vs_",as.character(input$selectFC[2]),"/DEG_pattern.txt")
            summary_box1 <- paste0("Divisive clustering_",as.character(input$selectFC[1]),"_vs_",as.character(input$selectFC[2]),"/Divisive_boxplot.pdf")
            fs <- c(fs, DEG_pattern,summary_box1)
            clusters <- multi_pattern2()$df
            clusters$cluster <- paste0("Group",clusters$cluster)
            write.table(clusters, DEG_pattern, quote = F, row.names = F, sep = "\t")
            clusters <- multi_pattern2()$df
            data <- as.data.frame(multi_deg_norm_count())
            clusters$cluster <- paste0("Group",clusters$cluster)
            clusters <- merge(clusters, data, by=0)
            cluster_names <- unique(clusters$cluster)
            print(head(clusters))
            for(name in cluster_names){
              clusterCount <- dplyr::filter(clusters, cluster == name)
              rownames(clusterCount) <- clusterCount$Row.names
              clusterCount <- clusterCount[,-1:-3]
              DEG_pattern_count <- paste0("Divisive clustering_",as.character(input$selectFC[1]),"_vs_",as.character(input$selectFC[2]),"/group_list/",name,".txt")
              fs <- c(fs,DEG_pattern_count)
              write.table(clusterCount, DEG_pattern_count, quote = F, row.names = T, col.names=NA, sep = "\t")
            }
            clusters <- multi_pattern2()$df
            clusterNumber <- length(unique(clusters$cluster))
            pdf_height <- pdf_h(clusterNumber)+2
            pdf_width <- pdf_w(clusterNumber)+2
            pdf(summary_box1, height = pdf_height, width = pdf_width)
            print(multi_boxplot_reactive()+
                    theme(axis.text.x= element_text(size = 8),
                          axis.text.y= element_text(size = 8),
                          title = element_text(size = 8),text = element_text(size = 8)))
            dev.off()
            if(!is.null(input$multi_pattern1_count_rows_selected)){
              box1 <- paste0("Divisive clustering_",as.character(input$selectFC[1]),"_vs_",as.character(input$selectFC[2]),"/GOI_boxplot.pdf")
              fs <- c(fs,box1)
              data <- multi_GOIbox()
              rowlist <- rownames(data)
              pdf_height <- pdf_h(rowlist)
              pdf_width <- pdf_w(rowlist)
              pdf(box1, height = pdf_height, width = pdf_width)
              print(GOIboxplot(data = data))
              dev.off()
            }
            if(!is.null(input$Gene_set7) && !is.null(input$multi_whichGroup1_1)){
              enrich1 <- paste0("Divisive clustering_",as.character(input$selectFC[1]),"_vs_",as.character(input$selectFC[2]),"/dotplot_",input$Gene_set7,".pdf")
              enrichtxt <- paste0("Divisive clustering_",as.character(input$selectFC[1]),"_vs_",as.character(input$selectFC[2]),"/enrichment_",input$Gene_set7,".txt")
              if(!is.null(input$multi_whichGroup1_1)){
                fs <- c(fs,enrich1,enrichtxt)
                p1 <- multi_enrich_H()
                pdf(enrich1, height = 6, width = 8)
                print(p1)
                dev.off()
              }
              p <- multi_enrich2()
              if(input$multi_whichGroup1_2 != "not selected"){
                cnet1 <- paste0("Divisive clustering_",as.character(input$selectFC[1]),"_vs_",as.character(input$selectFC[2]),"/cnet_",input$Gene_set7,".pdf")
                fs <- c(fs,cnet1)
                pdf(cnet1, height = 6, width = 6)
                print(p)
                dev.off()
                write.table(multi_enrich_div_table(), enrichtxt, row.names = F, sep = "\t", quote = F)
              }
            }
          }
        }
        print("finish divisive clustering")
        if(length(input$selectFC2) == 2){
          if(!is.null(multi_kmeans_box())){
            dir.create(paste0("kmeans clustering_",as.character(input$selectFC2[1]),"_vs_",as.character(input$selectFC2[2]),"/"),showWarnings = FALSE)
            kmeans_pattern <- paste0("kmeans clustering_",as.character(input$selectFC2[1]),"_vs_",as.character(input$selectFC2[2]),"/kmeans_pattern.txt")
            summary_kmeansbox1 <- paste0("kmeans clustering_",as.character(input$selectFC2[1]),"_vs_",as.character(input$selectFC2[2]),"/kmeans_boxplot.pdf")
            kmeans_heat<- paste0("kmeans clustering_",as.character(input$selectFC2[1]),"_vs_",as.character(input$selectFC2[2]),"/kmeans_heatmap.pdf")
            fs <- c(fs, kmeans_pattern,summary_kmeansbox1,kmeans_heat)
            write.table(multi_kmeans_cluster(), kmeans_pattern, row.names = T, col.names=NA, sep = "\t", quote = F)
            clusters <- multi_kmeans_cluster()
            clusterNumber <- length(unique(clusters$Cluster))
            pdf_height <- pdf_h(clusterNumber)+2
            pdf_width <- pdf_w(clusterNumber)+2
            pdf(summary_kmeansbox1, height = pdf_height, width = pdf_width)
            print(multi_kmeans_box()+
                    theme(axis.text.x= element_text(size = 8),
                          axis.text.y= element_text(size = 8),
                          title = element_text(size = 8),text = element_text(size = 8)))
            dev.off()
            if(is.null(input$multi_kmeans_count_table_rows_selected)) ht <- multi_kmeans() else ht <- multi_kmeans_GOI()
            pdf(kmeans_heat, height = 10, width = 7)
            set.seed(123)
            draw(ht)
            dev.off()
            if(!is.null(input$multi_kmeans_count_table_rows_selected)){
              box2 <- paste0("kmeans clustering_",as.character(input$selectFC2[1]),"_vs_",as.character(input$selectFC2[2]),"/GOI_boxplot.pdf")
              fs <- c(fs,box2)
              data <- multi_kmeans_GOIbox()
              rowlist <- rownames(data)
              pdf_height <- pdf_h(rowlist)
              pdf_width <- pdf_w(rowlist)
              pdf(box2, height = pdf_height, width = pdf_width)
              print(GOIboxplot(data = data))
              dev.off()
            }
            dir.create(paste0("kmeans clustering_",as.character(input$selectFC2[1]),"_vs_",as.character(input$selectFC2[2]),"/cluster_list/"),showWarnings = FALSE)
            clusters <- unique(multi_kmeans_cluster()$Cluster)
            for(name in clusters){
              clusterCount <- dplyr::filter(multi_kmeans_cluster(), Cluster == name)
              clusterCount <- clusterCount[,-1]
              kmeans_pattern_count <- paste0("kmeans clustering_",as.character(input$selectFC2[1]),"_vs_",as.character(input$selectFC2[2]),"/cluster_list/",name,".txt")
              fs <- c(fs, kmeans_pattern_count)
              write.table(clusterCount, kmeans_pattern_count, row.names = T, col.names=NA, sep = "\t", quote = F)
            }
            if(!is.null(input$Gene_set8) && !is.null(input$multi_whichGroup2_1)){
              enrich2 <- paste0("kmeans clustering_",as.character(input$selectFC2[1]),"_vs_",as.character(input$selectFC2[2]),"/dotplot_",input$Gene_set8,".pdf")
              enrichtxt2 <- paste0("kmeans clustering_",as.character(input$selectFC2[1]),"_vs_",as.character(input$selectFC2[2]),"/enrichment_",input$Gene_set8,".txt")
              if(!is.null(input$multi_whichGroup2_1)){
                fs <- c(fs,enrich2,enrichtxt2)
                p1 <- multi_enrich_H2()
                pdf(enrich2, height = 6, width = 8)
                print(p1)
                dev.off()
              }
              p <- multi_enrich12()
              if(input$multi_whichGroup2_2 != "not selected"){
                cnet2 <- paste0("kmeans clustering_",as.character(input$selectFC2[1]),"_vs_",as.character(input$selectFC2[2]),"/cnet_",input$Gene_set8,".pdf")
                fs <- c(fs,cnet2)
                pdf(cnet2, height = 6, width = 6)
                print(p)
                dev.off()
                write.table(multi_enrich_k_table(), enrichtxt2, row.names = F, sep = "\t", quote = F)
              }
            }
          }
        }
        print("finish kmeans clustering")
        if(!is.null(input$selectEnrich_pair)){
          if(length(input$selectEnrich_pair) == 2 && input$Species6 != "not selected"){
            if(!is.null(multi_enrich1_H())){
              dir.create("GSEA/",showWarnings = FALSE)
              multiEnrich_table <- paste0("GSEA/GSEA_",input$Gene_set6,".txt")
              multiEnrich <- paste0("GSEA/GSEA_",input$Gene_set6,".pdf")
              fs <- c(fs, multiEnrich,multiEnrich_table)
              write.table(multi_GSEA_table(), multiEnrich_table, row.names = F, sep = "\t", quote = F)
              p1 <- multi_enrich1_H()
              pdf(multiEnrich, height = 5, width = 7)
              print(p1)
              dev.off()
            }}}
        print("GSEA")
        report <- paste0(format(Sys.time(), "%Y%m%d_"),"MultiDEG_report",".docx")
        fs <- c(fs,report)
        rmarkdown::render("multi_report.Rmd", output_format = "word_document", output_file = report,
                          params = list(raw_count = multi_d_row_count_matrix(),
                                        multi_norm_count_matrix = multi_norm_count_matrix(),
                                        input = input,
                                        multi_metadata = multi_metadata(),
                                        multi_umap_plot = multi_umap_plot(),
                                        multi_boxplot_reactive = multi_boxplot_reactive(),
                                        multi_enrich_div_table = multi_enrich_div_table(),
                                        multi_kmeans_box = multi_kmeans_box(),
                                        multi_enrich_k_table = multi_enrich_k_table(),
                                        multi_GSEA_table = multi_GSEA_table(),
                                        multi_pattern1 = multi_pattern1(),
                                        multi_deg_count1 = multi_deg_count1()), 
                          envir = new.env(parent = globalenv()),intermediates_dir = tempdir(),encoding="utf-8"
        )
        zip(zipfile=fname, files=fs)
      })
    },
    contentType = "application/zip"
  )
  
  output$download_multi_enrichment_cnet2 = downloadHandler(
    filename = function(){
      paste(download_multi_overview_dir(), paste0(input$Gene_set8,"_kmeans_cnet.pdf"), sep="_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p <- multi_enrich12()
        if(input$multi_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$multi_pdf_height
        if(input$multi_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$multi_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(p)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  multi_enrich_div_table <- reactive({
    if(input$Species6 != "Xenopus laevis" && input$Ortholog6 != "Arabidopsis thaliana" && input$Species6 != "Arabidopsis thaliana"){
      return(enrich_for_table(data = multi_enrich_viewer2(), H_t2g = multi_Hallmark_set2(), Gene_set = input$Gene_set7))
    }else return(as.data.frame(multi_enrich_viewer2()))
  })
  
  multi_enrich_k_table <- reactive({
    if(input$Species6 != "Xenopus laevis" && input$Ortholog6 != "Arabidopsis thaliana" && input$Species6 != "Arabidopsis thaliana"){
      return(enrich_for_table(data = multi_enrich_viewer12(), H_t2g = multi_Hallmark_set3(), Gene_set = input$Gene_set8))
    }else return(as.data.frame(multi_enrich_viewer12()))
  })
  
  output$multi_enrichment_result <- DT::renderDataTable({
    multi_enrich_div_table()
  })
  
  output$multi_enrichment_result2 <- DT::renderDataTable({
    multi_enrich_k_table()
  })
  
  output$download_multi_enrichment_table = downloadHandler(
    filename = function() {
      paste(download_multi_overview_dir(), paste0(input$Gene_set7,"_divisive_enrichment.txt"), sep="_")
    },
    content = function(file){write.table(multi_enrich_div_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  output$download_multi_enrichment_table2 = downloadHandler(
    filename = function() {
      paste(download_multi_overview_dir(), paste0(input$Gene_set8,"_kmeans_enrichment.txt"), sep="_")
    },
    content = function(file){write.table(multi_enrich_k_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  #Multi PCA ------------------------------------------------------------------------------
  multi_pca_plot <- reactive({
    data <- multi_deg_norm_count()
    if(is.null(data)){
      return(NULL)
    }
    PCAplot(data,legend = input$PCA_legend_multi)
  })

  multi_pca_data <- reactive({
    row_count <- multi_d_row_count_matrix()
    deg_norm_count <- multi_deg_norm_count()
    if(is.null(row_count) || is.null(deg_norm_count)){
      return(NULL)
    }
    PCAdata(row_count = row_count, deg_norm_count = deg_norm_count)
  })
  
  output$download_multi_PCA = downloadHandler(
    filename = function(){
      paste0(download_multi_overview_dir(), "_PCA-MDS-dendrogram.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$multi_pdf_height == 0){
          pdf_height <- 3.5
        }else pdf_height <- input$multi_pdf_height
        if(input$multi_pdf_width == 0){
          pdf_width <- 9
        }else pdf_width <- input$multi_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(multi_pca_plot())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$multi_PCA <- renderPlot({
    p <- multi_pca_plot()
    if(is.null(p)){
      return(NULL)
    }else{
      print(p)
    }
  })
  
  output$multi_PCA_data <- DT::renderDataTable({
    multi_pca_data()
  })
  
  output$download_multi_PCA_table = downloadHandler(
    filename = function() {
      paste0(download_multi_overview_dir(), "_PCA_table.txt")
    },
    content = function(file){write.table(multi_pca_data(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  
  output$multi_umap_n <- renderUI({
    info <- umap_neighbor_info(multi_d_row_count_matrix())
    if(!info$valid){
      return(NULL)
    }
    sliderInput("multi_n_neighbors", "n_neighbors", min = 2,
                max = info$max, step = 1,
                value = info$value)
  })
  
  multi_umap_plot <- reactive({
    data <- multi_d_row_count_matrix()
    if(is.null(input$multi_n_neighbors)){
      return(NULL)
    }else{
      if(is.null(data)){
        return(NULL)
      }else{
        p<- try(umap_plot(data = data, n_neighbors = input$multi_n_neighbors,lab=input$multi_umap_label), silent = TRUE)
        return(p)
      }
    }
  })
  
  output$multi_umap_error <- renderText({
    data <- multi_d_row_count_matrix()
    err <- umap_error_message(data, input$multi_n_neighbors)
    if(!is.null(err)){
      return(err)
    }
    p <- multi_umap_plot()
    if(inherits(p, "try-error")){
      return(conditionMessage(attr(p, "condition")))
    }
    NULL
  })
  
  
  output$multi_umap <- renderPlot({
    p <- multi_umap_plot()
    withProgress(message = "umap",{
      if(inherits(p, "try-error") || is.null(p)) return(NULL)
      print(p)
    })
  })
  
  output$download_multi_umap = downloadHandler(
    filename = function(){
      paste0(download_multi_overview_dir(), paste0(input$multi_n_neighbors,"_neighbors_umap.pdf"))
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$multi_pdf_height == 0){
          pdf_height <- 3.5
        }else pdf_height <- input$multi_pdf_height
        if(input$multi_pdf_width == 0){
          pdf_width <- 4.7
        }else pdf_width <- input$multi_pdf_width
        plot_obj <- multi_umap_plot()
        if(inherits(plot_obj, "try-error")){
          stop(conditionMessage(attr(plot_obj, "condition")), call. = FALSE)
        }
        if(is.null(plot_obj)){
          stop("umap: plot is not available", call. = FALSE)
        }
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_obj)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  # 3 conditions ------------------------------------------------------------------------------
  org2 <- reactive({
    return(org(Species = input$Species2,Ortholog = input$Ortholog2))
  })
  ortholog2 <- reactive({
    return(no_org_ID(count = row_count_matrix2(),Species = input$Species2,Ortholog = input$Ortholog2,Biomart_archive=input$Biomart_archive2))
  })
  isoform2 <- reactive({
    return(isoform_ID(count = row_count_matrix2(),Species = input$Species2,Ortholog = input$Ortholog2,Biomart_archive=input$Biomart_archive2,RNA_type=input$Level_cond3))
  })
  gene_type2 <- reactive({
    return(gene_type(my.symbols=rownames(row_count_matrix2()),org=org2(),Species=input$Species2,RNA_type=input$Level_cond3))
  })
  org_code2 <- reactive({
    return(org_code(Species = input$Species2, Ortholog= input$Ortholog2))
  })
  observeEvent(pre_d_row_count_matrix2(),({
    updateSelectizeInput(session,inputId = "sample_order_cond3","Select samples:",
                         choices = colnames(pre_d_row_count_matrix2()),selected = colnames(pre_d_row_count_matrix2()))
  }))
  observeEvent(d_row_count_matrix2(), ({
    updateCollapse(session,id =  "input_collapse_panel2", open="D_row_count_matrix_panel2")
  }))
  d_row_count_matrix2 <- reactive({
    count <- pre_d_row_count_matrix2()
    order <- input$sample_order_cond3
    data <- try(safe_ordered_matrix(count, order), silent = TRUE)
    if(inherits(data, "try-error")) validate("")
    return(data)
  })
  row_count_matrix2 <- reactive({
    if (input$data_file_type2 == "Row3"){
      tmp <- input$file4$datapath
      if(is.null(input$file4) && input$goButton2 > 0 )  tmp = example_data_path("data/example4.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example4.txt")
      return(read_df(tmp = tmp))
    }
    if (input$data_file_type2 == "Row4"){
      tmp <- input$file5$datapath
      if(is.null(input$file5) && input$goButton2 > 0 )  tmp = example_data_path("data/example2.csv", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example2.csv")
      return(read_df(tmp = tmp))
    }
    if (input$data_file_type2 == "RowRecode_cond3"){
      tmp <- input$file_recode_cond3$datapath
      if(is.null(input$file_recode_cond3) && input$goButton2 > 0 )  {
        return(read_df(tmp = example_data_path("data/example4.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example4.txt")))
      }
      if(!is.null(tmp)){
        load(tmp)
        return(rawcount)
      }
    }
  })
  metadata2 <- reactive({
    if (input$data_file_type2 == "Row3"){
      return(NULL)
    }
    if (input$data_file_type2 == "Row4"){
      tmp <- input$file6$datapath
      if(is.null(input$file6) && input$goButton2 > 0 )  tmp = example_data_path("data/example6.csv", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example6.csv")
      df <- read_df(tmp = tmp)
      if(!is.null(df)) rownames(df) <- gsub("-",".",rownames(df))
      return(df)
    }
    if (input$data_file_type2 == "RowRecode_cond3"){
      tmp <- input$file_recode_cond3$datapath
      if(!is.null(tmp)){
        load(tmp)
        if(!is.null(metadata)) return(metadata) else return(NULL)
      }
    }
  })
  norm_count_matrix2 <- reactive({
    tmp <- input$norm_file2$datapath
    if (input$data_file_type2 == "RowRecode_cond3" && is.null(tmp)){
      recode <- input$file_recode_cond3$datapath
      if(!is.null(recode)){
        load(recode)
        return(norm_count_matrix) 
      }
    }
    df <- read_df(tmp = tmp)
    if(!is.null(df)){
      if(!is.null(input$sample_order_cond3)){
      if (input$data_file_type2 == "Row4"){
        meta <- anno_rep_meta(metadata2())
        if (!is.null(df) && !is.null(meta)){
          row_t <- t(df)
          meta <- data.frame(characteristics = meta[,1], row.names = rownames(meta))
          colname <- colnames(meta)
          data <- merge(meta, row_t, by=0, sort = F)
          if(dim(data)[1] == 0) {
            rownames(meta) <- gsub("\\.","-",rownames(meta))
            data <- merge(meta, row_t, by=0, sort = F)
            if(dim(data)[1] == 0) {
              rownames(row_t) <- gsub("\\.","-",rownames(row_t))
              data <- merge(meta, row_t, by=0, sort = F)
              validate("Error: failed to merge count data with metadata. Please check row names of matadata.")
            }
          }
          rownames(data) <- data$characteristics
          data2 <- data[, - which(colnames(data) %in% c("Row.names", colname))]
          data2_t <- t(data2)
          df <- apply(data2_t, 2, as.numeric)
          rownames(df) <- rownames(data2_t)
          df <- as.data.frame(df)
          if(input$Species2 != "not selected"){
            if(gene_type2() != "SYMBOL"){
              rownames(df) < gsub("\\..*","", rownames(df))
            }
          }
        }else return(NULL)
      }else df <- anno_rep(df)
      order <- input$sample_order_cond3
      df <- try(df[,order])
      if(length(df) == 1){
        if(inherits(df, "try-error")) validate("")
      }
      if(input$Species2 != "not selected"){
        if(gene_type2() != "SYMBOL"){
          rownames(df) < gsub("\\..*","", rownames(df))
        }
      }
      }
    }else df <- NULL
    return(df)
  })
  
  pre_d_row_count_matrix2 <- reactive({
    row <- row_count_matrix2()
    if (input$data_file_type2 == "Row3"){
      if(is.null(row)) {
        return(NULL)
      }else{
        return(anno_rep(row))
      }
    }
    if (input$data_file_type2 == "Row4"){
      meta <- anno_rep_meta(metadata2())
      if (is.null(row) || is.null(meta)){
        return(NULL)
      } else {
        row_t <- t(row)
        meta <- data.frame(characteristics = meta[,1], row.names = rownames(meta))
        colname <- colnames(meta)
        data <- merge(meta, row_t, by=0, sort = F)
        if(dim(data)[1] == 0) {
          rownames(meta) <- gsub("\\.","-",rownames(meta))
          data <- merge(meta, row_t, by=0, sort = F)
          if(dim(data)[1] == 0) {
            rownames(row_t) <- gsub("\\.","-",rownames(row_t))
            data <- merge(meta, row_t, by=0, sort = F)
            validate("Error: failed to merge count data with metadata. Please check row names of matadata.")
          }
        }
        rownames(data) <- data$characteristics
        data2 <- data[, - which(colnames(data) %in% c("Row.names", colname))]
        data2_t <- t(data2)
        data3 <- apply(data2_t, 2, as.numeric)
        rownames(data3) <- rownames(data2_t)
        if(input$Species2 != "not selected"){
          if(gene_type2() != "SYMBOL"){
            rownames(data3) <- gsub("\\..*","", rownames(data3))
          }
        }
        return(data3)
      }
    }
    if (input$data_file_type2 == "RowRecode_cond3"){
      tmp <- input$file_recode_cond3$datapath
      if(!is.null(tmp)){
        load(tmp)
        return(d_rawcount) 
      }else if(input$goButton2 >0) return(row_count_matrix2())
    }
  })
  
  observeEvent(input$goButton2,({
    updateSelectInput(session,inputId = "Species2","Species",species_list, selected = "Mus musculus")
  }))
  
  observeEvent(input$file7, ({
    updateCollapse(session,id =  "input_collapse_panel3", open="Norm_count_matrix_panel")
  }))
  observeEvent(input$file8, ({
    updateCollapse(session,id =  "input_collapse_panel3", open="norm_Metadata_panel")
  }))
  output$Row_count_matrix2 <- DT::renderDataTable({
    row_count_matrix2()
  })
  output$Metadata2 <- DT::renderDataTable({
    metadata2()
  })
  output$D_Row_count_matrix2 <- DT::renderDataTable({
    d_row_count_matrix2()
  })
  
  output$download_cond3_d_row_count = downloadHandler(
    filename = function() {
      if (input$data_file_type2 == "Row3"){
        paste(gsub("\\..+$", "", input$file4), "defined_row_count.txt", sep = "-")
      }else{
        paste(gsub("\\..+$", "", input$file5), paste0(gsub("\\..+$", "", input$file6), ".txt"), sep = "-")
      }
    },
    content = function(file){write.table(d_row_count_matrix2(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  
  gene_ID <- reactive({
    res <- d_row_count_matrix2()
    if(is.null(res)){
      return(NULL)
    }else{
      if(input$Species2 != "not selected") {
        if(gene_type2() != "SYMBOL"){
          gene_IDs <- ensembl2symbol(gene_type=gene_type2(),data = res,Species=input$Species2,
                                     Ortholog=ortholog2(),Isoform=isoform2(),org = org2(),merge=FALSE,rowname=FALSE)
          return(gene_IDs)
        }
      }else{ return(NULL) }
    }
  })
  
  
  
  # 3 conditions DEG ------------------------------------------------------------------------------
  cond3_group_info <- reactive({
    count <- d_row_count_matrix2()
    if(is.null(count)){
      return(NULL)
    }
    collist <- gsub("\\_.+$", "", colnames(count))
    if(length(unique(collist)) != 3){
      return(NULL)
    }
    vec <- c()
    for (i in 1:length(unique(collist))) {
      num <- length(collist[collist == unique(collist)[i]])
      vec <- c(vec, num)
    }
    conditions <- as.factor(rep(paste("C", 1:length(unique(collist)), sep=""), times = vec))
    patterns <- GetPatterns(conditions)
    rownames(patterns) <- tolower(rownames(patterns))
    eename <- rownames(patterns)[which(rowSums(patterns) == length(unique(collist)))]
    list(
      count = data.matrix(count),
      collist = collist,
      conditions = conditions,
      patterns = patterns,
      eename = eename
    )
  })

  cond3_multiout_cache_key <- reactive({
    if (input$data_file_type2 == "RowRecode_cond3"){
      return(NULL)
    }
    info <- cond3_group_info()
    if(is.null(info)){
      return(NULL)
    }
    digest::digest(list(
      EBSeq_mode = input$EBSeq_mode,
      count = info$count,
      conditions = as.character(info$conditions),
      patterns = info$patterns
    ))
  })

  cond3_multiout_cache <- reactiveValues(
    key = NULL,
    value = NULL
  )

  MultiOut <- reactive({
    if (input$data_file_type2 == "RowRecode_cond3"){
      tmp <- input$file_recode_cond3$datapath
      if(!is.null(tmp)){
        load(tmp)
      }else if(input$goButton2 >0) {
        demo_dds_path <- file.path(tempdir(),"dds.rds")
        validate(need(file.exists(demo_dds_path), "Example Recode.Rdata is not available in this deployment. Please upload a Recode.Rdata file."))
        dds <- readRDS(demo_dds_path)
      }
      return(dds)
    }else{
      if(input$EBSeq_mode == TRUE) time <- "a few" else time <- "5 - 10"
      withProgress(message = paste0("EBSeq multiple comparison test takes ",time," minutes"),{
        info <- cond3_group_info()
        cache_key <- cond3_multiout_cache_key()
        if(!is.null(info)){
          if(!is.null(cache_key) && identical(isolate(cond3_multiout_cache$key), cache_key) && !is.null(isolate(cond3_multiout_cache$value))){
            return(isolate(cond3_multiout_cache$value))
          }
          ngvector <- NULL
          Sizes <- MedianNorm(info$count)
          stopifnot(length(info$eename) == 1)
          MultiOut <- NULL
          MultiOut <- EBMultiTest(Data = info$count, NgVector = ngvector, fast=input$EBSeq_mode,
                                  Conditions = info$conditions, AllParti = info$patterns,
                                  sizeFactors = Sizes, maxround = 5)
          stopifnot(!is.null(MultiOut))
          cond3_multiout_cache$key <- cache_key
          cond3_multiout_cache$value <- MultiOut
          return(MultiOut)
        }
        incProgress(1)
      })
    }
  })
  
  cond3_multi_summary <- reactive({
    info <- cond3_group_info()
    if(is.null(info)){
      return(NULL)
    }
    multi_out <- MultiOut()
    MultiPP <- GetMultiPP(multi_out)
    PP <- as.data.frame(MultiPP$PP)
    colnames(PP) <- tolower(colnames(PP))
    pos <- which(names(PP) == info$eename)
    probs <- rowSums(PP[,-pos])
    results <- cbind(PP, MultiPP$MAP[rownames(PP)], probs)
    colnames(results) <- c(colnames(PP), "MAP", "PPDE")
    ord <- order(results[,"PPDE"], decreasing = TRUE)
    results <- results[ord,]
    list(
      results = as.data.frame(results),
      pattern = as.data.frame(MultiPP$Pattern),
      ord = ord,
      multi_out = multi_out,
      collist = info$collist
    )
  })

  deg_result2_raw <- reactive({
    summary <- cond3_multi_summary()
    if(is.null(summary)){
      return(NULL)
    }
    summary$results
  })

  deg_result2 <- reactive({
    res <- deg_result2_raw()
    if(is.null(res)){
      return(NULL)
    }
    if(input$Species2 != "not selected" && gene_type2() != "SYMBOL"){
      gene_IDs  <- gene_ID()
      res$Row.names <- rownames(res)
      data2 <- merge(res, gene_IDs, by="Row.names")
      rownames(data2) <- data2$Row.names
      res <- data2[,-1]
    }
    return(res)
  })
  deg_result2_pattern <- reactive({
    summary <- cond3_multi_summary()
    if(is.null(summary)){
      return(NULL)
    }
    return(summary$pattern)
  })
  deg_result2_condmean <- reactive({
    summary <- cond3_multi_summary()
    if(is.null(summary)){
      return(NULL)
    }
    if(input$EBSeq_mode == TRUE) {
      MultiFC <- GetMultiFC(summary$multi_out)
      res <- MultiFC$CondMeans[summary$ord,]
      rownames(res) <- rownames(summary$results)
    }else {
      res <- v1_GetMultiFC(summary$multi_out,collist=summary$collist,count=d_row_count_matrix2(),
                           EBSeq_mode = input$EBSeq_mode)
      res <- res[rownames(summary$results),]
    }
    return(as.data.frame(res))
  })
  
  deg_norm_count2_raw <- reactive({
    count <- d_row_count_matrix2()
    if(is.null(count)){
      return(NULL)
    }else{
      collist <- gsub("\\_.+$", "", colnames(count))
      if(length(unique(collist)) == 3){
        if(is.null(norm_count_matrix2())){
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
          normalized_counts <- as.data.frame(GetNormalizedMat(count, Sizes))
        }else{
          normalized_counts <- norm_count_matrix2()
        }
        return(normalized_counts)
      }else return(NULL)
    }
  })

  deg_norm_count2 <- reactive({
    normalized_counts <- deg_norm_count2_raw()
    if(is.null(normalized_counts)){
      return(NULL)
    }
    if(input$Species2 != "not selected" && gene_type2() != "SYMBOL"){
      gene_IDs  <- gene_ID()
      normalized_counts$Row.names <- rownames(normalized_counts)
      data2 <- merge(normalized_counts, gene_IDs, by="Row.names")
      data2$Unique_ID <- paste(data2$SYMBOL,data2$Row.names, sep = "\n- ")
      rownames(data2) <- data2$Row.names
      normalized_counts <- data2[,-1]
    }
    return(normalized_counts)
  })
  
  output$DEG_result2_1 <- DT::renderDataTable({
    data_3degcount2_1()
  })
  output$DEG_result2_2 <- DT::renderDataTable({
    data_3degcount2_2()
  })
  output$DEG_result2_3 <- DT::renderDataTable({
    data_3degcount2_3()
  })
  output$Normalized_Count_matrix2 <- DT::renderDataTable({
    deg_norm_count2()
  })
  
  download_cond3_dir <-reactive({
    if (input$data_file_type2 == "Row3"){
      dir_name <- paste(gsub("\\..+$", "", input$file4), "EBSeq" , sep ="-")
    }else{
      dir_name <- paste(gsub("\\..+$", "", input$file5), "EBSeq" , sep ="-")
    }
    dir_name <- paste0(dir_name, paste0("_fc", input$fc2))
    dir_name <- paste0(dir_name, paste0("_fdr", input$fdr2))
    dir_name <- paste0(dir_name, paste0("_basemean", input$basemean2))
    return(dir_name)
  })
  
  output$download_cond3_norm_count = downloadHandler(
    filename = function() {paste(download_cond3_dir(),"normalized_count.txt", sep = "-")},
    content = function(file){write.table(deg_norm_count2(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  
  output$download_3cond_DEG_table1 = downloadHandler(
    filename = function() {paste(download_cond3_dir(),"DEG_result1.txt", sep = "-")},
    content = function(file){write.table(data_3degcount2_1(), file, row.names = F, sep = "\t", quote = F)}
  )
  output$download_3cond_DEG_table2 = downloadHandler(
    filename = function() {paste(download_cond3_dir(),"DEG_result2.txt", sep = "-")},
    content = function(file){write.table(data_3degcount2_2(), file, row.names = F, sep = "\t", quote = F)}
  )
  output$download_3cond_DEG_table3 = downloadHandler(
    filename = function() {paste(download_cond3_dir(),"DEG_result3.txt", sep = "-")},
    content = function(file){write.table(data_3degcount2_3(), file, row.names = F, sep = "\t", quote = F)}
  )

  cond3_result_table <- reactive({
    data_3degcount1(
      gene_type = gene_type2(),
      data = deg_norm_count2(),
      result_Condm = deg_result2_condmean(),
      result_FDR = deg_result2(),
      specific = 1,
      result_list = TRUE
    ) %>% as.data.frame()
  })
  
  output$cond3_result <- DT::renderDataTable({
    cond3_result_table()
  })
  output$download_cond3_result = downloadHandler(
    filename = function() {paste(download_cond3_dir(),"DEG_result_ALL.txt", sep = "-")},
    content = function(file){
      table <- cond3_result_table()
      if(gene_type2() != "SYMBOL"){
        if(length(grep("SYMBOL", colnames(table))) != 0){
          table$Unique_ID <- gsub("\n"," ",table$Unique_ID)
        }
      }
      write.table(table,file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  
  
  #3conditions DEG vis------------------------
  output$not_cond3 <- renderText({
    count <- d_row_count_matrix2()
    if(!is.null(count)){
      collist <- gsub("\\_.+$", "", colnames(count))
      if(length(unique(collist)) != 3) print("Uploaded count data is in an inappropriate format. Please refer to the RNAseqChef manual for guidance and make the necessary corrections.")
    }
  })
  output$not_cond3_select <- renderText({
    count <- d_row_count_matrix2()
    if(!is.null(count)){
      collist <- gsub("\\_.+$", "", colnames(count))
      if(length(unique(collist)) != 3) print(paste0("Your count file contains data for ", length(unique(collist))," conditions. Please select samples to narrow it down to 3 conditions."))
    }
  })
  #3conditions DEG_1------------------------
  cond3_degcount1_list <- reactive({
    lapply(1:3, function(i) {
      data_3degcount1(
        gene_type = gene_type2(),
        data = deg_norm_count2(),
        result_Condm = deg_result2_condmean(),
        result_FDR = deg_result2(),
        specific = i,
        fc = input$fc2,
        fdr = input$fdr2,
        basemean = input$basemean2
      )
    })
  })

  cond3_degcount2_list <- reactive({
    data_list <- cond3_degcount1_list()
    lapply(data_list, function(data3) {
      data_3degcount2(
        gene_type = gene_type2(),
        data3 = data3,
        Species = input$Species2,
        Ortholog = ortholog2(),
        Isoform = isoform2(),
        org = org2()
      )
    })
  })

  data_3degcount1_1 <- reactive({
    cond3_degcount1_list()[[1]]
  })
  
  data_3degcount2_1 <- reactive({
    cond3_degcount2_list()[[1]]
  })
  
  #3conditions scatter + heatmap_1
  
  cond3_scatter1_plot <- reactive({
    p <- cond3_scatter_plot(gene_type=gene_type2(),data = deg_norm_count2(), data4 = data_3degcount2_1(),
                            result_Condm = deg_result2_condmean(), result_FDR = deg_result2(), specific = 1,
                            fc = input$fc2, fdr = input$fdr2, basemean = input$basemean2, Species = input$Species2)
    return(p)
  })
  
  output$scatter_1 <- renderPlot({
    if(is.null(deg_norm_count2())){
      return(NULL)
    }else{
      withProgress(message = "scatter plot",{
        suppressWarnings(print(cond3_scatter1_plot()))
      })
    }
  })
  
  output$download_3cond_scatter1 = downloadHandler(
    filename = function(){
      paste0(download_cond3_dir(), "_DEG_overview.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$cond3_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$cond3_pdf_height
        if(input$cond3_pdf_width == 0){
          pdf_width <- 10
        }else pdf_width <- input$cond3_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(cond3_scatter1_plot()) 
        print(cond3_scatter2_plot())
        print(cond3_scatter3_plot())
        dev.off()
      })
    }
  )
  #3conditions DEG_2------------------------
  data_3degcount1_2 <- reactive({
    cond3_degcount1_list()[[2]]
  })
  
  data_3degcount2_2 <- reactive({
    cond3_degcount2_list()[[2]]
  })
  
  #3conditions scatter + heatmap_2
  cond3_scatter2_plot <- reactive({
    p <- cond3_scatter_plot(gene_type=gene_type2(),data = deg_norm_count2(), data4 = data_3degcount2_2(),
                            result_Condm = deg_result2_condmean(), result_FDR = deg_result2(), specific = 2,
                            fc = input$fc2, fdr = input$fdr2, basemean = input$basemean2, Species = input$Species2)
    return(p)
  })
  
  output$scatter_2 <- renderPlot({
    if(is.null(deg_norm_count2())){
      return(NULL)
    }else{
      print(cond3_scatter2_plot())
    }
  })
  
  
  #3conditions DEG_3------------------------
  data_3degcount1_3 <- reactive({
    cond3_degcount1_list()[[3]]
  })
  
  data_3degcount2_3 <- reactive({
    cond3_degcount2_list()[[3]]
  })
  
  #3conditions scatter + heatmap_3
  cond3_scatter3_plot <- reactive({
    p <- cond3_scatter_plot(gene_type=gene_type2(),data = deg_norm_count2(), data4 = data_3degcount2_3(),
                            result_Condm = deg_result2_condmean(), result_FDR = deg_result2(), specific = 3,
                            fc = input$fc2, fdr = input$fdr2, basemean = input$basemean2, Species = input$Species2)
    return(p)
  })
  
  output$scatter_3 <- renderPlot({
    if(is.null(deg_norm_count2())){
      return(NULL)
    }else{
      print(cond3_scatter3_plot())
    }
  })
  
  
  output$download_3cond_scatter = downloadHandler(
    filename = function(){
      if(input$cond3_GOI_color_type == "default") paste0(download_cond3_dir(), "_GOI_scatter.pdf") else paste0(download_cond3_dir(), "_GOI_scatter_",input$cond3_GOI_color_pathway2,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$cond3_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$cond3_pdf_height
        if(input$cond3_pdf_width == 0){
          pdf_width <- 10
        }else pdf_width <- input$cond3_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(cond3_scatter_plot(gene_type=gene_type2(),data = deg_norm_count2(), data4 = data_3degcount2_1(),
                                 result_Condm = deg_result2_condmean(), result_FDR = deg_result2(), 
                                 fc = input$fc2, fdr = input$fdr2, basemean = input$basemean2,
                                 y_axis = input$cond3_scatter_yrange,x_axis = input$cond3_scatter_xrange,heatmap = FALSE,id_cut=input$cond3_uniqueID_cut,
                                 specific = cond3_specific_group2(), GOI = input$GOI2, Species = input$Species2,brush=brush_info_cond3(),
                                 GOI_color_type=input$cond3_GOI_color_type,cond3_pathway_color_gene=cond3_pathway_color_gene()))
        dev.off()
      })
    }
  )
  
  
  #3conditions PCA--------------
  cond3_pca_plot <- reactive({
    data <- deg_norm_count2_raw()
    if(is.null(data)){
      return(NULL)
    }
    PCAplot(data = data, legend = input$PCA_legend_cond3)
  })

  cond3_pca_table <- reactive({
    data <- deg_norm_count2_raw()
    if(is.null(data)){
      return(NULL)
    }
    PCAdata(row_count = data, deg_norm_count = data)
  })

  output$PCA2 <- renderPlot({
    plot_obj <- cond3_pca_plot()
    if(is.null(plot_obj)){
      return(NULL)
    }else{
      print(plot_obj)
    }
  })
  output$PCA3_data <- DT::renderDataTable({
    cond3_pca_table()
  })
  
  output$download_3cond_PCA = downloadHandler(
    filename = function(){
      paste0(download_cond3_dir(), "_PCA-MDS-dendrogram.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$cond3_pdf_height == 0){
          pdf_height <- 3.5
        }else pdf_height <- input$cond3_pdf_height
        if(input$cond3_pdf_width == 0){
          pdf_width <- 9
        }else pdf_width <- input$cond3_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(cond3_pca_plot())
        dev.off()
      })
    }
  )
  
  output$download_cond3_pca_table = downloadHandler(
    filename = function() {
      paste0(download_cond3_dir(), "PCA_table.txt")
    },
    content = function(file){write.table(cond3_pca_table(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  #3conditions GOI------------------------------------------------------
  GOI_list2 <- reactive({
    count <- deg_norm_count2()
    if(is.null(d_row_count_matrix2())){
      return(NULL)
    }else{
      if(is.null(input$cond3_GOI_color_type)) validate("")
      if(input$cond3_GOI_color_type != "default" && !is.null(cond3_pathway_color_gene())){
        count <- count %>% dplyr::filter(row.names(.) %in% rownames(cond3_pathway_color_gene()))
      }
      if(gene_type2() != "SYMBOL"){
        if(input$Species2 != "not selected"){
          GOI <- count$Unique_ID
        }else GOI <- rownames(count)
      }else{
        if(input$Species2 != "not selected"){
          GOI <- rownames(count)
        }else GOI <- rownames(count)
      }
      return(normalize_goi_choices(GOI))
    }
  })

  cond3_specific_group <- reactive({
    if(is.null(d_row_count_matrix2())){
      return(NULL)
    }else{
      return(unique(gsub("\\_.+$", "", colnames(d_row_count_matrix2()))))
    }
  })
  
  output$cond3_GOI_pair<- renderUI({
    if(is.null(d_row_count_matrix2())){
      return(NULL)
    }else{
      withProgress(message = "Preparing GOI pair",{
        selectInput("cond3_GOIpair", "Group of interest:", c(cond3_specific_group()),multiple = F)
      })
    }
  })
  
  observeEvent(GOI_list2(), {
    choices <- GOI_list2()
    if(is.null(choices)){
      choices <- character(0)
    }
    selected <- intersect(isolate(input$GOI2), choices)
    updateSelectizeInput(session, "GOI2", choices = choices, selected = selected,
                         options = goi_selectize_options, server = TRUE)
  }, ignoreNULL = FALSE)
  output$cond3_GOI_color_type <- renderUI({
    if(input$Species2 == "not selected"){
      radioButtons("cond3_GOI_color_type","Filter", c("All genes" = "default"),selected="default")
    }else{
      radioButtons("cond3_GOI_color_type","Filter", c("All genes" = "default","Pathway of interest"="pathway"),selected="default")
    }
  })
  observeEvent(input$Species2, ({
    if(input$Species2 == "not selected"){
      updateSelectInput(session, "cond3_GOI_color_pathway1","Select a gene set for gene extraction","")
    }else  if(input$Species2 != "Xenopus laevis" && input$Ortholog2 != "Arabidopsis thaliana" && input$Species2 != "Arabidopsis thaliana"){
      updateSelectInput(session, "cond3_GOI_color_pathway1","Select a gene set for gene extraction",gene_set_list) 
    }else {
      updateSelectInput(session, "cond3_GOI_color_pathway1","Select a gene set for gene extraction",c("KEGG", "GO biological process", 
                                                                                                "GO cellular component","GO molecular function")) 
    }
  }))
  cond3_pathway_color_gene_list <- reactive({
    genes <- GeneList_for_enrichment(Species = input$Species2, Ortholog = input$Ortholog2,
                                     Gene_set=input$cond3_GOI_color_pathway1, org = org2(), 
                                     Biomart_archive=input$Biomart_archive,gene_type=gene_type2())
    return(genes)
  })
  cond3_GOI_color_pathway_list <- reactive({
    list <-cond3_pathway_color_gene_list()
    if(!is.null(input$cond3_GOI_color_pathway1)){
      if(input$Species2 == "Xenopus laevis" || input$Ortholog2 == "Arabidopsis thaliana" || input$Species2 == "Arabidopsis thaliana"){
        print(head(list))
        if(input$cond3_GOI_color_pathway1 == "GO biological process") list <- list %>% dplyr::filter(Ontology == "BP")
        if(input$cond3_GOI_color_pathway1 == "GO cellular component") list <- list %>% dplyr::filter(Ontology == "CC")
        if(input$cond3_GOI_color_pathway1 == "GO molecular function") list <- list %>% dplyr::filter(Ontology == "MF") 
      }}
    list <- unique(list$gs_name)
    return(list)
  })
  observeEvent(input$cond3_GOI_color_pathway1, ({
    if(is.null(input$cond3_GOI_color_pathway1)){
      updateSelectInput(session, "cond3_GOI_color_pathway2","","")
    }else{
      updateSelectInput(session, "cond3_GOI_color_pathway2","",cond3_GOI_color_pathway_list())
    }
  }))
  output$GOIreset_cond3 <- renderUI({
    actionButton("GOIreset_cond3", "GOI reset")
  })
  observeEvent(input$GOIreset_cond3, {
    choices <- GOI_list2()
    if(is.null(choices)){
      choices <- character(0)
    }
    updateSelectizeInput(session, "GOI2", choices = choices,
                         selected = character(0),
                         options = goi_selectize_options, server = TRUE)
  })
  
  output$cond3_xrange <- renderUI({
    data <- range_for_GOIscatter()
    if(!is.null(data)){
      min <- floor(min(data$FC_x))
      max <- ceiling(max(data$FC_x))
      sliderInput("cond3_scatter_xrange","X_axis range:",min = min-1,
                  max=max+1, step = 0.5,
                  value = c(min, max))
    }
  })
  output$cond3_yrange <- renderUI({
    data <- range_for_GOIscatter()
    if(!is.null(data)){
      min <- floor(min(data$FC_y))
      max <- ceiling(max(data$FC_y))
      sliderInput("cond3_scatter_yrange","Y_axis range:",min = min-1,
                  max=max+1, step = 0.5,
                  value = c(min, max))
    }
  })
  cond3_pathway_color_gene <- reactive({
    ##extract pathway genes
    if(is.null(input$cond3_GOI_color_pathway1) || is.null(input$cond3_GOI_color_pathway2) || 
       is.null(gene_type2()) || is.null(org2()) || input$cond3_GOI_color_pathway1 == "" || 
       input$cond3_GOI_color_pathway2 == "" || !input$cond3_GOI_color_pathway2 %in% cond3_GOI_color_pathway_list()) validate("")
    genes <- cond3_pathway_color_gene_list()
    genes <- try(dplyr::filter(genes, gs_name == input$cond3_GOI_color_pathway2))
    if(length(genes) == 1) if(class(genes)=="try-error") validate("")
    
    my.symbols <- as.character(genes$entrez_gene)
    if(gene_type2() == "non-model organism"){
      gene_IDs <-  try(dplyr::filter(ortholog2(), ENTREZID %in% my.symbols))
      if(length(gene_IDs) == 1) if(class(gene_IDs)=="try-error") validate("")
      df <- data.frame(gene = gene_IDs$ENSEMBL, row.names = gene_IDs$ENSEMBL)
    }else if(gene_type2() == "isoform"){
      gene_IDs <-  try(dplyr::filter(isoform2(), ENTREZID %in% my.symbols))
      if(length(gene_IDs) == 1) if(class(gene_IDs)=="try-error") validate("")
      df <- data.frame(gene = gene_IDs$Transcript_ID, row.names = gene_IDs$Transcript_ID)
    }else{
      if(gene_type2() == "SYMBOL") columns <- c("ENTREZID", sgd_symbol_column(org2())) else columns <- c("ENTREZID","ENSEMBL")
      if(org2()$packageName == "org.At.tair.db") {
        if(gene_type2() == "SYMBOL"){
          gene_IDs <- AnnotationDbi::select(org2(), keys = my.symbols,
                                            keytype = "TAIR",
                                            columns = c("TAIR","SYMBOL"))
          colnames(gene_IDs) <- c("TAIR","GeneID")
        }else{
          gene_IDs <- data.frame(TAIR = my.symbols, GeneID = my.symbols)
        }
      } else {
        gene_IDs <- AnnotationDbi::select(org2(), keys = my.symbols,
                                          keytype = "ENTREZID",
                                          columns = columns)
        colnames(gene_IDs) <- c("entrezid","GeneID")
        gene_IDs <- na.omit(gene_IDs) 
      }
      gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
      df <- data.frame(gene = gene_IDs$GeneID, row.names = gene_IDs$GeneID)
      if(dim(df)[1] == 0) validate("No filtered genes.")
    }
    return(df)
  })
  output$cond3_uniqueID_cut <- renderUI({
    if(gene_type2() != "SYMBOL" && input$Species2 != "not selected") 
      radioButtons("cond3_uniqueID_cut","Short unique ID",c("ON"=TRUE,"OFF"=FALSE),selected = TRUE)
  })
  
  cond3_specific_group2 <- reactive({
    if(!is.null(input$cond3_GOIpair)){
      return(which(cond3_specific_group() == input$cond3_GOIpair))
    }
  })
  range_for_GOIscatter <- reactive({
    if(!is.null(cond3_specific_group2())){
      return(cond3_scatter_range(gene_type=gene_type2(),data = deg_norm_count2(), data4 = data_3degcount2_1(),
                                 result_Condm = deg_result2_condmean(), result_FDR = deg_result2(), 
                                 fc = input$fc2, fdr = input$fdr2, basemean = input$basemean2,specific = cond3_specific_group2()))
    }
  })
  output$cond3_GOIscatter <- renderPlot({
    if(!is.null(input$cond3_GOIpair) && !is.null(input$cond3_scatter_yrange) && !is.null(input$cond3_scatter_xrange)){
      cond3_scatter_plot(gene_type=gene_type2(),data = deg_norm_count2(), data4 = data_3degcount2_1(),
                         result_Condm = deg_result2_condmean(), result_FDR = deg_result2(), 
                         fc = input$fc2, fdr = input$fdr2, basemean = input$basemean2,
                         y_axis = input$cond3_scatter_yrange,x_axis = input$cond3_scatter_xrange,heatmap = FALSE,id_cut=input$cond3_uniqueID_cut,
                         specific = cond3_specific_group2(), GOI = input$GOI2, Species = input$Species2,brush=brush_info_cond3(),
                         GOI_color_type=input$cond3_GOI_color_type,cond3_pathway_color_gene=cond3_pathway_color_gene())
    }
  })
  
  brush_info_cond3 <- reactive({
    if(!is.null(input$cond3_GOIpair) && !is.null(input$cond3_scatter_yrange) && !is.null(input$cond3_scatter_xrange)){
      gene_type=gene_type2()
      data = deg_norm_count2()
      data4 = data_3degcount2_1()
      result_Condm = deg_result2_condmean() 
      result_FDR = deg_result2()
      fc = input$fc2
      fdr = input$fdr2
      basemean = input$basemean2
      Species = input$Species2
      specific = cond3_specific_group2()
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
      if(gene_type != "SYMBOL"){
        if(length(grep("SYMBOL", colnames(data))) != 0){
          data <- data[, - which(colnames(data) == "SYMBOL")]
        }
      }
      if(specific == 1) {
        specific = collist[1]
        FC_xlab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[2], ")"))
        FC_ylab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[3], ")"))
        result_Condm$FC_x <- log2((result_Condm$C1 + 0.01)/(result_Condm$C2 + 0.01))
        result_Condm$FC_y <- log2((result_Condm$C1 + 0.01)/(result_Condm$C3 + 0.01))
        Pattern1 <- "Pattern4"
        Pattern2 <- "Pattern5"
      }
      if(specific == 2) {
        specific = collist[2]
        FC_xlab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[1], ")"))
        FC_ylab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[3], ")"))
        result_Condm$FC_x <- log2((result_Condm$C2 + 0.01)/(result_Condm$C1 + 0.01))
        result_Condm$FC_y <- log2((result_Condm$C2 + 0.01)/(result_Condm$C3 + 0.01))
        Pattern1 <- "Pattern3"
        Pattern2 <- "Pattern5"
      }
      if(specific == 3) {
        specific = collist[3]
        FC_xlab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[1], ")"))
        FC_ylab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[2], ")"))
        result_Condm$FC_x <- log2((result_Condm$C3 + 0.01)/(result_Condm$C1 + 0.01))
        result_Condm$FC_y <- log2((result_Condm$C3 + 0.01)/(result_Condm$C2 + 0.01))
        Pattern1 <- "Pattern2"
        Pattern2 <- "Pattern5"
      }
      result_FDR$FDR <- 1 - result_FDR$PPDE
      result <- merge(result_Condm, result_FDR, by=0)
      data$Row.names <- rownames(data)
      data2 <- merge(data, result, by="Row.names")
      result <- dplyr::filter(data2, apply(data2[,2:(Cond_1 + Cond_2 + Cond_3)],1,mean) > basemean)
      sig <- rep(3, nrow(result))
      sig[which(result$FDR <= fdr & result$FC_x < log2(1/fc) & result$FC_y < log2(1/fc) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 2
      sig[which(result$FDR <= fdr & result$FC_x > log2(fc) & result$FC_y > log2(fc) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 1
      data3 <- data.frame(Row.names = result$Row.names, FC_x = result$FC_x,
                          FC_y = result$FC_y, padj = result$FDR, sig = sig, FC_xy = result$FC_x * result$FC_y)
      if((sum(sig == 1) >= 1) && (sum(sig == 2) >= 1)){
        new.levels <- c( paste0(paste0(specific,"_high: "), sum(sig == 1)), paste0(paste0(specific,"_low: "), sum(sig == 2)), "NS" )
        col = c("red","blue", "darkgray")}
      if((sum(sig == 1) >= 1) && (sum(sig == 2) == 0)){
        new.levels <- c(paste0(paste0(specific,"_high: "), sum(sig == 1)), "NS" )
        col = c("red", "darkgray")}
      if((sum(sig == 1) == 0) && (sum(sig == 2) >= 1)){
        new.levels <- c(paste0(paste0(specific,"_low: "), sum(sig == 2)), "NS" )
        col = c("blue", "darkgray")}
      if((sum(sig == 1) == 0) && (sum(sig == 2) == 0)){
        new.levels <- c("NS")
        col = "darkgray"}
      data3$sig <- factor(data3$sig, labels = new.levels)
      if(gene_type != "SYMBOL"){
        if(Species != "not selected"){
          data3 <- merge(data3, data, by="Row.names")
        }
      }
      return(brushedPoints(data3, input$plot1_brush_cond3,xvar = "FC_x",yvar="FC_y"))
    }
  })
  
  
  cond3_GOIcount <- reactive({
    count <- deg_norm_count2()
    if(gene_type2() != "SYMBOL"){
      if(input$Species2 != "not selected"){
        if(!is.null(brush_info_cond3())){
          brush <- brush_info_cond3()
          if(dim(brush)[1] != 0) {
            if(input$cond3_GOI_color_type != "default" && !is.null(cond3_pathway_color_gene())){
              brush <- brush %>% dplyr::filter(Row.names %in% rownames(cond3_pathway_color_gene()))
            }
            Unique_ID <- brush$Unique_ID 
            }else Unique_ID <- input$GOI2
        }
        label_data <- as.data.frame(Unique_ID, row.names = Unique_ID)
        data <- merge(count, label_data, by="Unique_ID")
        if(input$cond3_uniqueID_cut) {
          id_list <- gsub("\\\n.+$", "", data$Unique_ID)
          dup_list <- unique(id_list[duplicated(id_list)])
          for(i in 1:length(data$Unique_ID)){
            if(gsub("\\\n.+$", "", data$Unique_ID[i]) == "NA") {
              data$Unique_ID[i] <- gsub(".+\\s", "", data$Unique_ID[i])
            }else if(! gsub("\\\n.+$", "", data$Unique_ID[i]) %in% dup_list) {
              data$Unique_ID[i] <- gsub("\\\n.+$", "", data$Unique_ID[i])
            }
          }
        }
        rownames(data) <- data$Unique_ID
        data <- data[, - which(colnames(data) == "SYMBOL")]
        data <- data[,-1]
      }else{
        if(!is.null(brush_info_cond3())){
          brush <- brush_info_cond3()
          if(dim(brush)[1] != 0) {
            if(input$cond3_GOI_color_type != "default" && !is.null(cond3_pathway_color_gene())){
              brush <- brush %>% dplyr::filter(Row.names %in% rownames(cond3_pathway_color_gene()))
            }
            Row.names <- brush$Row.names 
          }else Row.names <- input$GOI2
        }
        count$Row.names <- rownames(count)
        label_data <- as.data.frame(Row.names, row.names = Row.names)
        data <- merge(count, label_data, by="Row.names")
        rownames(data) <- data$Row.names
        data <- data[,-1]
      }
    }else{
      if(!is.null(brush_info_cond3())){
        brush <- brush_info_cond3()
        if(dim(brush)[1] != 0) {
          if(input$cond3_GOI_color_type != "default" && !is.null(cond3_pathway_color_gene())){
            brush <- brush %>% dplyr::filter(Row.names %in% rownames(cond3_pathway_color_gene()))
          }
          Row.names <- brush$Row.names 
        }else Row.names <- input$GOI2
      }
      count$Row.names <- rownames(count)
      label_data <- as.data.frame(Row.names, row.names = Row.names)
      data <- merge(count, label_data, by="Row.names")
      rownames(data) <- data$Row.names
      data <- data[,-1]
    }
    return(data)
  })
  
  cond3_GOIheat <- reactive({
    data <- cond3_GOIcount()
    if(is.null(data)){
      ht <- NULL
    }else{
      data.z <- genefilter::genescale(data, axis=1, method="Z")
      data.z <- na.omit(data.z)
      ht <- GOIheatmap(data.z)
    }
    return(ht)
  })
  
  output$cond3_GOIheatmap <- renderPlot({
    if(is.null(d_row_count_matrix2())){
      return(NULL)
    }else{
      if(!is.null(brush_info_cond3())){
        if(!is.null(input$GOI2) || dim(brush_info_cond3())[1] != 0){
          withProgress(message = "Boxplot",{
            suppressWarnings(print(cond3_GOIheat()))
            incProgress(1)
          })
        }}
    }
  })
  
  cond3_GOIbox <- reactive({
    count <- deg_norm_count2()
    data <- cond3_GOIcount()
    if(is.null(data)){
      p <- NULL
    }else{
      p <- GOIboxplot(data = data) + scale_fill_manual(values=c("gray", "#4dc4ff", "#ff8082"))
    }
    return(p)
  })
  
  
  output$cond3_GOIboxplot <- renderPlot({
    if(is.null(d_row_count_matrix2())){
      return(NULL)
    }else{
      if(!is.null(brush_info_cond3())){
        if(!is.null(input$GOI2) || dim(brush_info_cond3())[1] != 0){
          withProgress(message = "Boxplot",{
            suppressWarnings(print(cond3_GOIbox()))
            incProgress(1)
          })
        }}
    }
  })
  
  output$download_3cond_GOIbox = downloadHandler(
    filename = function(){
      if(input$cond3_GOI_color_type == "default") paste0(download_cond3_dir(), "_GOIboxplot.pdf") else paste0(download_cond3_dir(), "_GOIboxplot_",input$cond3_GOI_color_pathway2,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- cond3_GOIcount()
        rowlist <- rownames(data)
        if(input$cond3_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)
        }else pdf_height <- input$cond3_pdf_height
        if(input$cond3_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)
        }else pdf_width <- input$cond3_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(cond3_GOIbox())
        dev.off()
        
        incProgress(1)
      })
    }
  )
  output$download_3cond_GOIheat = downloadHandler(
    filename = function(){
      if(input$cond3_GOI_color_type == "default") paste0(download_cond3_dir(), "_GOIheatmap.pdf") else paste0(download_cond3_dir(), "_GOIheatmap_",input$cond3_GOI_color_pathway2,".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- cond3_GOIcount()
        rowlist <- rownames(data)
        if(input$cond3_pdf_height == 0){
          pdf_height <- 10
        }else pdf_height <- input$cond3_pdf_height
        if(input$cond3_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$cond3_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(cond3_GOIheat())
        dev.off()
        incProgress(1)
      })
    }
  )
  #3conditions enrichment ------------------------------------------------------------------------------
  output$Spe2 <- renderText({
    if(input$Species2 == "not selected") print("Please select 'Species'")
  })
  #3conditions enrichment_1 ------------------------------------------------------------------------------
  enrich3_1 <- reactive({
    if(input$Species2 != "Xenopus laevis" && input$Ortholog2 !=  "Arabidopsis thaliana" && input$Species2 !=  "Arabidopsis thaliana"){
      return(keggEnrichment1(data3 = data_3degcount1_1(),data4 = data_3degcount2_1(),
                             Species = input$Species2, Gene_set = input$Gene_set2,
                             org = org2(),H_t2g = Hallmark_cond3()))
    }else return(keggEnrichment1_xenopus(data3 = data_3degcount1_1(),data4 = data_3degcount2_1(),
                                         Species = input$Species2, Gene_set = input$Gene_set2,
                                         org = org2(),org_code = org_code2()))
  })
  
  enrichment3_1_1 <- reactive({
    return(enrichment3_1(data3 = data_3degcount1_1(),data4 = data_3degcount2_1(),
                         cnet_list2 = enrich3_1()))
  })
  
  keggEnrichment2_1 <- reactive({
    return(keggEnrichment2(data3 = data_3degcount1_1(), data4 = data_3degcount2_1(), 
                           cnet_list2 = enrich3_1()))
  })
  
  output$keggenrichment2_1 <- renderPlot({
    if(!is.null(input$Gene_set2) && input$Species2 != "not selected"){
      withProgress(message = "Enrichment analysis",{
        p <- keggEnrichment2_1()
        if(!is.null(p)){
          print(p) 
        }
        incProgress(1)
      })
    }
  })
  
  #3conditions enrichment_2 ------------------------------------------------------------------------------
  Hallmark_cond3 <- reactive({
    return(GeneList_for_enrichment(Species = input$Species2, Ortholog=input$Ortholog2,Biomart_archive=input$Biomart_archive2, Gene_set = input$Gene_set2, org = org2()))
  })
  
  enrich3_2 <- reactive({
    if(input$Species2 != "Xenopus laevis" && input$Ortholog2 !=  "Arabidopsis thaliana" && input$Species2 !=  "Arabidopsis thaliana"){
      return(keggEnrichment1(data3 = data_3degcount1_2(),data4 = data_3degcount2_2(),
                             Species = input$Species2, Gene_set = input$Gene_set2,
                             org = org2(),H_t2g = Hallmark_cond3()))
    }else  return(keggEnrichment1_xenopus(data3 = data_3degcount1_2(),data4 = data_3degcount2_2(),
                                          Species = input$Species2, Gene_set = input$Gene_set2,
                                          org = org2(),org_code = org_code2()))
  })
  enrichment3_2_1 <- reactive({
    return(enrichment3_1(data3 = data_3degcount1_2(),data4 = data_3degcount2_2(),
                         cnet_list2 = enrich3_2()))
  })
  keggEnrichment2_2 <- reactive({
    return(keggEnrichment2(data3 = data_3degcount1_2(), data4 = data_3degcount2_2(), 
                           cnet_list2 = enrich3_2()))
  })
  
  output$keggenrichment2_2 <- renderPlot({
    if(!is.null(input$Gene_set2) && input$Species2 != "not selected"){
      p <- keggEnrichment2_2()
      if(!is.null(p)){
        print(p) 
      }
    }
  })
  #3conditions enrichment_3 ------------------------------------------------------------------------------
  enrich3_3 <- reactive({
    if(input$Species2 != "Xenopus laevis" && input$Ortholog2 !=  "Arabidopsis thaliana" && input$Species2 !=  "Arabidopsis thaliana"){
      return(keggEnrichment1(data3 = data_3degcount1_3(),data4 = data_3degcount2_3(),
                             Species = input$Species2, Gene_set = input$Gene_set2,
                             org = org2(),H_t2g = Hallmark_cond3()))
    }else return(keggEnrichment1_xenopus(data3 = data_3degcount1_3(),data4 = data_3degcount2_3(),
                                         Species = input$Species2, Gene_set = input$Gene_set2,
                                         org = org2(),org_code = org_code2()))
  })
  enrichment3_3_1 <- reactive({
    data <- enrichment3_1(data3 = data_3degcount1_3(),data4 = data_3degcount2_3(),
                          cnet_list2 = enrich3_3())
    return(data)
  })
  keggEnrichment2_3 <- reactive({
    return(keggEnrichment2(data3 = data_3degcount1_3(), data4 = data_3degcount2_3(), 
                           cnet_list2 = enrich3_3()))
  })
  
  output$keggenrichment2_3 <- renderPlot({
    if(!is.null(input$Gene_set2) && input$Species2 != "not selected"){
      p <- keggEnrichment2_3()
      if(!is.null(p)){
        print(p) 
      }
    }
  })
  
  output$Gene_set2 <- renderUI({
    if(input$Species2 != "Xenopus laevis" && input$Ortholog2 !=  "Arabidopsis thaliana" && input$Species2 !=  "Arabidopsis thaliana"){
      selectInput('Gene_set2', 'Gene Set', gene_set_list)
    }else selectInput('Gene_set2', 'Gene Set', c("KEGG", "GO biological process", 
                                                 "GO cellular component","GO molecular function"))
  })
  
  cond3_enrich_table1 <- reactive({
    if(input$Species2 != "Xenopus laevis" && input$Ortholog2 !=  "Arabidopsis thaliana" && input$Species2 !=  "Arabidopsis thaliana"){
      return(enrich_for_table(data = as.data.frame(enrichment3_1_1()), H_t2g = Hallmark_cond3(), Gene_set = input$Gene_set2))
    }else return(as.data.frame(enrichment3_1_1()))
  })
  cond3_enrich_table2 <- reactive({
    if(input$Species2 != "Xenopus laevis" && input$Ortholog2 !=  "Arabidopsis thaliana" && input$Species2 !=  "Arabidopsis thaliana"){
      return(enrich_for_table(data = as.data.frame(enrichment3_2_1()), H_t2g = Hallmark_cond3(), Gene_set = input$Gene_set2))
    }else return(as.data.frame(enrichment3_2_1()))
  })
  cond3_enrich_table3 <- reactive({
    if(input$Species2 != "Xenopus laevis" && input$Ortholog2 !=  "Arabidopsis thaliana" && input$Species2 !=  "Arabidopsis thaliana"){
      return(enrich_for_table(data = as.data.frame(enrichment3_3_1()), H_t2g = Hallmark_cond3(), Gene_set = input$Gene_set2))
    }else return(as.data.frame(enrichment3_3_1()))
  })
  
  
  output$enrichment3_result_1 <- DT::renderDataTable({
    cond3_enrich_table1()
  })
  output$enrichment3_result_2 <- DT::renderDataTable({
    cond3_enrich_table2()
  })
  output$enrichment3_result_3 <- DT::renderDataTable({
    cond3_enrich_table3()
  })
  
  output$download_cond3_enrichment_table1 = downloadHandler(
    filename = function() {
      paste(download_cond3_dir(), paste0(input$Gene_set2,"-enrichment1.txt"), sep="_")
    },
    content = function(file){write.table(cond3_enrich_table1(), file, row.names = F, sep = "\t", quote = F)}
  )
  output$download_cond3_enrichment_table2 = downloadHandler(
    filename = function() {
      paste(download_cond3_dir(), paste0(input$Gene_set2,"-enrichment2.txt"), sep="_")
    },
    content = function(file){write.table(cond3_enrich_table2(), file, row.names = F, sep = "\t", quote = F)}
  )
  output$download_cond3_enrichment_table3 = downloadHandler(
    filename = function() {
      paste(download_cond3_dir(), paste0(input$Gene_set2,"-enrichment3.txt"), sep="_")
    },
    content = function(file){write.table(cond3_enrich_table3(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  output$download_cond3_enrichment = downloadHandler(
    filename = function(){
      paste(download_cond3_dir(), paste0(input$Gene_set2,"-enrichment.pdf"), sep="_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- keggEnrichment2_1()
        p2 <- keggEnrichment2_2()
        p3 <- keggEnrichment2_3()
        if(input$cond3_pdf_height == 0){
          pdf_height <- 12
        }else pdf_height <- input$cond3_pdf_height
        if(input$cond3_pdf_width == 0){
          pdf_width <- 15
        }else pdf_width <- input$cond3_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(p1, p2, p3, nrow =3))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_3cond_report = downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%m%d_"),download_cond3_dir(), "_3conditionsDEG",".zip")
    },
    content = function(fname){
      withProgress(message = "Preparing download, please wait",{
        fs <- c()
        setwd(tempdir())
        print(fs)
        dir.create("DEG_result/",showWarnings = FALSE)
        dir.create("Clustering/",showWarnings = FALSE)
        DEG <- "DEG_result/DEG_result.txt"
        result1 <- "DEG_result/DEG_signature1.txt"
        result2 <- "DEG_result/DEG_signature2.txt"
        result3 <- "DEG_result/DEG_signature3.txt"
        count <- "DEG_result/normalized_count.txt"
        PCA <- "Clustering/clustering.pdf"
        PCA_table <- "Clustering/pca.txt"
        scatter <- "DEG_result/scatter_plot.pdf"
        fs <- c(DEG, result1,result2,result3,count,PCA,PCA_table,scatter)
        table <- data_3degcount1(gene_type=gene_type2(),data = deg_norm_count2(),result_Condm = deg_result2_condmean(),
                                 result_FDR = deg_result2(), specific = 1,result_list=TRUE)
        if(gene_type2() != "SYMBOL"){
          if(length(grep("SYMBOL", colnames(table))) != 0){
            table$Unique_ID <- gsub("\n"," ",table$Unique_ID)
          }
        }
        write.table(table, DEG, row.names = T, col.names=NA, sep = "\t", quote = F)
        write.table(deg_norm_count2(), count, row.names = T, col.names=NA, sep = "\t", quote = F)
        write.table(data_3degcount2_1(), result1, row.names = F, sep = "\t", quote = F)
        write.table(data_3degcount2_2(), result2, row.names = F, sep = "\t", quote = F)
        write.table(data_3degcount2_3(), result3, row.names = F, sep = "\t", quote = F)
        write.table(PCAdata(row_count = deg_norm_count2(), deg_norm_count = deg_norm_count2()), PCA_table, row.names = T, col.names=NA, sep = "\t", quote = F)
        pdf(PCA, height = 3.5, width = 9)
        print(PCAplot(data = deg_norm_count2(),legend = input$PCA_legend_cond3))
        dev.off()
        pdf(scatter, height = 6, width = 10)
        print(cond3_scatter1_plot()) 
        print(cond3_scatter2_plot())
        print(cond3_scatter3_plot())
        dev.off()
        if(!is.null(input$cond3_GOIpair) && !is.null(input$cond3_scatter_yrange) && !is.null(input$cond3_scatter_xrange)){
          dir.create("GOI_profiling/",showWarnings = FALSE)
          if(input$cond3_GOI_color_type == "default") goiscatter <- "GOI_profiling/scatter_plot.pdf" else goiscatter <- paste0("GOI_profiling/scatter_plot_",input$cond3_GOI_color_pathway2,".pdf") 
          fs <- c(fs,goiscatter)
          pdf(goiscatter, height = 6, width = 10)
          print(cond3_scatter_plot(gene_type=gene_type2(),data = deg_norm_count2(), data4 = data_3degcount2_1(),
                                   result_Condm = deg_result2_condmean(), result_FDR = deg_result2(), 
                                   fc = input$fc2, fdr = input$fdr2, basemean = input$basemean2,
                                   y_axis = input$cond3_scatter_yrange,x_axis = input$cond3_scatter_xrange,heatmap = FALSE,id_cut=input$cond3_uniqueID_cut,
                                   specific = cond3_specific_group2(), GOI = input$GOI2, Species = input$Species2,brush=brush_info_cond3(),
                                   GOI_color_type=input$cond3_GOI_color_type,cond3_pathway_color_gene=cond3_pathway_color_gene()))
          dev.off()
          if(!is.null(brush_info_cond3())){
            if(!is.null(input$GOI2) || dim(brush_info_cond3())[1] != 0){
              if(input$cond3_GOI_color_type == "default") {
                boxplot <- "GOI_profiling/boxplot.pdf"
                heat <- "GOI_profiling/heatmap.pdf"
              }else {
                boxplot <- paste0("GOI_profiling/boxplot_",input$cond3_GOI_color_pathway2,".pdf") 
                heat <- paste0("GOI_profiling/heatmap_",input$cond3_GOI_color_pathway2,".pdf") 
              }
              fs <- c(fs,boxplot,heat)
              data <- cond3_GOIcount()
              rowlist <- rownames(data)
              pdf_height <- pdf_h(rowlist)
              pdf_width <- pdf_w(rowlist)
              pdf(boxplot, height = pdf_height, width = pdf_width)
              print(cond3_GOIbox())
              dev.off()
              pdf(heat, height = 10, width = 7)
              print(cond3_GOIheat())
              dev.off()
            }}
        }
        if(input$Species2 != "not selected" && !is.null(input$Gene_set2)){
          dir.create("Enrichment_analysis/",showWarnings = FALSE)
          enrichplot <- paste0("Enrichment_analysis/Enrichment_analysis_",input$Gene_set2,".pdf")
          enrich_table1 <- paste0("Enrichment_analysis/Enrichment_analysis_",input$Gene_set2,"_1.txt")
          enrich_table2 <- paste0("Enrichment_analysis/Enrichment_analysis_",input$Gene_set2,"_2.txt")
          enrich_table3 <- paste0("Enrichment_analysis/Enrichment_analysis_",input$Gene_set2,"_3.txt")
          fs <- c(fs,enrichplot,enrich_table1,enrich_table2,enrich_table3)
          p1 <- keggEnrichment2_1()
          p2 <- keggEnrichment2_2()
          p3 <- keggEnrichment2_3()
          pdf(enrichplot, height = 12, width = 15)
          print(plot_grid(p1, p2, p3, nrow =3))
          dev.off()
          write.table(cond3_enrich_table1(), enrich_table1, row.names = F, sep = "\t", quote = F)
          write.table(cond3_enrich_table2(), enrich_table2, row.names = F, sep = "\t", quote = F)
          write.table(cond3_enrich_table3(), enrich_table3, row.names = F, sep = "\t", quote = F)
        }
        recode <- "Recode.Rdata"
        rawcount <- row_count_matrix2()
        metadata <- metadata2()
        d_rawcount <- d_row_count_matrix2()
        norm_count_matrix = norm_count_matrix2()
        dds <- MultiOut()
        fs <- c(fs,recode)
        save(rawcount,metadata,d_rawcount,dds,norm_count_matrix,file=recode)
        report <- paste0(format(Sys.time(), "%Y%m%d_"),"3conditions_report",".docx")
        fs <- c(fs,report)
        rmarkdown::render("3conditions_report.Rmd", output_format = "word_document", output_file = report,
                          params = list(raw_count = d_row_count_matrix2(),
                                        input = input,
                                        metadata2 = metadata2(),
                                        norm_count_matrix2 = norm_count_matrix2(),
                                        deg_norm_count2 = deg_norm_count2(),
                                        cond3_enrich_table1 = cond3_enrich_table1(),
                                        cond3_enrich_table2 = cond3_enrich_table2(),
                                        cond3_enrich_table3 = cond3_enrich_table3()), 
                          envir = new.env(parent = globalenv()),intermediates_dir = tempdir(),encoding="utf-8"
        )
        zip(zipfile=fname, files=fs)
      })
    },
    contentType = "application/zip"
  )
  #normalized count analysis
  #norm_count_input-------------------
  org3 <- reactive({
    if (is.null(input$Species3) || identical(input$Species3, "not selected")) {
      return(NULL)
    }
    return(org(Species = input$Species3,Ortholog = input$Ortholog3))
  })
  ortholog3 <- reactive({
    count <- norm_count_input()
    if (is.null(count) || is.null(input$Species3) || identical(input$Species3, "not selected")) {
      return(NULL)
    }
    return(no_org_ID(count = norm_count_input(),Species = input$Species3,Ortholog = input$Ortholog3,Biomart_archive=input$Biomart_archive3))
  })
  isoform3 <- reactive({
    count <- norm_count_input()
    if (is.null(count) || is.null(input$Species3) || identical(input$Species3, "not selected")) {
      return(NULL)
    }
    return(isoform_ID(count = norm_count_input(),Species = input$Species3,Ortholog = input$Ortholog3,Biomart_archive=input$Biomart_archive3,RNA_type=input$Level_norm))
  })
  gene_type3 <- reactive({
    count <- norm_count_input()
    if (is.null(count) || is.null(input$Species3) || identical(input$Species3, "not selected")) {
      return("SYMBOL")
    }
    gene_ids <- rownames(count)
    if (is.null(gene_ids) || !length(gene_ids)) {
      return("SYMBOL")
    }
    org_value <- org3()
    if (is.null(org_value)) {
      return("SYMBOL")
    }
    return(gene_type(my.symbols = gene_ids, org = org_value, Species = input$Species3, RNA_type = input$Level_norm))
  })
  org_code3 <- reactive({
    return(org_code(Species = input$Species3, Ortholog= input$Ortholog3))
  })
  observeEvent(pre_d_norm_count_matrix(),({
    updateSelectizeInput(session,inputId = "sample_order_norm","Select samples:",
                         choices = colnames(pre_d_norm_count_matrix()),selected = colnames(pre_d_norm_count_matrix()))
  }))
  observeEvent(d_norm_count_matrix(), ({
    updateCollapse(session,id =  "norm_input_collapse_panel", open="D_norm_count_matrix_panel")
  }))
  
  observeEvent(input$Species3, ({
    if(input$Species3 == "not selected"){
      updateSelectInput(session, "gene_set_forFilter","Select a gene set for gene extraction","")
    }else  if(input$Species3 != "Xenopus laevis" && input$Ortholog3 != "Arabidopsis thaliana" && input$Species3 != "Arabidopsis thaliana"){
        updateSelectInput(session, "gene_set_forFilter","Select a gene set for gene extraction",gene_set_list) 
    }else {
      updateSelectInput(session, "gene_set_forFilter","Select a gene set for gene extraction",c("KEGG", "GO biological process", 
                                                        "GO cellular component","GO molecular function")) 
    }
  }))
  norm_gene_pathway_filter_list <- reactive({
    list <-norm_gene_list()
    if (is.null(list)) {
      return(character(0))
    }
    if(!is.null(input$gene_set_forFilter)){
      if(input$Species3 == "Xenopus laevis" || input$Ortholog3 == "Arabidopsis thaliana" || input$Species3 == "Arabidopsis thaliana"){
        if(input$gene_set_forFilter == "GO biological process") list <- list %>% dplyr::filter(Ontology == "BP")
        if(input$gene_set_forFilter == "GO cellular component") list <- list %>% dplyr::filter(Ontology == "CC")
        if(input$gene_set_forFilter == "GO molecular function") list <- list %>% dplyr::filter(Ontology == "MF") 
      }}
    list <- unique(list$gs_name)
    return(list)
  })
  observeEvent(input$gene_set_forFilter, ({
    if(is.null(input$gene_set_forFilter)){
      updateSelectInput(session, "gene_set_forFilter2","","")
    }else{
      updateSelectInput(session, "gene_set_forFilter2","",norm_gene_pathway_filter_list())
    }
  }))
  
  d_norm_count_matrix <- reactive({
    count <- pre_d_norm_count_matrix()
    order <- input$sample_order_norm
    data <- try(safe_ordered_matrix(count, order), silent = TRUE)
    if(inherits(data, "try-error")) validate("")
    return(data)
  })
  norm_count_input <- reactive({
    withProgress(message = "Importing normalized count matrix, please wait",{
      if (input$data_file_type3 == "Row5"){
        tmp <- input$file7$datapath
        if(is.null(input$file7) && input$goButton3 > 0 )  tmp = example_data_path("data/example8.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example8.txt")
        return(read_df(tmp = tmp))
      }else{
        tmp <- input$file8$datapath
        if(is.null(input$file8) && input$goButton3 > 0 )  tmp = example_data_path("data/example2.csv", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example2.csv")
        return(read_df(tmp = tmp))
      }
      incProgress(1)
    })
  })
  norm_metadata <- reactive({
    if (input$data_file_type3 == "Row5"){
      return(NULL)
    }else{
      tmp <- input$file9$datapath
      if(is.null(input$file9) && input$goButton3 > 0 )  tmp <- example_data_path("data/example9.csv", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example9.csv")
      df <- read_df(tmp = tmp)
      if(!is.null(df)) rownames(df) <- gsub("-",".",rownames(df))
      return(df)
    }
  })
  norm_metadata_is_single_column <- reactive({
    meta <- norm_metadata()
    if (is.null(meta)) {
      return(FALSE)
    }
    ncol(as.data.frame(meta, stringsAsFactors = FALSE)) == 1
  })
  norm_has_matrix_data <- function(x) {
    !is.null(x) && NROW(x) > 0 && NCOL(x) > 0
  }
  norm_gene_list <- reactive({
    if (is.null(norm_count_input()) || is.null(input$Species3) || identical(input$Species3, "not selected") ||
        is.null(input$gene_set_forFilter) || identical(input$gene_set_forFilter, "")) {
      return(NULL)
    }
    genes <- GeneList_for_enrichment(Species = input$Species3, Ortholog = input$Ortholog3,
                                     Gene_set=input$gene_set_forFilter, org = org3(), 
                                     Biomart_archive=input$Biomart_archive3,gene_type=gene_type3())
    return(genes)
  })
  gene_list <- reactive({
    df <- NULL
    if(input$norm_filter == "custom"){
    data <- input$file10$datapath
    if(is.null(input$file10) && input$goButton3 == 0) return(NULL)
    if(is.null(input$file10) && input$goButton3 > 0 )  data = example_data_path("data/enrich_example.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/enrich_example.txt")
    df <- read_df(tmp = data)
    gene <-c(rownames(df))
    df <- as.data.frame(gene, row.names = gene)
    }else if(input$norm_filter == "Gene set" && input$Species3 != "not selected"){
      print("c")
      if(is.null(input$gene_set_forFilter) || is.null(input$gene_set_forFilter2) || 
         input$gene_set_forFilter == "" || input$gene_set_forFilter2 == "" ||
         !input$gene_set_forFilter2 %in% norm_gene_pathway_filter_list()) validate("")
      genes <- norm_gene_list()
      if(is.null(genes)) validate("")
      genes <- try(dplyr::filter(genes, gs_name == input$gene_set_forFilter2))
      if(inherits(genes, "try-error")) validate("")
      print("e")
      print(gene_type3())
      print(head(isoform3()))
      my.symbols <- as.character(genes$entrez_gene)
      if(gene_type3() == "non-model organism"){
        gene_IDs <-  try(dplyr::filter(ortholog3(), ENTREZID %in% my.symbols))
        if(inherits(gene_IDs, "try-error")) validate("")
        df <- data.frame(gene = gene_IDs$ENSEMBL, row.names = gene_IDs$ENSEMBL)
      }else if(gene_type3() == "isoform"){
        gene_IDs <-  try(dplyr::filter(isoform3(), ENTREZID %in% my.symbols))
        if(inherits(gene_IDs, "try-error")) validate("")
        df <- data.frame(gene = gene_IDs$Transcript_ID, row.names = gene_IDs$Transcript_ID)
      }else{
      if(gene_type3() == "SYMBOL") columns <- c("ENTREZID", sgd_symbol_column(org3())) else columns <- c("ENTREZID","ENSEMBL")
      if(org3()$packageName == "org.At.tair.db") {
        if(gene_type3() == "SYMBOL"){
          gene_IDs <- AnnotationDbi::select(org3(), keys = my.symbols,
                                            keytype = "TAIR",
                                            columns = c("TAIR","SYMBOL"))
          colnames(gene_IDs) <- c("TAIR","GeneID")
        }else{
          gene_IDs <- data.frame(TAIR = my.symbols, GeneID = my.symbols)
        }
      } else {
      gene_IDs <- AnnotationDbi::select(org3(), keys = my.symbols,
                                        keytype = "ENTREZID",
                                        columns = columns)
      print(head(gene_IDs))
      colnames(gene_IDs) <- c("entrezid","GeneID")
      }
      gene_IDs <- na.omit(gene_IDs)
      gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
      df <- data.frame(gene = gene_IDs$GeneID, row.names = gene_IDs$GeneID)
      if(dim(df)[1] == 0) validate("No filtered genes.")
      print("d")
      }
    }
    return(df)
  })
  pre_d_norm_count_matrix <- reactive({
    withProgress(message = "Creating defined count matrix, please wait", {
      
      row <- norm_count_input()
      gene_list <- gene_list()
      
      if (input$data_file_type3 == "Row5") {
        
        if (is.null(row)) {
          return(NULL)
        } else {
          print("a1")
          if (input$norm_filter != "not selected" && input$Species3 != "not selected") {
            print("b1")
            if (is.null(input$gene_set_forFilter) || is.null(input$gene_set_forFilter2) ||
                input$gene_set_forFilter == "" || input$gene_set_forFilter2 == "" ||
                !input$gene_set_forFilter2 %in% norm_gene_pathway_filter_list()) validate("")
            print("c1")
            if (!is.null(gene_list)) {
              row <- merge(row, gene_list, by = 0)
              print(head(row))
              if (dim(row)[1] == 0) validate("No filtered genes.")
              rownames(row) <- row$Row.names
              row <- row[, -1, drop = FALSE]
              row <- row[, -which(colnames(row) == "gene"), drop = FALSE]
            }
          }
          return(anno_rep(row))
        }
        
      } else {
        
        meta_raw <- norm_metadata()
        
        if (is.null(row) || is.null(meta_raw)) {
          return(NULL)
        } else {
          
          row_t <- t(row)
          meta_raw <- as.data.frame(meta_raw, stringsAsFactors = FALSE)
          
          if (ncol(meta_raw) == 1) {
            # 今まで通りの処理
            meta <- anno_rep_meta(meta_raw)
            meta <- as.data.frame(meta, stringsAsFactors = FALSE)
            meta <- data.frame(
              characteristics = meta[[1]],
              row.names = rownames(meta),
              stringsAsFactors = FALSE
            )
            
            colname <- colnames(meta)
            data <- merge(meta, row_t, by = 0, sort = FALSE)
            
            if (dim(data)[1] == 0) {
              rownames(meta) <- gsub("\\.", "-", rownames(meta))
              data <- merge(meta, row_t, by = 0, sort = FALSE)
              if (dim(data)[1] == 0) {
                rownames(row_t) <- gsub("\\.", "-", rownames(row_t))
                data <- merge(meta, row_t, by = 0, sort = FALSE)
                validate("Error: failed to merge count data with metadata. Please check row names of metadata.")
              }
            }
            
            rownames(data) <- data$characteristics
            data2 <- data[, -which(colnames(data) %in% c("Row.names", colname)), drop = FALSE]
            
          } else {
            # 2列以上なら multi_d_row_count_matrix と同様の処理
            meta <- meta_raw
            colname <- colnames(meta)
            data <- merge(meta, row_t, by = 0, sort = FALSE)
            
            if (dim(data)[1] == 0) {
              rownames(meta) <- gsub("\\.", "-", rownames(meta))
              data <- merge(meta, row_t, by = 0, sort = FALSE)
              if (dim(data)[1] == 0) {
                rownames(row_t) <- gsub("\\.", "-", rownames(row_t))
                data <- merge(meta, row_t, by = 0, sort = FALSE)
                validate("Error: failed to merge count data with metadata. Please check row names of metadata.")
              }
            }
            
            rownames(data) <- data[, 1]
            data2 <- data[, -which(colnames(data) %in% c("Row.names", colname)), drop = FALSE]
            data2 <- data2[, 1:length(rownames(row)), drop = FALSE]
          }
          
          data2_t <- t(data2)
          data3 <- apply(data2_t, 2, as.numeric)
          data3 <- as.matrix(data3)
          rownames(data3) <- rownames(data2_t)
          if (!is.null(colnames(data2_t))) {
            colnames(data3) <- colnames(data2_t)
          }
          
          row <- data3
          
          if (input$norm_filter != "not selected" && input$Species3 != "not selected") {
            if (!is.null(gene_list)) {
              row <- merge(row, gene_list, by = 0)
              if (dim(row)[1] == 0) validate("No filtered genes.")
              rownames(row) <- row$Row.names
              row <- row[, -1, drop = FALSE]
              row <- row[, -which(colnames(row) == "gene"), drop = FALSE]
            }
          }
          
          return(row)
        }
      }
      
      incProgress(1)
    })
  })
  
  d_norm_count_matrix_cutofff <- reactive({
    data <- d_norm_count_matrix()
    if (is.null(data)) {
      return(NULL)
    }
    
    withProgress(message = "Creating defined count matrix, please wait", {
      if (input$basemean3 != 0) {
        keep <- rowMeans(data, na.rm = TRUE) > input$basemean3
        data <- data[keep, , drop = FALSE]
      }
      incProgress(1)
      data
    })
  })
  
  gene_ID_norm <- reactive({
      if(input$Species3 != "not selected"){
        gene_type_value <- gene_type3()
        if(!identical(gene_type_value, "SYMBOL")){
          res <- d_norm_count_matrix()
          if(!is.null(res)){
          gene_IDs <- ensembl2symbol(gene_type=gene_type_value,data = res,Species=input$Species3,
                                     Ortholog=ortholog3(),Isoform=isoform3(),org = org3(),merge=FALSE)
          return(gene_IDs)
          }
        }
      }
      return(NULL)
  })
  
  observeEvent(input$goButton3,({
    updateSelectInput(session,inputId = "Species3","Species",species_list, selected = "Homo sapiens")
  }))
  
  observeEvent(input$file7, ({
    updateCollapse(session,id =  "norm_input_collapse_panel", open="Norm_count_panel")
  }))
  observeEvent(input$file8, ({
    updateCollapse(session,id =  "norm_input_collapse_panel", open="Norm_count_panel")
  }))
  observeEvent(input$file10, ({
    updateCollapse(session,id =  "norm_input_collapse_panel", open="gene_list_panel")
  }))
  observeEvent(input$file9, ({
    updateCollapse(session,id =  "norm_input_collapse_panel", open="norm_Metadata_panel")
  }))
  output$norm_count_input1 <- DT::renderDataTable({
    norm_count_input()
  })
  output$norm_Metadata <- DT::renderDataTable({
    norm_metadata()
  })
  output$d_norm_count <- DT::renderDataTable({
    d_norm_count_matrix()
  })
  output$Gene_list <- DT::renderDataTable({
    gene_list()
  })
  download_name_stem <- function(file_input, fallback) {
    if (is.null(file_input)) {
      return(fallback)
    }
    if (is.data.frame(file_input) || is.list(file_input)) {
      file_name <- file_input$name
    } else {
      file_name <- file_input
    }
    file_name <- file_name[!is.na(file_name) & nzchar(file_name)]
    if (!length(file_name)) {
      return(fallback)
    }
    tools::file_path_sans_ext(basename(file_name[[1]]))
  }
  output$download_d_norm_count = downloadHandler(
    filename = function() {
      if (identical(input$data_file_type3, "Row5")) {
        paste0(download_name_stem(input$file7, "normalized_count"), ".txt")
      } else {
        paste0(
          download_name_stem(input$file8, "normalized_count"),
          "-",
          download_name_stem(input$file9, "metadata"),
          ".txt"
        )
      }
    },
    content = function(file){write.table(d_norm_count_matrix(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  
  download_norm_dir <-reactive({
    if (identical(input$data_file_type3, "Row5")) {
      dir_name <- paste0(download_name_stem(input$file7, "normalized_count"), "_")
    } else {
      dir_name <- paste0(
        download_name_stem(input$file8, "normalized_count"),
        "-",
        download_name_stem(input$file9, "metadata"),
        "_",
        input$basemean3,
        "_"
      )
    }
    return(dir_name)
  })
  
  #norm clustering--------------------
  output$download_norm_PCA = downloadHandler(
    filename = function(){
      paste0(download_norm_dir(), "PCA-MDS-dendrogram.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$norm_pdf_height == 0){
          pdf_height <- 3.5
        }else pdf_height <- input$norm_pdf_height
        if(input$norm_pdf_width == 0){
          pdf_width <- 9
        }else pdf_width <- input$norm_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(PCAplot(data = d_norm_count_matrix_cutofff(),legend = input$PCA_legend_norm))
        dev.off()
      })
    }
  )
  
  output$norm_PCA <- renderPlot({
    withProgress(message = "Clustering",{
      if(is.null(d_norm_count_matrix_cutofff())){
        return(NULL)
      }else{
        print(PCAplot(data = d_norm_count_matrix_cutofff(),legend = input$PCA_legend_norm))
      }
    })
  })
  
  output$norm_PCA_table <- DT::renderDataTable({
    PCAdata(row_count = d_norm_count_matrix_cutofff(),deg_norm_count = d_norm_count_matrix_cutofff())
  })
  output$d_norm_count_cutoff <- DT::renderDataTable({
    d_norm_count_cutoff_uniqueID()
  })
  
  output$download_norm_pca_table = downloadHandler(
    filename = function() {
      paste0(download_norm_dir(), "PCA_table.txt")
    },
    content = function(file){write.table(PCAdata(row_count = d_norm_count_matrix_cutofff(), deg_norm_count = d_norm_count_matrix_cutofff()), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  
  output$norm_umap_n <- renderUI({
    info <- umap_neighbor_info(d_norm_count_matrix_cutofff())
    if(!info$valid){
      return(NULL)
    }
    sliderInput("norm_n_neighbors", "n_neighbors", min = 2,
                max = info$max, step = 1,
                value = info$value)
  })
  
  norm_umap_plot <- reactive({
    data <- d_norm_count_matrix_cutofff()
    if(is.null(input$norm_n_neighbors)){
      return(NULL)
    }else{
      if(is.null(data)){
        return(NULL)
      }else{
        p<- try(umap_plot(data = data, n_neighbors = input$norm_n_neighbors,lab=input$norm_umap_label), silent = TRUE)
        return(p)
      }
    }
  })
  
  output$norm_umap_error <- renderText({
    data <- d_norm_count_matrix_cutofff()
    err <- umap_error_message(data, input$norm_n_neighbors)
    if(!is.null(err)){
      return(err)
    }
    p <- norm_umap_plot()
    if(inherits(p, "try-error")){
      return(conditionMessage(attr(p, "condition")))
    }
    NULL
  })
  
  
  output$norm_umap <- renderPlot({
    p <- norm_umap_plot()
    withProgress(message = "umap",{
      if(inherits(p, "try-error") || is.null(p)) return(NULL)
      print(p)
    })
  })
  
  output$download_norm_umap = downloadHandler(
    filename = function(){
      paste0(download_norm_dir(), paste0(input$norm_n_neighbors,"_neighbors_umap.pdf"))
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$norm_pdf_height == 0){
          pdf_height <- 3.5
        }else pdf_height <- input$norm_pdf_height
        if(input$norm_pdf_width == 0){
          pdf_width <- 4.7
        }else pdf_width <- input$norm_pdf_width
        plot_obj <- norm_umap_plot()
        if(inherits(plot_obj, "try-error")){
          stop(conditionMessage(attr(plot_obj, "condition")), call. = FALSE)
        }
        if(is.null(plot_obj)){
          stop("umap: plot is not available", call. = FALSE)
        }
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_obj)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  #norm GOI------------------------------------------------------
  output$filtered_regionGOI <- renderText({
    if(is.null(GOI_list3())){
      return(NULL)
    }else{ 
      if(length(input$selectFC_normGOI) == 2){
        print(paste0("The number of genes after the filtration (basemean > ", input$basemean3,", |log2(", 
                     input$selectFC_normGOI[1],"/", input$selectFC_normGOI[2],")| > ", log2(input$fc3),"): ", 
                     length(GOI_list3())))
      }else{
        print(paste0("The number of genes after the filtration (basemean > ", input$basemean3,"): ", 
                     length(GOI_list3())))
      }
    }
  })
  output$selectFC_normGOI <- renderUI({
    if(is.null(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      selectizeInput("selectFC_normGOI", "Option: select a pair for fold change cut-off", 
                     choices =c(unique(unique(gsub("\\_.*","", colnames(d_norm_count_matrix_cutofff()))))),
                     selected = "", options = list(maxItems = 2)
                     )
    }
  })
  
  d_norm_count_cutoff_uniqueID <- reactive({
    count <- d_norm_count_matrix_cutofff()
    if(!norm_has_matrix_data(count)){
      return(NULL)
    }
    count <- as.data.frame(count)
    count <- dplyr::filter(count, apply(count,1,mean) > input$basemean3)
    if(!norm_has_matrix_data(count)){
      return(NULL)
    }
    if(!is.null(gene_ID_norm())){
        gene_IDs  <- gene_ID_norm()
        data2 <- merge(count, gene_IDs, by= 0)
        rownames(data2) <- data2[,1]
        data2 <- data2[, - which(colnames(data2) == "Row.names.y")]
        data2$Unique_ID <- paste(data2$SYMBOL,data2$Row.names, sep = "\n- ")
        count <- data2[,-1]
    }
    return(count)
  })
  d_norm_count_uniqueID <- reactive({
    count <- d_norm_count_matrix()
    if(!norm_has_matrix_data(count)){
      return(NULL)
    }
    count <- as.data.frame(count)
    if(!is.null(gene_ID_norm())){
        gene_IDs  <- gene_ID_norm()
        data2 <- merge(count, gene_IDs, by= 0)
        rownames(data2) <- data2[,1]
        data2 <- data2[, - which(colnames(data2) == "Row.names.y")]
        data2$Unique_ID <- paste(data2$SYMBOL,data2$Row.names, sep = "\n- ")
        count <- data2[,-1]
    }
    return(count)
  })
  GOI_list3 <- reactive({
    count <- preGOI_list3()
    if(is.null(count)){
      return(NULL)
    }else{
      if(dim(count)[1] == 0){
        validate(paste0("The number of genes after the filtration: 0",".\n",
                        "Please adjust cut-off conditions."))
      }
      if(gene_type3() != "SYMBOL"){
        if(input$Species3 != "not selected"){
          GOI <- count$Unique_ID
        }else GOI <- rownames(count)
      }else{
        if(input$Species3 != "not selected"){
          GOI <- rownames(count)
        }else GOI <- rownames(count)
      }
      return(normalize_goi_choices(GOI))
    }
  })
  preGOI_list3 <- reactive({
    if(input$normGOI_filter_on == "ON"){
      print(input$selectFC_normGOI)
      if(length(input$selectFC_normGOI) == 2){
        data <- d_norm_count_cutoff_uniqueID()
        if(dim(data)[1] != 0){
          cond1 <- input$selectFC_normGOI[1]
          cond2 <- input$selectFC_normGOI[2]
          cond1_num <- data %>% dplyr::select(.,starts_with(cond1)) %>% colnames() %>% length()
          cond2_num <- data %>% dplyr::select(.,starts_with(cond2)) %>% colnames() %>% length()
          cond1_ave <- data %>% dplyr::select(starts_with(cond1)) %>% rowSums(na.rm=TRUE)/cond1_num
          cond2_ave <- data %>% dplyr::select(starts_with(cond2)) %>% rowSums(na.rm=TRUE)/cond2_num
          Log2FoldChange <- log((cond1_ave + 0.01)/(cond2_ave + 0.01),2)
          data$Log2FoldChange <- Log2FoldChange
          data2 <- data %>% dplyr::filter(abs(Log2FoldChange) > log2(input$fc3))
          data2 <- data2[, - which(colnames(data2) == "Log2FoldChange")]
        }else data2 <- NULL
      }else data2 <- d_norm_count_uniqueID()
    }else data2 <- d_norm_count_uniqueID()
    return(data2)
  })
  
  observeEvent(GOI_list3(), {
    choices <- GOI_list3()
    if(is.null(choices)){
      choices <- character(0)
    }
    selected <- intersect(isolate(input$GOI3), choices)
    updateSelectizeInput(session, "GOI3", choices = choices, selected = selected,
                         options = goi_selectize_options, server = TRUE)
  }, ignoreNULL = FALSE)

  observeEvent(input$GOIreset_norm, {
    choices <- GOI_list3()
    if(is.null(choices)){
      choices <- character(0)
    }
    updateSelectizeInput(session, "GOI3", choices = choices,
                         selected = character(0),
                         options = goi_selectize_options, server = TRUE)
  })
  output$GOIreset_norm <- renderUI({
    actionButton("GOIreset_norm", "GOI reset")
  })
  output$norm_uniqueID_cut <- renderUI({
    if(input$Species3 != "not selected" && !identical(gene_type3(), "SYMBOL")) 
      radioButtons("norm_uniqueID_cut","Short unique ID",c("ON"=TRUE,"OFF"=FALSE),selected = TRUE)
  })
  GOI3_INPUT <- reactive({
    if(is.null(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      if(input$GOI3_type == "ALL") return(GOI_list3())
      if(input$GOI3_type == "custom") return(input$GOI3)
    }
  })
  norm_GOIcount <- reactive({
    count <- d_norm_count_uniqueID() 
    if(gene_type3() != "SYMBOL"){
      if(input$Species3 != "not selected"){
        Unique_ID <- GOI3_INPUT()
        label_data <- as.data.frame(Unique_ID, row.names = Unique_ID)
        data <- merge(count, label_data, by="Unique_ID")
        if(input$norm_uniqueID_cut) {
          id_list <- gsub("\\\n.+$", "", data$Unique_ID)
          dup_list <- unique(id_list[duplicated(id_list)])
          for(i in 1:length(data$Unique_ID)){
            if(gsub("\\\n.+$", "", data$Unique_ID[i]) == "NA") {
              data$Unique_ID[i] <- gsub(".+\\s", "", data$Unique_ID[i])
            }else if(! gsub("\\\n.+$", "", data$Unique_ID[i]) %in% dup_list) {
              data$Unique_ID[i] <- gsub("\\\n.+$", "", data$Unique_ID[i])
            }
          }
        }
        rownames(data) <- data$Unique_ID
        data <- data[, - which(colnames(data) == "SYMBOL")]
        data <- data[, - which(colnames(data) == "Unique_ID")]
      }else{
        Row.names <- GOI3_INPUT()
        count$Row.names <- rownames(count)
        label_data <- as.data.frame(Row.names, row.names = Row.names)
        data <- merge(count, label_data, by="Row.names")
        rownames(data) <- data$Row.names
        data <- data[, - which(colnames(data) == "Row.names")]
      }
    }else{
      Row.names <- GOI3_INPUT()
      count$Row.names <- rownames(count)
      label_data <- as.data.frame(Row.names, row.names = Row.names)
      data <- merge(count, label_data, by="Row.names")
      rownames(data) <- data$Row.names
      data <- data[, - which(colnames(data) == "Row.names")]
    }
    return(data)
  })
  
  norm_GOIheat <- reactive({
    data <- norm_GOIcount()
    if(is.null(data)){
      ht <- NULL
    }else{
      data.z <- genefilter::genescale(data, axis=1, method="Z")
      data.z <- na.omit(data.z)
      ht <- GOIheatmap(data.z, type = input$GOI3_type, GOI = input$GOI3)
    }
    return(ht)
  })
  
  output$norm_GOIheatmap <- renderPlot({
    if(is.null(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      if(!is.null(GOI3_INPUT())){
        withProgress(message = "heatmap",{
          suppressWarnings(print(norm_GOIheat()))
          incProgress(1)
        })
      }
    }
  })
  output$statistics <- renderUI({
    if(!is.null(d_norm_count_matrix_cutofff()) && !is.null(GOI3_INPUT())){
      data <- norm_GOIcount()
      if(!is.null(data)){
        collist <- gsub("\\_.+$", "", colnames(data))
        collist <- unique(collist)
        if(length(collist) == 2){
          selectInput("statistics","statistics",choices = c("not_selected","Welch's t-test","Wilcoxon test"),selected="not_selected",multiple = F)
        }else{
          selectInput("statistics","statistics",choices = c("not_selected","TukeyHSD","Dunnet's test","Dunn's test"),selected="not_selected", multiple = F)
        }
      }}
  })
  output$Color_rev_norm <- renderUI({
    if(input$PlotType != "Errorplot"){
    if(is.null(input$Color_norm)) validate("")
    if(input$Color_norm != "default") {
    radioButtons('Color_rev_norm', 'Color reverse', 
                 choiceNames = c("OFF","ON"),choiceValues = c("OFF","ON"),
                 selected = "OFF")
    }}
  })
  output$Color_norm <- renderUI({
    if(input$PlotType != "Errorplot"){
    radioButtons('Color_norm', 'Color palette', 
                 choiceNames = GOI_color_palette, 
                 choiceValues = GOI_color_palette,
                 selected = "default")
    }
  })
  output$Color_design_norm <- renderUI({
    if(input$PlotType != "Errorplot"){
    radioButtons('Color_design_norm', 'graph design', 
                 choiceNames = c("new","<=v1.1.0"), 
                 choiceValues = c("new","<=v1.1.0"),
                 selected = "new")
    }
  })
  statistical_analysis_goi <- reactive({
    data <- norm_GOIcount()
    if(is.null(data) || is.null(input$statistics)){
      p <- NULL
    }else{
      p <- GOIboxplot(data = data,statistical_test =input$statistics,color_design=input$Color_design_norm,ylab=input$norm_ylab,
                      plottype=input$PlotType,color = input$Color_norm,rev=input$Color_rev_norm,ymin=input$norm_ymin)
    }
    return(p)
  })
  
  norm_GOIbox <- reactive({
    if(!is.null(statistical_analysis_goi())){
      if(input$statistics == "not_selected"){
        return(statistical_analysis_goi())
      }else return(statistical_analysis_goi()[["plot"]])
    }
  })
  
  
  output$norm_GOIboxplot <- renderPlot({
    if(is.null(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      if(!is.null(GOI3_INPUT())){
        if(length(rownames(norm_GOIcount())) >200){
          validate("Unable to display more than 200 genes. Please adjust the threshold to narrow down the number of genes to less than 200, or utilize the 'Custom' mode.")
        }
        withProgress(message = "Boxplot",{
          suppressWarnings(print(norm_GOIbox()))
          incProgress(1)
        })
      }
    }
  })
  
  norm_GOIbox_statistic <- reactive({
    if(!is.null(statistical_analysis_goi())){
      if(input$statistics != "not_selected"){
        data <- as.data.frame(statistical_analysis_goi()[["statistical_test"]])
        colnames(data)[1] <- "gene"
        if(input$statistics == "TukeyHSD" || input$statistics == "Dunnet's test"){
          data <- data[, - which(colnames(data) == "term")]
        }else{
          data <- data[, - which(colnames(data) == ".y.")]
          data <- data[, - which(colnames(data) == "n1")]
          data <- data[, - which(colnames(data) == "n2")]
        }
        data <- data[, - which(colnames(data) == "y.position")]
        data <- data[, - which(colnames(data) == "groups")]
        data <- data[, - which(colnames(data) == "xmin")]
        data <- data[, - which(colnames(data) == "xmax")]
        return(data)
      }else return(NULL)}
  })
  
  output$statistical_table <- DT::renderDataTable({
    if(is.null(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      if(!is.null(GOI3_INPUT())){
        if(length(rownames(norm_GOIcount())) >200){
          validate("Cannot display more than 200 genes.")
        }
        norm_GOIbox_statistic()
      }}
  })
  
  output$download_statisics = downloadHandler(
    filename = function(){
      paste0(download_norm_dir(), "GOIboxplot_",input$statistics,".txt")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        write.table(norm_GOIbox_statistic(),file, row.names = F, col.names=TRUE, sep = "\t", quote = F)
        incProgress(1)
      })
    }
  )
  
  output$download_norm_GOIbox = downloadHandler(
    filename = function(){
      paste0(download_norm_dir(), "GOIboxplot.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- norm_GOIcount()
        rowlist <- rownames(data)
        if(input$norm_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)
        }else pdf_height <- input$norm_pdf_height
        if(input$norm_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)
        }else pdf_width <- input$norm_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(norm_GOIbox())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_norm_GOIheat = downloadHandler(
    filename = function(){
      paste0(download_norm_dir(), "GOIheatmap.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- norm_GOIcount()
        rowlist <- rownames(data)
        if(input$norm_pdf_height == 0){
          pdf_height <- 10
        }else pdf_height <- input$norm_pdf_height
        if(input$norm_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$norm_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(norm_GOIheat())
        dev.off()
        incProgress(1)
      })
    }
  )
  #norm corrplot---------------------
  updateCounter_corr <- reactiveValues(i = 0)
  norm_corr_screen_cache <- reactiveValues(
    key = NULL,
    value = NULL
  )
  
  observe({
    input$corr_start
    isolate({
      updateCounter_corr$i <- updateCounter_corr$i + 1
    })
  })
  
  
  #Restart
  observeEvent(input$GOI_x, {
    isolate(updateCounter_corr$i == 0)
    updateCounter_corr$i <- 0
  }) 
  observeEvent(input$corr_statistics, {
    isolate(updateCounter_corr$i == 0)
    updateCounter_corr$i <- 0
  }) 
  GOI_list3_corr <- reactive({
    count <- d_norm_count_cutoff_uniqueID()
    if(is.null(count)){
      return(NULL)
    }else{
      if(gene_type3() != "SYMBOL"){
        if(input$Species3 != "not selected"){
          GOI <- count$Unique_ID
        }else GOI <- rownames(count)
      }else{
        if(input$Species3 != "not selected"){
          GOI <- rownames(count)
        }else GOI <- rownames(count)
      }
      return(normalize_goi_choices(GOI))
    }
  })
  observeEvent(GOI_list3_corr(), {
    goi_choices <- GOI_list3_corr()
    if(is.null(goi_choices)){
      goi_choices <- character(0)
    }
    goi_choices <- c("", goi_choices)
    selected_x <- isolate(input$GOI_x)
    if(is.null(selected_x) || !selected_x %in% goi_choices){
      selected_x <- ""
    }
    selected_y <- isolate(input$GOI_y)
    if(is.null(selected_y) || !selected_y %in% goi_choices){
      selected_y <- ""
    }
    color_choices <- unique(c("", "sample_name", goi_choices))
    selected_color <- isolate(input$corr_color)
    if(is.null(selected_color) || !selected_color %in% color_choices){
      selected_color <- ""
    }
    selected_color_selected <- isolate(input$corr_color_selected)
    if(is.null(selected_color_selected) || !selected_color_selected %in% color_choices){
      selected_color_selected <- ""
    }
    updateSelectizeInput(session, "GOI_x", choices = goi_choices, selected = selected_x, server = TRUE)
    updateSelectizeInput(session, "GOI_y", choices = goi_choices, selected = selected_y, server = TRUE)
    updateSelectizeInput(session, "corr_color", choices = color_choices, selected = selected_color, server = TRUE)
    updateSelectizeInput(session, "corr_color_selected", choices = color_choices, selected = selected_color_selected, server = TRUE)
  }, ignoreNULL = FALSE)
  output$norm_corr_cutoff <- renderUI({
    if(!is.null(pre_pre_norm_GOI_corrplot())){
    radioButtons("norm_corr_cutoff", "cut-off (padj or pval)", c('padj'="padj",'pval'="pval"), selected = "padj")
    }
  })
  norm_GOIcount_corr <- reactive({
    if(!is.null(input$GOI_x)){
      if(input$GOI_x != ""){
        data <- d_norm_count_cutoff_uniqueID()
        if(gene_type3() != "SYMBOL"){
          if(input$Species3 != "not selected"){
            rownames(data) <- data$Unique_ID
            data <- data[, - which(colnames(data) == "SYMBOL")]
            data <- data[, - which(colnames(data) == "Unique_ID")]
          }
        }
    return(data)
      }
    }
  })
  
  norm_GOI_corrplot <- reactive({
    data <- norm_GOIcount_corr()
    if(is.null(data) || is.null(input$GOI_x)){
      p <- NULL
    }else{
      if(input$corr_mode == "corr_mode2"){
      }else{
        if(!is.null(pre_norm_GOI_corrplot())){
          df2 <- pre_norm_GOI_corrplot()
            if(dim(brush_info_corr())[1]!=0){
              label_data <- brush_info_corr()$prey
            }else if(!is.null(input$statistical_table_corrplot_rows_selected)){
              label_data <- df2[input$statistical_table_corrplot_rows_selected,]$prey
            }else label_data <- NULL
          if(!is.null(label_data)) {
            Color <- c("blue","green","darkgray","red")
            if(length(df2$prey[df2$color == "negative_correlation"]) == 0) Color <- c("green","darkgray","red")
            for(name in label_data){
              df2$color[df2$prey == name] <- "GOI"
            }
            df2$color <- factor(df2$color, levels = c("negative_correlation","GOI","NS", "positive_correlation"))
          }else{
            Color <- c("blue","darkgray","red")
            df2$color <- factor(df2$color, levels = c("negative_correlation","NS", "positive_correlation"))
            if(length(df2$prey[df2$color == "negative_correlation"]) == 0) Color <- c("darkgray","red")
          }
          if(input$corr_statistics == "spearman") ylab <- paste0("Spearman's correlation") else ylab <- paste0("Pearson's correlation")
          p <- ggplot(df2,aes(x=rank,y=corr_score,col=color))+
            ggrastr::geom_point_rast(aes(color = color),size = 0.4)  +
            theme_bw()+ scale_color_manual(values = Color)+
            theme(legend.position = "top" , legend.title = element_blank(),
                  axis.text.x= ggplot2::element_text(size = 12),
                  axis.text.y= ggplot2::element_text(size = 12),
                  text = ggplot2::element_text(size = 12),
                  title = ggplot2::element_text(size = 12))  + ylab(ylab)
          if(!is.null(label_data)) {
            p <- p + geom_point(data=dplyr::filter(df2, color == "GOI"),color="green", size=1)
            p <- p + ggrepel::geom_text_repel(data = dplyr::filter(df2, color == "GOI"), mapping = aes(label = prey),
                                               color = "black",segment.color = "black",
                                               box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"),
                                               force = 1, fontface = "bold.italic",
                                               bg.color = grDevices::adjustcolor("white", alpha.f = 0.6), bg.r = 0.15)+ 
              guides(
                fill = guide_legend(
                  override.aes = aes(label = "")
                )
              )
          }
        }
      }
    }
    return(p)
  })
  norm_GOI_corrplot_pair <- reactive({
    data <- norm_GOIcount_corr()
    if(is.null(data) || is.null(input$GOI_x) || is.null(input$corr_color) || is.null(corr_table_pair())){
      p <- NULL
    }else{
      if(input$corr_mode == "corr_mode2"){
        data <- as.data.frame(t(data))
        p <- corr_plot_pair(data = data, corr_color = input$corr_color,GOI_x = input$GOI_x,GOI_y = input$GOI_y)
        if(corr_table_pair()$pvalue < 0.0001) pvalue <- "p < 0.0001" else pvalue <- paste0("p = ", round(corr_table_pair()$pvalue, digits = 5))
        p <- p +  labs(caption = paste0("r = ", round(corr_table_pair()$corr_score, digits = 4), ", ", pvalue))
      }
    }
    return(p)
  })
  norm_corr_screen_key <- reactive({
    data <- norm_GOIcount_corr()
    if(is.null(data) || is.null(input$GOI_x) || input$GOI_x == ""){
      return(NULL)
    }
    digest::digest(list(
      GOI_x = input$GOI_x,
      corr_statistics = input$corr_statistics,
      data = data
    ))
  })
  pre_pre_norm_GOI_corrplot <- reactive({
    data <- norm_GOIcount_corr()
    if(is.null(data) || input$corr_start == 0  || updateCounter_corr$i == 0){
      df2 <- NULL
    }else{
      cache_key <- norm_corr_screen_key()
      if(!is.null(cache_key) &&
         identical(isolate(norm_corr_screen_cache$key), cache_key) &&
         !is.null(isolate(norm_corr_screen_cache$value))){
        return(isolate(norm_corr_screen_cache$value))
      }
      withProgress(message = "Correlation analysis takes a few minutes",{
          row_ids <- rownames(data)
          bait <- input$GOI_x
          bait_values <- as.numeric(data[bait, , drop = TRUE])
          result_mat <- matrix(NA_real_, nrow = length(row_ids), ncol = 3)
          progress_step <- max(1L, floor(length(row_ids) / 40L))
          last_progress <- 0
          for(idx in seq_along(row_ids)){
            corr <- suppressWarnings(cor.test(
              x = bait_values,
              y = as.numeric(data[row_ids[idx], , drop = TRUE]),
              method = "spearman"
            ))
            result_mat[idx, 1] <- unname(corr$statistic)
            result_mat[idx, 2] <- unname(corr$estimate)
            result_mat[idx, 3] <- corr$p.value
            if(idx %% progress_step == 0L || idx == length(row_ids)){
              current_progress <- idx / length(row_ids)
              incProgress(current_progress - last_progress)
              last_progress <- current_progress
            }
          }
          df2 <- data.frame(
            prey = row_ids,
            bait = bait,
            statistics = result_mat[, 1],
            corr_score = result_mat[, 2],
            pvalue = result_mat[, 3],
            method = input$corr_statistics,
            stringsAsFactors = FALSE
          )
          padj <- p.adjust(df2$pvalue,method="BH")
          df2$padj <- padj
          df2 <- na.omit(df2)
          df2 <- df2%>% dplyr::arrange(-corr_score, padj) %>%
            dplyr::mutate(rank = row_number())
          rownames(df2) <- df2$rank
          if(!is.null(cache_key)){
            norm_corr_screen_cache$key <- cache_key
            norm_corr_screen_cache$value <- df2
          }
      })
    }
    return(df2)
  })
  output$corr_fdr <- renderUI({
    if(!is.null(pre_pre_norm_GOI_corrplot())){
      numericInput("corr_fdr","padj (or pval) cut-off",max = 1, value = 0.05,step=0.001)
    }
  })

  pre_norm_GOI_corrplot   <- reactive({
    df2 <- pre_pre_norm_GOI_corrplot()
    if(!is.null(df2) && !is.null(input$norm_corr_cutoff)){
    df2$color <- "NS"
    if(input$norm_corr_cutoff == "padj"){
      df2$color[df2$padj < input$corr_fdr & df2$corr_score > 0] <- "positive_correlation"
      df2$color[df2$padj < input$corr_fdr & df2$corr_score < 0] <- "negative_correlation"
    }else{
      df2$color[df2$pvalue < input$corr_fdr & df2$corr_score > 0] <- "positive_correlation"
      df2$color[df2$pvalue < input$corr_fdr & df2$corr_score < 0] <- "negative_correlation"
    }
    df2 <- df2 %>% dplyr::select(prey,color,everything())
    return(df2)
    }
  })
  norm_GOI_corrplot_selected <- reactive({
    data <- norm_GOIcount_corr()
    if(is.null(data) || is.null(input$GOI_x) || is.null(input$norm_corr_selected_list) || 
       is.null(corr_table()) || is.null(input$corr_color_selected)){
      p <- NULL
    }else{
      if(input$corr_mode == "corr_mode1"){
        GOI_y <- input$norm_corr_selected_list
        if(GOI_y != ""){
        data <- as.data.frame(t(data))
        p <- corr_plot_pair(data = data, corr_color = input$corr_color_selected,GOI_x = input$GOI_x,GOI_y = GOI_y)
        if(corr_table()$pvalue[corr_table()$prey == GOI_y] < 0.0001) pvalue <- "p < 0.0001" else pvalue <- paste0("p = ", round(corr_table()$pvalue[corr_table()$prey == GOI_y], digits = 5))
        if(corr_table()$padj[corr_table()$prey == GOI_y] < 0.0001) padj <- "padj < 0.0001" else padj <- paste0("padj = ", round(corr_table()$padj[corr_table()$prey == GOI_y], digits = 5))
        p <- p +  labs(caption = paste0("r = ", round(corr_table()$corr_score[corr_table()$prey == GOI_y], digits = 4), ", ", pvalue,", ",padj))
        }
      }
    }
    return(p)
  })
  norm_GOI_corrplot_selected_for_download <- reactive({
    data <- norm_GOIcount_corr()
    if(is.null(data) || is.null(input$GOI_x) || is.null(input$norm_corr_selected_list)){
      df <- NULL
    }else{
      if(input$corr_mode == "corr_mode1"){
        df2 <- pre_norm_GOI_corrplot()
        if(dim(brush_info_corr())[1]!=0){
          GOI_y <- brush_info_corr()$prey
        }else GOI_y <- df2[input$statistical_table_corrplot_rows_selected,]$prey
        data <- as.data.frame(t(data))
        df <- list()
        processNum <- length(GOI_y)
        withProgress(message = paste0("Prepare correlation plots of all selected genes (",length(GOI_y),")"),{
        for(y in GOI_y){
          p <- corr_plot_pair(data = data, corr_color = input$corr_color_selected,GOI_x = input$GOI_x,GOI_y = y)
          if(corr_table()$pvalue[corr_table()$prey == y] < 0.0001) pvalue <- "p < 0.0001" else pvalue <- paste0("p = ", round(corr_table()$pvalue[corr_table()$prey == y], digits = 5))
          if(corr_table()$padj[corr_table()$prey == y] < 0.0001) padj <- "padj < 0.0001" else padj <- paste0("padj = ", round(corr_table()$padj[corr_table()$prey == y], digits = 5))
          p <- p +  labs(caption = paste0("r = ", round(corr_table()$corr_score[corr_table()$prey == y], digits = 4), ", ", pvalue,", ",padj))
          df[[y]] <- p
          incProgress(1/processNum, message = paste0("Prepare ",y))
        }
        })
      }
    }
    return(df)
  })
  
  brush_info_corr <- reactive({
    data <- as.data.frame(pre_norm_GOI_corrplot())
    if(dim(data)[1] != 0) return(brushedPoints(data, input$plot1_brush_corr,xvar = "rank",yvar="corr_score"))
  })
  
  output$norm_corrplot <- renderPlot({
    if(is.null(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      if(input$corr_mode == "corr_mode2"){
        if(!is.null(norm_GOI_corrplot_pair()) && !is.null(input$GOI_x) && !is.null(input$GOI_y) &&
           input$GOI_x != "" && input$GOI_y != ""){
          norm_GOI_corrplot_pair()
        }
      }else{
        if(!is.null(norm_GOI_corrplot()) && !is.null(input$GOI_x) && input$GOI_x != ""){
          norm_GOI_corrplot()
        }
        
      }
    }
  })
  output$norm_corrplot_selected <- renderPlot({
    if(is.null(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      if(input$corr_mode == "corr_mode1"){
        if(!is.null(norm_GOI_corrplot_selected()) && 
           !is.null(input$GOI_x) && input$GOI_x != ""){
          norm_GOI_corrplot_selected()
        }
      }
    }
  })
  corr_table <- reactive({
    data <- norm_GOIcount_corr()
    if(input$corr_mode == "corr_mode2"){
    }else{
      if(!is.null(pre_norm_GOI_corrplot())){
        df <- pre_norm_GOI_corrplot()
        return(df)
      }
    }
  })
  corr_table_forOutput <- reactive({
    df <- corr_table()
        if(gene_type3() != "SYMBOL"){
          if(input$Species3 != "not selected"){
            colnames(df)[1] <- "prey_UniqueID"
            df$prey <- gsub(".+\\-","",df$prey_UniqueID)
            df$prey <- gsub(" ","",df$prey)
            df <- df %>% dplyr::select(prey, color, everything())
            df$prey_UniqueID <- gsub("\n"," ",df$prey_UniqueID)
            df$bait <- gsub("\n"," ",df$bait)
          }
        }
        return(df)
  })
  corr_table_pair <- reactive({
    data <- norm_GOIcount_corr()
    if(input$corr_mode == "corr_mode2"){
      if(is.null(data) || is.null(input$GOI_x) || is.null(input$GOI_y)){
        corr <- NULL
      }else{
        if(input$GOI_x != "" && input$GOI_y != ""){
        corr <- cor.test(x=as.numeric(data[input$GOI_x,]),y=as.numeric(data[input$GOI_y,]),method=input$corr_statistics)
        df <- data.frame(x_axis = input$GOI_x, y_axis = input$GOI_y, statistics=corr$statistic,corr_score = corr$estimate,
                         pvalue = corr$p.value, method = corr$method,alternative=corr$alternative)
        return(df)
        }
      }
    }
  })
  corr_table_pair_forOutput <- reactive({
    df <- corr_table_pair()
    if(gene_type3() != "SYMBOL"){
      if(input$Species3 != "not selected"){
        colnames(df)[1] <- "x_axis_UniqueID"
        df$x_axis <- gsub(".+\\-","",df$x_axis_UniqueID)
        df$x_axis <- gsub(" ","",df$x_axis)
        df <- df %>% dplyr::select(x_axis, x_axis_UniqueID, everything())
        df$x_axis_UniqueID <- gsub("\n"," ",df$x_axis_UniqueID)
        df$y_axis <- gsub("\n"," ",df$y_axis)
      }
    }
    print(df)
          return(df)
  })
  output$statistical_table_corrplot <- renderDataTable({
    if(input$corr_mode == "corr_mode2"){
      if(!is.null(norm_GOI_corrplot_pair()) && !is.null(input$GOI_x) && !is.null(input$GOI_y) &&
         input$GOI_x != "" && input$GOI_y != ""){
        corr_table_pair_forOutput()
      }
    }else{
      if(!is.null(corr_table()) && !is.null(input$GOI_x) && input$GOI_x != ""){
        corr_table_forOutput()
      }
    }
  })
  corr_pair <- reactive({
    if(input$corr_mode == "corr_mode2"){
      if(!is.null(norm_GOI_corrplot_pair()) && !is.null(input$GOI_x) && !is.null(input$GOI_y) &&
         input$GOI_x != "" && input$GOI_y != ""){
        return(paste0(input$GOI_x,"-",input$GOI_y))
      }
    }else{
      if(!is.null(corr_table()) && !is.null(input$GOI_x) && 
         input$corr_start > 0 && updateCounter_corr$i > 0 && input$GOI_x != ""){
        return(paste0(input$GOI_x,"-screen"))
      }
    }
  })

  output$norm_corr_selected_list <- renderUI({
    if(input$corr_mode == "corr_mode1"){
      if(!is.null(corr_table()) && !is.null(input$GOI_x) && input$GOI_x != ""){
        if(dim(brush_info_corr())[1]!=0){
          GOI_y <- brush_info_corr()$prey
        }else{
          GOI_y <- corr_table()[input$statistical_table_corrplot_rows_selected,]$prey
        }
        selectInput("norm_corr_selected_list","select GOI", GOI_y, multiple = F)
      }
    }
  })
  output$download_norm_corr = downloadHandler(
    filename = function(){
      paste0(download_norm_dir(), "correlation_plot_",corr_pair(),".pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$norm_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$norm_pdf_height
        if(input$norm_pdf_width == 0){
          pdf_width <- 5
        }else pdf_width <- input$norm_pdf_width
        if(input$corr_mode == "corr_mode2"){
          p <- norm_GOI_corrplot_pair()
        }else{
          p <- norm_GOI_corrplot()
        }
        pdf(file, height = pdf_height, width = pdf_width)
        print(p)
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_norm_corr_selected = downloadHandler(
    filename = function(){
      paste0(download_norm_dir(), "correlation_plot_",input$GOI_x,"-selected.pdf")
    },
    content = function(file) {
      p <- norm_GOI_corrplot_selected_for_download()
      GOI_y <- names(p)
      processNum <- length(GOI_y)
      withProgress(message = paste0("Download correlation plots of all selected genes (",length(GOI_y),")"),{
        if(input$norm_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$norm_pdf_height
        if(input$norm_pdf_width == 0){
          pdf_width <- 5
        }else pdf_width <- input$norm_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        for(y in GOI_y){
          print(p[[y]]) 
          incProgress(1/processNum, message = paste0("Download ",y))
        }
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_statisics_corrplot = downloadHandler(
    filename = function(){
      paste0(download_norm_dir(), "correlation_",corr_pair(),".txt")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$corr_mode == "corr_mode2"){
          table <- corr_table_pair_forOutput()
        }else{
          table <- corr_table_forOutput()
        }
        write.table(table,file, row.names = F, col.names=TRUE, sep = "\t", quote = F)
        incProgress(1)
      })
    }
  )
  #norm kmeans------------------------------------------------------
  updateCounter_kmeans <- reactiveValues(i = 0)
  norm_kmeans_running <- reactiveVal(FALSE)
  norm_kmeans_status <- reactiveVal("")
  set_norm_kmeans_status <- function(message = NULL) {
    if (is.null(message) || !nzchar(message)) {
      norm_kmeans_status("")
      removeNotification(id = "norm_kmeans_progress")
      return(invisible(NULL))
    }
    message <- paste(as.character(message), collapse = "\n")
    current <- isolate(norm_kmeans_status())
    if (!identical(current, message)) {
      norm_kmeans_status(message)
      showNotification(
        message,
        id = "norm_kmeans_progress",
        duration = NULL,
        closeButton = FALSE,
        type = "message"
      )
    }
    invisible(message)
  }
  reset_norm_kmeans_state <- function() {
    updateCounter_kmeans$i <- 0
    norm_kmeans_running(FALSE)
    set_norm_kmeans_status(NULL)
    freezeReactiveValue(input, "kmeans_order")
    updateSelectInput(session, "kmeans_order",
                      choices = character(0), selected = character(0))
    freezeReactiveValue(input, "norm_select_kmean")
    updateSelectInput(session, "norm_select_kmean",
                      choices = character(0), selected = character(0))
  }
  build_norm_kmeans_order <- function(data.z, cl) {
    if (is.null(cl) || !length(cl)) {
      return(character(0))
    }
    cluster_ids <- sort(unique(cl))
    if (length(cluster_ids) <= 1) {
      return(as.character(cluster_ids))
    }
    cl_df <- data.frame(cluster = cl, stringsAsFactors = FALSE)
    data2 <- merge(cl_df, data.z, by = 0)
    rownames(data2) <- data2[, 1]
    data2 <- data2[, -1, drop = FALSE]
    df <- data.frame(matrix(nrow = 0, ncol = ncol(data2) - 1))
    if (ncol(data2) > 1) {
      colnames(df) <- colnames(data2[, -1, drop = FALSE])
    }
    for (i in cluster_ids) {
      data3 <- data2 %>% dplyr::filter(cluster == i)
      if (nrow(data3) == 0) {
        next
      }
      data4 <- colSums(data3[, -1, drop = FALSE], na.rm = TRUE)
      df <- rbind(df, data4)
    }
    if (nrow(df) <= 1) {
      return(as.character(cluster_ids))
    }
    as.character(cluster_ids[hclust(dist(df), "average")$order])
  }
  
  observeEvent(input$kmeans_start, {
    norm_kmeans_running(TRUE)
    set_norm_kmeans_status("k-means clustering is running.\nPreparing clusters and heatmap...")
    updateCounter_kmeans$i <- updateCounter_kmeans$i + 1
  }, ignoreInit = TRUE)
  
  output$norm_kmeans_progress_text <- renderText({
    norm_kmeans_status()
  })
  
  
  #Restart
  observeEvent(input$norm_kmeans_number, {
    reset_norm_kmeans_state()
  }, ignoreInit = TRUE) 
  observeEvent(input$kmeans_cv, {
    reset_norm_kmeans_state()
  }, ignoreInit = TRUE) 
  observeEvent(input$selectFC_norm, {
    reset_norm_kmeans_state()
  }, ignoreInit = TRUE) 
  observeEvent(input$fc3, {
    reset_norm_kmeans_state()
  }, ignoreInit = TRUE) 
  observeEvent(input$basemean3, {
    reset_norm_kmeans_state()
  }, ignoreInit = TRUE) 
  observeEvent(list(input$data_file_type3, norm_count_input(), norm_metadata()), {
    reset_norm_kmeans_state()
  }, ignoreInit = TRUE)
  output$selectFC_norm <- renderUI({
    if(is.null(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      selectizeInput("selectFC_norm", "Option 1: select a pair for fold change cut-off", c(unique(unique(gsub("\\_.*","", colnames(d_norm_count_matrix_cutofff()))))),
                     selected = "", multiple = TRUE, 
                     options = list(maxItems = 2))
    }
  })
  norm_select <- reactive({
    if (input$data_file_type3 == "Row5" || norm_metadata_is_single_column()) {
      pair <- unique(gsub("\\_.*","", colnames(d_norm_count_matrix_cutofff())))
    }else{
      count <- d_norm_count_matrix_cutofff_fc2()
      meta <- norm_metadata()
      if (is.null(count) || is.null(meta)) {
        return(character(0))
      }
      collist <- gsub("\\_.+$", "", colnames(count))
      collist <- gsub(" ", ".", collist)
      meta <- as.data.frame(meta, stringsAsFactors = FALSE)
      if (ncol(meta) < 2) {
        return(character(0))
      }
      meta <- data.frame(condition=factor(meta[,1]), type=factor(meta[,2]))
      pair <- paste0(meta$condition,meta$type)
      pair <- gsub("\\+",".",pair)
    }
    return(unique(pair))
  })
  output$selectFC_norm2 <- renderUI({
    if(is.null(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      selectizeInput("selectFC_norm2", "Option 2: select a pair for fold change cut-off", norm_select(),
                     selected = "", multiple = TRUE, 
                     options = list(maxItems = 2))
    }
  })
  output$filtered_region <- renderText({
    if(is.null(d_norm_count_matrix_cutofff_fc())){
      return(NULL)
    }else{ 
      if(length(input$selectFC_norm) == 2 && length(input$selectFC_norm2) == 2){
      print(paste0("The number of genes after the filtration (basemean > ", input$basemean3,", |log2(", 
                   input$selectFC_norm[1],"/", input$selectFC_norm[2],")| > ", log2(input$fc3),", |log2(",
                   input$selectFC_norm2[1],"/", input$selectFC_norm2[2],")| > ", log2(input$fc3),"): ", 
                   length(rownames(d_norm_count_matrix_cutofff_fc()))))
      }else if(length(input$selectFC_norm) == 2 && length(input$selectFC_norm2) != 2){
        print(paste0("The number of genes after the filtration (basemean > ", input$basemean3,", |log2(", 
                     input$selectFC_norm[1],"/", input$selectFC_norm[2],")| > ", log2(input$fc3),"): ", 
                     length(rownames(d_norm_count_matrix_cutofff_fc()))))
      }else if(length(input$selectFC_norm) != 2 && length(input$selectFC_norm2) == 2){
        print(paste0("The number of genes after the filtration (basemean > ", input$basemean3,", |log2(",
                     input$selectFC_norm2[1],"/", input$selectFC_norm2[2],")| > ", log2(input$fc3),"): ", 
                     length(rownames(d_norm_count_matrix_cutofff_fc()))))
      }else if(length(input$selectFC_norm) != 2 && length(input$selectFC_norm2) != 2){
        print(paste0("The number of genes after the filtration (basemean > ", input$basemean3,"): ", 
                     length(rownames(d_norm_count_matrix_cutofff_fc()))))
      }
    }
  })
  d_norm_count_matrix_cutofff_fc2 <- reactive({
    data <- d_norm_count_cutoff_uniqueID()
    if(!norm_has_matrix_data(data)){
      return(NULL)
    }
    if(!is.null(gene_ID_norm())){
        rownames(data) <- data$Unique_ID
        data <- data[, - which(colnames(data) == "SYMBOL")]
        data <- data[, - which(colnames(data) == "Unique_ID")]
    }
    if(length(input$selectFC_norm) == 2){
      if(dim(data)[1] != 0){
        cond1 <- input$selectFC_norm[1]
        cond2 <- input$selectFC_norm[2]
        cond1_num <- data %>% dplyr::select(.,starts_with(cond1)) %>% colnames() %>% length()
        cond2_num <- data %>% dplyr::select(.,starts_with(cond2)) %>% colnames() %>% length()
        cond1_ave <- data %>% dplyr::select(starts_with(cond1)) %>% rowSums(na.rm=TRUE)/cond1_num
        cond2_ave <- data %>% dplyr::select(starts_with(cond2)) %>% rowSums(na.rm=TRUE)/cond2_num
        Log2FoldChange <- log((cond1_ave + 0.01)/(cond2_ave + 0.01),2)
        data$Log2FoldChange <- Log2FoldChange
        print(dplyr::filter(data, abs(Log2FoldChange) <= log2(input$fc3)))
        data2 <- data %>% dplyr::filter(abs(Log2FoldChange) > log2(input$fc3))
        data2 <- data2[, - which(colnames(data2) == "Log2FoldChange")]
      }else data2 <- NULL
    }else data2 <- data
    return(data2)
  })
  d_norm_count_matrix_cutofff_fc <- reactive({
    data <- d_norm_count_matrix_cutofff_fc2()
    if(!norm_has_matrix_data(data)){
      return(NULL)
    }
    if(length(input$selectFC_norm2) == 2){
      if(dim(data)[1] != 0){
        cond1 <- input$selectFC_norm2[1]
        cond2 <- input$selectFC_norm2[2]
        cond1_num <- data %>% dplyr::select(.,starts_with(cond1)) %>% colnames() %>% length()
        cond2_num <- data %>% dplyr::select(.,starts_with(cond2)) %>% colnames() %>% length()
        cond1_ave <- data %>% dplyr::select(starts_with(cond1)) %>% rowSums(na.rm=TRUE)/cond1_num
        cond2_ave <- data %>% dplyr::select(starts_with(cond2)) %>% rowSums(na.rm=TRUE)/cond2_num
        Log2FoldChange <- log((cond1_ave + 0.01)/(cond2_ave + 0.01),2)
        data$Log2FoldChange <- Log2FoldChange
        data2 <- data %>% dplyr::filter(abs(Log2FoldChange) > log2(input$fc3))
        data2 <- data2[, - which(colnames(data2) == "Log2FoldChange")]
      }else data2 <- NULL
    }else data2 <- data
    return(data2)
  })
  output$norm_kmeans_num <- renderUI({
    if(!norm_has_matrix_data(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      withProgress(message = "Preparing kmeans clustering",{
        sliderInput("norm_kmeans_number", "k-means number", min = 1,
                    max=20, step = 1,
                    value = 1)
      })
    }
  })
  
  output$kmeans_cv <- renderUI({
    sliderInput("kmeans_cv", "Most variable genes:", min = 0,
                max=8000, step = 100,
                value = 2000)
  })
  
  
  norm_count_matrix_cutoff2 <- reactive({
    data <- d_norm_count_matrix_cutofff_fc()
    if(!norm_has_matrix_data(data) || is.null(input$kmeans_cv)){
      return(NULL)
    }else{
      top_n <- min(nrow(data), input$kmeans_cv)
      if (top_n < 1) {
        return(NULL)
      }
      data2 <- data[order(apply(data,1,mad), decreasing = T)[seq_len(top_n)], , drop = FALSE]
      if (!norm_has_matrix_data(data2)) {
        return(NULL)
      }
      return(data2)
    }
  })
  
  norm_data_z <- reactive({
    data <- norm_count_matrix_cutoff2()
    if(!norm_has_matrix_data(data)){
      return(NULL)
    }else{
      data.z <- genefilter::genescale(data, axis = 1, method = "Z")
      data.z <- na.omit(data.z)
      if(!norm_has_matrix_data(data.z)){
        return(NULL)
      }
      return(data.z)
    }
  })
  
  norm_kmeans_run <- eventReactive(input$kmeans_start, {
    data <- norm_count_matrix_cutoff2()
    data.z <- norm_data_z()
    if(!norm_has_matrix_data(data) || is.null(data.z) || nrow(data.z) < 2 || ncol(data.z) < 1 || is.null(input$norm_kmeans_number) ||
       input$norm_kmeans_number < 1){
      norm_kmeans_running(FALSE)
      set_norm_kmeans_status(NULL)
      return(NULL)
    }else{
      tryCatch({
        withProgress(message = "k-means clustering", value = 0,{
          cluster_n <- min(input$norm_kmeans_number, nrow(data.z))
          if (cluster_n < 1) {
            norm_kmeans_running(FALSE)
            set_norm_kmeans_status(NULL)
            return(NULL)
          }
          set_norm_kmeans_status("k-means clustering is running.\nComputing consensus clusters...")
          incProgress(0.15, detail = "Computing consensus clusters")
          set.seed(123)
          cl = consensus_kmeans(data.z, cluster_n, 100)
          names(cl) <- rownames(data.z)
          if(length(unique(cl)) != cluster_n && cluster_n > 1){
            set.seed(123)
            cl = consensus_kmeans(data.z, cluster_n-1, 100)
            names(cl) <- rownames(data.z)
            if(length(unique(cl)) != cluster_n-1 && cluster_n > 2){
              set.seed(123)
              cl = consensus_kmeans(data.z, cluster_n-2, 100)
              names(cl) <- rownames(data.z)
            }
          }
          set_norm_kmeans_status("k-means clustering is running.\nPreparing cluster order...")
          incProgress(0.55, detail = "Preparing cluster order")
          order <- build_norm_kmeans_order(data.z, cl)
          incProgress(0.30, detail = "Finalizing k-means inputs")
        })
        return(list(
          data = data,
          data.z = data.z,
          cl = cl,
          order = order
        ))
      }, error = function(e) {
        norm_kmeans_running(FALSE)
        set_norm_kmeans_status(NULL)
        stop(e)
      })
    }
  }, ignoreInit = TRUE)
  pre_norm_kmeans <- reactive({
    run <- norm_kmeans_run()
    if (is.null(run)) {
      return(NULL)
    }
    run$cl
  })
  pre_norm_kmeans_order <- reactive({
    run <- norm_kmeans_run()
    if (is.null(run)) {
      return(character(0))
    }
    run$order
  })
  
  observeEvent(pre_norm_kmeans_order(), {
    order <- as.character(pre_norm_kmeans_order())
    selected <- intersect(isolate(input$kmeans_order), order)
    if (!length(selected)) {
      selected <- order
    }
    freezeReactiveValue(input, "kmeans_order")
    updateSelectInput(session, "kmeans_order",
                      choices = order, selected = selected)
  }, ignoreInit = TRUE)
  
  norm_kmeans <- reactive({
    run <- norm_kmeans_run()
    if(is.null(run) || is.null(pre_norm_kmeans()) || updateCounter_kmeans$i == 0){
      return(NULL)
    }else{
      cluster_order <- as.character(input$kmeans_order)
      if (!length(cluster_order)) {
        cluster_order <- as.character(pre_norm_kmeans_order())
      }
      if(length(cluster_order) == length(unique(pre_norm_kmeans()))){
        set.seed(123)
        ht <- Heatmap(run$data.z, name = "z-score",
                      column_order = colnames(run$data.z),
                      clustering_method_columns = 'ward.D2',
                      cluster_row_slices = F, split = factor(pre_norm_kmeans(),levels = cluster_order),
                      show_row_names = F,column_names_side = "top",use_raster = TRUE)
      }else validate("Select all clusters from 'Order of clusters on heatmap'")
      return(ht)
    }
  })
  norm_kmeans_GOI <- reactive({
    ht <- norm_kmeans()
    run <- norm_kmeans_run()
    if (is.null(run)) {
      return(NULL)
    }
    data.z <- run$data.z
    if(is.null(ht) || is.null(pre_norm_kmeans()) || 
       updateCounter_kmeans$i == 0){
      return(NULL)
    }else{
      if(!is.null(input$norm_kmeans_count_table_rows_selected)){
        clusters <- norm_kmeans_cluster()
        if(!norm_has_matrix_data(clusters)){
          return(ht)
        }
        data <- clusters[input$norm_kmeans_count_table_rows_selected, , drop = FALSE]
        if(!norm_has_matrix_data(data)){
          return(ht)
        }
        lab <- rownames(data)
        if(input$Species3 != "not selected"){
          if(gene_type3() != "SYMBOL"){
            lab <- data$Unique_ID
          }
        }
        indexes <- which(rownames(data.z) %in% lab)
        labels <- rownames(data.z)[indexes]
        set.seed(123)
        ht <- ht + rowAnnotation(
          link = anno_mark(at = indexes, labels = labels,which="row",link_width = unit(1, "cm"),
                           labels_gp = gpar(fontface = "italic")),
          width = unit(1, "cm") + max_text_width(labels))
      }
      return(ht)
    }
  })
  
  pre_norm_kmeans_cluster <- reactive({
    ht <- norm_kmeans()
    run <- norm_kmeans_run()
    if (is.null(run)) {
      return(NULL)
    }
    data.z <- run$data.z
    data <- run$data
    if(is.null(ht) || !norm_has_matrix_data(data.z) || !norm_has_matrix_data(data)){
      return(NULL)
    }else{
      r.dend <- suppressWarnings(row_dend(ht))
      rcl.list <- suppressWarnings(row_order(ht))
      lapply(rcl.list, function(x) length(x))
      Cluster <- NULL
      expected_clusters <- length(unique(pre_norm_kmeans()))
      if(expected_clusters > 0){
        if(length(lapply(rcl.list, function(x) length(x))) != expected_clusters){
          return(NULL)
        }else{
          out <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
          cluster_order <- as.character(input$kmeans_order)
          if (!length(cluster_order)) {
            cluster_order <- as.character(pre_norm_kmeans_order())
          }
          for (i in cluster_order){
            genes <- safe_cluster_genes(data.z, suppressWarnings(row_order(ht)[[i]]))
            if (!length(genes)) next
            clu <- data.frame(GeneID = genes, Cluster = paste("cluster", i, sep=""), stringsAsFactors = FALSE)
            out <- rbind(out, clu)
          }
          if (!nrow(out)) return(NULL)
          colnames(out) <- c("GeneID", "Cluster")
          out <- as.data.frame(out)
          rownames(out) <- out$GeneID
          clusterCount <- merge(out, data, by=0)
          rownames(clusterCount) <- clusterCount$GeneID
          if(input$Species3 != "not selected"){
            if(gene_type3() != "SYMBOL"){
              rownames(clusterCount) <- gsub(".+\\ ", "", clusterCount$GeneID)
            }}
          clusterCount <- clusterCount[,-1:-2]
          return(clusterCount)
        }
      }else return(NULL)
    }
  })
  norm_kmeans_cluster <- reactive({
    count <- pre_norm_kmeans_cluster()
    if(!norm_has_matrix_data(count)){
      return(NULL)
    }
    if(!is.null(gene_ID_norm())){
        gene_IDs  <- gene_ID_norm()
        data2 <- merge(count, gene_IDs, by= 0)
        rownames(data2) <- data2[,1]
        data2 <- data2[, - which(colnames(data2) == "Row.names.y")]
        data2$Unique_ID <- paste(data2$SYMBOL,data2$Row.names, sep = "\n- ")
        count <- data2[,-1]
      }
    return(count)
  })
  
  observeEvent(norm_kmeans_cluster(), {
    clusters <- norm_kmeans_cluster()
    choices <- if (is.null(clusters)) character(0) else unique(clusters$Cluster)
    selected <- intersect(isolate(input$norm_select_kmean), choices)
    if (!length(selected)) {
      selected <- choices
    }
    freezeReactiveValue(input, "norm_select_kmean")
    updateSelectInput(session, "norm_select_kmean",
                      choices = choices, selected = selected)
  }, ignoreInit = TRUE)
  
  ###----
  norm_kmeans_box <- reactive({
    res <- norm_kmeans_cluster()
    ma  <- as.data.frame(norm_count_matrix_cutoff2())
    meta <- norm_metadata()
    meta_has_single_column <- FALSE
    
    if (is.null(ma) || is.null(res)) return(NULL)
    
    withProgress(message = "Boxplot", {
      if (input$data_file_type3 == "Row5") {
        collist <- gsub("\\_.+$", "", colnames(ma))
        meta <- data.frame(
          condition = factor(collist),
          row.names = colnames(ma),
          stringsAsFactors = FALSE
        )
      } else {
        if (is.null(meta)) {
          return(NULL)
        }
        meta <- as.data.frame(meta, stringsAsFactors = FALSE)
        if (ncol(meta) == 0) {
          return(NULL)
        }
        meta_has_single_column <- ncol(meta) == 1
        
        if (meta_has_single_column) {
          collist <- gsub("\\_.+$", "", colnames(ma))
          meta <- data.frame(
            condition = factor(collist),
            row.names = colnames(ma),
            stringsAsFactors = FALSE
          )
        } else {
          meta <- data.frame(
            condition = factor(meta[[1]]),
            type      = factor(meta[[2]]),
            row.names = colnames(ma),
            stringsAsFactors = FALSE
          )
          meta$type <- factor(meta$type, levels = unique(meta$type), ordered = TRUE)
        }
      }
      
      meta$condition <- factor(meta$condition, levels = unique(meta$condition), ordered = TRUE)
      
      res$Cluster <- gsub("cluster", "", res$Cluster)
      res <- data.frame(genes = rownames(res), cluster = res$Cluster, stringsAsFactors = FALSE)
      
      table <- rownames_to_column(as.data.frame(ma), "genes") %>%
        gather("sample", "expression", -genes) %>%
        right_join(distinct(res[, c("genes", "cluster")]), by = "genes") %>%
        left_join(rownames_to_column(as.data.frame(meta), "sample"), by = "sample") %>%
        as.data.frame()
      
      table$cluster <- as.integer(table$cluster)
      table <- na.omit(table)
      
      if (input$data_file_type3 == "Row5" || meta_has_single_column) {
        p <- try(
          degPlotCluster(table, time = "condition", process = TRUE, min_genes = 0) +
            scale_color_brewer(palette = "Set1", direction = -1) +
            theme_bw(base_size = 15) +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 90, hjust = 1))
        )
      } else {
        p <- try(
          degPlotCluster(table, time = "condition", color = "type", process = TRUE, min_genes = 0) +
            scale_color_brewer(palette = "Set1", direction = -1) +
            theme_bw(base_size = 15) +
            theme(legend.position = "top",
                  axis.text.x = element_text(angle = 90, hjust = 1))
        )
      }
      print(p)
      if (inherits(p, "try-error")) {
        stop(as.character(p))
      }
      
      return(p)
    })
  })
  
  output$norm_kmeans_boxplot<- renderPlot({
    if(is.null(norm_kmeans_box())|| length(input$selectFC_norm2) != 2){
      return(NULL)
    }else{
      print(norm_kmeans_box())
    }
  })
  
  output$download_norm_kmeans_boxplot = downloadHandler(
    filename = function(){
      kmeans_label <- input$norm_kmeans_number
      if (is.null(kmeans_label) || !length(kmeans_label) || is.na(kmeans_label[[1]])) {
        kmeans_label <- "kmeans"
      } else {
        kmeans_label <- as.character(kmeans_label[[1]])
      }
      paste0(download_norm_dir(), kmeans_label, "_kmeans_boxplot.pdf")
    },
    contentType = "application/pdf",
    content = function(file) {
      withProgress(message = "Preparing download",{
        plot_obj <- isolate(norm_kmeans_box())
        if (is.null(plot_obj)) {
          stop("k-means boxplot is not available for download.")
        }
        cluster_values <- tryCatch(unique(plot_obj$data$cluster), error = function(e) integer())
        cluster_count <- length(cluster_values)
        if(input$norm_pdf_height == 0){
          pdf_height <- if (cluster_count > 0) pdf_h(seq_len(cluster_count)) + 2 else 10
        }else pdf_height <- input$norm_pdf_height
        if(input$norm_pdf_width == 0){
          pdf_width <- if (cluster_count > 0) pdf_w(seq_len(cluster_count)) + 2 else 7
        }else pdf_width <- input$norm_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        on.exit(dev.off(), add = TRUE)
        print(plot_obj)
        incProgress(1)
      })
    }
  )
  outputOptions(output, "download_norm_kmeans_boxplot", suspendWhenHidden = FALSE)
  ###---
  
  norm_kmeans_pattern_extract <- reactive({
    clusters <- norm_kmeans_cluster()
    if(is.null(clusters)){
      return(NULL)
    }else{
      if(is.null(input$norm_select_kmean)){
        return(NULL)
      }else{
        cluster_name <- input$norm_select_kmean
        clusterCount <- dplyr::filter(clusters, Cluster %in% cluster_name)
        clusterCount <- clusterCount[,-1]
        if(input$Species3 != "not selected"){
          if(gene_type3() != "SYMBOL"){
            clusterCount <- clusterCount[, - which(colnames(clusterCount) == "Unique_ID")]
          }
        }
        return(clusterCount)
      }
    }
  })
  
  output$norm_kmeans_extract_table <- renderDT({
    norm_kmeans_pattern_extract()
  })
  
  output$norm_kmeans_heatmap <- renderPlot({
    if(is.null(input$norm_kmeans_count_table_rows_selected)) ht <- norm_kmeans() else ht <- norm_kmeans_GOI()
    if(is.null(ht)){
      return(NULL)
    }else{
      if (isTRUE(norm_kmeans_running())) {
        set_norm_kmeans_status("k-means clustering is running.\nRendering heatmap...")
      }
      withProgress(message = "Draw heatmap", value = 0,{
        set.seed(123)
        incProgress(0.2, detail = "Rendering heatmap")
        draw(ht)
        incProgress(0.8)
      })
      if (isTRUE(norm_kmeans_running())) {
        norm_kmeans_running(FALSE)
        set_norm_kmeans_status(NULL)
      }
    }
  })
  
  output$norm_kmeans_count_table <- renderDT({
    norm_kmeans_cluster()
  })
  
  
  norm_kmeans_GOIbox <- reactive({
    if(!is.null(input$norm_kmeans_count_table_rows_selected)){
      clusters <- norm_kmeans_cluster()
      if(!norm_has_matrix_data(clusters)){
        return(NULL)
      }
      data <- clusters[input$norm_kmeans_count_table_rows_selected, , drop = FALSE]
      if(!norm_has_matrix_data(data)){
        return(NULL)
      }
      if(input$Species3 != "not selected"){
        if(gene_type3() != "SYMBOL"){
          rownames(data) <- paste(data$SYMBOL,rownames(data), sep = "\n- ")
          data <- data[, - which(colnames(data) == "SYMBOL")]
          data <- data[, - which(colnames(data) == "Unique_ID")]
        }
      }
      data <- data[, - which(colnames(data) == "Cluster"), drop = FALSE]
      if(!norm_has_matrix_data(data)){
        return(NULL)
      }
      return(data)
    }
  })
  
  output$norm_kmeans_box <- renderPlot({
    if(!is.null(input$norm_kmeans_count_table_rows_selected)){
      data <- norm_kmeans_GOIbox()
      if(!norm_has_matrix_data(data)){
        return(NULL)
      }
      GOIboxplot(data = data)
    }
  })
  
  output$download_norm_kmeans_box = downloadHandler(
    filename = function(){
      paste0(download_norm_dir(), "_boxplot.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- norm_kmeans_GOIbox()
        rowlist <- rownames(data)
        if(input$norm_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)
        }else pdf_height <- input$norm_pdf_height
        if(input$norm_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)
        }else pdf_width <- input$norm_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(GOIboxplot(data = data))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  observeEvent(input$norm_kmeans_count_table_rows_selected, ({
    updateCollapse(session,id =  "norm_kmeans_collapse_panel", open="kmeans_box_panel")
  }))
  observeEvent(norm_count_input(), ({
    updateCollapse(session,id =  "norm_input_collapse_panel", open="Norm_count_panel")
  }))
  observeEvent(norm_metadata(), ({
    updateCollapse(session,id =  "norm_input_collapse_panel", open="norm_Metadata_panel")
  }))
  
  
  output$download_norm_kmeans_cluster = downloadHandler(
    filename = function() {
      paste0(download_norm_dir(), "kmeans_count_table.txt")
    },
    content = function(file){write.table(norm_kmeans_cluster(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  output$download_norm_kmeans_heatmap = downloadHandler(
    filename = function() {
      kmeans_label <- input$norm_kmeans_number
      if (is.null(kmeans_label) || !length(kmeans_label) || is.na(kmeans_label[[1]])) {
        kmeans_label <- "kmeans"
      } else {
        kmeans_label <- as.character(kmeans_label[[1]])
      }
      paste0(download_norm_dir(), kmeans_label, "_kmeans_heatmap.pdf")
    },
    content = function(file){
      withProgress(message = "Preparing download",{
        if(input$norm_pdf_height == 0){
          pdf_height <- 10
        }else pdf_height <- input$norm_pdf_height
        if(input$norm_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$norm_pdf_width
        if(is.null(input$norm_kmeans_count_table_rows_selected)) ht <- norm_kmeans() else ht <- norm_kmeans_GOI()
        pdf(file, height = pdf_height, width = pdf_width)
        print(ht)
        dev.off()
        incProgress(1)
      })
    }
  )
  outputOptions(output, "download_norm_kmeans_heatmap", suspendWhenHidden = FALSE)
  
  output$download_norm_kmeans_extract_count = downloadHandler(
    filename = function() {
      paste0(download_norm_dir(),"kmeans_selected_table.txt")
    },
    content = function(file){write.table(norm_kmeans_pattern_extract(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  
  #venn diagram ------------------
  
  output$select_file2 <- renderUI({
    if(is.null(overlap_list())){
      return(NULL)
    }else{
      selectInput("selectfile", "gene_list", choices = c("not selected", overlap_list()), selected = "not selected",multiple = F)
    }
  })
  
  files_table <- reactive({
    upload = list()
    name = c()
    if(is.null(input$files)){
      if(input$goButton_venn > 0 ){
        df <- list()
        df[["day0"]] <-  c(rownames(read.table(example_data_path("data/example11.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example11.txt"), header = T, row.names = 1, sep="\t")))
        df[["day1"]] <- c(rownames(read.table(example_data_path("data/example12.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example12.txt"), header = T, row.names = 1, sep="\t")))
        df[["day5"]] <- c(rownames(read.table(example_data_path("data/example13.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example13.txt"), header = T, row.names = 1, sep="\t")))
        return(df)
      }
      return(NULL)
    }else{
      for(nr in 1:length(input$files[, 1])){
        df <- read_df(input$files[[nr, 'datapath']])
        name <- c(name, gsub(paste0("\\.",tools::file_ext(input$files[[nr, 'datapath']]),"$"), "", input$files[nr,]$name))
        upload[[nr]] <- c(rownames(df))
      }
      names(upload) <- name
      return(upload)
    }
  })
  
  observeEvent(input$goButton_venn,({
    updateSelectInput(session,inputId = "Species7","Species",species_list, selected = "Mus musculus")
  }))
  
  output$eulerr_label <- renderUI({
    if(input$venn_type == "eulerr") radioButtons("eulerr_label","label",c("ON"="ON","OFF"="OFF"),"ON")
  })
  output$venn <- renderPlot({
    req(files_table())
    
    gene_list <- files_table()
    names(gene_list) <- vapply(names(gene_list), function(x) {
      x <- gsub("_", " ", x)
      paste(strwrap(x, width = 15), collapse = "\n")
    }, character(1))
    
    if (input$venn_type == "default" || is.null(input$eulerr_label)) {
      venn::venn(
        gene_list,
        ilabels = "counts",
        zcolor = "style",
        opacity = 0,
        ilcs = 1.5,
        sncs = 1.5
      )
    } else {
      if (input$eulerr_label == "ON") {
        label <- list(cex = 2)
      } else {
        label <- NULL
      }
      
      plot(
        euler(gene_list, shape = "ellipse"),
        labels = label,
        quantities = list(type = "counts", cex = 2),
        edges = list(col = as.vector(seq_along(names(gene_list))), lex = 2),
        fills = list(fill = rep("white", length(names(gene_list)))),
        legend = list(side = "right", cex = 2)
      )
    }
  })
  
  overlap_list <- reactive({
    req(files_table())
    ints <- venn::extractInfo(files_table(), what = "intersections", use.names = TRUE)
    names(ints)
  })
  
  overlap_table2 <- reactive({
    req(files_table())
    ints <- venn::extractInfo(files_table(), what = "intersections", use.names = TRUE)
    
    do.call(rbind, lapply(names(ints), function(nm) {
      data.frame(
        Gene = unlist(ints[[nm]]),
        Group = nm,
        stringsAsFactors = FALSE
      )
    }))
  })
  
  
  output$venn_result <- renderDataTable({
    overlap_table2()
  })
  
  output$download_vennplot = downloadHandler(
    filename ="venn_diagram.pdf",
    content = function(file){
      if(is.null(files_table())){
        return(NULL)
      }else{
        withProgress(message = "Preparing download",{
          gene_list <- files_table()
          for(i in 1:length(names(gene_list))){
            names(gene_list)[i] <- gsub("_", " ", names(gene_list)[i])
            names(gene_list)[i] <- paste(strwrap(names(gene_list)[i], width = 15),collapse = "\n")
          }
          if(input$venn_pdf_height == 0){
            pdf_height <- 3
          }else pdf_height <- input$venn_pdf_height
          if(input$venn_pdf_width == 0){
            pdf_width <- 3
          }else pdf_width <- input$venn_pdf_width
          pdf(file, height = pdf_height, width = pdf_width)
          if(input$venn_type == "default" || is.null(input$eulerr_label)) print(venn::venn(gene_list, ilabels = TRUE, zcolor = "style", opacity = 0, ilcs = 1.5, sncs = 1.5)) else{
            if(input$eulerr_label =="ON") label=list(cex=0.8) else label=NULL
            print(plot(euler(gene_list, shape = "ellipse"), 
                 labels = label,quantities = list(type="counts",cex=0.8),
                 edges = list(col=as.vector(seq(1,length(names(gene_list)))),lex = 0.8),
                 fills = list(fill=rep("white",length(names(gene_list)))),legend = list(side = "right",cex=0.8)) )
          }
          dev.off()
          incProgress(1)
        })
      }
    }
  )
  
  
  countfiles_integrated <- reactive({
    if(is.null(input$countfiles)){
      if(input$goButton_venn > 0){
        df <- list()
        df["day0"] <- list(read.table(example_data_path("data/day0.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/day0.txt"), header = T, row.names = 1))
        df["day1"] <- list(read.table(example_data_path("data/day1.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/day1.txt"), header = T, row.names = 1))
        df["day5"] <- list(read.table(example_data_path("data/day5.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/day5.txt"), header = T, row.names = 1))
        return(df)
      }
      return(NULL)
    }else{
      upload = list()
      name = c()
      for(nr in 1:length(input$countfiles[, 1])){
        df <- read_df(input$countfiles[[nr, 'datapath']])
        file_name <- gsub(paste0("\\.",tools::file_ext(input$countfiles[[nr, 'datapath']]),"$"), "", input$countfiles[nr,]$name)
        name <- c(name, file_name)
        upload[[nr]] <- list(df)
      }
      names(upload) <- name
      return(upload)
    }
  })
  
  integrated_count <- reactive({
    if(!is.null(input$selectfile)){
      files <- countfiles_integrated()
      if(is.null(files)){
        return(NULL)
      }else{
        if(input$selectfile != "not selected"){
          cluster_file <- overlap_table2()
          rownames(cluster_file) <- cluster_file$Gene
          gene_set <- dplyr::filter(cluster_file, Group == input$selectfile)
        }
        matrix_list <- list()
        for (name in names(files)) {
          matrix <- as.data.frame(files[name])
          matrix_2 <- matrix
          matrix_3 <- merge(matrix, matrix_2, by = 0)[,-2:-(1 + length(colnames(matrix)))]
          matrix_list[name] <- list(matrix_3)
        }
        base <- Reduce(function(x, y) merge(x, y, by = "Row.names"), matrix_list)
        rownames(base) <- base$Row.names
        base <- data.matrix(base[,-1])
        if(input$selectfile != "not selected"){
          base <- merge(gene_set, base, by = 0)
          if(length(colnames(gene_set)) != 0){
            base <- base[,-2:-(1 + length(colnames(gene_set)))]
          }
          rownames(base) <- base$Row.names
          base <- data.matrix(base[,-1])
        }
        colnames(base) <- gsub("\\.y$", "", colnames(base))
        if(length(names(files)) == 1) {
          colnames(base) <- gsub(names(files),"",colnames(base))
          colnames(base) <- gsub("^\\.", "", colnames(base))
        }
        return(base)
      }
    }else return(NULL)
  })
  
  integrated_count_z <- reactive({
    if(!is.null(input$selectfile)){
      files <- countfiles_integrated()
      if(is.null(files)){
        return(NULL)
      }else{
        if(input$selectfile != "not selected"){
          cluster_file <- overlap_table2()
          rownames(cluster_file) <- cluster_file$Gene
          gene_set <- dplyr::filter(cluster_file, Group == input$selectfile)
        }
        matrix_list <- list()
        matrix_z_list <- list()
        for (name in names(files)) {
          matrix <- as.data.frame(files[name])
          matrix_2 <- matrix
          if(input$pre_zscoring == "TRUE"){
            matrix_z <- genefilter::genescale(matrix_2, axis = 1, method = "Z")
            matrix_z <- na.omit(matrix_z)
            matrix_z <- merge(matrix, matrix_z, by = 0)[,-2:-(1 + length(colnames(matrix)))]
            matrix_z_list[name] <- list(matrix_z)
          }
          matrix_3 <- merge(matrix, matrix_2, by = 0)[,-2:-(1 + length(colnames(matrix)))]
          matrix_list[name] <- list(matrix_3)
        }
        base <- Reduce(function(x, y) merge(x, y, by = "Row.names"), matrix_list)
        rownames(base) <- base$Row.names
        base <- data.matrix(base[,-1])
        if(input$selectfile != "not selected"){
          base <- merge(gene_set, base, by = 0)
          if(length(colnames(gene_set)) != 0){
            base <- base[,-2:-(1 + length(colnames(gene_set)))]
          }
          rownames(base) <- base$Row.names
          base <- data.matrix(base[,-1])
        }
        if(input$pre_zscoring == "TRUE"){
          base_z <- Reduce(function(x, y) merge(x, y, by = "Row.names"), matrix_z_list)
          rownames(base_z) <- base_z$Row.names
          base_z <- data.matrix(base_z[,-1])
          if(input$selectfile != "not selected"){
            base_z <- merge(gene_set, base_z, by = 0)
            if(length(colnames(gene_set)) != 0){
              base_z <- base_z[,-2:-(1 + length(colnames(gene_set)))]
            }
            rownames(base_z) <- base_z$Row.names
            base_z <- data.matrix(base_z[,-1])
          }
        }
        if(input$pre_zscoring != "TRUE"){
          base <- genefilter::genescale(base, axis = 1, method = "Z")
          base_z <- na.omit(base)
        }
        colnames(base_z) <- gsub("\\.y$", "", colnames(base_z))
        if(length(names(files)) == 1) {
          colnames(base_z) <- gsub(names(files),"",colnames(base_z))
          colnames(base_z) <- gsub("^\\.", "", colnames(base_z))
        }
        return(base_z)
      }
    }else return(NULL)
  })
  
  
  integrated_heatmap <- reactive({
    if(!is.null(input$selectfile)){
      base_z <- integrated_count_z()
      if(input$selectfile == "not selected" || is.null(base_z)){
        return(NULL)
      }else{
        withProgress(message = "Integrated heatmap",{
          cond <- gsub("\\_.+$", "", colnames(base_z))
          cond <- gsub(".+\\.", "", cond)
          cond <- factor(cond, levels = unique(cond), ordered = TRUE)
          cond_color <- rainbow_hcl(length(unique(cond)),c=100)
          names(cond_color) <- unique(cond)
          if(length(rownames(base_z)) <= 50){
            ht <- Heatmap(base_z, name = "z-score",
                          clustering_method_columns = 'ward.D2',
                          cluster_row_slices = T, show_row_names = T,
                          top_annotation = HeatmapAnnotation(condition = cond, col = list(condition = cond_color)),
                          column_names_side = "top",use_raster = TRUE,
                          row_names_gp = gpar(fontface = "italic"))
          }else{
            ht <- Heatmap(base_z, name = "z-score",
                          clustering_method_columns = 'ward.D2',
                          cluster_row_slices = T, show_row_names = F,use_raster = TRUE,
                          top_annotation = HeatmapAnnotation(condition = cond, col = list(condition = cond_color)),
                          column_names_side = "top")
          }
          incProgress(1)
          return(draw(ht))
        })
      }
    }
  })
  
  output$GOIbox_venn <- renderPlot({
    if(!is.null(input$integrated_count_table_rows_selected)){
      data <- as.data.frame(integrated_count())
      GOIboxplot(data = data[input$integrated_count_table_rows_selected,])
    }
  })
  
  output$intheatmap <- renderPlot({
    print(integrated_heatmap())
  })
  output$integrated_count_table <- renderDT({
    integrated_count()
  })
  output$integrated_count_z_table <- renderDataTable({
    integrated_count_z()
  })
  
  observeEvent(input$integrated_count_table_rows_selected, ({
    updateCollapse(session,id =  "int_result_collapse_panel", open="integrated_count_box_panel")
  }))
  
  output$download_venn_result = downloadHandler(
    filename ="venn_result.txt",
    content = function(file){write.table(overlap_table2(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  output$download_integrated_count_table = downloadHandler(
    filename = function(){
      paste(input$selectfile, "integrated_count_table.txt", sep="_")
    },
    content = function(file){write.table(integrated_count(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  output$download_integrated_z_count_table = downloadHandler(
    filename = function(){
      paste(input$selectfile, "integrated_zscored_count_table.txt", sep="_")
    },
    content = function(file){write.table(integrated_count_z(), file, row.names = T, col.names=NA, sep = "\t", quote = F)}
  )
  output$download_integrated_heatmap = downloadHandler(
    filename = function(){
      paste(input$selectfile, "integrated_heatmap.pdf", sep="_")
    },
    content = function(file){
      withProgress(message = "Preparing download",{
        if(input$venn_pdf_height == 0){
          pdf_height <- 8
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 8
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(integrated_heatmap())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_GOIbox_venn = downloadHandler(
    filename = function(){
      paste(input$selectfile, "boxplot.pdf", sep="_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- as.data.frame(integrated_count())
        rowlist <- rownames(data[input$integrated_count_table_rows_selected,])
        if(input$venn_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(GOIboxplot(data = data[input$integrated_count_table_rows_selected,]))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  #venn enrichment---------
  output$venn_whichGroup1 <- renderUI({
    clusters <- overlap_table2()
    if(is.null(clusters)){
      return(NULL)
    }else{
      selectInput("venn_whichGroup1", "gene list", choices = c(sort(unique(clusters$Group))),multiple = TRUE)
    }
  })
  
  pre_venn_enrich_input1 <- reactive({
    if(is.null(input$venn_whichGroup1)){
      return(NULL)
    }else{
      clusters <- overlap_table2()
      cluster_name <- input$venn_whichGroup1
      clusterCount <- data.frame(GeneID = NA, Group=NA)
      for(name in cluster_name){
        clusterCount2 <- dplyr::filter(clusters, Group == name)
        clusterCount2 <- data.frame(GeneID = clusterCount2$Gene, Group = clusterCount2$Group)
        clusterCount <- rbind(clusterCount, clusterCount2) 
      }
      clusterCount <- na.omit(clusterCount)
      return(clusterCount)
    }
  })
  
  venn_enrich_input1 <- debounce(pre_venn_enrich_input1,1000)
  
  output$Gene_set9 <- renderUI({
    if(input$Species7 != "Xenopus laevis" && input$Ortholog7 != "Arabidopsis thaliana" && input$Species7 != "Arabidopsis thaliana"){
      selectInput('Gene_set9', 'Gene Set', gene_set_list)
    }else selectInput('Gene_set9', 'Gene Set', c("KEGG", "GO biological process", 
                                                 "GO cellular component","GO molecular function"))
  })
  
  output$venn_Spe <- renderText({
    if(input$Species7 == "not selected") print("Please select 'Species'")
  })
  
  org7 <- reactive({
    return(org(Species = input$Species7,Ortholog = input$Ortholog7))
  })
  ortholog7 <- reactive({
    return(no_org_ID(gene_list = overlap_table2(),Species = input$Species7,Ortholog = input$Ortholog7,Biomart_archive=input$Biomart_archive7))
  })
  gene_type7 <- reactive({
    data <- overlap_table2()
    data <- data.frame(row.names = unique(data[,1]))
    return(gene_type(my.symbols=rownames(data),org=org7(),Species=input$Species7,RNA_type=input$Level_venn))
  })
  org_code7 <- reactive({
    return(org_code(Species = input$Species7, Ortholog= input$Ortholog7))
  })
  
  venn_Hallmark_set <- reactive({
    return(GeneList_for_enrichment(Species = input$Species7, Ortholog=input$Ortholog7,Biomart_archive=input$Biomart_archive7, Gene_set = input$Gene_set9, org = org7()))
  })
  
  venn_enrich_base1 <- reactive({
    enrich_viewer_forMulti1(gene_type=gene_type7(),df = venn_enrich_input1(), Species = input$Species7,
                            Ortholog=ortholog7(),Isoform=isoform7(), org = org7())
  })

  venn_enrich_base2 <- reactive({
    enrich_viewer_forMulti1(gene_type=gene_type7(),df = venn_enrich_input2(), Species = input$Species7,
                            Ortholog=ortholog7(),Isoform=isoform7(), org = org7())
  })

  venn_enrich_viewer2 <- reactive({
    if(input$Species7 != "Xenopus laevis" && input$Ortholog7 != "Arabidopsis thaliana" && input$Species7 != "Arabidopsis thaliana"){
      return(enrich_viewer_forMulti2(gene_type=gene_type7(),df = venn_enrich_input1(), Species = input$Species7,Ortholog = ortholog7(),Isoform=isoform7(), org = org7(),
                                     org_code = org_code7(),H_t2g = venn_Hallmark_set(),Gene_set = input$Gene_set9))
    }else return(enrich_viewer_forMulti2_xenopus(gene_type=gene_type7(),df = venn_enrich_input1(), Species = input$Species7,Ortholog = ortholog7(),Isoform=isoform7(), org = org7(),
                                                 org_code = org_code7(),Gene_set = input$Gene_set9))
  })
  venn_enrich_h <- reactive({
    if(input$Species7 != "Xenopus laevis" && input$Ortholog7 != "Arabidopsis thaliana" && input$Species7 != "Arabidopsis thaliana"){
      return(enrich_gene_list(data = venn_enrich_base1(),
                              Gene_set = input$Gene_set9, org = org7(), H_t2g = venn_Hallmark_set()))
    }else return(enrich_gene_list_xenopus(data = venn_enrich_base1(),
                                          Gene_set = input$Gene_set9, org = org7(),org_code = org_code7()))
  })
  venn_enrich_H <- reactive({
    return(enrich_genelist(data = venn_enrich_base1(),
                           enrich_gene_list = venn_enrich_h(),section = "venn",group_order = input$venn_whichGroup1))
  })
  
  output$venn_enrichment1 <- renderPlot({
    dotplot_for_output(data = venn_enrich_viewer2(),
                       plot_genelist = venn_enrich_H(), Gene_set = input$Gene_set9, 
                       Species = input$Species7)
  })
  
  output$download_venn_cluster_enrichment = downloadHandler(
    filename = function(){
      paste(input$venn_whichGroup1, paste0(input$Gene_set9,"_enrichment.pdf"), sep="_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- multi_enrich_H()
        if(input$venn_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 8
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        dotplot_for_output(data = venn_enrich_viewer2(),
                           plot_genelist = venn_enrich_H(), Gene_set = input$Gene_set9, 
                           Species = input$Species7)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  venn_enrich_table <- reactive({
    if(input$Species7 != "Xenopus laevis" && input$Ortholog7 != "Arabidopsis thaliana" && input$Species7 != "Arabidopsis thaliana"){
      return(enrich_for_table(data = as.data.frame(venn_enrich_viewer2()), H_t2g = venn_Hallmark_set(), Gene_set = input$Gene_set9))
    }else return(as.data.frame(venn_enrich_viewer2()))
  })
  
  output$venn_enrichment_result <- DT::renderDataTable({
    venn_enrich_table()
  })
  
  output$download_venn_enrichment_table = downloadHandler(
    filename = function() {
      paste(input$venn_whichGroup1, paste0(input$Gene_set9,"_enrichment.txt"), sep="_")
    },
    content = function(file){write.table(venn_enrich_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  output$venn_whichGroup2 <- renderUI({
    clusters <- overlap_table2()
    if(is.null(clusters)){
      return(NULL)
    }else{
      selectInput("venn_whichGroup2", "gene list", choices = c("not selected",unique(clusters$Group)),selected = "not selected",multiple = FALSE)
    }
  })
  
  venn_enrich_input2 <- reactive({
    if(is.null(input$venn_whichGroup2) || input$venn_whichGroup2 == "not selected"){
      return(NULL)
    }else{
      clusters <- overlap_table2()
      cluster_name <- input$venn_whichGroup2
      clusterCount <- data.frame(GeneID = NA, Group=NA)
      for(name in cluster_name){
        clusterCount2 <- dplyr::filter(clusters, Group == name)
        clusterCount2 <- data.frame(GeneID = clusterCount2$Gene, Group = clusterCount2$Group)
        clusterCount <- rbind(clusterCount, clusterCount2) 
      }
      clusterCount <- na.omit(clusterCount)
      return(clusterCount)
    }
  })
  
  venn_enrich2 <- reactive({
    cnet_global(data = venn_enrich_base2(), 
                group = input$venn_whichGroup2, enrich_gene_list = venn_enrich_h())
  })
  output$venn_enrichment2 <- renderPlot({
    cnet_for_output(data = venn_enrich_input2(), plot_data = venn_enrich2(), 
                    Gene_set = input$Gene_set9, Species = input$Species7)
  })
  
  output$download_venn_enrichment_cnet = downloadHandler(
    filename = function(){
      paste(input$venn_whichGroup2, paste0(input$Gene_set9,"_cnet.pdf"), sep="_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p <- venn_enrich2()
        if(input$venn_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$venn_pdf_height
        if(input$venn_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$venn_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(p)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  # enrichment viewer ------------------------------------------------------------------------------
  output$Spe3 <- renderText({
    if(input$Species4 == "not selected") print("Please select 'Species'")
  })
  Hallmark_enrich <- reactive({
    return(GeneList_for_enrichment(Species = input$Species4, Ortholog=input$Ortholog4,gene_type=gene_type4(),
                                   Biomart_archive=input$Biomart_archive4, Gene_set = input$Gene_set3, org = org4(), Custom_gene_list = Custom_input()))
  })
  org4 <- reactive({
    return(org(Species = input$Species4,Ortholog = input$Ortholog4))
  })
  isoform4 <- reactive({
    return(isoform_ID(gene_list = enrich_input(),Species = input$Species4,Ortholog = input$Ortholog4,Biomart_archive=input$Biomart_archive4,RNA_type=input$Level_enrich))
  })
  ortholog4 <- reactive({
    return(no_org_ID(gene_list = enrich_input(),Species = input$Species4,Ortholog = input$Ortholog4,Biomart_archive=input$Biomart_archive4))
  })
  gene_type4 <- reactive({
    data <- enrich_input()
    data <- data.frame(row.names = unique(data[,1]))
    return(gene_type(my.symbols=rownames(data),org=org4(),Species=input$Species4,RNA_type=input$Level_enrich))
  })
  org_code4 <- reactive({
    return(org_code(Species = input$Species4, Ortholog= input$Ortholog4))
  })
  
  pre_enrich_input <- reactive({
    upload = list()
    name = c()
    tmp <- NULL
    if(!is.null(input$enrich_data_file[, 1])){
      for(nr in 1:length(input$enrich_data_file[, 1])){
        if(tools::file_ext(input$enrich_data_file[[nr, 'datapath']]) == "xlsx") {
          df <- try(as.data.frame(readxl::read_xlsx(input$enrich_data_file[[nr, 'datapath']])))
        }
        if(tools::file_ext(input$enrich_data_file[[nr, 'datapath']]) == "csv") df <- read.csv(input$enrich_data_file[[nr, 'datapath']], header=TRUE, sep = ",",quote = "")
        if(tools::file_ext(input$enrich_data_file[[nr, 'datapath']]) == "txt" || 
           tools::file_ext(input$enrich_data_file[[nr, 'datapath']]) == "tsv") df <- read.table(input$enrich_data_file[[nr, 'datapath']], header=TRUE, sep = "\t",quote = "")
        file_name <- gsub(paste0("\\.",tools::file_ext(input$enrich_data_file[[nr, 'datapath']]),"$"), "", input$enrich_data_file[nr,]$name)
        name <- c(name, file_name)
        upload[nr] <- list(df)
      }
      names(upload) <- name
      if(length(names(upload)) == 1){
        tmp <- upload[[name]]
        rownames(tmp) = gsub("\"", "", rownames(tmp))
        if(str_detect(colnames(tmp)[1], "^X\\.")){
          colnames(tmp) = sub("^X\\.", "", colnames(tmp))
        }
        if(rownames(tmp)[1] == 1){
          if(dim(tmp)[2] >= 2){
          tmp <- data.frame(Gene = tmp[,1], Group = tmp[,2])
          }else{
            tmp <- data.frame(Gene = tmp[,1], 
                              Group = gsub(paste0("\\.",tools::file_ext(input$enrich_data_file[[1, 'datapath']]),"$"), "", input$enrich_data_file[1,]$name))
          }
        }else{
          tmp <- data.frame(Gene = rownames(tmp), Group = tmp[,1])
        }
      }else{
        df2 <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
        for(file in names(upload)){
          df <- upload[[file]]
          if(rownames(df)[1] == 1){
            df[,1] = gsub("\"", "", df[,1])
            df <- data.frame(Gene = df[,1], Group = file)
          }else{
            rownames(df) = gsub("\"", "", rownames(df))
            df <- data.frame(Gene = rownames(df), Group = file)
          }
          df2 <- rbind(df2,df)
        }
        tmp <- df2
      }
    }else{
      if(input$goButton4 > 0 )  tmp = read_gene_list(example_data_path("data/enrich_example.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/enrich_example.txt"))
    }
    return(tmp)
  })
  output$pre_enrich_input_choice <- renderUI({
    if(!is.null(pre_enrich_input())){
      list <- unique(pre_enrich_input()[,2])
      if(length(list) > 20 && length(input$enrich_data_file[, 1]) == 1){ 
        selectInput("pre_enrich_input_choice","Group name",c("File name","Second column"),selected = "File name",multiple = F)
      }else return(NULL)
    }
  })
  output$enrich_input_choice <- renderUI({
    if(!is.null(pre_enrich_input())){
      list <- unique(pre_enrich_input()[,2])
      if(!is.null(input$pre_enrich_input_choice)){
        if(input$pre_enrich_input_choice == "File name"){
          file_name <- gsub(paste0("\\.",tools::file_ext(input$enrich_data_file[[1, 'datapath']]),"$"), "", input$enrich_data_file[1,]$name)
          list <- file_name
        }
      }
      list <- sort(list)
      selectInput("enrich_input_choice","Order of groups",list, selected = list, multiple = TRUE)
    }
  })
  enrich_input <- reactive({
    list <- pre_enrich_input()
    if(!is.null(list) && !is.null(input$enrich_input_choice)){
      if(!is.null(input$pre_enrich_input_choice)){
        if(input$pre_enrich_input_choice == "File name"){
          file_name <- gsub(paste0("\\.",tools::file_ext(input$enrich_data_file[[1, 'datapath']]),"$"), "", input$enrich_data_file[1,]$name)
          list$Group <- file_name
        }
      }
    group <- factor(input$enrich_input_choice, levels = input$enrich_input_choice)
    list <- list %>% dplyr::filter(Group %in% group)
    return(list)
    }
  })
  
  
  output$enrichment_input <- DT::renderDataTable({
    as.data.frame(enrich_input())
  })
  
  
  enrich_viewer1 <- reactive({
    print(head(isoform4()))
    return(gene_list_convert_for_enrichment(gene_type=gene_type4(),data= enrich_input(), org=org4(),Species = input$Species4,Ortholog=ortholog4(),Isoform=isoform4()))
  })
  
  enrich_viewer2 <- reactive({
    req(input$Gene_set3)
    
    data3 <- enrich_viewer1()
    if (is.null(data3) || nrow(data3) == 0) {
      return(NULL)
    }
    
    withProgress(message = "enrichment analysis", {
      df <- data.frame(matrix(rep(NA, 10), nrow = 1))[numeric(0), ]
      colnames(df) <- c(
        "ID", "Description", "GeneRatio", "BgRatio", "pvalue",
        "p.adjust", "qvalue", "geneID", "Count", "Group"
      )
      
      if (input$Species4 != "Xenopus laevis" &&
          input$Ortholog4 != "Arabidopsis thaliana" &&
          input$Species4 != "Arabidopsis thaliana") {
        H_t2g <- Hallmark_enrich()
        H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
      }
      
      for (name in unique(data3$Group)) {
        ids <- unique(na.omit(data3$ENTREZID[data3$Group == name]))
        if (length(ids) == 0) next
        
        em <- NULL
        
        if (input$Species4 != "Xenopus laevis" &&
            input$Ortholog4 != "Arabidopsis thaliana" &&
            input$Species4 != "Arabidopsis thaliana") {
          
          em <- clusterProfiler::enricher(
            gene = ids,
            TERM2GENE = H_t2g2,
            pvalueCutoff = 0.5
          )
          
        } else {
          if (input$Gene_set3 == "KEGG") {
            em <- clusterProfiler::enrichKEGG(
              gene = ids,
              organism = org_code4(),
              pvalueCutoff = 0.05,
              keyType = "ncbi-geneid"
            )
          }
          
          if (input$Gene_set3 == "GO biological process") {
            em <- clusterProfiler::enrichGO(
              gene = ids,
              OrgDb = org4(),
              ont = "BP",
              pvalueCutoff = 0.05
            )
          }
          
          if (input$Gene_set3 == "GO cellular component") {
            em <- clusterProfiler::enrichGO(
              gene = ids,
              OrgDb = org4(),
              ont = "CC",
              pvalueCutoff = 0.05
            )
          }
          
          if (input$Gene_set3 == "GO molecular function") {
            em <- clusterProfiler::enrichGO(
              gene = ids,
              OrgDb = org4(),
              ont = "MF",
              pvalueCutoff = 0.05
            )
          }
        }
        
        if (is.null(em)) next
        
        em_df <- as.data.frame(em)
        if (nrow(em_df) == 0) next
        
        if ("ENTREZID" %in% AnnotationDbi::keytypes(org4())) {
          em_df <- as.data.frame(clusterProfiler::setReadable(em, org4(), keyType = "ENTREZID"))
        }
        
        em_df$Group <- name
        df <- rbind(df, em_df)
      }
      
      if (nrow(df) == 0) {
        return(NULL)
      }
      
      df$Description <- gsub("HALLMARK_", "", df$Description)
      df$GeneRatio <- DOSE::parse_ratio(df$GeneRatio)
      df
    })
  })
  
  # enrichment plot ------------------------------------------------------------------------------
  enrich_venn <- reactive({
    if(input$Species4 != "Xenopus laevis" && input$Ortholog4 != "Arabidopsis thaliana" && input$Species4 != "Arabidopsis thaliana"){
      return(enrich_gene_list(data = enrich_viewer1(), Gene_set = input$Gene_set3,
                              org = org4(), H_t2g = Hallmark_enrich()))
    }else return(enrich_gene_list_xenopus(data = enrich_viewer1(), Gene_set = input$Gene_set3,
                                          org = org4(), org_code = org_code4()))
  })
  enrich_H <- reactive({
    return(enrich_genelist(data = enrich_viewer1(), enrich_gene_list = enrich_venn(), group_order=input$enrich_input_choice,
                           showCategory = input$enrich_showCategory,section = "enrichmentviewer"))
  })
  
  output$enrichment3 <- renderPlot({
    dotplot_for_output(data = enrich_viewer2(),
                       plot_genelist = enrich_H(), Gene_set = input$Gene_set3, 
                       Species = input$Species4)
  })
  
  enrichGroup <- reactive({
    data3 <- enrich_viewer1()
    if(is.null(data3)){
      return(NULL)
    }else{
      group <- unique(data3$Group)
      return(group)
    }
  })
  
  output$whichGroup <- renderUI({
    if(is.null(enrichGroup())){
      return(NULL)
    }else{
      selectInput("which_group", "Group", choices = c(enrichGroup()), multiple = FALSE)
    }
  })
  
  enrich2 <- reactive({
    cnet_global(data = enrich_viewer1(), group = input$which_group, enrich_gene_list = enrich_venn(),
                showCategory = input$enrich_showCategory)
  })
  
  output$enrichment4 <- renderPlot({
    cnet_for_output(data = enrich_input(), plot_data = enrich2(), 
                    Gene_set = input$Gene_set3, Species = input$Species4)
  })
  
  output$download_enrichment = downloadHandler(
    filename = function(){
      paste(input$enrich_data_file, paste0(input$Gene_set3,".pdf"), sep ="-")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- enrich_H()
        if(input$enrich_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 6.5
        }else pdf_width <- input$enrich_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(p1))
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_enrichment_cnet = downloadHandler(
    filename = function(){
      paste(input$enrich_data_file, paste(input$Gene_set3,paste0(input$which_group,".pdf"), sep = "_"), sep ="-")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- enrich2()
        if(input$enrich_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$enrich_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(p1))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$Gene_set3 <- renderUI({
    if(input$Species4 != "Xenopus laevis" && input$Ortholog4 != "Arabidopsis thaliana" && input$Species4 != "Arabidopsis thaliana"){
      selectInput('Gene_set3', 'Gene Set', gene_set_list)
    }else selectInput('Gene_set3', 'Gene Set', c("KEGG", "GO biological process", 
                                                 "GO cellular component","GO molecular function"))
  })
  Custom_input <- reactive({
    tmp <- input$custom_input$datapath
    data <- read_gene_list(tmp)
    df <- gene_list_convert_for_enrichment(gene_type=gene_type4(),data= data, org=org4(),Species = input$Species4,Ortholog=ortholog4(),Isoform=isoform4())
    return(df)
  })
  output$Custom_input <- renderUI({
    if(is.null(input$Gene_set3)){
      return(NULL)
    }else{
      if(input$Gene_set3 == "Custom gene set"){
        fileInput("custom_input",
                  "Select a file (txt, csv, xlsx)",
                  accept = c("txt", "csv","xlsx"),
                  multiple = FALSE,
                  width = "80%")
      }else return(NULL)
    }
  })
  
  enrich_viewer_table <- reactive({
    if(input$Species4 != "Xenopus laevis" && input$Ortholog4 != "Arabidopsis thaliana" && input$Species4 != "Arabidopsis thaliana"){
      return(enrich_for_table(data = as.data.frame(enrich_viewer2()), H_t2g = Hallmark_enrich(), Gene_set = input$Gene_set3))
    }else return(as.data.frame(enrich_viewer2()))
  })
  
  output$enrichment_result <- DT::renderDataTable({
    enrich_viewer_table()
  })
  
  output$download_enrichment_table = downloadHandler(
    filename = function() {
      paste(input$enrich_data_file, paste0(input$Gene_set3,"_table.txt"), sep ="-")
    },
    content = function(file){write.table(enrich_viewer_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  
  observeEvent(input$goButton4,({
    updateSelectInput(session,inputId = "Species4","Species",species_list, selected = "Mus musculus")
  }))
  
  #motif----------------
  output$motif_Spe <- renderText({
    if(input$Species4 == "not selected") {
      print("Please select 'Species'")
    }else{
      if(input$Species4 != "Homo sapiens" && input$Species4 != "Mus musculus"){
        print("'Homo sapiens' or 'Mus musculus'")
      }
    }
  })
  output$promoter_upstream <- renderUI({
    numericInput("promoter_upstream", "upstream", min   = 0, max   = Inf, value = 400)
  })
  output$promoter_downstream <- renderUI({
    numericInput("promoter_downstream", "downstream", min   = 0, max   = Inf, value = 100)
  })
  output$promoter_padj <- renderUI({
    numericInput("promoter_padj", "pvalue", min   = 0, max   = 0.05, value = 0.05)
  })
  
  promoter <- reactive({
    if(input$motifButton > 0){
      x <- getTargetSeq(Species = input$Species4, upstream = input$promoter_upstream,
                        downstream = input$promoter_downstream)
      return(x)
    }
  })
  updateCounter <- reactiveValues(i = 0)
  
  observe({
    input$motifButton
    isolate({
      updateCounter$i <- updateCounter$i + 1
    })
  })
  
  
  #Restart
  defaultvalues <- observeEvent(enrich_motif(), {
    isolate(updateCounter$i == 0)
    updateCounter$i <- 0
  }) 
  enrich_motif <- reactive({
    if(updateCounter$i > 0 && input$motifButton > 0){
      return(MotifAnalysis(data= enrich_input(), org=org4(),Species = input$Species4, x = promoter()))
    }
  })
  output$motif_plot <- renderPlot({
    if(input$motifButton > 0 && !is.null(enrich_motif())){
      Motifplot(df2 = enrich_motif(), showCategory = input$enrich_showCategory, padj = input$promoter_padj,data= enrich_input(), group_order=input$enrich_input_choice)
    }
  })
  motif_table <- reactive({
    if(input$motifButton > 0 && !is.null(enrich_motif())){
      df2 <- enrich_motif()
      df <- data.frame(matrix(rep(NA, 11), nrow=1))[numeric(0), ]
      for(name in names(df2)){
        res <- df2[[name]]
        res <- dplyr::filter(res, X1 > -log10(input$promoter_padj))
        res <- res %>% dplyr::arrange(-X1.1)
        df <- rbind(df, res)
      }
      colnames(df) <- c("motif.id", "motif.name","motif.percentGC", "negLog10P", "negLog10Padj", "log2enr",
                        "pearsonResid", "expForegroundWgtWithHits", "sumForegroundWgtWithHits", "sumBackgroundWgtWithHits",
                        "Group")
      df$Group <- gsub("\n"," ",df$Group)
      df$padj <- 10^(-df$negLog10Padj)
      return(df)
    }
  })
  output$motif_warning <- renderText({
    if(input$motifButton > 0 && updateCounter$i > 0){
      if(length(motif_table()$motif.id) == 0){
        print("Cannot detect any motifs.")
      }
    }else{
      return(NULL)
    }
  })
  output$motif_result <- DT::renderDT({
    if(input$motifButton > 0 && !is.null(enrich_motif())){
      motif_table()
    }
  })
  
  promoter_motif_region <- reactive({
    target_motif <- motif_table()[input$motif_result_rows_selected,]
    if(input$motifButton > 0 && !is.null(input$motif_result_rows_selected)){
      res <- MotifRegion(data= enrich_input(), target_motif = target_motif,Species = input$Species4, x = promoter())  
      return(res)
    }
  })
  
  output$promoter_motif_region_table <- renderDataTable({
    target_motif <- motif_table()[input$motif_result_rows_selected,]
    if(input$motifButton > 0 && !is.null(input$motif_result_rows_selected)){
      promoter_motif_region()
    }
  })
  
  output$download_motif_table = downloadHandler(
    filename = function() {
      paste(input$enrich_data_file, 
            paste("Motif_upstream",input$promoter_upstream,"_downstream",input$promoter_downstream,"_table.txt",sep = ""), sep ="-")
    },
    content = function(file){write.table(motif_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  output$download_motif_plot = downloadHandler(
    filename = function() {
      paste(input$enrich_data_file, 
            paste("Motif_upstream",input$promoter_upstream,"_downstream",input$promoter_downstream,".pdf",sep = ""), sep ="-")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- Motifplot(df2 = enrich_motif(), showCategory = input$enrich_showCategory, padj = input$promoter_padj,data= enrich_input(), group_order=input$enrich_input_choice)
        if(input$enrich_pdf_height == 0){
          pdf_height <- 6
        }else pdf_height <- input$enrich_pdf_height
        if(input$enrich_pdf_width == 0){
          pdf_width <- 6
        }else pdf_width <- input$enrich_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(p1)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$download_promoter_motif_region = downloadHandler(
    filename = function() {
      paste(input$enrich_data_file, 
            paste("Motif_upstream",input$promoter_upstream,"_downstream",input$promoter_downstream,"_promoter_region.txt",sep = ""), sep ="-")
    },
    content = function(file){write.table(promoter_motif_region(), file, row.names = F, sep = "\t", quote = F)}
  )
  observeEvent(motif_table()[input$motif_result_rows_selected,], ({
    updateCollapse(session,id =  "Promoter_motif_collapse_panel", open="Promoter_motif_region_panel")
  }))
  
  
  #volcano navi------------------------------------------------------
  org5 <- reactive({
    return(org(Species = input$Species5,Ortholog = input$Ortholog5))
  })
  ortholog5 <- reactive({
    return(no_org_ID(count = degresult(),Species = input$Species5,Ortholog = input$Ortholog5,Biomart_archive=input$Biomart_archive5))
  })
  isoform5 <- reactive({
    return(isoform_ID(count = degresult(),Species = input$Species5,Ortholog = input$Ortholog5,Biomart_archive=input$Biomart_archive5,RNA_type=input$Level_volcano))
  })
  gene_type5 <- reactive({
    return(gene_type(my.symbols=rownames(degresult()),org=org5(),Species=input$Species5,RNA_type=input$Level_volcano))
  })
  org_code5 <- reactive({
    return(org_code(Species = input$Species5, Ortholog= input$Ortholog5))
  })
  
  degresult <- reactive({
    withProgress(message = "Importing normalized count matrix, please wait",{
      tmp <- input$deg_file1$datapath
      if(is.null(input$deg_file1) && input$goButton5 > 0 )  tmp = example_data_path("data/DEGexample.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/DEGexample.txt")
      if(!is.null(tmp)){
        df <- read_df(tmp = tmp)
        if (is.null(df)) {
          return(NULL)
        }
        numeric_cols <- intersect(c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), colnames(df))
        for (col_name in numeric_cols) {
          df[[col_name]] <- suppressWarnings(as.numeric(df[[col_name]]))
        }
        if ("padj" %in% colnames(df)) {
          df <- df %>% dplyr::filter(!is.na(padj))
        }
        if ("log2FoldChange" %in% colnames(df)) {
          df <- df %>% dplyr::filter(!is.na(log2FoldChange))
        }
        return(df)
      }
      incProgress(1)
    })
  })
  norm_count_input_for_deg <- reactive({
    withProgress(message = "Importing normalized count matrix, please wait",{
      tmp <- input$deg_file2$datapath
      if(is.null(input$deg_file2) && input$goButton5 > 0 )  tmp = example_data_path("data/day0.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/day0.txt")
      return(read_df(tmp = tmp))
      incProgress(1)
    })
  })
  norm_count_combined_DEG <- reactive({
    count <- norm_count_input_for_deg()
    result <- degresult()
    if(is.null(result)) {
      return(NULL)
    }else{
      sign <- 1
      if(input$volcano_inputType == "reverseON") sign <- -1
      result <- data.frame(row.names = rownames(result),
                           log2FoldChange = sign * result$log2FoldChange, padj = result$padj)
      if(is.null(count)){
        return(result)
      }else{
        data <- merge(result, count, by=0)
        rownames(data) <- data$Row.names
        data <- data[,-1]
        return(data)
      }
    }
    
  })
  
  gene_ID_DEG <- reactive({
    res <- norm_count_combined_DEG()
    if(is.null(res)){
      return(NULL)
    }else{
      if(input$Species5 != "not selected"){
        if(gene_type5() != "SYMBOL"){
          gene_IDs <- ensembl2symbol(gene_type=gene_type5(),data = res,Species=input$Species5,
                                     Ortholog=ortholog5(),Isoform=isoform5(),org = org5(),merge=FALSE)
          return(gene_IDs)
        }
      }else{ return(NULL) }
    }
  })
  
  observeEvent(input$goButton5,({
    updateSelectInput(session,inputId = "Species5","Species",species_list, selected = "Mus musculus")
  }))
  
  observeEvent(input$deg_file1, ({
    updateCollapse(session,id =  "DEG_input_collapse_panel", open="DEG_panel")
  }))
  observeEvent(input$deg_file2, ({
    updateCollapse(session,id =  "DEG_input_collapse_panel", open="normalized_count_matrix_panel")
  }))
  output$DEG_input1 <- DT::renderDataTable({
    degresult()
  })
  output$norm_count_forDEG <- DT::renderDataTable({
    norm_count_input_for_deg()
  })
  
  download_DEG_dir <-reactive({
    dir_name <- paste(
      download_name_stem(input$deg_file1, if (input$goButton5 > 0) "DEGexample" else "DEG"),
      paste(input$fc4, input$fdr4,sep="_"),
      sep = "_"
    )
    return(dir_name)
  })
  output$deg_uniqueID_cut <- renderUI({
    if(gene_type5() != "SYMBOL" && input$Species5 != "not selected") 
      radioButtons("deg_uniqueID_cut","Short unique ID",c("ON"=TRUE,"OFF"=FALSE),selected = TRUE)
  })
  DEG_uniqueID <- reactive({
    count <- norm_count_combined_DEG()
    if(input$Species5 != "not selected"){
      if(gene_type5() != "SYMBOL"){
        gene_IDs  <- gene_ID_DEG()
        data2 <- merge(count, gene_IDs, by= 0)
        rownames(data2) <- data2[,1]
        data2 <- data2[, - which(colnames(data2) == "Row.names.y")]
        data2$Unique_ID <- paste(data2$SYMBOL,data2$Row.names, sep = "\n- ")
        count <- data2[,-1]
      }
    }
    return(count)
  })
  output$deg_GOI_color_type <- renderUI({
    if(input$Species5 == "not selected"){
      radioButtons("deg_GOI_color_type","Color type", c("All genes" = "default"),selected="default")
    }else{
      radioButtons("deg_GOI_color_type","Filter", c("All genes" = "default","Pathway of interest"="pathway"),selected="default")
    }
  })
  observeEvent(input$Species5, ({
    if(input$Species5 == "not selected"){
      updateSelectInput(session, "deg_GOI_color_pathway1","Select a gene set for gene extraction","")
    }else  if(input$Species5 != "Xenopus laevis" && input$Ortholog5 != "Arabidopsis thaliana" && input$Species5 != "Arabidopsis thaliana"){
      updateSelectInput(session, "deg_GOI_color_pathway1","Select a gene set for gene extraction",gene_set_list) 
    }else {
      updateSelectInput(session, "deg_GOI_color_pathway1","Select a gene set for gene extraction",c("KEGG", "GO biological process", 
                                                                                                    "GO cellular component","GO molecular function")) 
    }
  }))
  deg_GOI_color_pathway_list <- reactive({
      list <-deg_pathway_color_gene_list()
      if(!is.null(input$deg_GOI_color_pathway1)){
        if(input$Species5 == "Xenopus laevis" || input$Ortholog5 == "Arabidopsis thaliana" || input$Species5 == "Arabidopsis thaliana"){
          if(input$deg_GOI_color_pathway1 == "GO biological process") list <- list %>% dplyr::filter(Ontology == "BP")
          if(input$deg_GOI_color_pathway1 == "GO cellular component") list <- list %>% dplyr::filter(Ontology == "CC")
          if(input$deg_GOI_color_pathway1 == "GO molecular function") list <- list %>% dplyr::filter(Ontology == "MF") 
        }}
      list <- unique(list$gs_name)
      return(list)
  })
  observeEvent(input$deg_GOI_color_pathway1, ({
    if(is.null(input$deg_GOI_color_pathway1)){
      updateSelectInput(session, "deg_GOI_color_pathway2","","")
    }else{
      updateSelectInput(session, "deg_GOI_color_pathway2","",deg_GOI_color_pathway_list())
    }
  }))
  deg_pathway_color_gene_list <- reactive({
    genes <- GeneList_for_enrichment(Species = input$Species5, Ortholog = input$Ortholog5,
                                     Gene_set=input$deg_GOI_color_pathway1, org = org5(), 
                                     Biomart_archive=input$Biomart_archive5,gene_type=gene_type5())
    return(genes)
  })
  deg_pathway_color_gene <- reactive({
    ##extract pathway genes
    if(is.null(input$deg_GOI_color_pathway1) || is.null(input$deg_GOI_color_pathway2) ||
       input$deg_GOI_color_pathway1 == "" || input$deg_GOI_color_pathway2 == "" || 
       !input$deg_GOI_color_pathway2 %in% deg_GOI_color_pathway_list()) validate("")
    genes <- deg_pathway_color_gene_list()
    genes <- try(dplyr::filter(genes, gs_name == input$deg_GOI_color_pathway2))
    if(length(genes) == 1) if(class(genes)=="try-error") validate("")
    
    my.symbols <- as.character(genes$entrez_gene)
    if(gene_type5() == "non-model organism"){
      gene_IDs <-  try(dplyr::filter(ortholog5(), ENTREZID %in% my.symbols))
      if(length(gene_IDs) == 1) if(class(gene_IDs)=="try-error") validate("")
      df <- data.frame(gene = gene_IDs$ENSEMBL, row.names = gene_IDs$ENSEMBL)
    }else if(gene_type5() == "isoform"){
      gene_IDs <-  try(dplyr::filter(isoform5(), ENTREZID %in% my.symbols))
      if(length(gene_IDs) == 1) if(class(gene_IDs)=="try-error") validate("")
      df <- data.frame(gene = gene_IDs$Transcript_ID, row.names = gene_IDs$Transcript_ID)
    }else{
      if(gene_type5() == "SYMBOL") columns <- c("ENTREZID", sgd_symbol_column(org5())) else columns <- c("ENTREZID","ENSEMBL")
      if(org5()$packageName == "org.At.tair.db") {
        if(gene_type5() == "SYMBOL"){
          gene_IDs <- AnnotationDbi::select(org5(), keys = my.symbols,
                                            keytype = "TAIR",
                                            columns = c("TAIR","SYMBOL"))
          colnames(gene_IDs) <- c("TAIR","GeneID")
        }else{
          gene_IDs <- data.frame(TAIR = my.symbols, GeneID = my.symbols)
        }
      } else {
      gene_IDs <- AnnotationDbi::select(org5(), keys = my.symbols,
                                        keytype = "ENTREZID",
                                        columns = columns)
      colnames(gene_IDs) <- c("entrezid","GeneID")
      }
      gene_IDs <- na.omit(gene_IDs)
      gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
      df <- data.frame(gene = gene_IDs$GeneID, row.names = gene_IDs$GeneID)
      if(dim(df)[1] == 0) validate("No filtered genes.")
    }
    return(df)
  })
  GOI_DEG <- reactive({
    count <- DEG_uniqueID()
    if(is.null(count)){
      return(NULL)
    }else{
      if(is.null(input$deg_GOI_color_type)) validate("")
      if(input$deg_GOI_color_type != "default" && !is.null(deg_pathway_color_gene())){
        count <- count %>% dplyr::filter(row.names(.) %in% rownames(deg_pathway_color_gene()))
      }
      if(gene_type5() != "SYMBOL"){
        if(input$Species5 != "not selected"){
          GOI <- count$Unique_ID
        }else GOI <- rownames(count)
      }else{
        if(input$Species5 != "not selected"){
          GOI <- rownames(count)
        }else GOI <- rownames(count)
      }
      return(normalize_goi_choices(GOI))
    }
  })
  
  observeEvent(GOI_DEG(), {
    choices <- GOI_DEG()
    if(is.null(choices)){
      choices <- character(0)
    }
    selected <- intersect(isolate(input$degGOI), choices)
    updateSelectizeInput(session, "degGOI", choices = choices, selected = selected,
                         options = goi_selectize_options, server = TRUE)
  }, ignoreNULL = FALSE)

  observeEvent(input$GOIreset_deg, {
    choices <- GOI_DEG()
    if(is.null(choices)){
      choices <- character(0)
    }
    updateSelectizeInput(session, "degGOI", choices = choices,
                         selected = character(0),
                         options = goi_selectize_options, server = TRUE)
  })
  output$GOIreset_deg <- renderUI({
    actionButton("GOIreset_deg", "GOI reset")
  })
  
  DEG_uniqueID_base <- reactive({
    data <- DEG_uniqueID()
    if(is.null(data)){
      return(NULL)
    }
    data <- as.data.frame(data)
    if(!nrow(data) || !all(c("log2FoldChange", "padj") %in% colnames(data))){
      return(NULL)
    }
    data
  })

  DEG_uniqueID_plot <- reactive({
    data <- DEG_uniqueID_base()
    if(is.null(data)){
      return(NULL)
    }
    data$padj[data$padj == 0] <- 10^(-300)
    data <- na.omit(data)
    if(!nrow(data)){
      return(NULL)
    }
    data$Row.names <- rownames(data)
    data
  })
  
  output$deg_volcano_x <- renderUI({
    data <- DEG_uniqueID_plot()
    if(is.null(data)){
      return(NULL)
    }
    min <- floor(min(data$log2FoldChange))
    max <- ceiling(max(data$log2FoldChange))
    sliderInput("deg_xrange","X_axis range:",min = min-1,
                max=max+1, step = 0.1,
                value = c(min, max))
  })
  output$deg_volcano_y <- renderUI({
    data <- DEG_uniqueID_plot()
    if(is.null(data)){
      return(NULL)
    }
    max <- ceiling(max(-log10(data$padj)))
    sliderInput("deg_yrange","Y_axis range:",min = 0, max= max+1, step = 1,
                value = max)
  })
  DEGlist_volcano <- reactive({
    data <- DEG_uniqueID_plot()
    if(is.null(data) || is.null(input$deg_xrange) || is.null(input$deg_yrange)){
      return(NULL)
    }
    data$color <- "NS"
    data$color[data$log2FoldChange < -log2(input$fc4) & data$padj < input$fdr4] <- "down"
    data$color[data$log2FoldChange > log2(input$fc4) & data$padj < input$fdr4] <- "up"
    data
  })
  output$uplist_volcano <- renderDataTable({
    data <- DEGlist_volcano()
    if(is.null(data)){
      return(NULL)
    }
    data %>% dplyr::filter(color == "up")
  })
  output$downlist_volcano <- renderDataTable({
    data <- DEGlist_volcano()
    if(is.null(data)){
      return(NULL)
    }
    data %>% dplyr::filter(color == "down")
  })
  deg_volcano <- reactive({
    data <- DEG_uniqueID_plot()
    if(is.null(data) || is.null(input$deg_xrange) || is.null(input$deg_yrange)){
      return(NULL)
    }
      brushed <- brush_info_volcano()
      if(!is.null(brushed) && nrow(brushed) != 0){
        if(gene_type5() != "SYMBOL"){
          if(input$Species5 != "not selected"){
            label_data <- brushed$Unique_ID
          }else{
            label_data <- rownames(brushed)
          }
        }else{
          label_data <- rownames(brushed)
        }
      }else{
        if(!is.null(input$degGOI)){
          label_data <- input$degGOI
        }else label_data <- NULL
      }
      data$padj[data$padj <= 10^(-300)] <- 10^(-300)
      if(input$deg_GOI_color_type == "default"){
      data$color <- "NS"
      data$color[data$log2FoldChange < -log2(input$fc4) & data$padj < input$fdr4] <- "down"
      data$color[data$log2FoldChange > log2(input$fc4) & data$padj < input$fdr4] <- "up"
      if(!is.null(label_data)) {
        Color <- c(down = "blue", GOI = "green", NS = "darkgray", up = "red")
        for(name in label_data){
          if(gene_type5() != "SYMBOL"){
            if(input$Species5 != "not selected"){
              data$color[data$Unique_ID == name] <- "GOI"
            }else{
              data$color[data$Row.names == name] <- "GOI"
            }
          }else{
            data$color[data$Row.names == name] <- "GOI"
          }
        }
        legend_levels <- c("down", "GOI", "NS", "up")
        data$color <- factor(data$color, levels = legend_levels)
      }else{
        Color <- c(down = "blue", NS = "darkgray", up = "red")
        legend_levels <- c("down", "NS", "up")
        data$color <- factor(data$color, levels = legend_levels)
      }
      }else{
        df <- deg_pathway_color_gene()
        data$color <- "others"
        Color <- c(others = "lightgray", down = "blue", NS = "gray1", up = "red")
        for(name in rownames(df)){
          data$color[data$Row.names == name] <- "NS"
          data$color[data$Row.names == name & data$log2FoldChange < -log2(input$fc4) & data$padj < input$fdr4] <- "down"
          data$color[data$Row.names == name & data$log2FoldChange > log2(input$fc4) & data$padj < input$fdr4] <- "up"
        }
        data$color <- factor(data$color, levels = c("others", "NS","down","up"))
        
        if(!is.null(label_data)) {
          Color <- c(others = "lightgray", down = "blue", GOI = "green", NS = "gray1", up = "red")
          for(name2 in label_data){
            if(gsub(".+\\s", "", name2) %in% rownames(df)){
              if(gene_type5() != "SYMBOL"){
                if(input$Species5 != "not selected"){
                  data <- data %>% dplyr::mutate(color=if_else(Unique_ID==name2, "GOI", color))
                }else{
                  data <- data %>% dplyr::mutate(color=if_else(Row.names==name2, "GOI", color))
                }
              }else{
                data <- data %>% dplyr::mutate(color=if_else(Row.names==name2, "GOI", color))
              }
            }
          }
          legend_levels <- c("others", "down", "GOI", "NS", "up")
          data$color <- factor(data$color, levels = legend_levels)
        }else{
          Color <- c(others = "lightgray", down = "blue", NS = "gray1", up = "red")
          legend_levels <- c("others", "down", "NS", "up")
          data$color <- factor(data$color, levels = legend_levels)
        }
      }
      legend_breaks <- legend_levels[legend_levels %in% unique(as.character(stats::na.omit(data$color)))]
      data$minusLog10padj<--log10(data$padj)
      v <- ggplot(data, aes(x = log2FoldChange, y = minusLog10padj)) + ggrastr::geom_point_rast(aes(color = color),size = 0.4)
      v <- v  + geom_vline(xintercept = c(-log2(input$fc4), log2(input$fc4)), linetype = c(2, 2), color = c("black", "black")) +
        geom_hline(yintercept = c(-log10(input$fdr4)), linetype = 2, color = c("black"))
      v <- v +theme_bw()+ scale_color_manual(values = Color, breaks = legend_breaks)+
        theme(legend.position = "top" , legend.title = element_blank(),
              axis.text.x= ggplot2::element_text(size = 12),
              axis.text.y= ggplot2::element_text(size = 12),
              text = ggplot2::element_text(size = 12),
              title = ggplot2::element_text(size = 12)) +
        xlab("log2 fold change") + ylab("-log10(padj)") +
        xlim(input$deg_xrange)+
        ylim(c(0, input$deg_yrange))
      if(input$deg_GOI_color_type == "pathway") {
        v <- v + geom_point(data=dplyr::filter(data, color == "NS"),color="gray1", size=1)
        v <- v + geom_point(data=dplyr::filter(data, color == "up"),color="red", size=1)
        v <- v + geom_point(data=dplyr::filter(data, color == "down"),color="blue", size=1)
      }
      if(!is.null(label_data)) {
        if(gene_type5() != "SYMBOL"){
          if(input$Species5 != "not selected"){
            if(input$deg_uniqueID_cut) {
              id_list <- gsub("\\\n.+$", "", data$Unique_ID)
              dup_list <- unique(id_list[duplicated(id_list)])
              for(i in 1:length(data$Unique_ID)){
                if(gsub("\\\n.+$", "", data$Unique_ID[i]) == "NA") {
                  data$Unique_ID[i] <- gsub(".+\\s", "", data$Unique_ID[i])
                }else if(! gsub("\\\n.+$", "", data$Unique_ID[i]) %in% dup_list) {
                  data$Unique_ID[i] <- gsub("\\\n.+$", "", data$Unique_ID[i])
                }
              }
            }
            v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
            v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Unique_ID),
                                              box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"),
                                              force = 1, fontface = "bold.italic",
                                              bg.color = grDevices::adjustcolor("white", alpha.f = 0.6), bg.r = 0.15)
          }else{
            v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
            v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Row.names),
                                              box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"),
                                              force = 1, fontface = "bold.italic",
                                              bg.color = grDevices::adjustcolor("white", alpha.f = 0.6), bg.r = 0.15)
          }
        }else{
          v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
          v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Row.names),
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"),
                                            force = 1, fontface = "bold.italic",
                                            bg.color = grDevices::adjustcolor("white", alpha.f = 0.6), bg.r = 0.15)
        }
      }
      return(v)
  })
  brush_info_volcano <- reactive({
    data <- DEG_uniqueID_plot()
    if(is.null(data)){
      return(NULL)
    }
    data$minusLog10padj<--log10(data$padj)
    brushedPoints(data, input$plot1_brush_volcano,xvar = "log2FoldChange",yvar="minusLog10padj")
  })
  
  output$deg_volcano1 <- renderPlot({
    if(!is.null(input$deg_xrange)){
      if(is.null(DEG_uniqueID())){
        return(NULL)
      }else{
        deg_volcano()
      }
    }
  })
  
  output$download_volcano_navi = downloadHandler(
    filename = function(){
      paste0(download_DEG_dir(), "_volcano.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$volcano_pdf_height == 0){
          pdf_height <- 5
        }else pdf_height <- input$volcano_pdf_height
        if(input$volcano_pdf_width == 0){
          pdf_width <- 5
        }else pdf_width <- input$volcano_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(deg_volcano())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  
  deg_GOIcount <- reactive({
    count <- norm_count_input_for_deg()
    if(is.null(count)){
      return(NULL)
    }else{
      count <- DEG_uniqueID()
      if(is.null(count)){
        return(NULL)
      }
      brushed <- brush_info_volcano()
      if(is.null(brushed) || dim(brushed)[1] == 0){
        if(is.null(input$degGOI) || length(input$degGOI) == 0){
          return(NULL)
        }
        if(gene_type5() != "SYMBOL"){
          if(input$Species5 != "not selected"){
            Unique_ID <- input$degGOI
            label_data <- as.data.frame(Unique_ID, row.names = Unique_ID)
            data <- merge(count, label_data, by="Unique_ID")
            if(input$deg_uniqueID_cut) {
              id_list <- gsub("\\\n.+$", "", data$Unique_ID)
              dup_list <- unique(id_list[duplicated(id_list)])
              for(i in 1:length(data$Unique_ID)){
                if(gsub("\\\n.+$", "", data$Unique_ID[i]) == "NA") {
                  data$Unique_ID[i] <- gsub(".+\\s", "", data$Unique_ID[i])
                }else if(! gsub("\\\n.+$", "", data$Unique_ID[i]) %in% dup_list) {
                  data$Unique_ID[i] <- gsub("\\\n.+$", "", data$Unique_ID[i])
                }
              }
            }
            rownames(data) <- data$Unique_ID
            data <- data[, - which(colnames(data) == "SYMBOL")]
            data <- data[, - which(colnames(data) == "Unique_ID")]
          }else{
            Row.names <- input$degGOI
            count$Row.names <- rownames(count)
            label_data <- as.data.frame(Row.names, row.names = Row.names)
            data <- merge(count, label_data, by="Row.names")
            rownames(data) <- data$Row.names
            data <- data[, - which(colnames(data) == "Row.names")]
          }
        }else{
          Row.names <- input$degGOI
          count$Row.names <- rownames(count)
          label_data <- as.data.frame(Row.names, row.names = Row.names)
          data <- merge(count, label_data, by="Row.names")
          rownames(data) <- data$Row.names
          data <- data[, - which(colnames(data) == "Row.names")]
        }
      }else{
        data <- brushed
        if(input$deg_GOI_color_type != "default" && !is.null(deg_pathway_color_gene())){
          data <- data %>% dplyr::filter(row.names(.) %in% rownames(deg_pathway_color_gene()))
        }
        data <- data[, - which(colnames(data) == "minusLog10padj")]
        if(gene_type5() != "SYMBOL"){
          if(input$Species5 != "not selected"){
            if(input$deg_uniqueID_cut) {
              id_list <- gsub("\\\n.+$", "", data$Unique_ID)
              dup_list <- unique(id_list[duplicated(id_list)])
              for(i in 1:length(data$Unique_ID)){
                if(gsub("\\\n.+$", "", data$Unique_ID[i]) == "NA") {
                  data$Unique_ID[i] <- gsub(".+\\s", "", data$Unique_ID[i])
                }else if(! gsub("\\\n.+$", "", data$Unique_ID[i]) %in% dup_list) {
                  data$Unique_ID[i] <- gsub("\\\n.+$", "", data$Unique_ID[i])
                }
              }
            }
            rownames(data) <- data$Unique_ID
            data <- data[, - which(colnames(data) == "SYMBOL")]
            data <- data[, - which(colnames(data) == "Unique_ID")]
          }
        }
      }
      data <- data[, - which(colnames(data) == "log2FoldChange")]
      data <- data[, - which(colnames(data) == "padj")]
      numeric_cols <- vapply(data, is.numeric, logical(1))
      data <- data[, numeric_cols, drop = FALSE]
      if (ncol(data) == 0) {
        return(NULL)
      }
      return(data)
    }
  })
  
  DEG_GOIheat <- reactive({
    data <- deg_GOIcount()
    if(is.null(data)){
      ht <- NULL
    }else{
      data.z <- genefilter::genescale(as.matrix(data), axis=1, method="Z")
      data.z <- na.omit(data.z)
      ht <- GOIheatmap(data.z)
    }
    return(ht)
  })
  
  output$deg_GOIheatmap <- renderPlot({
    if(is.null(DEG_uniqueID())){
      return(NULL)
    }else{
      brushed <- brush_info_volcano()
      if(!is.null(input$degGOI) || (!is.null(brushed) && dim(brushed)[1] != 0)){
        withProgress(message = "heatmap",{
          suppressWarnings(print(DEG_GOIheat()))
          incProgress(1)
        })
      }
    }
  })
  
  deg_GOIbox <- reactive({
    data <- deg_GOIcount()
    if(is.null(data)){
      p <- NULL
    }else{
      if(dim(data)[1] == 0) validate("")
      p <- GOIboxplot(data = data)
    }
    return(p)
  })
  
  
  output$deg_GOIbox <- renderPlot({
    if(is.null(DEG_uniqueID())){
      return(NULL)
    }else{
      brushed <- brush_info_volcano()
      if(!is.null(input$degGOI) || (!is.null(brushed) && dim(brushed)[1] != 0)){
        withProgress(message = "Boxplot",{
          suppressWarnings(print(deg_GOIbox()))
          incProgress(1)
        })
      }
    }
  })
  
  output$download_deg_GOIbox = downloadHandler(
    filename = function(){
      paste0(download_DEG_dir(), "GOIboxplot.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- deg_GOIcount()
        rowlist <- rownames(data)
        if(input$volcano_pdf_height == 0){
          pdf_height <- pdf_h(rowlist)
        }else pdf_height <- input$volcano_pdf_height
        if(input$volcano_pdf_width == 0){
          pdf_width <- pdf_w(rowlist)
        }else pdf_width <- input$volcano_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(deg_GOIbox())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_deg_heatmap = downloadHandler(
    filename = function(){
      paste0(download_DEG_dir(), "GOIheatmap.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- deg_GOIcount()
        rowlist <- rownames(data)
        if(input$volcano_pdf_height == 0){
          pdf_height <- 10
        }else pdf_height <- input$volcano_pdf_height
        if(input$volcano_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$volcano_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(DEG_GOIheat())
        dev.off()
        incProgress(1)
      })
    }
  )
  output$download_uplist_volcano = downloadHandler(
    filename = function() {
      paste0(download_DEG_dir(), "_up.txt")
    },
    content = function(file){
      table <- DEGlist_volcano()
      validate(need(!is.null(table), "No DEG data is available."))
      table <- table %>% dplyr::filter(color == "up")
      write.table(table, file, row.names = T,col.names = NA, sep = "\t", quote = F)
      }
  )
  output$download_downlist_volcano = downloadHandler(
    filename = function() {
      paste0(download_DEG_dir(), "_down.txt")
    },
    content = function(file){
      table <- DEGlist_volcano()
      validate(need(!is.null(table), "No DEG data is available."))
      table <- table %>% dplyr::filter(color == "down")
      write.table(table, file, row.names = T,col.names = NA, sep = "\t", quote = F)
    }
  )
  #msigdbr-----------------
  msigdbr_list <- reactive({
    if(input$msigdbr_Species == ""){
      return(NULL)
    }else{
      withProgress(message = "Prepare gene sets",{
        msigdbr_list <- msigdbr(species = input$msigdbr_Species) %>%
          as.data.frame()
        msigdbr_list$gs_name <- gsub("_"," ",msigdbr_list$gs_name)
        return(msigdbr_list)
      })
    }
  })
  
  msig_list <- reactive({
    if(input$msigdbr_Species == "" || input$msigdbr_gene_set == ""){
      return(NULL)
    }else{
      data <- msigdbr_list() %>% 
        dplyr::filter(gs_name == input$msigdbr_gene_set)
      return(data)
    }
  }) 
  
  output$msigdbr_geneset <- renderDataTable({
    msig_list()
  })
  
  observe({
    if(is.null(msigdbr_list())){
      return(NULL)
    }else{
      withProgress(message = "Prepare a list of gene sets",{
        msigdbr_gsname <- unique(msigdbr_list()$gs_name)
        updateSelectizeInput(session, inputId = "msigdbr_gene_set",choices = c("", msigdbr_gsname)) 
      })
    }
  })
  
  output$download_msigdbr_list = downloadHandler(
    filename = function(){
      paste0(input$msigdbr_gene_set, ".txt")
    },
    content = function(file) {write.table(msig_list(), file, quote = F, row.names = F, sep = "\t")})
  
  #Dorothea-----------------
  Dorothea_list <- reactive({
    if(input$dorothea_Species == "not selected" || is.null(Dorothea_type())){
      return(NULL)
    }else{
      withProgress(message = "Prepare gene sets",{
        dorothea_list <- dorothea(species = input$dorothea_Species, type = Dorothea_type(), confidence = input$dorothea_confidence)
        dorothea_list <- dorothea_list[,-2]
        colnames(dorothea_list)[1] <- "TF"
        return(dorothea_list)
      })
    }
  })
  
  Dorothea_type <- reactive({
    if(input$dorothea_TFtype == "activator") return("DoRothEA regulon (activator)")
    if(input$dorothea_TFtype == "repressor") return("DoRothEA regulon (repressor)")
    if(input$dorothea_TFtype == "not selected") return(NULL)
  })
  
  Dorothea_TF_list <- reactive({
    if(input$dorothea_Species == "not selected" || input$dorothea_target_set == ""){
      return(NULL)
    }else{
      data <- Dorothea_list() %>% 
        dplyr::filter(target == input$dorothea_target_set) %>%
        as.data.frame()
      return(data)
    }
  }) 
  
  Dorothea_target_list <- reactive({
    if(input$dorothea_Species == "not selected" || input$dorothea_tf_set == ""){
      return(NULL)
    }else{
      data <- Dorothea_list() %>% 
        dplyr::filter(TF == input$dorothea_tf_set) %>%
        as.data.frame()
      return(data)
    }
  }) 
  
  output$dorothea_tf_list <- renderDataTable({
    Dorothea_TF_list()
  })
  output$dorothea_target_list <- renderDataTable({
    Dorothea_target_list()
  })
  
  observe({
    if(is.null(Dorothea_list())){
      return(NULL)
    }else{
      withProgress(message = "Prepare a list of DoRothEA gene set",{
        dorothea_target <- unique(Dorothea_list()$target)
        updateSelectizeInput(session, inputId = "dorothea_target_set",choices = c("", dorothea_target)) 
        dorothea_tf <- unique(Dorothea_list()$TF)
        updateSelectizeInput(session, inputId = "dorothea_tf_set",choices = c("", dorothea_tf)) 
      })
    }
  })
  
  
  output$download_dorothea_tf_list = downloadHandler(
    filename = function(){
      paste("DoRothEA_TF-search_", input$dorothea_target_set, ".txt", sep = "")
    },
    content = function(file) {write.table(Dorothea_TF_list(), file, quote = F, row.names = F, sep = "\t")})
  
  output$download_dorothea_target_list = downloadHandler(
    filename = function(){
      paste("DoRothEA_Target-search_", input$dorothea_tf_set, ".txt", sep = "")
    },
    content = function(file) {write.table(Dorothea_target_list(), file, quote = F, row.names = F, sep = "\t")})
  
  observeEvent(input$dorothea_tf_set, ({
    updateCollapse(session,id =  "dorothea_collapse_panel", open="dorothea_target_panel")
  }))
  observeEvent(input$dorothea_target_set, ({
    updateCollapse(session,id =  "dorothea_collapse_panel", open="dorothea_tf_panel")
  }))
  
  output$reference_tab <- renderUI({
    req(identical(input$navBar, "reference"))
    includeHTML("www/reference.html")
  })
  
  output$change_log_tab <- renderUI({
    req(identical(input$navBar, "log"))
    includeHTML("www/change-log.html")
  })
  
  
  output$sessionInfo <- renderPrint({
    capture.output(sessionInfo())
  })
  
  ##ensembl2symbol----
  output$Spe_ens <- renderText({
    if(input$Species_ens == "not selected") print("Please select 'Species'")
  })
  org_ens <- reactive({
    return(org(Species = input$Species_ens,Ortholog = input$Ortholog_ens))
  })
  ortholog_ens <- reactive({
    return(no_org_ID(gene_list = data.frame(rownames(input_ens())),Species = input$Species_ens,Ortholog = input$Ortholog_ens,Biomart_archive=input$Biomart_archive_ens))
  })
  gene_type_ens <- reactive({
    data <- input_ens()
    return(gene_type(my.symbols=rownames(data),org=org_ens(),Species=input$Species_ens))
  })
  
  input_ens <- reactive({
    upload = list()
    tmp <- input$data_file_ens$datapath
    if(is.null(input$data_file_ens) && input$goButton_ens > 0 )  tmp = example_data_path("data/sample_ens.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/refs/tags/v1.1.4/data/sample_ens.txt")
    if(!is.null(tmp)){
    if(tools::file_ext(tmp) == "xlsx") {
      df <- try(readxl::read_xlsx(tmp))
      df <- as.data.frame(df)
    }
    if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",",quote = "")
    if(tools::file_ext(tmp) == "txt" || tools::file_ext(tmp) == "tsv") df <- read.table(tmp, header=TRUE, sep = "\t",quote = "")
    if(rownames(df)[1] == 1){
      if(dim(df)[2]==1) df <- data.frame(row.names = df[,1]) else{
        df <- data.frame(row.names = df[,1], gene=df[,2])
      }
    }
    rownames(df) = gsub("\"", "", rownames(df))
    if(dim(df)[2] != 0) {
      if(str_detect(colnames(df)[1], "^X\\.")) colnames(df) = sub("^X\\.", "", colnames(df))
    }else df <- data.frame(geneID = rownames(df),row.names = rownames(df))
    if(sum(!str_detect(rownames(df),"\\.")) == 0) rownames(df) <- gsub("\\..+$", "",rownames(df))
    return(df)
    }
  })
  output$input_ens <- DT::renderDataTable({
    input_ens()
  })
  id_convert_ens <- reactive({
    if(input$Species_ens != "not selected") {
      if(gene_type_ens() == "SYMBOL") validate("The gene names must be ENSEMBL ID.")
      data <- ensembl2symbol(gene_type=gene_type_ens(),
                             data = input_ens(), input$Species_ens,
                             Ortholog=ortholog_ens(),org = org_ens())
      data$Unique_ID <- paste0(data$SYMBOL,"\n- ",rownames(data))
      id_list <- gsub("\\\n.+$", "", data$Unique_ID)
      dup_list <- unique(id_list[duplicated(id_list)])
      for(i in 1:dim(data)[1]){
        if(! gsub("\\\n.+$", "", data$Unique_ID[i]) %in% dup_list) {
          data$Unique_ID[i] <- data$SYMBOL[i]
        }else if(gsub("\\\n.+$", "", data$Unique_ID[i]) == "NA") {
          data$Unique_ID[i] <- rownames(data)[i]
        }else data$Unique_ID[i] <- gsub("\\\n"," ", data$Unique_ID[i])
      }
      return(data)
    }
  })
  output$input_ens2symbol <- DT::renderDataTable({
    id_convert_ens()
  })
  output$download_ens2symbol = downloadHandler(
    filename = function(){
      paste0(sub("\\..+$", "", input$data_file_ens), "_to_symbol.txt")
    },
    content = function(file) {write.table(id_convert_ens(), file, quote = F, row.names = T,col.names = NA, sep = "\t")})
  
  
})
