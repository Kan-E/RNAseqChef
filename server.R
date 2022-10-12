popoverTempate <- 
  '<div class="popover popover-lg" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content"></div></div>'
shinyServer(function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  # pair-wise ------------------------------------------------------------------------------
  org1 <- reactive({
    return(org(Species = input$Species))
  })
  org_code1 <- reactive({
    return(org_code(Species = input$Species))
  })
  
  row_count_matrix <- reactive({
    withProgress(message = "Importing row count matrix, please wait",{
      if (input$data_file_type == "Row1"){
        tmp <- input$file3$datapath
        if(is.null(input$file3) && input$goButton > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example1.txt"
        return(read_df(tmp = tmp))
      }
      if (input$data_file_type == "Row2"){
        tmp <- input$file1$datapath
        if(is.null(input$file1) && input$goButton > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example2.csv"
        return(read_df(tmp = tmp))
      }
    })
  })
  metadata <- reactive({
    if (input$data_file_type != "Row2"){
      return(NULL)
    }else{
      tmp <- input$file2$datapath
      if(is.null(input$file2) && input$goButton > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example3.csv"
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
        upload <- anno_rep(read_df(input$norm_file1[[1, 'datapath']]))
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
  d_row_count_matrix <- reactive({
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
  
  
  # pair-wise DEG ------------------------------------------------------------------------------
  dds <- reactive({
    count <- d_row_count_matrix()
    file_name <- gsub("\\..+$", "", input$file1)
    collist <- gsub("\\_.+$", "", colnames(count))
    if (input$DEG_method == "DESeq2") {
      withProgress(message = "DESeq2",{
        group <- data.frame(con = factor(collist))
        dds<- DESeqDataSetFromMatrix(countData = round(count),colData = group, design = ~ con)
        dds$con <- factor(dds$con, levels = unique(collist))
        dds <- DESeq(dds)
        incProgress(1)
      })
    }
    if (input$DEG_method == "edgeR") {
      withProgress(message = "edgeR",{
        group <- factor(collist)
        dds <- DGEList(counts = count, group = group)
        keep <- filterByExpr(dds)
        dds = dds[keep, , keep.lib.sizes=FALSE]
        dds <- calcNormFactors(dds)
        dds <- estimateCommonDisp(dds)
        dds <- estimateTagwiseDisp(dds)
        incProgress(1)
      })
    }
    return(dds)
  })
  
  deg_result <- reactive({
    if(is.null(d_row_count_matrix())){
      return(NULL)
    }else{
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
          NormMat <- GetNormalizedMat(count, Sizes)
          EBOut <- NULL
          EBOut <- EBTest(Data = count, NgVector = ngvector, Conditions = conditions, sizeFactors = Sizes, maxround = 5)
          stopifnot(!is.null(EBOut))
          PP <- as.data.frame(GetPPMat(EBOut))
          fc_res <- PostFC(EBOut)
          results <- cbind(PP, fc_res$PostFC, fc_res$RealFC,unlist(EBOut$C1Mean)[rownames(PP)], unlist(EBOut$C2Mean)[rownames(PP)])
          colnames(results) <- c("PPEE", "PPDE", "PostFC", "RealFC","C1Mean","C2Mean")
          res <- results[order(results[,"PPDE"], decreasing = TRUE),]
          incProgress(1)
        })
      }
      res <- as.data.frame(res)
      if(input$Species != "not selected"){
        if(str_detect(rownames(count)[1], "ENS")){
          my.symbols <- gsub("\\..*","", rownames(res))
          gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL")
          res$Row.names <- gsub("\\..*","", rownames(res))
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data2 <- merge(res, gene_IDs, by="Row.names")
          rownames(data2) <- data2$Row.names
          res <- data2[,-1]
        }
      }
      return(res)
    }
  })
  
  deg_norm_count <- reactive({
    if(is.null(d_row_count_matrix())){
      return(NULL)
    }else{
      if(!is.null(norm_count_matrix())){
        return(norm_count_matrix())
      }else {
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
            normalized_counts <- as.data.frame(normalized_counts)
            my.symbols <- gsub("\\..*","", rownames(normalized_counts))
            gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                            keytype = "ENSEMBL",
                                            columns = c("ENSEMBL","SYMBOL"))
            colnames(gene_IDs) <- c("Row.names","SYMBOL")
            normalized_counts$Row.names <- gsub("\\..*","", rownames(normalized_counts))
            gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
            data2 <- merge(normalized_counts, gene_IDs, by="Row.names")
            rownames(data2) <- data2$Row.names
            normalized_counts <- data2[,-1]
          }
        }
        return(normalized_counts)
      }
    }
  })
  
  observeEvent(input$goButton,({
    updateSelectInput(session,inputId = "Species","Species",species_list, selected = "Mus musculus")
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
  output$DEG_result <- DT::renderDataTable({
    deg_result()
  })
  output$download_pair_d_row_count = downloadHandler(
    filename = function() {
      if (input$data_file_type == "Row1"){
        paste(gsub("\\..+$", "", input$file3), paste0(input$FDR_method,".txt"), sep ="-")
      }else{
        paste(gsub("\\..+$", "", input$file1), paste0(gsub("\\..+$", "", input$file2),".txt"), sep ="-")
      }},
    content = function(file){write.table(d_row_count_matrix(), file, row.names = T, sep = "\t", quote = F)}
  )
  output$download_pair_DEG_result = downloadHandler(
    filename = function() {
      if (input$data_file_type == "Row1"){
        paste(paste(gsub("\\..+$", "", input$file3), input$DEG_method, sep = "-"), paste0(input$FDR_method,".txt"), sep ="-")
      }else{
        paste(paste(gsub("\\..+$", "", input$file1), input$DEG_method, sep = "-"), paste0(input$FDR_method,".txt"), sep ="-")
      }},
    content = function(file){write.table(deg_result(), file, row.names = T, sep = "\t", quote = F)}
  )
  output$download_pair_norm_count = downloadHandler(
    filename = function() {
      if (input$data_file_type == "Row1"){
        paste(paste(gsub("\\..+$", "", input$file3), input$DEG_method, sep = "-"), "normalized.txt", sep ="-")
      }else{
        paste(paste(gsub("\\..+$", "", input$file1), input$DEG_method, sep = "-"), "normalized.txt", sep ="-")
      }},
    content = function(file){write.table(deg_norm_count(), file, row.names = T, sep = "\t", quote = F)}
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
  gene_ID_pair <- reactive({
    res <- d_row_count_matrix()
    if(is.null(res)){
      return(NULL)
    }else{
      if(input$Species != "not selected"){
        if(str_detect(rownames(res)[1], "ENS")){
          my.symbols <- gsub("\\..*","", rownames(res))
          gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          rownames(gene_IDs) <- gene_IDs$Row.names
          return(gene_IDs)
        }
      }else{ return(NULL) }
    }
  })
  
  
  data_degcount <- reactive({
    data <- deg_result()
    count <- deg_norm_count()
    if(is.null(count) || is.null(data)){
      return(NULL)
    }else{
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
          gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data <- merge(data, gene_IDs, by="Row.names")
          data$Unique_ID <- paste(data$SYMBOL,data$Row.names, sep = " - ")
          genenames <- as.vector(data$SYMBOL)
        }else{
          genenames=NULL
        }
      }else{
        if(input$Species != "not selected"){
          my.symbols <- data$Row.names
          gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                          keytype = "SYMBOL",
                                          columns = c("SYMBOL", "ENTREZID"))
          colnames(gene_IDs) <- c("Row.names", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data <- merge(data, gene_IDs, by="Row.names")
        }
        genenames <- as.vector(data$Row.names)
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
        if(str_detect(rownames(count)[1], "ENS")){
          up_all <- merge(up_all, gene_ID_pair(), by=0)
          rownames(up_all) <- up_all$Row.names
          up_all <- up_all[,-1]
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
        if(str_detect(rownames(count)[1], "ENS")){
          down_all <- merge(down_all, gene_ID_pair(), by=0)
          rownames(down_all) <- down_all$Row.names
          down_all <- down_all[,-1]
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
    content = function(file) {write.table(data_degcount_up(), file, quote = F, row.names = T, sep = "\t")})
  
  output$download_pair_deg_count_down = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), "_DEG_count_down.txt")
    },
    content = function(file) {write.table(data_degcount_down(), file, quote = F, row.names = T, sep = "\t")})
  
  # pair-wise MA ------------------------------------------------------------------------------
  output$MA <- renderPlot({
    withProgress(message = "MA plot and heatmap",{
      if(is.null(d_row_count_matrix())){
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
      data.z <- genescale(data2[,8:(7 + Cond_1 + Cond_2)], axis=1, method="Z")
      ht <- as.grob(GOIheatmap(data.z,show_row_names = FALSE))
    }
    p <- plot_grid(m1, ht, rel_widths = c(2, 1))
    return(p)
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
  GOI_list <- reactive({
      data <- data_degcount()
      count <- deg_norm_count()
      if(is.null(d_row_count_matrix())){
        return(NULL)
      }else{
        if(str_detect(rownames(count)[1], "ENS")){
          if(input$Species != "not selected"){
            GOI <- data$Unique_ID
          }else GOI <- data$Row.names
        }else{
          if(input$Species != "not selected"){
            GOI <- data$Row.names
          }else GOI <- data$Row.names
        }
        return(GOI)
      }
  })
  
  output$GOI <- renderUI({
    if(is.null(deg_norm_count())){
      return(NULL)
    }else{
      withProgress(message = "Preparing GOI list (about 10 sec)",{
      selectizeInput("GOI", "genes of interest (GOI)", c(GOI_list()),multiple = TRUE, 
                     options = list(delimiter = " ", create = T))
      })
    }
  })
  output$GOIreset_pair <- renderUI({
      actionButton("GOIreset_pair", "GOI reset")
  })
  observeEvent(input$GOIreset_pair, {
    withProgress(message = "Preparing GOI list (about 10 sec)",{
    updateSelectizeInput(session, "GOI", choices = c(GOI_list()), 
                         selected = character(0),
                         options = list(create=TRUE, 'plugins' = list('remove_button'), persist = FALSE))
  })
  })
  

  
  output$volcano_x <- renderUI({
    sliderInput("xrange","X_axis range:",min = -100,
                max=100, step = 1,
                value = c(-10, 10))
  })
  output$volcano_y <- renderUI({
    sliderInput("yrange","Y_axis range:",min = 0, max= 300, step = 1,
                value = 100)
  })
  
  pair_volcano <- reactive({
    if(!is.null(input$xrange)){
      data <- as.data.frame(data_degcount())
      count <- deg_norm_count()
      if(!is.null(input$GOI)){
        label_data <- input$GOI
      }else label_data <- NULL
      data$color <- "NS"
      data$color[data$log2FoldChange < -log2(input$fc) & data$padj < input$fdr] <- "down"
      data$color[data$log2FoldChange > log2(input$fc) & data$padj < input$fdr] <- "up"
      data$padj[data$padj == 0] <- 10^(-300)
      if(!is.null(label_data)) {
        Color <- c("blue","green","darkgray","red")
        for(name in label_data){
          if(str_detect(rownames(count)[1], "ENS")){
            if(input$Species != "not selected"){
              data$color[data$Unique_ID == name] <- "GOI"
            }else{
              data$color[data$Row.names == name] <- "GOI"
            }
          }else{
            data$color[data$Row.names == name] <- "GOI"
          }
        }
      }else{
        Color <- c("blue","darkgray","red")
      }
      
      v <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) + geom_point(aes(color = color),size = 0.4)
      v <- v  + geom_vline(xintercept = c(-log2(input$fc), log2(input$fc)), linetype = c(2, 2), color = c("black", "black")) +
        geom_hline(yintercept = c(-log10(input$fdr)), linetype = 2, color = c("black"))
      v <- v +theme_bw()+ scale_color_manual(values = Color)+
        theme(legend.position = "top" , legend.title = element_blank(),
              axis.text.x= ggplot2::element_text(size = 12),
              axis.text.y= ggplot2::element_text(size = 12),
              text = ggplot2::element_text(size = 12),
              title = ggplot2::element_text(size = 12)) +
        xlab("log2 fold change") + ylab("-log10(padj)") +
        xlim(input$xrange)+
        ylim(c(0, input$yrange))
      if(!is.null(label_data)) {
        if(str_detect(rownames(count)[1], "ENS")){
          if(input$Species != "not selected"){
            v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
            v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Unique_ID),
                                              box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), fontface = "bold.italic")
          }else{
            v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
            v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Row.names),
                                              box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), fontface = "bold.italic")
          }
        }else{
          v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
          v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Row.names),
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), fontface = "bold.italic")
        }
      }
      return(v)
    }else return(NULL)
  })
  
  output$volcano1 <- renderPlot({
    if(!is.null(input$xrange)){
      if(is.null(d_row_count_matrix())){
        return(NULL)
      }else{
        withProgress(message = "volcano plot",{
          print(pair_volcano())
          incProgress(1)
        })
      }
    }
  })
  
  output$download_pair_volcano = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), "_volcano.pdf")
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
    
    if(is.null(data2)){
      ht <- NULL
    }else{
      data.z <- genescale(data2[,8:(7 + Cond_1 + Cond_2)], axis=1, method="Z")
      ht <- GOIheatmap(data.z)
    }
    return(ht)
  })
  
  output$GOIheatmap <- renderPlot({
    if(is.null(d_row_count_matrix())){
      return(NULL)
    }else{
      if(!is.null(input$GOI)){
        withProgress(message = "Heatmap",{
          suppressWarnings(print(pair_GOIheatmap()))
          incProgress(1)
        })
      }
    }
  })
  
  output$download_pair_GOIheatmap = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), "_GOIheatmap.pdf")
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
  
  
  
  pair_GOIbox <- reactive({
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
    if(is.null(data2)){
      p <- NULL
    }else{
      data3 <- data2[,8:(7 + Cond_1 + Cond_2)]
      p <- GOIboxplot(data = data3) + scale_fill_manual(values=c("gray", "#ff8082"))
    }
    return(p)
  })
  
  output$GOIbox <- renderPlot({
    if(is.null(d_row_count_matrix())){
      return(NULL)
    }else{
      if(!is.null(input$GOI)){
        withProgress(message = "Boxplot",{
          suppressWarnings(print(pair_GOIbox()))
          incProgress(1)
        })
      }
    }
  })
  
  output$download_pair_GOIbox = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), "_GOIboxplot.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        data <- data_degcount()
        count <- deg_norm_count()
        if(str_detect(rownames(data)[1], "ENS")){
          if(length(grep("SYMBOL", colnames(data))) != 0){
            count <- count[, - which(colnames(count) == "SYMBOL")]
          }
        }
        if(str_detect(rownames(count)[1], "ENS")){
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
  
  # pair-wise PCA ------------------------------------------------------------------------------
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
        print(PCAplot(data = deg_norm_count()))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$PCA <- renderPlot({
    if(is.null(d_row_count_matrix())){
      return(NULL)
    }else{
      print(PCAplot(data = deg_norm_count()))
    }
  })
  
  output$pair_PCA_data <- DT::renderDataTable({
    PCAdata(row_count = d_row_count_matrix(), deg_norm_count = deg_norm_count())
  })
  
  output$download_pair_PCA_table = downloadHandler(
    filename = function() {
      paste0(download_pair_overview_dir(), "_PCA_table.txt")
    },
    content = function(file){write.table(PCAdata(row_count = d_row_count_matrix(), deg_norm_count = deg_norm_count()), file, row.names = T, sep = "\t", quote = F)}
  )
  # pair-wise enrichment ------------------------------------------------------------------------------
  output$Spe1 <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  Hallmark_set <- reactive({
    return(GeneList_for_enrichment(Species = input$Species, Gene_set = input$Gene_set, org = org1()))
  })
  
  enrichment_1_1 <- reactive({
    df <- enrichment_enricher()
    if(!is.null(input$Gene_set) && input$Species != "not selected" && !is.null(df)){
      data <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
      colnames(data) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
      for(name in names(df)){
        if(!is.null(df[[name]])) {
          group1 <- as.data.frame(df[[name]])
        }else group1 <- NULL
        group1$Group <- name
        data <- rbind(data, group1)
      }
            if(length(data$Description) != 0) {
              data["Description"] <- lapply(data["Description"], gsub, pattern="HALLMARK_", replacement = "")
              data$GeneRatio <- parse_ratio(data$GeneRatio)
              return(data)
            }else return(NULL)
    }else{return(NULL)}
  })
  
  enrichment_enricher <- reactive({
    data3 <- data_degcount2()
    if(!is.null(input$Gene_set) && input$Species != "not selected" && !is.null(data3)){
      withProgress(message = "enrichment analysis",{
        H_t2g <- Hallmark_set()
        H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
        em_up <- try(enricher(dplyr::filter(data3, group == "Up")$ENTREZID, TERM2GENE=H_t2g2, pvalueCutoff = 0.05))
        em_down <- try(enricher(dplyr::filter(data3, group == "Down")$ENTREZID, TERM2GENE=H_t2g2, pvalueCutoff = 0.05))
        df <- list()
        df[["Up"]] <- em_up
        df[["Down"]] <- em_down
        for(name in names(df)){
          if (length(as.data.frame(df[[name]])$ID) == 0) {
            df[[name]] <- NULL
          } else{
            df[[name]] <- setReadable(df[[name]], org1(), 'ENTREZID')
          }
        }
        incProgress(1)
        return(df)
      })
    }else{return(NULL)}
  })
  
  enrichment_1_gsea <- reactive({
    data <- data_degcount()
    data3 <- data_degcount2()
    count <- deg_norm_count()
    if(!is.null(input$Gene_set) && input$Species != "not selected" &&
       !is.null(data) && !is.null(data3) && !is.null(count)){
      data <- na.omit(data)
      geneList <- data$log2FoldChange
      names(geneList) = as.character(data$ENTREZID)
      geneList <- sort(geneList, decreasing = TRUE)
        withProgress(message = "GSEA",{
          H_t2g <- Hallmark_set()
          H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
          em3 <- GSEA(geneList, TERM2GENE = H_t2g2,pvalueCutoff = 0.05,exponent = 1, eps = 0, pAdjustMethod = "BH",
                      minGSSize = 50, maxGSSize = 500,by = "fgsea",verbose = F)
          if (length(as.data.frame(em3)$ID) == 0) {
            em4 <- NA
          } else{
            em4 <- setReadable(em3, org1(), 'ENTREZID')
          }
          return(em4)
          incProgress(1)
        })
    }else return(NULL)
  })
  
  
  # pair-wise enrichment plot ------------------------------------------------------------------------------
   pair_enrich1_H <- reactive({
    if(!is.null(input$Gene_set) && input$Species != "not selected"){
        count <- deg_norm_count()
        df <- enrichment_enricher()
        data3 <- data_degcount2()
        H_t2g <- Hallmark_set()
        H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
        if (is.null(df[["Up"]]) && is.null(df[["Down"]]))  {
          p1 <- NULL
        } else{
          data <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          colnames(data) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
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
            data$GeneRatio <- parse_ratio(data$GeneRatio)
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
                              scale_size(range=c(1, 6))+ theme_dose(font.size=12)+ylab(NULL)+xlab(NULL) + 
                              scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
            }}else p1 <- NULL
        }
        em3 <- enrichment_1_gsea()
        if (length(as.data.frame(em3)$ID) == 0) {
          p4 <- NULL
        } else{
          if (length(as.data.frame(em3)$ID) >= 5){
            p4 <- as.grob(gseaplot2(em3, 1:5, pvalue_table = F,base_size = 14))
          }else{
            p4 <- as.grob(gseaplot2(em3, 1:length(em3$ID), pvalue_table = F,base_size = 14))
          }
        }
        p <- plot_grid(p1, p4, nrow = 1)
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
      data <- data_degcount()
      data3 <- data_degcount2()
      count <- deg_norm_count()
      df <- enrichment_enricher()
      upgene <- data3[data3$log2FoldChange > log(input$fc, 2),]
      downgene <- data3[data3$log2FoldChange < log(1/input$fc, 2),]
      p <- list()
      for(name in names(df)){
        if(length(as.data.frame(df[[name]])$ID) == 0){
          cnet1 <- NULL
        } else {
          cnet1 <- setReadable(df[[name]], org1(), 'ENTREZID')
        }
        if (length(as.data.frame(cnet1)$ID) == 0) {
          p2 <- NULL
        } else{
          if(name == "Up") genes <- upgene
          if(name == "Down") genes <- downgene
          geneList <- genes$log2FoldChange
          names(geneList) = as.character(genes$ENTREZID)
          p2 <- try(as.grob(cnetplot(cnet1, foldChange=geneList,
                                     cex_label_gene = 0.7, cex_label_category = 0.75,
                                     cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none")))
          if(length(class(p2)) == 1){
            if(class(p2) == "try-error") p2 <- NULL
          }
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
    selectInput('Gene_set', 'Gene Set', gene_set_list)
  })
  
  pair_enrich_table <- reactive({
    return(enrich_for_table(data = as.data.frame(enrichment_1_1()), H_t2g = Hallmark_set(), Gene_set = input$Gene_set))
  })
  
  output$pair_enrichment_result <- DT::renderDataTable({
    pair_enrich_table()
  })
  
  pair_gsea_table <- reactive({
    data <- as.data.frame(enrichment_1_gsea())
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
        df["day0"] <- list(read.table("https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/iwat-raw/day0.txt",header = T, row.names = 1))
        df["day1"] <- list(read.table("https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/iwat-raw/day1.txt",header = T, row.names = 1))
        df["day5"] <- list(read.table("https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/iwat-raw/day5.txt",header = T, row.names = 1))
        return(df)
      }
      return(NULL)
    }else{
      for(nr in 1:length(input$file11[, 1])){
        df <- read_df(input$file11[[nr, 'datapath']])
        upload[gsub("\\..+$", "", input$file11[nr,]$name)] <- list(df)
      }
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
          keep <- filterByExpr(dds)
          dds = dds[keep, , keep.lib.sizes=FALSE]
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
  
  
  deg_result_batch <- reactive({
    count_files <- batch_files()
    count_list <- list()
    file_num <- length(names(count_files))
    for (name in names(count_files)) {
      if(name != "combined"){
      count <- count_files[[name]]
      collist <- gsub("\\_.+$", "", colnames(count))
      if (input$DEG_method == "DESeq2") {
          dds <- dds_batch()[[name]]
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
      if (input$DEG_method == "edgeR") {
          group <- factor(collist)
          dds <- dds_batch()[[name]]
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
          NormMat <- GetNormalizedMat(count, Sizes)
          EBOut <- NULL
          EBOut <- EBTest(Data = count, NgVector = ngvector, Conditions = conditions, sizeFactors = Sizes, maxround = 5)
          stopifnot(!is.null(EBOut))
          PP <- as.data.frame(GetPPMat(EBOut))
          fc_res <- PostFC(EBOut)
          results <- cbind(PP, fc_res$PostFC, fc_res$RealFC,unlist(EBOut$C1Mean)[rownames(PP)], unlist(EBOut$C2Mean)[rownames(PP)])
          colnames(results) <- c("PPEE", "PPDE", "PostFC", "RealFC","C1Mean","C2Mean")
          res <- results[order(results[,"PPDE"], decreasing = TRUE),]
          incProgress(1)
        })
      }
      res <- as.data.frame(res)
      if(input$Species != "not selected"){
        if(str_detect(rownames(count)[1], "ENS")){
          my.symbols <- gsub("\\..*","", rownames(res))
          gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL")
          res$Row.names <- gsub("\\..*","", rownames(res))
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data2 <- merge(res, gene_IDs, by="Row.names")
          rownames(data2) <- data2$Row.names
          res <- data2[,-1]
        }
      }
      count_list[name] <- list(res)
      }
    }
    return(count_list)
  })
  
  
  
  deg_norm_count_batch <- reactive({
    if(!is.null(norm_count_matrix())){
      files <- norm_count_matrix()
      files["combined"] <- list(batch_count_combined())
      return(files)
    }else{
      count_files <- batch_files()
      count_files["combined"] <- list(batch_count_combined())
      count_list <- list()
      for (name in names(count_files)) {
        count <- count_files[[name]]
        collist <- gsub("\\_.+$", "", colnames(count))
        group <- data.frame(con = factor(collist))
        if (input$DEG_method == "DESeq2") {
          if(name != "combined"){
          dds <- dds_batch()[[name]]
          contrast <- c("con", unique(collist))
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
            keep <- filterByExpr(dds)
            dds = dds[keep, , keep.lib.sizes=FALSE]
            dds <- calcNormFactors(dds, method = "TMM")
            normalized_counts <- cpm(dds)
          }
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
            normalized_counts <- as.data.frame(normalized_counts)
            my.symbols <- gsub("\\..*","", rownames(normalized_counts))
            gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                            keytype = "ENSEMBL",
                                            columns = c("ENSEMBL","SYMBOL"))
            colnames(gene_IDs) <- c("Row.names","SYMBOL")
            normalized_counts$Row.names <- gsub("\\..*","", rownames(normalized_counts))
            gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
            data2 <- merge(normalized_counts, gene_IDs, by="Row.names")
            rownames(data2) <- data2$Row.names
            normalized_counts <- data2[,-1]
          }
        }
        count_list[name] <- list(normalized_counts)
      }
      return(count_list)
    }
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
        base <- matrix_list[[1]]
        int_matrix <- lapply(matrix_list[-1], function(i) base <<- merge(base, i, by = "Row.names"))
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
          if(str_detect(rownames(res)[1], "ENS")){
            my.symbols <- gsub("\\..*","", rownames(res))
            gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                            keytype = "ENSEMBL",
                                            columns = c("ENSEMBL","SYMBOL"))
            colnames(gene_IDs) <- c("Row.names","SYMBOL")
            gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
            rownames(gene_IDs) <- gene_IDs$Row.names
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
            gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                            keytype = "ENSEMBL",
                                            columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
            colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
            gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
            data <- merge(data, gene_IDs, by="Row.names")
            data$Unique_ID <- paste(data$SYMBOL,data$Row.names, sep = " - ")
            genenames <- as.vector(data$SYMBOL)
          }else{
            genenames=NULL
          }
        }else{
          if(input$Species != "not selected"){
            my.symbols <- data$Row.names
            gene_IDs<-AnnotationDbi::select(org1(),keys = my.symbols,
                                            keytype = "SYMBOL",
                                            columns = c("SYMBOL", "ENTREZID"))
            colnames(gene_IDs) <- c("Row.names", "ENTREZID")
            gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
            data <- merge(data, gene_IDs, by="Row.names")
          }
          genenames <- as.vector(data$Row.names)
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
        data2 <- as.data.frame(data_degcount2_batch()[name])
        colnames(data2) <- sub(paste0(name,"."), "", colnames(data2))
        up_all <- dplyr::filter(data2, log2FoldChange > 0)
        rownames(up_all) <- up_all$Row.names
        up_all <- up_all[,8:(7 + Cond_1 + Cond_2)]
        if(input$Species != "not selected"){
          if(str_detect(rownames(count)[1], "ENS")){
            up_all <- merge(up_all, gene_ID_pair(), by=0)
            rownames(up_all) <- up_all$Row.names
            up_all <- up_all[,-1]
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
        data2 <- as.data.frame(data_degcount2_batch()[name])
        colnames(data2) <- sub(paste0(name,"."), "", colnames(data2))
        down_all <- dplyr::filter(data2, log2FoldChange < 0)
        rownames(down_all) <- down_all$Row.names
        down_all <- down_all[,8:(7 + Cond_1 + Cond_2)]
        if(input$Species != "not selected"){
          if(str_detect(rownames(count)[1], "ENS")){
            down_all <- merge(down_all, gene_ID_pair(), by=0)
            rownames(down_all) <- down_all$Row.names
            down_all <- down_all[,-1]
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
        m1 <- as.grob(ggmaplot(data, fdr = input$fdr, fc = input$fc, size = 0.1,
                               palette = c("#B31B21", "#1465AC", "darkgray"),
                               genenames = genenames,
                               legend = "top", top = 20,
                               font.label = c("bold", 6),font.legend = "bold",
                               font.main = "bold",xlab = xlab,
                               ggtheme = ggplot2::theme_minimal(),
                               select.top.method = "fc"))
        data2 <- as.data.frame(data_degcount2_batch()[name])
        colnames(data2) <- sub(paste0(name,"."), "", colnames(data2))
        if(is.null(data2)){
          ht <- NULL
        }else{
          data.z <- genescale(data2[,8:(7 + Cond_1 + Cond_2)], axis=1, method="Z")
          ht <- as.grob(Heatmap(data.z, name = "z-score",column_order = colnames(data.z),
                                clustering_method_columns = 'ward.D2',
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
      p2 <- PCAplot(data = data)
      pca_list[[name]] <- p2
    }
    return(pca_list)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("pair-wise_batch", "zip", sep=".")
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
        deg_r <- deg_result_batch()
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
            write.table(deg, DEG, quote = F, row.names = T, sep = "\t")
            write.table(up, UP, quote = F, row.names = T, sep = "\t")
            write.table(down, DOWN, quote = F, row.names = T, sep = "\t")
            pdf(MA, height = 4, width = 7)
            print(ma_r[[name]])
            dev.off() 
            }
            write.table(norm, norm_count, quote = F, row.names = T, sep = "\t")
            pdf(PCA, height = 3.5, width = 9)
            print(pca_r[[name]])
            dev.off()
        }
        zip(zipfile=fname, files=fs)
      })
    },
    contentType = "application/zip"
  )
  
  #multi DEG--------------------------
  org6 <- reactive({
    return(org(Species = input$Species6))
  })
  org_code6 <- reactive({
    return(org_code(Species = input$Species6))
  })
  
  multi_row_count_matrix <- reactive({
    withProgress(message = "Importing row count matrix, please wait",{
      if (input$multi_data_file_type == "Row1"){
        tmp <- input$multi_file1$datapath
        if(is.null(input$multi_file1) && input$goButton6 > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example4.txt"
        return(read_df(tmp = tmp))
      }
      if (input$multi_data_file_type == "Row2"){
        tmp <- input$multi_file2$datapath
        if(is.null(input$multi_file2) && input$goButton6 > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example8.txt"
        return(read_df(tmp = tmp))
      }
    })
  })
  multi_metadata <- reactive({
    if (input$multi_data_file_type != "Row2"){
      return(NULL)
    }else{
      tmp <- input$multi_file3$datapath
      if(is.null(input$multi_file3) && input$goButton6 > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example5.csv"
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
        rownames(data) <- data[,1]
        data2 <- data[, - which(colnames(data) %in% c("Row.names", colname))]
        data2 <- data2[,1:length(rownames(row))]
        data2_t <- t(data2)
        data3 <- apply(data2_t, 2, as.numeric)
        rownames(data3) <- rownames(data2_t)
        if(input$Species6 != "not selected"){
          if(str_detect(rownames(data3)[1], "ENS")){
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
  multi_dds <- reactive({
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
  
  multi_deg_result <- reactive({
    if(is.null(multi_d_row_count_matrix())){
      return(NULL)
    }else{
      withProgress(message = "Prepare a DEG result",{
        count <- multi_d_row_count_matrix()
        meta <- multi_metadata()
        dds <- multi_dds()
        df <- list()
        for(i in 1:choose(n=length(unique(dds$meta)),k=2)){
          res <- results(dds, contrast = c("meta", as.character(unique(dds$meta)[combn(x=length(unique(dds$meta)),m=2)[1,i]]),as.character(unique(dds$meta)[combn(x=length(unique(dds$meta)),m=2)[2,i]])))
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
          ihw_res <- ihw(pvalue ~ baseMean,  data=as.data.frame(res), alpha = 0.1)
          res$padj <- IHW::as.data.frame(ihw_res)$adj_pvalue
        }
        if(input$FDR_method6 == "Qvalue") {
          res <- results(dds)
          qvalue <- qvalue::qvalue(res$pvalue)
          res$padj <- qvalue$qvalues
        }
        res <- data.frame(gene=rownames(res), padj=res$padj)
        fc_matrix <- merge(fc_matrix, res, by="gene")
        rownames(fc_matrix)<- fc_matrix$gene
        fc_matrix <- fc_matrix[,-1]
        
        res <- fc_matrix
        if(input$Species6 != "not selected"){
          if(str_detect(rownames(count)[1], "ENS")){
            my.symbols <- gsub("\\..*","", rownames(res))
            gene_IDs<-AnnotationDbi::select(org6(),keys = my.symbols,
                                            keytype = "ENSEMBL",
                                            columns = c("ENSEMBL","SYMBOL"))
            colnames(gene_IDs) <- c("Row.names","SYMBOL")
            res$Row.names <- gsub("\\..*","", rownames(res))
            gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
            data2 <- merge(res, gene_IDs, by="Row.names")
            rownames(data2) <- data2$Row.names
            res <- data2[,-1]
          }
        }
        return(res)
      })
    }
  })
  
  multi_deg_norm_count <- reactive({
    if(is.null(multi_d_row_count_matrix())){
      return(NULL)
    }else{
      if(!is.null(multi_norm_count_matrix())){
        return(multi_norm_count_matrix())
      }else {
        count <- multi_d_row_count_matrix()
        dds <- multi_dds()
        normalized_counts <- counts(dds, normalized=TRUE)
        
        if(input$Species6 != "not selected"){
          if(str_detect(rownames(count)[1], "ENS")){
            normalized_counts <- as.data.frame(normalized_counts)
            my.symbols <- gsub("\\..*","", rownames(normalized_counts))
            gene_IDs<-AnnotationDbi::select(org6(),keys = my.symbols,
                                            keytype = "ENSEMBL",
                                            columns = c("ENSEMBL","SYMBOL"))
            colnames(gene_IDs) <- c("Row.names","SYMBOL")
            normalized_counts$Row.names <- gsub("\\..*","", rownames(normalized_counts))
            gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
            data2 <- merge(normalized_counts, gene_IDs, by="Row.names")
            rownames(data2) <- data2$Row.names
            normalized_counts <- data2[,-1]
          }
        }
        return(normalized_counts)
      }
    }
  })
  
  observeEvent(input$goButton6,({
    updateSelectInput(session,inputId = "Species6","Species",species_list, selected = "Mus musculus")
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
    content = function(file){write.table(multi_deg_result(), file, row.names = T, sep = "\t", quote = F)}
  )
  output$download_multi_norm_count = downloadHandler(
    filename = function() {
      if (input$multi_data_file_type == "Row1"){
        paste(gsub("\\..+$", "", input$multi_file1), "normalized.txt", sep ="-")
      }else{
        paste(gsub("\\..+$", "", input$multi_file2), "normalized.txt", sep ="-")
      }},
    content = function(file){write.table(multi_deg_norm_count(), file, row.names = T, sep = "\t", quote = F)}
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
  output$multi_DEG_total1 <- renderText({
    if(is.null(multi_pattern1())){
      return(NULL)
    }else{ 
      print(paste0("The DEG number after the filtration: ", length(multi_pattern1()$gene)))
    }
  })
  
  multi_gene_ID_pair <- reactive({
    res <- multi_d_row_count_matrix()
    if(is.null(res)){
      return(NULL)
    }else{
      if(input$Species6 != "not selected"){
        if(str_detect(rownames(res)[1], "ENS")){
          my.symbols <- gsub("\\..*","", rownames(res))
          gene_IDs<-AnnotationDbi::select(org6(),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          rownames(gene_IDs) <- gene_IDs$Row.names
          return(gene_IDs)
        }
      }else{ return(NULL) }
    }
  })
  
  multi_select <- reactive({
    dds <- multi_dds()
    return(unique(dds$meta))
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
    count <- multi_d_row_count_matrix()
    meta <- multi_metadata()
    dds <- multi_dds()
    if(is.null(dds) || length(input$selectFC) != 2){
      return(NULL)
    }else{
      withProgress(message = "Fold Change, FDR, and base mean cut-off",{
        if (input$multi_data_file_type == "Row1"){
          collist <- gsub("\\_.+$", "", colnames(count))
          meta <- data.frame(condition = factor(collist))
        }else meta <- data.frame(condition=factor(meta[,1]), type=factor(meta[,2]))
        res <- results(dds, contrast = c("meta", as.character(input$selectFC[1]),as.character(input$selectFC[2])))
        
        data <- as.data.frame(multi_deg_norm_count())
        collist <- gsub("\\_.+$", "", colnames(data))
        if(input$Species6 != "not selected"){
          if(str_detect(rownames(data)[1], "ENS")) data <- data[, - which(colnames(data) == "SYMBOL")]
        }
        data <- dplyr::filter(data, apply(data,1,mean) > input$basemean6)
        res <- as.data.frame(res)
        data2 <- merge(res,data, by=0)
        res <- data2[,1:7]
        rownames(res) <- res$Row.names
        
        sig_res_LRT <- res %>% 
          data.frame() %>%
          rownames_to_column(var="gene") %>%
          as_tibble() %>%
          filter(padj < input$fdr6) %>%
          filter(abs(log2FoldChange) > log2(input$fc6))
        return(sig_res_LRT)
      })
    }
  })
  
  multi_pattern1_2 <- reactive({
    sig_res_LRT <- multi_pattern1()
    dds <- multi_dds()
    if(length(sig_res_LRT$gene) == 0){
      return(NULL)
    }else{
      withProgress(message = "Select most significant genes",{
        clustering_sig_genes <- sig_res_LRT %>%
          arrange(padj) %>%
          head(n=input$topP)
        rld <- rlogTransformation(dds)
        rld_mat <- assay(rld)
        cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
        return(cluster_rlog)
      })
    }
  })
  
  
  multi_pattern2 <- reactive({
    count <- multi_d_row_count_matrix()
    meta <- multi_metadata()
    cluster_rlog <- multi_pattern1_2()
    if(is.null(cluster_rlog)|| length(input$selectFC) != 2){
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
          p <- degPlotCluster(table, time = "condition", process = TRUE)+ 
            scale_color_brewer(palette = "Set1", direction=-1)+
            theme_bw(base_size = 15)+ theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
        }
        else{
          p <- degPlotCluster(table, time = "condition", color=colnames(meta)[2], process = TRUE)+ 
            scale_color_brewer(palette = "Set1", direction=-1)+
            theme_bw(base_size = 15)+ theme(legend.position = "top")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
        }
        return(p)
      })
    }
  })
  
  output$multi_boxplot <- renderPlot({
    if(is.null(multi_boxplot_reactive())|| length(input$selectFC) != 2){
      return(NULL)
    }else{
      print(multi_boxplot_reactive())
    }
  })
  
  output$multi_select_file1 <- renderUI({
    clusters <- multi_pattern2()$df
    if(is.null(clusters)){
      return(NULL)
    }else{
      clusters$cluster <- paste0("Group",clusters$cluster)
      selectInput("multi_selectfile1", "cluster_list", choices = c(unique(clusters$cluster)), multiple = F)
    }
  })
  
  multi_pattern_extract <- reactive({
    data <- as.data.frame(multi_deg_norm_count())
    clusters <- multi_pattern2()$df
    if(is.null(data) || is.null(clusters)){
      return(NULL)
    }else{
      if(input$multi_selectfile1 == "not selected" || is.null(input$multi_selectfile1)){
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
        if(str_detect(rownames(data)[1], "ENS")){
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
      paste0(download_multi_overview_dir(), "_DEG_count_up.txt")
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
    content = function(file) {write.table(multi_pattern_extract(), file, quote = F, row.names = T, sep = "\t")})
  
  multi_GOIbox <- reactive({
    if(!is.null(input$multi_pattern1_count_rows_selected)){
    data <- multi_pattern_extract()[input$multi_pattern1_count_rows_selected,]
    if(input$Species6 != "not selected"){
      if(str_detect(rownames(data)[1], "ENS")){
        rownames(data) <- paste(data$SYMBOL,rownames(data), sep = " - ")
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
        if(is.null(multi_boxplot_reactive)|| length(input$selectFC) != 2){
          return(NULL)
        }else{
          clusters <- multi_pattern2()$df
          clusterNumber <- length(unique(clusters$cluster))
          print(clusterNumber)
          if(input$multi_pdf_height == 0){
            pdf_height <- pdf_h(clusterNumber)
          }else pdf_height <- input$multi_pdf_height
          if(input$multi_pdf_width == 0){
            pdf_width <- pdf_w(clusterNumber)
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
  output$selectFC2 <- renderUI({
    if(is.null(multi_deg_result())){
      return(NULL)
    }else{
      selectizeInput("selectFC2", "Select a pair for fold change cut-off", c(multi_select()),
                     selected = "", multiple = TRUE, 
                     options = list(maxItems = 2))
    }
  })
  output$topP2 <- renderUI({
    if(is.null(multi_deg_result())){
      return(NULL)
    }else{
      sliderInput("topP2", "Most significant genes", min = 1,
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
    if(is.null(multi_deg_count())){
      return(NULL)
    }else{ 
      print(paste0("The DEG number after the filtration: ", length(rownames(multi_deg_count1()))))
    }
  })
  
  multi_deg_count1 <- reactive({
    data <- as.data.frame(multi_deg_norm_count())
    count <- multi_d_row_count_matrix()
    meta <- multi_metadata()
    dds <- multi_dds()
    if(is.null(dds) || length(input$selectFC2) != 2){
      return(NULL)
    }else{
      if (input$multi_data_file_type == "Row1"){
        collist <- gsub("\\_.+$", "", colnames(count))
        meta <- data.frame(condition = factor(collist))
      }else meta <- data.frame(condition=factor(meta[,1]), type=factor(meta[,2]))
      res <- results(dds, contrast = c("meta", as.character(input$selectFC2[1]),as.character(input$selectFC2[2])))
      
      collist <- gsub("\\_.+$", "", colnames(data))
      if(input$Species6 != "not selected"){
        if(str_detect(rownames(data)[1], "ENS")) data <- data[, - which(colnames(data) == "SYMBOL")]
      }
      data <- dplyr::filter(data, apply(data,1,mean) > input$basemean6)
      res <- as.data.frame(res)
      data2 <- merge(res,data, by=0)
      res <- data2[,1:7]
      rownames(res) <- res$Row.names
      
      sig_res_LRT <- res %>% 
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble() %>%
        filter(padj < input$fdr6) %>%
        filter(abs(log2FoldChange) > log2(input$fc6))
      return(sig_res_LRT)
    }
  })
  
  multi_deg_count <- reactive({
    data <- as.data.frame(multi_deg_norm_count())
    sig_res_LRT <- multi_deg_count1()
    if(is.null(dds) || length(input$selectFC2) != 2){
      return(NULL)
    }else{
      sig_res_LRT <- sig_res_LRT %>%
        arrange(padj) %>%
        head(n=input$topP2)
      sig_res_LRT <- as.data.frame(sig_res_LRT)
      rownames(sig_res_LRT) <- sig_res_LRT$gene
      sig_res_LRT <- sig_res_LRT[,-1]
      
      if(input$Species6 != "not selected"){
        if(str_detect(rownames(data)[1], "ENS")) data <- data[, - which(colnames(data) == "SYMBOL")]
      }
      collist <- gsub("\\_.+$", "", colnames(data))
      data <- dplyr::filter(data, apply(data,1,mean) > input$basemean6)
      data2 <- merge(data, sig_res_LRT,by=0)
      rownames(data2) <- data2$Row.names
      data2<- data2[,-1]
      data2 <- data2[,1:length(collist)]
      return(data2)
    }
  })
  
  multi_data_z <- reactive({
    data <- multi_deg_count()
    if(is.null(data)){
      return(NULL)
    }else{
      data.z <- genescale(data, axis = 1, method = "Z")
      data.z <- na.omit(data.z)
      return(data.z)
    }
  })
  
  multi_kmeans <- reactive({
    data.z <- multi_data_z()
    if(is.null(data.z)){
      return(NULL)
    }else{
      withProgress(message = "k-means clustering",{
        ht <- Heatmap(data.z, name = "z-score",
                      column_order = colnames(data.z),
                      clustering_method_columns = 'ward.D2',
                      row_km= input$multi_kmeans_number, cluster_row_slices = F, row_km_repeats = 100,
                      show_row_names = F,column_names_side = "top")
        ht <- draw(ht)
        return(ht)
      })
    }
  })
  
  multi_kmeans_cluster <- reactive({
    ht <- multi_kmeans()
    data.z <- multi_data_z()
    data <- multi_deg_count()
    if(is.null(ht) || is.null(data.z)){
      return(NULL)
    }else{
      r.dend <- row_dend(ht)
      rcl.list <- row_order(ht)
      lapply(rcl.list, function(x) length(x))
      Cluster <- NULL
      if(!is.null(input$multi_kmeans_number)){
        if(length(lapply(rcl.list, function(x) length(x))) != input$multi_kmeans_number){
          return(NULL)
        }else{
          for (i in 1:length(row_order(ht))){ if (i == 1) {
            clu <- t(t(row.names(data.z[row_order(ht)[[i]],])))
            out <- cbind(clu, paste("cluster", i, sep=""))
            colnames(out) <- c("GeneID", "Cluster")} else {
              clu <- t(t(row.names(data.z[row_order(ht)[[i]],])))
              clu <- cbind(clu, paste("cluster", i, sep=""))
              out <- rbind(out, clu)}}
          out <- as.data.frame(out)
          rownames(out) <- out$GeneID
          clusterCount <- merge(out, data, by=0)
          rownames(clusterCount) <- clusterCount$GeneID
          clusterCount <- clusterCount[,-1:-2]
          return(clusterCount)
        }
      }else return(NULL)
    }
  })
  
  output$multi_kmeans_heatmap <- renderPlot({
    ht <- multi_kmeans()
    if(is.null(ht)){
      return(NULL)
    }else{
      print(ht)
    }
  })
  
  multi_kmeans_box <- reactive({
    res <- multi_kmeans_cluster()
    ma <- as.data.frame(multi_deg_count())
    meta <- multi_metadata()
    if(is.null(ma) || is.null(res)){
      return(NULL)
    }else{
      withProgress(message = "Boxplot",{
        if (input$multi_data_file_type == "Row1"){
          collist <- gsub("\\_.+$", "", colnames(ma))
          meta <- data.frame(condition = factor(collist),row.names = colnames(ma))
        }else {
          meta <- data.frame(condition=factor(meta[,1]), type=factor(meta[,2]), row.names = colnames(ma))
          meta[,2] <- factor(meta[,2], levels = unique(meta[,2]),ordered = TRUE)
        }
        meta$condition <- factor(meta$condition, levels = unique(meta$condition),ordered = TRUE)  
        res$Cluster <- gsub("cluster","",res$Cluster)
        res<- data.frame(genes=rownames(res),cluster=res$Cluster)
        table <- rownames_to_column(as.data.frame(ma), "genes") %>%
          gather("sample", "expression", -genes) %>%
          right_join(distinct(res[,c("genes", "cluster")]),
                     by = "genes") %>%
          left_join(rownames_to_column(as.data.frame(meta), "sample"),
                    by = "sample") %>%
          as.data.frame()
        table$cluster = as.integer(table$cluster)
        table<-na.omit(table)
        if (input$multi_data_file_type == "Row1"){
          p <- degPlotCluster(table, time = "condition", process = TRUE)+ 
            scale_color_brewer(palette = "Set1", direction=-1)+
            theme_bw(base_size = 15)+ theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
        }
        else{
          p <- degPlotCluster(table, time = "condition", color=colnames(meta)[2], process = TRUE)+ 
            scale_color_brewer(palette = "Set1", direction=-1)+
            theme_bw(base_size = 15)+ theme(legend.position = "top")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
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
          print(clusterNumber)
          if(input$multi_pdf_height == 0){
            pdf_height <- pdf_h(clusterNumber)
          }else pdf_height <- input$multi_pdf_height
          if(input$multi_pdf_width == 0){
            pdf_width <- pdf_w(clusterNumber)
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
    data <- as.data.frame(multi_deg_norm_count())
    if(input$Species6 != "not selected"){
      if(str_detect(rownames(data)[1], "ENS")){
        data <- data.frame(SYMBOL=data$SYMBOL, row.names = rownames(data))
        clusters <-merge(clusters,data, by=0)
        colnames(clusters)[1] <- "genes"
      }
    }
    clusters
  })
  
  output$download_multi_kmeans_cluster = downloadHandler(
    filename = function() {
      paste0(download_multi_overview_dir(), "kmeans_count_table.txt")
    },
    content = function(file){write.table(multi_kmeans_cluster(), file, row.names = T, sep = "\t", quote = F)}
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
        pdf(file, height = pdf_height, width = pdf_width)
        print(multi_kmeans())
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
      selectInput("multi_selectfile2", "cluster_list", choices = c(unique(clusters$Cluster)), multiple = F)
    }
  })
  
  multi_kmeans_pattern_extract <- reactive({
    count <- multi_d_row_count_matrix()
    clusters <- multi_kmeans_cluster()
    if(is.null(clusters)){
      return(NULL)
    }else{
      if(is.null(input$multi_selectfile2)){
        return(NULL)
      }else{
        cluster_name <- input$multi_selectfile2
        clusterCount <- dplyr::filter(clusters, Cluster == cluster_name)
        clusterCount <- clusterCount[,-1]
        return(clusterCount)
      }
    }
  })
  
  multi_kmeans_extract<- reactive({
    if(is.null(input$multi_selectfile2)){
      return(NULL)
    }else{
      clusters <- multi_kmeans_pattern_extract()
      data <- as.data.frame(multi_deg_norm_count())
      if(input$Species6 != "not selected"){
        if(str_detect(rownames(data)[1], "ENS")){
          data <- data.frame(SYMBOL=data$SYMBOL, row.names = rownames(data))
          clusters <- merge(clusters,data, by=0)
          rownames(clusters) <- clusters$Row.names
          clusters <- clusters[,-1]
        }
      }
      return(clusters)
    }
  })
  
  output$multi_pattern2_count <- DT::renderDT({
    if(is.null(input$multi_selectfile2)){
      return(NULL)
    }else{
      clusters <- multi_kmeans_pattern_extract()
      data <- as.data.frame(multi_deg_norm_count())
      if(input$Species6 != "not selected"){
        if(str_detect(rownames(data)[1], "ENS")){
          data <- data.frame(SYMBOL=data$SYMBOL, row.names = rownames(data))
          clusters <- merge(clusters,data, by=0)
          colnames(clusters)[1] <- "genes"
        }
      }
      clusters
    }
  })
  
  
  multi_kmeans_GOIbox <- reactive({
    if(!is.null(input$multi_pattern2_count_rows_selected)){
      data <- multi_kmeans_extract()[input$multi_pattern2_count_rows_selected,]
      if(input$Species6 != "not selected"){
        if(str_detect(rownames(data)[1], "ENS")){
          rownames(data) <- paste(data$SYMBOL,rownames(data), sep = " - ")
          data <- data[, - which(colnames(data) == "SYMBOL")]
        }
      }
      return(data)
    }
  })
  
  output$multi_kmeans_GOIboxplot <- renderPlot({
    if(!is.null(input$multi_pattern2_count_rows_selected)){
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
  
  
  observeEvent(input$multi_pattern2_count_rows_selected, ({
    updateCollapse(session,id =  "multi_collapse_panel2", open="multi_deg_kmeans_boxplot_panel")
  }))
  
  
  output$download_deg_kmeans_pattern_count = downloadHandler(
    filename = function() {
      paste(paste(download_multi_overview_dir(),input$multi_selectfile2, sep = "_"), 
            "kmeans_count_table.txt", sep = "_")
    },
    content = function(file){write.table(multi_kmeans_extract(), file, row.names = T, sep = "\t", quote = F)}
  )
  
  #Multi DEG enrichment------------
  output$multi_Spe1 <- renderText({
    if(input$Species6 == "not selected") print("Please select 'Species'")
  })
  multi_Hallmark_set <- reactive({
    return(GeneList_for_enrichment(Species = input$Species6, Gene_set = input$Gene_set6, org = org6()))
  })
  
  output$selectEnrich_pair <- renderUI({
    if(is.null(multi_deg_result())){
      return(NULL)
    }else{
      selectizeInput("selectEnrich_pair", "Select a pair for GSEA", c(multi_select()),
                     selected = "", multiple = TRUE, 
                     options = list(maxItems = 2))
    }
  })
  
  multi_enrich_pairFC <- reactive({
    count <- multi_d_row_count_matrix()
    meta <- multi_metadata()
    dds <- multi_dds()
    if(is.null(dds) || length(input$selectEnrich_pair) != 2){
      return(NULL)
    }else{
      if (input$multi_data_file_type == "Row1"){
        collist <- gsub("\\_.+$", "", colnames(count))
        meta <- data.frame(condition = factor(collist))
      }else meta <- data.frame(condition=factor(meta[,1]), type=factor(meta[,2]))
      res <- results(dds, contrast = c("meta", as.character(input$selectEnrich_pair[1]),as.character(input$selectEnrich_pair[2])))
      sig_res_LRT <- res %>% 
        data.frame()%>%
        rownames_to_column(var="Row.names")
      sig_res_LRT$log2FoldChange <- -1 * sig_res_LRT$log2FoldChange 
      
      if(str_detect(rownames(count)[1], "ENS")){
        if(input$Species6 != "not selected"){
          my.symbols <- sig_res_LRT$Row.names
          gene_IDs<-AnnotationDbi::select(org6(),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data <- merge(sig_res_LRT, gene_IDs, by="Row.names")
          data$Unique_ID <- paste(data$SYMBOL,data$Row.names, sep = " - ")
        }
      }else{
        if(input$Species6 != "not selected"){
          my.symbols <- sig_res_LRT$Row.names
          gene_IDs<-AnnotationDbi::select(org6(),keys = my.symbols,
                                          keytype = "SYMBOL",
                                          columns = c("SYMBOL", "ENTREZID"))
          colnames(gene_IDs) <- c("Row.names", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data <- merge(sig_res_LRT, gene_IDs, by="Row.names")
        }
      }
      return(data)
    }
  })
  
  multi_enrichment_1_gsea <- reactive({
    if(!is.null(input$Gene_set6) && input$Species6 != "not selected"){
      data <- multi_enrich_pairFC()
      count <- multi_deg_norm_count()
      data <- na.omit(data)
      geneList <- data$log2FoldChange
      names(geneList) = as.character(data$ENTREZID)
      geneList <- sort(geneList, decreasing = TRUE)
        withProgress(message = "GSEA",{
          H_t2g <- multi_Hallmark_set()
          H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
          em3 <- GSEA(geneList, TERM2GENE = H_t2g2,pvalueCutoff = 0.05,exponent = 1, eps = 0, pAdjustMethod = "BH",
                      minGSSize = 50, maxGSSize = 500,by = "fgsea",verbose = F)
          if (length(as.data.frame(em3)$ID) == 0) {
            em4 <- NA
          } else{
            em4 <- setReadable(em3, org6(), 'ENTREZID')
          }
          return(em4)
          incProgress(1)
        })
    }else return(NULL)
  })
  
  
  # Multi DEG enrichment plot ------------------------------------------------------------------------------
   multi_enrich1_H <- reactive({
    if(!is.null(input$Gene_set6) && input$Species6 != "not selected"){
        em3 <- multi_enrichment_1_gsea()
        if (length(as.data.frame(em3)$ID) == 0) {
          p4 <- NULL
        } else{
          if (length(as.data.frame(em3)$ID) >= 5){
            p4 <- as.grob(gseaplot2(em3, 1:5, pvalue_table = F,base_size = 14))
          }else{
            p4 <- as.grob(gseaplot2(em3, 1:length(em3$ID), pvalue_table = F,base_size = 14))
          }
        }
        p <- plot_grid(p4, nrow = 1)
        return(p)
    }else return(NULL)
  })
  
  output$multi_enrichment1 <- renderPlot({
    dotplot_for_output(data = multi_enrich_pairFC(),
                       plot_genelist = multi_enrich1_H(), Gene_set = input$Gene_set6, 
                       Species = input$Species6)
  })
  
  output$Gene_set6 <- renderUI({
    selectInput('Gene_set6', 'Gene Set', gene_set_list)
  })
  
  multi_GSEA_table <- reactive({
    data <- as.data.frame(multi_enrichment_1_gsea())
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
      return(data3) 
      }
    }
  })
  
  output$multi_GSEA_result <- DT::renderDataTable({
    if(is.null(dds) || length(input$selectEnrich_pair) != 2){
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
            pdf_width <- pdf_w(clusterNumber)
          }else pdf_width <- 7
          pdf(file, height = pdf_height, width = pdf_width)
        print(plot_grid(p1, nrow =1))
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
  
  #multi DEG enrichment 2--------
  multi_Hallmark_set2 <- reactive({
    return(GeneList_for_enrichment(Species = input$Species6, Gene_set = input$Gene_set7, org = org6()))
  })
  
  multi_Hallmark_set3 <- reactive({
    return(GeneList_for_enrichment(Species = input$Species6, Gene_set = input$Gene_set8, org = org6()))
  })
  
  
  output$Gene_set7 <- renderUI({
    selectInput('Gene_set7', 'Gene Set', gene_set_list)
  })
  
  output$Gene_set8 <- renderUI({
    selectInput('Gene_set8', 'Gene Set', gene_set_list)
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
      selectInput("multi_whichGroup1_1", "cluster_list", choices = c(unique(clusters$cluster)), multiple = TRUE)
    }
  })
  
  output$multi_whichGroup2_1 <- renderUI({
    clusters <- multi_kmeans_cluster()
    if(is.null(clusters)){
      return(NULL)
    }else{
      selectInput("multi_whichGroup2_1", "cluster_list", choices = c(unique(clusters$Cluster)),multiple = TRUE)
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
  
  multi_enrich_viewer2 <- reactive({
    return(enrich_viewer_forMulti2(df = multi_enrich_input1(), Species = input$Species6, org = org6(),
                                   org_code = org_code6(),H_t2g = multi_Hallmark_set2(),Gene_set = input$Gene_set7))
  })
  
  multi_enrich_viewer12 <- reactive({
    return(enrich_viewer_forMulti2(df = multi_enrich_input2(), Species = input$Species6, org = org6(),
                                   org_code = org_code6(),H_t2g = multi_Hallmark_set3(),Gene_set = input$Gene_set8))
  })
  multi_enrich_h <- reactive({
    return(enrich_gene_list(data = enrich_viewer_forMulti1(df = multi_enrich_input1(), Species = input$Species6, org = org6()),
                           Gene_set = input$Gene_set7, org = org6(), H_t2g = multi_Hallmark_set2()))
  })
  multi_enrich_H <- reactive({
    return(enrich_genelist(data = enrich_viewer_forMulti1(df = multi_enrich_input1(), Species = input$Species6, org = org6()),
                           enrich_gene_list = multi_enrich_h()))
  })
  multi_enrich_h2 <- reactive({
    return(enrich_gene_list(data = enrich_viewer_forMulti1(df = multi_enrich_input2(), Species = input$Species6, org = org6()),
                           Gene_set = input$Gene_set8, org = org6(), H_t2g = multi_Hallmark_set3()))
  })
  multi_enrich_H2 <- reactive({
    return(enrich_genelist(data = enrich_viewer_forMulti1(df = multi_enrich_input2(), Species = input$Species6, org = org6()),
                            enrich_gene_list = multi_enrich_h2()))
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
    cnet_global(data = enrich_viewer_forMulti1(df = multi_enrich_input3(), Species = input$Species6, org = org6()), 
                group = input$multi_whichGroup1_2, enrich_gene_list = multi_enrich_h())
  })
  
  multi_enrich12 <- reactive({
    cnet_global(data = enrich_viewer_forMulti1(df = multi_enrich_input4(), Species = input$Species6, org = org6()), 
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
    return(enrich_for_table(data = multi_enrich_viewer2(), H_t2g = multi_Hallmark_set2(), Gene_set = input$Gene_set7))
  })
  
  multi_enrich_k_table <- reactive({
    return(enrich_for_table(data = multi_enrich_viewer12(), H_t2g = multi_Hallmark_set3(), Gene_set = input$Gene_set8))
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
    content = function(file){write.table(multi_enrich_div_table, file, row.names = F, sep = "\t", quote = F)}
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
    PCAplot(data)
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
    if(is.null(multi_d_row_count_matrix())){
      return(NULL)
    }else{
      print(multi_pca_plot())
    }
  })
  
  output$multi_PCA_data <- DT::renderDataTable({
    PCAdata(row_count = multi_d_row_count_matrix(), deg_norm_count = multi_deg_norm_count())
  })
  
  output$download_multi_PCA_table = downloadHandler(
    filename = function() {
      paste0(download_multi_overview_dir(), "_PCA_table.txt")
    },
    content = function(file){write.table(PCAdata(row_count = multi_d_row_count_matrix(), deg_norm_count = multi_deg_norm_count()), file, row.names = T, sep = "\t", quote = F)}
  )
  
  output$multi_umap_n <- renderUI({
    sliderInput("multi_n_neighbors", "n_neighbors", min = 2,
                max=100, step = 1,
                value = 15)
  })
  
  multi_umap_plot <- reactive({
    data <- multi_d_row_count_matrix()
    if(is.null(input$multi_n_neighbors)){
      return(NULL)
    }else{
      if(is.null(data)){
        return(NULL)
      }else{
        p<- try(umap_plot(data = data, n_neighbors = input$multi_n_neighbors))
        return(p)
      }
    }
  })
  
  output$multi_umap_error <- renderText({
    p <- multi_umap_plot()
    if(length(p) == 1){
      if (class(p) == "try-error") {
        print("umap: number of neighbors must be smaller than number of items")
      }
    }else return(NULL)
  })
  
  
  output$multi_umap <- renderPlot({
    p <- multi_umap_plot()
    withProgress(message = "umap",{
      if(length(p) == 1){
        if (class(p) == "try-error") {
          return(NULL)
        }
      }else{
        print(p)
      }
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
        pdf(file, height = pdf_height, width = pdf_width)
        print(multi_umap())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  
  # 3 conditions ------------------------------------------------------------------------------
  org2 <- reactive({
    return(org(Species = input$Species2))
  })
  org_code2 <- reactive({
    return(org_code(Species = input$Species2))
  })
  
  row_count_matrix2 <- reactive({
    if (input$data_file_type2 == "Row3"){
      tmp <- input$file4$datapath
      if(is.null(input$file4) && input$goButton2 > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example4.txt"
      return(read_df(tmp = tmp))
    }else{
      tmp <- input$file5$datapath
      if(is.null(input$file5) && input$goButton2 > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example2.csv"
      return(read_df(tmp = tmp))
    }
  })
  metadata2 <- reactive({
    if (input$data_file_type2 == "Row3"){
      return(NULL)
    }else{
      tmp <- input$file6$datapath
      if(is.null(input$file6) && input$goButton2 > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example6.csv"
      df <- read_df(tmp = tmp)
      if(!is.null(df)) rownames(df) <- gsub("-",".",rownames(df))
      return(df)
    }
  })
  norm_count_matrix2 <- reactive({
    tmp <- input$norm_file2$datapath
    df <- read_df(tmp = tmp)
    if(!is.null(df)){
      df <- anno_rep(df)
      if(input$Species2 != "not selected"){
        if(str_detect(rownames(df)[1], "ENS")){
          rownames(df) < gsub("\\..*","", rownames(df))
        }
      }
    }else df <- NULL
    return(df)
  })
  
  d_row_count_matrix2 <- reactive({
    row <- row_count_matrix2()
    if (input$data_file_type2 == "Row3"){
      if(is.null(row)) {
        return(NULL)
      }else{
        return(anno_rep(row))
      }
    }else{
      meta <- anno_rep_meta(metadata2())
      if (is.null(row) || is.null(meta)){
        return(NULL)
      } else {
        row_t <- t(row)
        meta <- data.frame(characteristics = meta[,1], row.names = rownames(meta))
        colname <- colnames(meta)
        data <- merge(meta, row_t, by=0, sort = F)
        rownames(data) <- data$characteristics
        data2 <- data[, - which(colnames(data) %in% c("Row.names", colname))]
        data2_t <- t(data2)
        data3 <- apply(data2_t, 2, as.numeric)
        rownames(data3) <- rownames(data2_t)
        if(input$Species2 != "not selected"){
          if(str_detect(rownames(data3)[1], "ENS")){
            rownames(data3) < gsub("\\..*","", rownames(data3))
          }
        }
        return(data3)
      }
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
    content = function(file){write.table(d_row_count_matrix2(), file, row.names = T, sep = "\t", quote = F)}
  )
  
  gene_ID <- reactive({
    res <- d_row_count_matrix2()
    if(is.null(res)){
      return(NULL)
    }else{
      if(input$Species2 != "not selected"){
        if(str_detect(rownames(res)[1], "ENS")){
          my.symbols <- gsub("\\..*","", rownames(res))
          gene_IDs<-AnnotationDbi::select(org2(),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL")
          res$Row.names <- gsub("\\..*","", rownames(res))
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          return(gene_IDs)
        }
      }else{ return(NULL) }
    }
  })
  

  
  # 3 conditions DEG ------------------------------------------------------------------------------
  MultiOut <- reactive({
    withProgress(message = "EBSeq multiple comparison test takes 5 - 10 minutes",{
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
      incProgress(1)
    })
  })
  
  deg_result2 <- reactive({
    count <- d_row_count_matrix2()
    if(is.null(count)){
      return(NULL)
    }else{
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
      res <- as.data.frame(results)
      
      if(input$Species2 != "not selected"){
        if(str_detect(rownames(count)[1], "ENS")){
          gene_IDs  <- gene_ID()
          res$Row.names <- rownames(res)
          data2 <- merge(res, gene_IDs, by="Row.names")
          rownames(data2) <- data2$Row.names
          res <- data2[,-1]
        }
      }
      return(res)
    }
  })
  deg_result2_pattern <- reactive({
    count <- d_row_count_matrix2()
    if(is.null(count)){
      return(NULL)
    }else{
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
    }
  })
  deg_result2_condmean <- reactive({
    count <- d_row_count_matrix2()
    if(is.null(count)){
      return(NULL)
    }else{
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
    }
  })
  
  deg_norm_count2 <- reactive({
    count <- d_row_count_matrix2()
      if(is.null(count)){
        return(NULL)
      }else{
        if(is.null(norm_count_matrix2())){
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
        normalized_counts <- as.data.frame(GetNormalizedMat(count, Sizes))
        }else{
          normalized_counts <- norm_count_matrix2()
        }
        if(input$Species2 != "not selected"){
          if(str_detect(rownames(count)[1], "ENS")){
            gene_IDs  <- gene_ID()
            normalized_counts$Row.names <- rownames(normalized_counts)
            data2 <- merge(normalized_counts, gene_IDs, by="Row.names")
            data2$Unique_ID <- paste(data2$SYMBOL,data2$Row.names, sep = " - ")
            rownames(data2) <- data2$Row.names
            normalized_counts <- data2[,-1]
          }
        }
        return(normalized_counts)
      }
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
    content = function(file){write.table(deg_norm_count2(), file, row.names = T, sep = "\t", quote = F)}
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
  
  output$cond3_result <- DT::renderDataTable({
    data_3degcount1(data = deg_norm_count2(),result_Condm = deg_result2_condmean(),
                    result_FDR = deg_result2(), specific = 1,result_list=TRUE) %>% as.data.frame()
  })
  output$download_cond3_result = downloadHandler(
    filename = function() {paste(download_cond3_dir(),"DEG_result_ALL.txt", sep = "-")},
    content = function(file){write.table(data_3degcount1(data = deg_norm_count2(),result_Condm = deg_result2_condmean(),
                                                         result_FDR = deg_result2(), specific = 1,result_list=TRUE), 
                                         file, row.names = T, sep = "\t", quote = F)}
  )
  
  
  #3conditions DEG vis------------------------
  
  #3conditions DEG_1------------------------
  data_3degcount1_1 <- reactive({
    data3 <- data_3degcount1(data = deg_norm_count2(),result_Condm = deg_result2_condmean(),
                             result_FDR = deg_result2(), specific = 1,fc = input$fc2, 
                             fdr = input$fdr2, basemean = input$basemean2)
    return(data3)
  })
  
  data_3degcount2_1 <- reactive({
    return(data_3degcount2(data3 = data_3degcount1_1(), Species = input$Species2, org = org2()))
  })
  
  #3conditions scatter + heatmap_1
  
  cond3_scatter1_plot <- reactive({
    p <- cond3_scatter_plot(data = deg_norm_count2(), data4 = data_3degcount2_1(),
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
    data3 <- data_3degcount1(data = deg_norm_count2(),result_Condm = deg_result2_condmean(),
                             result_FDR = deg_result2(), specific = 2,fc = input$fc2, 
                             fdr = input$fdr2, basemean = input$basemean2)
    return(data3)
  })
  
  data_3degcount2_2 <- reactive({
    return(data_3degcount2(data3 = data_3degcount1_2(), Species = input$Species2, org = org2()))
  })
  
  #3conditions scatter + heatmap_2
  cond3_scatter2_plot <- reactive({
    p <- cond3_scatter_plot(data = deg_norm_count2(), data4 = data_3degcount2_2(),
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
    data3 <- data_3degcount1(data = deg_norm_count2(),result_Condm = deg_result2_condmean(),
                             result_FDR = deg_result2(), specific = 3,fc = input$fc2, 
                             fdr = input$fdr2, basemean = input$basemean2)
    return(data3)
  })
  
  data_3degcount2_3 <- reactive({
    return(data_3degcount2(data3 = data_3degcount1_3(), Species = input$Species2, org = org2()))
  })
  
  #3conditions scatter + heatmap_3
  cond3_scatter3_plot <- reactive({
    p <- cond3_scatter_plot(data = deg_norm_count2(), data4 = data_3degcount2_3(),
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
      paste0(download_cond3_dir(), "_GOI_scatter.pdf")
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
        print(cond3_scatter_plot(data = deg_norm_count2(), data4 = data_3degcount2_1(),
                                 result_Condm = deg_result2_condmean(), result_FDR = deg_result2(), 
                                 fc = input$fc2, fdr = input$fdr2, basemean = input$basemean2,
                                 y_axis = input$cond3_scatter_yrange,x_axis = input$cond3_scatter_xrange,heatmap = FALSE,
                                 specific = cond3_specific_group2(), GOI = input$GOI2, Species = input$Species2))
        dev.off()
      })
    }
  )
  
  
  #3conditions PCA--------------
  output$PCA2 <- renderPlot({
    if(is.null(deg_norm_count2())){
      return(NULL)
    }else{
    print(PCAplot(data = deg_norm_count2()))
    }
  })
  output$PCA3_data <- DT::renderDataTable({
    if(is.null(deg_norm_count2())){
      return(NULL)
    }else{
    PCAdata(row_count = deg_norm_count2(), deg_norm_count = deg_norm_count2())
    }
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
        print(PCAplot(data = deg_norm_count2()))
        dev.off()
      })
    }
  )
  
  output$download_cond3_pca_table = downloadHandler(
    filename = function() {
      paste0(download_cond3_dir(), "PCA_table.txt")
    },
    content = function(file){write.table(PCAdata(row_count = deg_norm_count2(), deg_norm_count = deg_norm_count2()), file, row.names = T, sep = "\t", quote = F)}
  )
  #3conditions GOI------------------------------------------------------
  GOI_list2 <- reactive({
      count <- deg_norm_count2()
      if(is.null(d_row_count_matrix2())){
        return(NULL)
      }else{
        if(str_detect(rownames(count)[1], "ENS")){
          if(input$Species2 != "not selected"){
            GOI <- count$Unique_ID
          }else GOI <- rownames(count)
        }else{
          if(input$Species2 != "not selected"){
            GOI <- rownames(count)
          }else GOI <- rownames(count)
        }
        return(GOI)
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
  
  output$GOI2 <- renderUI({
    if(is.null(d_row_count_matrix2())){
      return(NULL)
    }else{
      withProgress(message = "Preparing GOI list (about 10 sec)",{
        selectizeInput("GOI2", "genes of interest (GOI)", c(GOI_list2()),multiple = TRUE, options = list(delimiter = " ", create = T))
      })
    }
  })
  output$GOIreset_cond3 <- renderUI({
    actionButton("GOIreset_cond3", "GOI reset")
  })
  observeEvent(input$GOIreset_cond3, {
    withProgress(message = "Preparing GOI list (about 10 sec)",{
    updateSelectizeInput(session, "GOI2", choices = c(GOI_list2()), 
                         selected = character(0),
                         options = list(create=TRUE, 'plugins' = list('remove_button'), persist = FALSE))
    })
  })
  
  output$cond3_xrange <- renderUI({
    sliderInput("cond3_scatter_xrange","X_axis range:",min = -100,
                max=100, step = 1,
                value = c(-10, 10))
  })
  output$cond3_yrange <- renderUI({
    sliderInput("cond3_scatter_yrange","Y_axis range:",min = -100, max= 100, step = 1,
                value = c(-10, 10))
  })
  
  cond3_specific_group2 <- reactive({
    if(!is.null(input$cond3_GOIpair)){
      return(which(cond3_specific_group() == input$cond3_GOIpair))
    }
  })
  
  output$cond3_GOIscatter <- renderPlot({
    if(!is.null(input$cond3_GOIpair) && !is.null(input$cond3_scatter_yrange) && !is.null(input$cond3_scatter_xrange)){
    cond3_scatter_plot(data = deg_norm_count2(), data4 = data_3degcount2_1(),
                       result_Condm = deg_result2_condmean(), result_FDR = deg_result2(), 
                       fc = input$fc2, fdr = input$fdr2, basemean = input$basemean2,
                       y_axis = input$cond3_scatter_yrange,x_axis = input$cond3_scatter_xrange,heatmap = FALSE,
                       specific = cond3_specific_group2(), GOI = input$GOI2, Species = input$Species2)
    }
  })
  
  
  cond3_GOIcount <- reactive({
    count <- deg_norm_count2()
    if(str_detect(rownames(count)[1], "ENS")){
      if(input$Species2 != "not selected"){
        Unique_ID <- input$GOI2
        label_data <- as.data.frame(Unique_ID, row.names = Unique_ID)
        data <- merge(count, label_data, by="Unique_ID")
        rownames(data) <- data$Unique_ID
        data <- data[, - which(colnames(data) == "SYMBOL")]
        data <- data[,-1]
      }else{
        Row.names <- input$GOI2
        count$Row.names <- rownames(count)
        label_data <- as.data.frame(Row.names, row.names = Row.names)
        data <- merge(count, label_data, by="Row.names")
        rownames(data) <- data$Row.names
        data <- data[,-1]
      }
    }else{
      Row.names <- input$GOI2
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
      data.z <- genescale(data, axis=1, method="Z")
      data.z <- na.omit(data.z)
      ht <- GOIheatmap(data.z)
    }
    return(ht)
  })
  
  output$cond3_GOIheatmap <- renderPlot({
    if(is.null(d_row_count_matrix2())){
      return(NULL)
    }else{
      if(!is.null(input$GOI2)){
        withProgress(message = "Boxplot",{
          suppressWarnings(print(cond3_GOIheat()))
          incProgress(1)
        })
      }
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
      if(!is.null(input$GOI2)){
        withProgress(message = "Boxplot",{
          suppressWarnings(print(cond3_GOIbox()))
          incProgress(1)
        })
      }
    }
  })
  
  output$download_3cond_GOIbox = downloadHandler(
    filename = function(){
      paste0(download_cond3_dir(), "_GOIboxplot.pdf")
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
      paste0(download_cond3_dir(), "_GOIheatmap.pdf")
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
    return(keggEnrichment1(data3 = data_3degcount1_1(),data4 = data_3degcount2_1(),
                         Species = input$Species2, Gene_set = input$Gene_set2,
                         org = org2(),H_t2g = Hallmark_cond3()))
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
    return(GeneList_for_enrichment(Species = input$Species2, Gene_set = input$Gene_set2, org = org2()))
  })
  
  enrich3_2 <- reactive({
    return(keggEnrichment1(data3 = data_3degcount1_2(),data4 = data_3degcount2_2(),
                           Species = input$Species2, Gene_set = input$Gene_set2,
                           org = org2(),H_t2g = Hallmark_cond3()))
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
    return(keggEnrichment1(data3 = data_3degcount1_3(),data4 = data_3degcount2_3(),
                           Species = input$Species2, Gene_set = input$Gene_set2,
                           org = org2(),H_t2g = Hallmark_cond3()))
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
    selectInput('Gene_set2', 'Gene Set', gene_set_list)
  })
  
  cond3_enrich_table1 <- reactive({
    return(enrich_for_table(data = as.data.frame(enrichment3_1_1()), H_t2g = Hallmark_cond3(), Gene_set = input$Gene_set2))
  })
  cond3_enrich_table2 <- reactive({
    return(enrich_for_table(data = as.data.frame(enrichment3_2_1()), H_t2g = Hallmark_cond3(), Gene_set = input$Gene_set2))
  })
  cond3_enrich_table3 <- reactive({
    return(enrich_for_table(data = as.data.frame(enrichment3_3_1()), H_t2g = Hallmark_cond3(), Gene_set = input$Gene_set2))
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
  
  #normalized count analysis
  #norm_count_input-------------------
  org3 <- reactive({
    return(org(Species = input$Species3))
  })
  org_code3 <- reactive({
    return(org_code(Species = input$Species3))
  })
  
  norm_count_input <- reactive({
    withProgress(message = "Importing normalized count matrix, please wait",{
      if (input$data_file_type3 == "Row5"){
        tmp <- input$file7$datapath
        if(is.null(input$file7) && input$goButton3 > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example4.txt"
        return(read_df(tmp = tmp))
      }else{
        tmp <- input$file8$datapath
        if(is.null(input$file8) && input$goButton3 > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example4.txt"
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
      if(is.null(input$file9) && input$goButton3 > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example9.csv"
      df <- read_df(tmp = tmp)
      if(!is.null(df)) rownames(df) <- gsub("-",".",rownames(df))
      return(df)
    }
  })
  gene_list <- reactive({
    data <- input$file10$datapath
    if(is.null(input$file10) && input$goButton3 == 0) return(NULL)
    if(is.null(input$file10) && input$goButton3 > 0 )  data = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/enrich_example.txt"
    df <- read_df(tmp = data)
    gene <-c(rownames(df))
    df <- as.data.frame(gene, row.names = gene)
    return(df)
  })
  d_norm_count_matrix <- reactive({
    withProgress(message = "Creating defined count matrix, please wait",{
      row <- norm_count_input()
      gene_list <- gene_list()
      if (input$data_file_type3 == "Row5"){
        if(is.null(row)) {
          return(NULL)
        }else{
          if(!is.null(gene_list)){
            row <- merge(row,gene_list, by=0)
            rownames(row) <- row$Row.names
            row <- row[,-1]
            row <- row[, - which(colnames(row) == "gene")]
          }
            return(anno_rep(row))
        }
      }else{
        meta <- anno_rep_meta(norm_metadata())
        if (is.null(row) || is.null(meta)){
          return(NULL)
        } else {
          row_t <- t(row)
          meta <- data.frame(characteristics = meta[,1], row.names = rownames(meta))
          colname <- colnames(meta)
          data <- merge(meta, row_t, by=0, sort = F)
          rownames(data) <- data$characteristics
          data2 <- data[, - which(colnames(data) %in% c("Row.names", colname))]
          data2_t <- t(data2)
          data3 <- apply(data2_t, 2, as.numeric)
          rownames(data3) <- rownames(data2_t)
          row <- data3
          if(!is.null(gene_list)){
            row <- merge(row,gene_list, by=0)
            rownames(row) <- row$Row.names
            row <- row[,-1]
            row <- row[, - which(colnames(row) == "gene")]
          }
          return(row)
        }
      }
      incProgress(1)
    })
  })
  
  d_norm_count_matrix_cutofff <- reactive({
    data <- d_norm_count_matrix()
    if(is.null(data)){
      return(NULL)
    }else{
      withProgress(message = "Creating defined count matrix, please wait",{
        if(input$basemean3 != 0){
          data <- dplyr::filter(data, apply(data,1,mean) > input$basemean3)
        }
        return(data)
        incProgress(1)
      })
    }
  })
  
  d_norm_count_matrix2 <- reactive({
    data <- d_norm_count_matrix()
    if(is.null(data)){
      return(NULL)
    }else{
      if(input$Species3 != "not selected"){
        if(str_detect(rownames(data)[1], "ENS")){
          gene_IDs <- gene_ID_norm()
          data <- merge(data, gene_IDs, by=0)
          data <- data[, - which(colnames(data) == "Row.names.y")]
          rownames(data) <- data[,1] 
          data <- data[,-1]
        }
      }
      return(data)
    }
  })
  
  gene_ID_norm <- reactive({
    res <- d_norm_count_matrix()
    if(is.null(res)){
      return(NULL)
    }else{
      if(input$Species3 != "not selected"){
        if(str_detect(rownames(res)[1], "ENS")){
          my.symbols <- rownames(res)
          gene_IDs<-AnnotationDbi::select(org3(),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          rownames(gene_IDs) <- gene_IDs$Row.names
          return(gene_IDs)
        }
      }else{ return(NULL) }
    }
  })
  
  observeEvent(input$goButton3,({
    updateSelectInput(session,inputId = "Species3","Species",species_list, selected = "Mus musculus")
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
    d_norm_count_matrix2()
  })
  output$Gene_list <- DT::renderDataTable({
    gene_list()
  })
  output$download_d_norm_count = downloadHandler(
    filename = function() {
      if (input$data_file_type3 == "Row5"){
        paste(gsub("\\..+$", "", input$file7), ".txt")
      }else{
        paste(gsub("\\..+$", "", input$file8), paste0(gsub("\\..+$", "", input$file9),".txt"), sep ="-")
      }},
    content = function(file){write.table(d_norm_count_matrix2(), file, row.names = T, sep = "\t", quote = F)}
  )
  
  download_norm_dir <-reactive({
    if (input$data_file_type3 == "Row5"){
      dir_name <- paste0(gsub("\\..+$", "", input$file7), "_")
    }else{
      dir_name <- paste0(paste(gsub("\\..+$", "", input$file8),paste(gsub("\\..+$", "", input$file9),input$basemean3,sep = "_"),sep = "-"), "_")
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
        print(PCAplot(data = d_norm_count_matrix_cutofff()))
        dev.off()
      })
    }
  )
  
  output$norm_PCA <- renderPlot({
    withProgress(message = "Clustering",{
      if(is.null(d_norm_count_matrix_cutofff())){
        return(NULL)
      }else{
        print(PCAplot(data = d_norm_count_matrix_cutofff()))
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
    content = function(file){write.table(PCAdata(row_count = d_norm_count_matrix_cutofff(), deg_norm_count = d_norm_count_matrix_cutofff()), file, row.names = T, sep = "\t", quote = F)}
  )
  
  output$norm_umap_n <- renderUI({
    sliderInput("norm_n_neighbors", "n_neighbors", min = 2,
                 max=100, step = 1,
                 value = 15)
  })
  
  norm_umap_plot <- reactive({
      data <- d_norm_count_matrix_cutofff()
      if(is.null(input$norm_n_neighbors)){
        return(NULL)
      }else{
        if(is.null(data)){
          return(NULL)
        }else{
          p<- try(umap_plot(data = data, n_neighbors = input$norm_n_neighbors))
          return(p)
        }
      }
  })
  
  output$norm_umap_error <- renderText({
    p <- norm_umap_plot()
    if(length(p) == 1){
    if (class(p) == "try-error") {
      print("umap: number of neighbors must be smaller than number of items")
      }
      }else return(NULL)
  })
  
  
  output$norm_umap <- renderPlot({
    p <- norm_umap_plot()
    withProgress(message = "umap",{
      if(length(p) == 1){
        if (class(p) == "try-error") {
          return(NULL)
        }
      }else{
        print(p)
      }
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
        pdf(file, height = pdf_height, width = pdf_width)
        print(norm_umap_plot())
        dev.off()
        incProgress(1)
      })
    }
  )
  
  #norm GOI------------------------------------------------------
  d_norm_count_cutoff_uniqueID <- reactive({
    count <- d_norm_count_matrix_cutofff()
    if(input$Species3 != "not selected"){
      if(str_detect(rownames(count)[1], "ENS")){
        gene_IDs  <- gene_ID_norm()
        data2 <- merge(count, gene_IDs, by= 0)
        rownames(data2) <- data2[,1]
        data2 <- data2[, - which(colnames(data2) == "Row.names.y")]
        data2$Unique_ID <- paste(data2$SYMBOL,data2$Row.names, sep = " - ")
        count <- data2[,-1]
      }
    }
    return(count)
  })
  
  GOI_list3 <- reactive({
      count <- d_norm_count_cutoff_uniqueID()
      if(is.null(count)){
        return(NULL)
      }else{
        if(str_detect(rownames(count)[1], "ENS")){
          if(input$Species3 != "not selected"){
            GOI <- count$Unique_ID
          }else GOI <- rownames(count)
        }else{
          if(input$Species3 != "not selected"){
            GOI <- rownames(count)
          }else GOI <- rownames(count)
        }
        return(GOI)
      }
  })
  
  output$GOI3 <- renderUI({
    if(is.null(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      withProgress(message = "Preparing GOI list (about 10 sec)",{
        selectizeInput("GOI3", "genes of interest (GOI)", c(GOI_list3()),multiple = TRUE, options = list(delimiter = " ", create = T))
      })
    }
  })
  output$GOIreset_norm <- renderUI({
    actionButton("GOIreset_norm", "GOI reset")
  })
  observeEvent(input$GOIreset_norm, {
    withProgress(message = "Preparing GOI list (about 10 sec)",{
    updateSelectizeInput(session, "GOI3", choices = c(GOI_list3()), 
                         selected = character(0),
                         options = list(create=TRUE, 'plugins' = list('remove_button'), persist = FALSE))
    })
  })
  
  norm_GOIcount <- reactive({
    count <- d_norm_count_cutoff_uniqueID()
    if(str_detect(rownames(count)[1], "ENS")){
      if(input$Species3 != "not selected"){
        Unique_ID <- input$GOI3
        label_data <- as.data.frame(Unique_ID, row.names = Unique_ID)
        data <- merge(count, label_data, by="Unique_ID")
        rownames(data) <- data$Unique_ID
        data <- data[, - which(colnames(data) == "SYMBOL")]
        data <- data[, - which(colnames(data) == "Unique_ID")]
      }else{
        Row.names <- input$GOI3
        count$Row.names <- rownames(count)
        label_data <- as.data.frame(Row.names, row.names = Row.names)
        data <- merge(count, label_data, by="Row.names")
        rownames(data) <- data$Row.names
        data <- data[, - which(colnames(data) == "Row.names")]
      }
    }else{
      Row.names <- input$GOI3
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
      data.z <- genescale(data, axis=1, method="Z")
      data.z <- na.omit(data.z)
      ht <- GOIheatmap(data.z)
    }
    return(ht)
  })
  
  output$norm_GOIheatmap <- renderPlot({
    if(is.null(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      if(!is.null(input$GOI3)){
        withProgress(message = "heatmap",{
          suppressWarnings(print(norm_GOIheat()))
          incProgress(1)
        })
      }
    }
  })
  
  norm_GOIbox <- reactive({
    data <- norm_GOIcount()
    if(is.null(data)){
      p <- NULL
    }else{
      p <- GOIboxplot(data = data)
    }
    return(p)
  })
  
  
  output$norm_GOIboxplot <- renderPlot({
    if(is.null(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      if(!is.null(input$GOI3)){
        withProgress(message = "Boxplot",{
          suppressWarnings(print(norm_GOIbox()))
          incProgress(1)
        })
      }
    }
  })
  
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
  #norm kmeans------------------------------------------------------
  output$norm_kmeans_num <- renderUI({
    if(is.null(d_norm_count_matrix_cutofff())){
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
    data <- d_norm_count_matrix_cutofff()
    if(is.null(data) || is.null(input$kmeans_cv)){
      return(NULL)
    }else{
      data2 <- data[order(apply(data,1,mad), decreasing = T)[1:input$kmeans_cv],]
      return(data2)
    }
  })
  
  norm_data_z <- reactive({
    data <- norm_count_matrix_cutoff2()
    if(is.null(data)){
      return(NULL)
    }else{
      data.z <- genescale(data, axis = 1, method = "Z")
      data.z <- na.omit(data.z)
      return(data.z)
    }
  })
  
  norm_kmeans <- reactive({
    data.z <- norm_data_z()
    if(is.null(data.z)){
      return(NULL)
    }else{
      withProgress(message = "k-means clustering",{
        ht <- Heatmap(data.z, name = "z-score",
                      column_order = colnames(data.z),
                      clustering_method_columns = 'ward.D2',
                      row_km= input$norm_kmeans_number, cluster_row_slices = F, row_km_repeats = 100,
                      show_row_names = F,column_names_side = "top")
        ht <- draw(ht)
        return(ht)
      })
    }
  })
  
  norm_kmeans_cluster <- reactive({
    ht <- norm_kmeans()
    data.z <- norm_data_z()
    data <- norm_count_matrix_cutoff2()
    if(is.null(ht) || is.null(data.z)){
      return(NULL)
    }else{
      r.dend <- row_dend(ht)
      rcl.list <- row_order(ht)
      lapply(rcl.list, function(x) length(x))
      Cluster <- NULL
      if(!is.null(input$norm_kmeans_number)){
        if(length(lapply(rcl.list, function(x) length(x))) != input$norm_kmeans_number){
          return(NULL)
        }else{
          for (i in 1:length(row_order(ht))){ if (i == 1) {
            clu <- t(t(row.names(data.z[row_order(ht)[[i]],])))
            out <- cbind(clu, paste("cluster", i, sep=""))
            colnames(out) <- c("GeneID", "Cluster")} else {
              clu <- t(t(row.names(data.z[row_order(ht)[[i]],])))
              clu <- cbind(clu, paste("cluster", i, sep=""))
              out <- rbind(out, clu)}}
          out <- as.data.frame(out)
          rownames(out) <- out$GeneID
          clusterCount <- merge(out, data, by=0)
          rownames(clusterCount) <- clusterCount$GeneID
          clusterCount <- clusterCount[,-1:-2]
          return(clusterCount)
        }
      }else return(NULL)
    }
  })
  
  output$norm_kmeans_heatmap <- renderPlot({
    ht <- norm_kmeans()
    if(is.null(ht)){
      return(NULL)
    }else{
      print(ht)
    }
  })
  
  output$norm_kmeans_count_table <- renderDT({
    norm_kmeans_cluster()
  })
  
  
  norm_kmeans_GOIbox <- reactive({
    if(!is.null(input$norm_kmeans_count_table_rows_selected)){
      data <- norm_kmeans_cluster()[input$norm_kmeans_count_table_rows_selected,]
      if(input$Species3 != "not selected"){
        if(str_detect(rownames(data)[1], "ENS")){
          rownames(data) <- paste(data$SYMBOL,rownames(data), sep = " - ")
          data <- data[, - which(colnames(data) == "SYMBOL")]
        }
      }
      data <- data[, - which(colnames(data) == "Cluster")]
      return(data)
    }
  })
  
  output$norm_kmeans_box <- renderPlot({
    if(!is.null(input$norm_kmeans_count_table_rows_selected)){
      GOIboxplot(data = norm_kmeans_GOIbox())
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
    content = function(file){write.table(norm_kmeans_cluster(), file, row.names = T, sep = "\t", quote = F)}
  )
  output$download_norm_kmeans_heatmap = downloadHandler(
    filename = function() {
      paste0(download_norm_dir(), paste(input$norm_kmeans_number,"kmeans_heatmap.pdf",sep = "_"))
    },
    content = function(file){
      withProgress(message = "Preparing download",{
        if(input$norm_pdf_height == 0){
          pdf_height <- 10
        }else pdf_height <- input$norm_pdf_height
        if(input$norm_pdf_width == 0){
          pdf_width <- 7
        }else pdf_width <- input$norm_pdf_width
        pdf(file, height = pdf_height, width = pdf_width)
        print(norm_kmeans())
        dev.off()
        incProgress(1)
      })
    }
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
        df[["day0"]] <-  c(rownames(read.table("https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example11.txt",header = T, row.names = 1,sep="\t")))
        df[["day1"]] <- c(rownames(read.table("https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example12.txt",header = T, row.names = 1,sep="\t")))
        df[["day5"]] <- c(rownames(read.table("https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/example13.txt",header = T, row.names = 1,sep="\t")))
        return(df)
      }
      return(NULL)
    }else{
      for(nr in 1:length(input$files[, 1])){
        df <- read_df(input$files[[nr, 'datapath']])
        name <- c(name, gsub("\\..+$", "", input$files[nr,]$name))
        upload[[nr]] <- c(rownames(df))
      }
      names(upload) <- name
      return(upload)
    }
  })
  
  observeEvent(input$goButton_venn,({
    updateSelectInput(session,inputId = "Species7","Species",species_list, selected = "Mus musculus")
  }))
  
  output$venn <- renderPlot({
    if(is.null(files_table())){
      return(NULL)
    }else{
      gene_list <- files_table()
      venn::venn(gene_list, ilabels = TRUE, zcolor = "style", opacity = 0, ilcs = 1.5, sncs = 1.5)
    }
  })
  
  overlap_list <- reactive({
    gene_list <- files_table()
    if(!is.null(gene_list)){
      data <- names(attr(venn(gene_list), "intersections"))
      return(data)
    }else return(NULL)
  })
  
  overlap_table2 <- reactive({
    df <- data.frame("Gene" = NA, "Group"=NA)
    gene_list <- files_table()
    if(is.null(gene_list)){
      return(NULL)
    }else{
      for (name in names(attr(venn(gene_list),"intersections"))){
        data <- as.data.frame(attr(venn(gene_list),"intersections")[name])
        data <- cbind(data, name)
        colnames(data) <- c("Gene", "Group")
        df <- rbind(df, data)
      }
      df <- na.omit(df)
      return(df)
    }
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
          if(input$venn_pdf_height == 0){
            pdf_height <- 3
          }else pdf_height <- input$venn_pdf_height
          if(input$venn_pdf_width == 0){
            pdf_width <- 3
          }else pdf_width <- input$venn_pdf_width
          pdf(file, height = pdf_height, width = pdf_width)
          print((venn::venn(gene_list, ilabels = TRUE, zcolor = "style", opacity = 0, ilcs = 1, sncs = 1 )))
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
        df["day0"] <- list(read.table("https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/day0.txt",header = T, row.names = 1))
        df["day1"] <- list(read.table("https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/day1.txt",header = T, row.names = 1))
        df["day5"] <- list(read.table("https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/day5.txt",header = T, row.names = 1))
        return(df)
      }
      return(NULL)
    }else{
      upload = list()
      name = c()
      for(nr in 1:length(input$countfiles[, 1])){
        df <- read_df(input$countfiles[[nr, 'datapath']])
        name <- c(name, gsub("\\..+$", "", input$countfiles[nr,]$name))
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
        base <- matrix_list[[1]]
        int_matrix <- lapply(matrix_list[-1], function(i) base <<- merge(base, i, by = "Row.names"))
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
            matrix_z <- genescale(matrix_2, axis = 1, method = "Z")
            matrix_z <- na.omit(matrix_z)
            matrix_z <- merge(matrix, matrix_z, by = 0)[,-2:-(1 + length(colnames(matrix)))]
            matrix_z_list[name] <- list(matrix_z)
          }
          matrix_3 <- merge(matrix, matrix_2, by = 0)[,-2:-(1 + length(colnames(matrix)))]
          matrix_list[name] <- list(matrix_3)
        }
        base <- matrix_list[[1]]
        int_matrix <- lapply(matrix_list[-1], function(i) base <<- merge(base, i, by = "Row.names"))
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
          base_z <- matrix_z_list[[1]]
          int_matrix <- lapply(matrix_z_list[-1], function(i) base_z <<- merge(base_z, i, by = "Row.names"))
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
          base <- genescale(base, axis = 1, method = "Z")
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
          if(length(rownames(base_z)) <= 50){
            ht <- Heatmap(base_z, name = "z-score",
                          clustering_method_columns = 'ward.D2',
                          cluster_row_slices = T, show_row_names = T,
                          top_annotation = HeatmapAnnotation(condition = cond),column_names_side = "top",
                          row_names_gp = gpar(fontface = "italic"))
          }else{
            ht <- Heatmap(base_z, name = "z-score",
                          clustering_method_columns = 'ward.D2',
                          cluster_row_slices = T, show_row_names = F,
                          top_annotation = HeatmapAnnotation(condition = cond),column_names_side = "top")
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
    content = function(file){write.table(integrated_count(), file, row.names = T, sep = "\t", quote = F)}
  )
  output$download_integrated_z_count_table = downloadHandler(
    filename = function(){
      paste(input$selectfile, "integrated_zscored_count_table.txt", sep="_")
    },
    content = function(file){write.table(integrated_count_z(), file, row.names = T, sep = "\t", quote = F)}
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
      selectInput("venn_whichGroup1", "gene list", choices = c(unique(clusters$Group)),multiple = TRUE)
    }
  })
  
  venn_enrich_input1 <- reactive({
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
  
  output$Gene_set9 <- renderUI({
    selectInput('Gene_set9', 'Gene Set', gene_set_list)
  })
  
  output$venn_Spe <- renderText({
    if(input$Species7 == "not selected") print("Please select 'Species'")
  })
  
  org7 <- reactive({
    return(org(Species = input$Species7))
  })
  org_code7 <- reactive({
    return(org_code(Species = input$Species7))
  })
  
  venn_Hallmark_set <- reactive({
    return(GeneList_for_enrichment(Species = input$Species7, Gene_set = input$Gene_set9, org = org7()))
  })
  
  venn_enrich_viewer2 <- reactive({
    return(enrich_viewer_forMulti2(df = venn_enrich_input1(), Species = input$Species7, org = org7(),
                                   org_code = org_code7(),H_t2g = venn_Hallmark_set(),Gene_set = input$Gene_set9))
  })
  venn_enrich_h <- reactive({
    return(enrich_gene_list(data = enrich_viewer_forMulti1(df = venn_enrich_input1(), Species = input$Species7, org = org7()),
                           Gene_set = input$Gene_set9, org = org7(), H_t2g = venn_Hallmark_set()))
  })
  venn_enrich_H <- reactive({
    return(enrich_genelist(data = enrich_viewer_forMulti1(df = venn_enrich_input1(), Species = input$Species7, org = org7()),
                           enrich_gene_list = venn_enrich_h()))
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
    return(enrich_for_table(data = as.data.frame(venn_enrich_viewer2()), H_t2g = venn_Hallmark_set(), Gene_set = input$Gene_set9))
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
    cnet_global(data = enrich_viewer_forMulti1(df = venn_enrich_input2(), Species = input$Species7, org = org7()), 
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
    return(GeneList_for_enrichment(Species = input$Species4, Gene_set = input$Gene_set3, org = org4(), Custom_gene_list = Custom_input()))
  })
  org4 <- reactive({
    return(org(Species = input$Species4))
  })
  org_code4 <- reactive({
    return(org_code(Species = input$Species4))
  })
  
  enrich_input <- reactive({
    tmp <- input$enrich_data_file$datapath
    if(is.null(input$enrich_data_file) && input$goButton4 > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/enrich_example.txt"
    return(read_gene_list(tmp))
  })
  
  Custom_input <- reactive({
    tmp <- input$custom_input$datapath
    data <- read_gene_list(tmp)
    df <- gene_list_convert_for_enrichment(data= data, Species = input$Species4)
    return(df)
  })
  
  output$enrichment_input <- DT::renderDataTable({
    as.data.frame(enrich_input())
  })
  
  
  enrich_viewer1 <- reactive({
    return(gene_list_convert_for_enrichment(data= enrich_input(), Species = input$Species4))
  })
  
  enrich_viewer2 <- reactive({
    if(!is.null(input$Gene_set3)){
      data3 <- enrich_viewer1()
      if(is.null(data3)){
        return(NULL)
      }else{
          withProgress(message = "enrichment analysis",{
            H_t2g <- Hallmark_enrich()
            if(is.null(H_t2g)){
              df <- NULL
            }else{
              H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
            df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
            colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
            for (name in unique(data3$Group)) {
              em <- enricher(data3$ENTREZID[data3$Group == name], TERM2GENE=H_t2g2, qvalueCutoff = 0.05)
              if (length(as.data.frame(em)$ID) != 0) {
                if(length(colnames(as.data.frame(em))) == 9){
                  cnet1 <- as.data.frame(setReadable(em, org4(), 'ENTREZID'))
                  cnet1$Group <- name
                  df <- rbind(df, cnet1)
                }
              }
            }
            }
            if(length(df$ID) !=0){
              df["Description"] <- lapply(df["Description"], gsub, pattern="HALLMARK_", replacement = "")
              df$GeneRatio <- parse_ratio(df$GeneRatio)
              return(df)
            }else return(NULL)
          })
      }
    }
  })
  
  # enrichment plot ------------------------------------------------------------------------------
  enrich_venn <- reactive({
    return(enrich_gene_list(data = enrich_viewer1(), Gene_set = input$Gene_set3,
                            org = org4(), H_t2g = Hallmark_enrich()))
  })
  enrich_H <- reactive({
    return(enrich_genelist(data = enrich_viewer1(), enrich_gene_list = enrich_venn(), 
                           showCategory = input$enrich_showCategory))
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
      paste(gsub("\\..+$", "", input$enrich_data_file), paste0(input$Gene_set3,".pdf"), sep ="-")
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
      paste(gsub("\\..+$", "", input$enrich_data_file), paste(input$Gene_set3,paste0(input$which_group,".pdf"), sep = "_"), sep ="-")
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
    selectInput('Gene_set3', 'Gene Set', c(gene_set_list, "Custom gene set"))
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
    return(enrich_for_table(data = as.data.frame(enrich_viewer2()), H_t2g = Hallmark_enrich(), Gene_set = input$Gene_set3))
  })
  
  output$enrichment_result <- DT::renderDataTable({
    enrich_viewer_table()
  })
  
  output$download_enrichment_table = downloadHandler(
    filename = function() {
      paste(gsub("\\..+$", "", input$enrich_data_file), paste0(input$Gene_set3,"_table.txt"), sep ="-")
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
  
  enrich_motif <- reactive({
    if(input$motifButton > 0){
      return(MotifAnalysis(data= enrich_input(), Species = input$Species4, x = promoter()))
    }
  })
  output$motif_plot <- renderPlot({
    if(input$motifButton > 0 && !is.null(enrich_motif())){
      Motifplot(df2 = enrich_motif(), showCategory = input$enrich_showCategory, padj = input$promoter_padj)
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
      df$padj <- 10^(-df$negLog10Padj)
      return(df)
    }
  })
  output$motif_warning <- renderText({
    if(input$motifButton > 0){
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
      paste(gsub("\\..+$", "", input$enrich_data_file), 
            paste("Motif_upstream",input$promoter_upstream,"_downstream",input$promoter_downstream,"_table.txt",sep = ""), sep ="-")
    },
    content = function(file){write.table(motif_table(), file, row.names = F, sep = "\t", quote = F)}
  )
  output$download_motif_plot = downloadHandler(
    filename = function() {
      paste(gsub("\\..+$", "", input$enrich_data_file), 
            paste("Motif_upstream",input$promoter_upstream,"_downstream",input$promoter_downstream,".pdf",sep = ""), sep ="-")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        p1 <- Motifplot(df2 = enrich_motif(), showCategory = input$enrich_showCategory, padj = input$promoter_padj)
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
      paste(gsub("\\..+$", "", input$enrich_data_file), 
            paste("Motif_upstream",input$promoter_upstream,"_downstream",input$promoter_downstream,"_promoter_region.txt",sep = ""), sep ="-")
    },
    content = function(file){write.table(promoter_motif_region(), file, row.names = F, sep = "\t", quote = F)}
  )
  observeEvent(motif_table()[input$motif_result_rows_selected,], ({
    updateCollapse(session,id =  "Promoter_motif_collapse_panel", open="Promoter_motif_region_panel")
  }))
  
  
  #volcano navi------------------------------------------------------
  org5 <- reactive({
    return(org(Species = input$Species5))
  })
  org_code5 <- reactive({
    return(org_code(Species = input$Species5))
  })
  
  degresult <- reactive({
    withProgress(message = "Importing normalized count matrix, please wait",{
      tmp <- input$deg_file1$datapath
      if(is.null(input$deg_file1) && input$goButton5 > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/DEGexample.txt"
      return(read_df(tmp = tmp))
      incProgress(1)
    })
  })
  norm_count_input_for_deg <- reactive({
    withProgress(message = "Importing normalized count matrix, please wait",{
      tmp <- input$deg_file2$datapath
      if(is.null(input$deg_file2) && input$goButton5 > 0 )  tmp = "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/day0.txt"
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
      result <- data.frame(row.names = rownames(result),
                           log2FoldChange = result$log2FoldChange, padj = result$padj)
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
        if(str_detect(rownames(res)[1], "ENS")){
          my.symbols <- rownames(res)
          gene_IDs<-AnnotationDbi::select(org5(),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          rownames(gene_IDs) <- gene_IDs$Row.names
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
    dir_name <- paste(paste(gsub("\\..+$", "", input$deg_file1),sep = "-"), paste(input$fc4, input$fdr4,sep="_"),sep="_")
    return(dir_name)
  })
  
  DEG_uniqueID <- reactive({
    count <- norm_count_combined_DEG()
    if(input$Species5 != "not selected"){
      if(str_detect(rownames(count)[1], "ENS")){
        gene_IDs  <- gene_ID_DEG()
        data2 <- merge(count, gene_IDs, by= 0)
        rownames(data2) <- data2[,1]
        data2 <- data2[, - which(colnames(data2) == "Row.names.y")]
        data2$Unique_ID <- paste(data2$SYMBOL,data2$Row.names, sep = " - ")
        count <- data2[,-1]
      }
    }
    return(count)
  })
  
  GOI_DEG <- reactive({
    withProgress(message = "Preparing GOI list (about 10 sec)",{
      count <- DEG_uniqueID()
      if(is.null(count)){
        return(NULL)
      }else{
        if(str_detect(rownames(count)[1], "ENS")){
          if(input$Species5 != "not selected"){
            GOI <- count$Unique_ID
          }else GOI <- rownames(count)
        }else{
          if(input$Species5 != "not selected"){
            GOI <- rownames(count)
          }else GOI <- rownames(count)
        }
        return(GOI)
      }
    })
  })
  
  output$degGOI <- renderUI({
    if(is.null(DEG_uniqueID())){
      return(NULL)
    }else{
      withProgress(message = "Preparing GOI list (about 10 sec)",{
        selectizeInput("degGOI", "genes of interest (GOI)", c(GOI_DEG()),multiple = TRUE, options = list(delimiter = " ", create = T))
      })
    }
  })
  output$GOIreset_deg <- renderUI({
    actionButton("GOIreset_deg", "GOI reset")
  })
  observeEvent(input$GOIreset_deg, {
    withProgress(message = "Preparing GOI list (about 10 sec)",{
      updateSelectizeInput(session, "degGOI", choices = c(GOI_DEG()), 
                           selected = character(0),
                           options = list(create=TRUE, 'plugins' = list('remove_button'), persist = FALSE))
    })
  })
  
  
  output$deg_volcano_x <- renderUI({
    sliderInput("deg_xrange","X_axis range:",min = -100,
                max=100, step = 1,
                value = c(-10, 10))
  })
  output$deg_volcano_y <- renderUI({
    sliderInput("deg_yrange","Y_axis range:",min = 0, max= 300, step = 1,
                value = 100)
  })
  
  deg_volcano <- reactive({
    if(!is.null(input$deg_xrange)){
      data <- as.data.frame(DEG_uniqueID())
      if(!is.null(input$degGOI)){
        label_data <- input$degGOI
      }else label_data <- NULL
      if(input$volcano_inputType == "reverseON"){
        data$log2FoldChange <- -1 * data$log2FoldChange
      }
      data$color <- "NS"
      data$Row.names <- rownames(data)
      data$color[data$log2FoldChange < -log2(input$fc4) & data$padj < input$fdr4] <- "down"
      data$color[data$log2FoldChange > log2(input$fc4) & data$padj < input$fdr4] <- "up"
      data$padj[data$padj == 0] <- 10^(-300)
      if(!is.null(label_data)) {
        Color <- c("blue","green","darkgray","red")
        for(name in label_data){
          if(str_detect(rownames(data)[1], "ENS")){
            if(input$Species5 != "not selected"){
              data$color[data$Unique_ID == name] <- "GOI"
            }else{
              data$color[data$Row.names == name] <- "GOI"
            }
          }else{
            data$color[data$Row.names == name] <- "GOI"
          }
        }
      }else{
        Color <- c("blue","darkgray","red")
      }
      
      v <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) + geom_point(aes(color = color),size = 0.4)
      v <- v  + geom_vline(xintercept = c(-log2(input$fc4), log2(input$fc4)), linetype = c(2, 2), color = c("black", "black")) +
        geom_hline(yintercept = c(-log10(input$fdr4)), linetype = 2, color = c("black"))
      v <- v +theme_bw()+ scale_color_manual(values = Color)+
        theme(legend.position = "top" , legend.title = element_blank(),
              axis.text.x= ggplot2::element_text(size = 12),
              axis.text.y= ggplot2::element_text(size = 12),
              text = ggplot2::element_text(size = 12),
              title = ggplot2::element_text(size = 12)) +
        xlab("log2 fold change") + ylab("-log10(padj)") +
        xlim(input$deg_xrange)+
        ylim(c(0, input$deg_yrange))
      if(!is.null(label_data)) {
        if(str_detect(rownames(data)[1], "ENS")){
          if(input$Species5 != "not selected"){
            v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
            v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Unique_ID),
                                              box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1, fontface = "bold.italic")
          }else{
            v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
            v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Row.names),
                                              box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1, fontface = "bold.italic")
          }
        }else{
          v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
          v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Row.names),
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1, fontface = "bold.italic")
        }
      }
      return(v)
    }else return(NULL)
  })
  
  output$deg_volcano1 <- renderPlot({
    if(!is.null(input$deg_xrange)){
      if(is.null(DEG_uniqueID())){
        return(NULL)
      }else{
        withProgress(message = "volcano plot",{
          print(deg_volcano())
          incProgress(1)
        })
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
      if(str_detect(rownames(count)[1], "ENS")){
        if(input$Species5 != "not selected"){
          Unique_ID <- input$degGOI
          label_data <- as.data.frame(Unique_ID, row.names = Unique_ID)
          data <- merge(count, label_data, by="Unique_ID")
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
      data <- data[, - which(colnames(data) == "log2FoldChange")]
      data <- data[, - which(colnames(data) == "padj")]
      return(data)
    }
  })
  
  DEG_GOIheat <- reactive({
    data <- deg_GOIcount()
    if(is.null(data)){
      ht <- NULL
    }else{
      data.z <- genescale(data, axis=1, method="Z")
      data.z <- na.omit(data.z)
      ht <- GOIheatmap(data.z)
    }
    return(ht)
  })
  
  output$deg_GOIheatmap <- renderPlot({
    if(is.null(DEG_uniqueID())){
      return(NULL)
    }else{
      if(!is.null(input$degGOI)){
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
      p <- GOIboxplot(data = data)
    }
    return(p)
  })
  
  
  output$deg_GOIbox <- renderPlot({
    if(is.null(DEG_uniqueID())){
      return(NULL)
    }else{
      if(!is.null(input$degGOI)){
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
  
  #msigdbr-----------------
  msigdbr_list <- reactive({
    if(input$msigdbr_Species == ""){
      return(NULL)
    }else{
      withProgress(message = "Prepare gene sets",{
        msigdbr_list <- msigdbr(species = input$msigdbr_Species)
        return(msigdbr_list)
      })
    }
  })
  
  msig_list <- reactive({
    if(input$msigdbr_Species == "" || input$msigdbr_gene_set == ""){
      return(NULL)
    }else{
      data <- msigdbr_list() %>% 
        dplyr::filter(gs_name == input$msigdbr_gene_set) %>%
        as.data.frame()
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
        dorothea_list <- dorothea(species = input$dorothea_Species, org = org(input$dorothea_Species), type = Dorothea_type(), confidence = input$dorothea_confidence)
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
  
})