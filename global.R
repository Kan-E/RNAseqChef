library(shiny)
library(DT)
library(gdata)
library(rstatix)
library(multcomp)
library(tidyverse)
library(tools)
library(ggpubr)
library(venn)
library(ggrepel)
library(ggdendro)
library(ggplotify)
library(gridExtra)
library(cowplot)
library(DESeq2)
library(EBSeq)
library(ggnewscale)
library(edgeR)
library(IHW)
library(qvalue)
library(DEGreport)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(org.Dm.eg.db)
library(org.Ce.eg.db)
library(org.Xl.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(msigdbr)
library(genefilter)
library(ComplexHeatmap)
library(shinyBS, verbose = FALSE)
library(plotly,verbose=FALSE)
library('shinyjs', verbose = FALSE)
library(BiocManager)
library(clusterProfiler.dplyr)
library(dorothea)
library(umap)
library(biomaRt)
options(repos = BiocManager::repositories())

msigdbr_species <- msigdbr_species()$species_name
gene_set_list <- c("MSigDB Hallmark", "KEGG", "Reactome", "PID (Pathway Interaction Database)",
                   "BioCarta","WikiPathways", "GO biological process", 
                   "GO cellular component","GO molecular function", "Human phenotype ontology", 
                   "DoRothEA regulon (activator)", "DoRothEA regulon (repressor)",
                   "Transcription factor targets", "miRNA target")
species_list <- c("not selected", "Homo sapiens", "Mus musculus", "Rattus norvegicus", "Xenopus laevis", 
                  "Drosophila melanogaster", "Caenorhabditis elegans")
read_df <- function(tmp){
  if(is.null(tmp)) {
    return(NULL)
  }else{
    if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
    if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1,quote = "")
    if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1,quote = "")
    rownames(df) = gsub("\"", "", rownames(df))
    if(length(colnames(df)) != 0){
    if(str_detect(colnames(df)[1], "^X\\.")){
    colnames(df) = str_sub(colnames(df), start = 3, end = -2) 
    }
    }
    return(df)
  }
}
read_gene_list <- function(tmp){
  if(is.null(tmp)) {
    return(NULL)
  }else{
    if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE)
    if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",",quote = "")
    if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t",quote = "")
    rownames(df) = gsub("\"", "", rownames(df))
    if(str_detect(colnames(df)[1], "^X\\.")){
      colnames(df) = str_sub(colnames(df), start = 3, end = -2) 
    }
    return(df)
  }
}
gene_list_convert_for_enrichment <- function(data, Species){
    if(is.null(data) || Species == "not selected"){
      return(NULL)
    }else{
      df <- data.frame(GeneID = data[,1], Group = data[,2])
      my.symbols <- df$GeneID
      if(str_detect(df$GeneID[1], "ENS")){
        gene_IDs<-AnnotationDbi::select(org(Species),keys = my.symbols,
                                        keytype = "ENSEMBL",
                                        columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
        colnames(gene_IDs) <- c("GeneID","SYMBOL", "ENTREZID")
      }else{
        gene_IDs <- AnnotationDbi::select(org(Species), keys = my.symbols,
                                          keytype = "SYMBOL",
                                          columns = c("ENTREZID", "SYMBOL"))
        colnames(gene_IDs) <- c("GeneID","ENTREZID")
      }
      gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
      data <- merge(df, gene_IDs, by="GeneID")
      return(data)
    }
}
dorothea <- function(species, confidence = "recommend",type, org){
  if(species == "Mus musculus"){
    net <- dorothea::dorothea_mm
  }else{
    net <- dorothea::dorothea_hs
  }
  if(confidence == "recommend"){
  net2 <- net %>% filter(confidence != "D") %>% filter(confidence != "E")
  }else net2 <- net
  if(type == "DoRothEA regulon (activator)") net2 <- net2%>% filter(mor == 1)
  if(type == "DoRothEA regulon (repressor)") net2 <- net2%>% filter(mor == -1)
  my.symbols <- gsub("\\..*","", net2$target)
  gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                  keytype = "SYMBOL",
                                  columns = c("SYMBOL", "ENTREZID"))
  colnames(gene_IDs) <- c("target", "ENTREZID")
  gene_IDs <- gene_IDs %>% distinct(target, .keep_all = T)
  gene_IDs <- na.omit(gene_IDs)
  net2 <- merge(net2, gene_IDs, by="target")
  net3 <- data.frame(gs_name = net2$tf, entrez_gene = net2$ENTREZID, target = net2$target, confidence = net2$confidence)
  net3 <- dplyr::arrange(net3, gs_name)
  if(species != "Mus musculus" && species != "Homo sapiens"){
    genes <- net3$entrez_gene
    switch (species,
            "Rattus norvegicus" = set <- "rnorvegicus_gene_ensembl",
            "Xenopus laevis" = set <- "xtropicalis_gene_ensembl",
            "Drosophila melanogaster" = set <- "dmelanogaster_gene_ensembl",
            "Caenorhabditis elegans" = set <- "celegans_gene_ensembl")
    convert = useMart("ensembl", dataset = set, host="https://dec2021.archive.ensembl.org")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
    genes = getLDS(attributes = c("entrezgene_id"), filters = "entrezgene_id",
                   values = genes ,mart = human,
                   attributesL = c("entrezgene_id"),
                   martL = convert, uniqueRows=T)
    colnames(genes) <- c("entrez_gene", "converted_entrez_gene")
    genes <- genes %>% distinct(converted_entrez_gene, .keep_all = T)
    merge <- merge(net3, genes, by = "entrez_gene") 
    net3 <- data.frame(gs_name = merge$gs_name, entrez_gene = merge$converted_entrez_gene, confidence = merge$confidence)
    net3 <- dplyr::arrange(net3, gs_name)
  }
  return(net3)
}
org <- function(Species){
  if(Species != "not selected"){
    switch (Species,
            "Mus musculus" = org <- org.Mm.eg.db,
            "Homo sapiens" = org <- org.Hs.eg.db,
            "Rattus norvegicus" = org <- org.Rn.eg.db,
            "Xenopus laevis" = org <- org.Xl.eg.db,
            "Drosophila melanogaster" = org <- org.Dm.eg.db,
            "Caenorhabditis elegans" = org <- org.Ce.eg.db)
    return(org)
  }
}
org_code <- function(Species){
  if(Species != "not selected"){
    switch (Species,
            "Mus musculus" = org_code <- "mmu",
            "Homo sapiens" = org_code <- "hsa",
            "Rattus norvegicus" = org_code <- "rno",
            "Xenopus laevis" = org_code <- "xla",
            "Drosophila melanogaster" = org_code <- "dme",
            "Caenorhabditis elegans" = org_code <- "cel")
    return(org_code)
  }
}
PCAdata <- function(row_count, deg_norm_count){
  if(is.null(row_count)){
    return(NULL)
  }else{
    data <- deg_norm_count
    if(length(grep("SYMBOL", colnames(data))) != 0){
      data <- data[, - which(colnames(data) == "SYMBOL")]
    }
    if(length(grep("Unique_ID", colnames(data))) != 0){
      data <- data[, - which(colnames(data) == "Unique_ID")]
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
    return(pca$rotation)
  } 
}
PCAplot <- function(data){
  if(length(grep("SYMBOL", colnames(data))) != 0){
    data <- data[, - which(colnames(data) == "SYMBOL")]
  }
  if(length(grep("Unique_ID", colnames(data))) != 0){
    data <- data[, - which(colnames(data) == "Unique_ID")]
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
    theme(legend.position = "none",
          axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),axis.text.y=element_blank(),
          panel.grid.minor.x=element_blank(),
          axis.title.y=element_blank(),panel.background=element_rect(fill="white"))+
    coord_flip()+ scale_y_reverse(expand=c(0.6, 0))
  p2 <- plot_grid(g1, g2, g3, nrow = 1)
  return(p2)
}
umap_plot <- function(data, n_neighbors){
  umap <- umap::umap(t(data),n_neighbors = n_neighbors, random_state = 123)
  data2 <- umap$layout %>% as.data.frame()
  label<- colnames(data)
  label<- gsub("\\_.+$", "", label)
  p<- ggplot(data2, mapping = aes(V1,V2, color = label, label = colnames(data)))+
    geom_point()+geom_text_repel()+ xlab("UMAP_1") + ylab("UMAP_2")+
    theme(panel.background =element_rect(fill=NA,color=NA),panel.border = element_rect(fill = NA),
          aspect.ratio=1)
  return(p)
}
pdf_h <- function(rowlist){
  if ((length(rowlist) > 81) && (length(rowlist) <= 200)) pdf_hsize <- 22.5
  if ((length(rowlist) > 64) && (length(rowlist) <= 81)) pdf_hsize <- 20.25
  if ((length(rowlist) > 49) && (length(rowlist) <= 64)) pdf_hsize <- 18
  if ((length(rowlist) > 36) && (length(rowlist) <= 49)) pdf_hsize <- 15.75
  if ((length(rowlist) > 25) && (length(rowlist) <= 36)) pdf_hsize <- 13.5
  if ((length(rowlist) > 16) && (length(rowlist) <= 25)) pdf_hsize <- 11.5
  if ((length(rowlist) > 12) && (length(rowlist) <= 16)) pdf_hsize <- 9
  if ((length(rowlist) > 9) && (length(rowlist) <= 12)) pdf_hsize <- 7.5
  if ((length(rowlist) > 6) && (length(rowlist) <= 9)) pdf_hsize <- 7.5
  if ((length(rowlist) > 4) && (length(rowlist) <= 6)) pdf_hsize <- 6
  if (length(rowlist) == 4) pdf_hsize <- 6
  if (length(rowlist) == 3) pdf_hsize <- 3
  if (length(rowlist) == 2) pdf_hsize <- 3
  if (length(rowlist) == 1) pdf_hsize <- 3
  if (length(rowlist) > 200) pdf_hsize <- 30
  return(pdf_hsize)
}
pdf_w <- function(rowlist){
  if ((length(rowlist) > 81) && (length(rowlist) <= 200)) pdf_wsize <- 22.5
  if ((length(rowlist) > 64) && (length(rowlist) <= 81)) pdf_wsize <- 20.25
  if ((length(rowlist) > 49) && (length(rowlist) <= 64)) pdf_wsize <- 18
  if ((length(rowlist) > 36) && (length(rowlist) <= 49)) pdf_wsize <- 15.75
  if ((length(rowlist) > 25) && (length(rowlist) <= 36)) pdf_wsize <- 13.5
  if ((length(rowlist) > 16) && (length(rowlist) <= 25)) pdf_wsize <- 11.5
  if ((length(rowlist) > 12) && (length(rowlist) <= 16)) pdf_wsize <- 9
  if ((length(rowlist) > 9) && (length(rowlist) <= 12)) pdf_wsize <- 9
  if ((length(rowlist) > 6) && (length(rowlist) <= 9)) pdf_wsize <- 6.75
  if ((length(rowlist) > 4) && (length(rowlist) <= 6)) pdf_wsize <- 9
  if (length(rowlist) == 4) pdf_wsize <- 6
  if (length(rowlist) == 3) pdf_wsize <- 9
  if (length(rowlist) == 2) pdf_wsize <- 6
  if (length(rowlist) == 1) pdf_wsize <- 3
  if (length(rowlist) > 200) pdf_wsize <- 30
  return(pdf_wsize)
}
pdfSize_for_GOI <- paste(strong("Heatmap:"),"height = 10, width = 7 <br>", 
                         strong("Boxplot:"),"<br>",
                         "Gene number = 1,","height = 3, width = 3 <br>",
                         "Gene number = 2,","height = 3, width = 6 <br>",
                         "Gene number = 3,","height = 3, width = 9 <br>",
                         "Gene number = 4,","height = 6, width = 6 <br>",
                         "Gene number = 5 ~ 6,","height = 6, width = 9 <br>",
                         "Gene number = 7 ~ 9,","height = 7.5, width = 6.75 <br>",
                         "Gene number = 10 ~ 12,","height = 7.5, width = 9 <br>",
                         "Gene number = 13 ~ 16,","height = 9, width = 9 <br>",
                         "Gene number = 17 ~ 25,","height = 11.5, width = 11.5 <br>",
                         "Gene number = 26 ~ 36,","height = 13.5, width = 13.5 <br>",
                         "Gene number = 37 ~ 49,","height = 15.75, width = 15.75 <br>",
                         "Gene number = 50 ~ 64,","height = 18, width = 18 <br>",
                         "Gene number = 65 ~ 81,","height = 20.5, width = 20.5 <br>",
                         "Gene number = 82 ~ 200,","height = 22.5, width = 22.5 <br>",
                         "Gene number > 200,", "height = 30, width = 30 <br>")

data_3degcount1 <- function(data,result_Condm, result_FDR, specific, fc, fdr, basemean,result_list=NULL){
  if(is.null(data)){
    return(NULL)
  }else{
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
    if(str_detect(rownames(data)[1], "ENS")){
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
      data3 <- data2[,- which(colnames(data2) == "Pattern1")]
      data3 <- data3[,- which(colnames(data3) == "Pattern2")]
      data3 <- data3[,- which(colnames(data3) == "Pattern3")]
      data3 <- data3[,- which(colnames(data3) == "Pattern4")]
      data3 <- data3[,- which(colnames(data3) == "Pattern5")]
      data3 <- data3[,- which(colnames(data3) == "C1")]
      data3 <- data3[,- which(colnames(data3) == "C2")]
      data3 <- data3[,- which(colnames(data3) == "C3")]
      data3 <- data3[,- which(colnames(data3) == "PPDE")]
      Pattern <- rep(3, nrow(data3))
      Pattern[which(data3$MAP == "Pattern1")] = paste(collist[1], "=", collist[2], "=", collist[3])
      Pattern[which(data3$MAP == "Pattern2" & data3$FC_y > 0)] = paste(collist[1], "=", collist[2], ">", collist[3])
      Pattern[which(data3$MAP == "Pattern2" & data3$FC_y < 0)] = paste(collist[1], "=", collist[2], "<", collist[3])
      Pattern[which(data3$MAP == "Pattern3" & data3$FC_x > 0)] = paste(collist[1], "=", collist[3], ">", collist[2])
      Pattern[which(data3$MAP == "Pattern3" & data3$FC_x < 0)] = paste(collist[1], "=", collist[3], "<", collist[2])
      Pattern[which(data3$MAP == "Pattern4" & data3$FC_x > 0 & (data3$FC_y > 0))] = paste(collist[1], ">", collist[2], "=", collist[3])
      Pattern[which(data3$MAP == "Pattern4" & data3$FC_x < 0 & (data3$FC_y < 0))] = paste(collist[1], "<", collist[2], "=", collist[3])
      Pattern[which(data3$MAP == "Pattern5" & data3$FC_x > 0 & (data3$FC_y > 0) & ((data3$FC_x - data3$FC_y) > 0))] = paste(collist[1], ">", collist[2], ">", collist[3])
      Pattern[which(data3$MAP == "Pattern5" & data3$FC_x > 0 & (data3$FC_y > 0) & ((data3$FC_x - data3$FC_y) < 0))] = paste(collist[1], ">", collist[3], ">", collist[2])
      Pattern[which(data3$MAP == "Pattern5" & data3$FC_x < 0 & (data3$FC_y < 0) & ((data3$FC_x - data3$FC_y) > 0))] = paste(collist[1], "<", collist[3], "<", collist[2])
      Pattern[which(data3$MAP == "Pattern5" & data3$FC_x < 0 & (data3$FC_y < 0) & ((data3$FC_x - data3$FC_y) < 0))] = paste(collist[1], "<", collist[2], "<", collist[3])
      Pattern[which(data3$MAP == "Pattern5" & data3$FC_x > 0 & (data3$FC_y < 0))] = paste(collist[2], "<", collist[1], "<", collist[3])
      Pattern[which(data3$MAP == "Pattern5" & data3$FC_x < 0 & (data3$FC_y > 0))] = paste(collist[3], "<", collist[1], "<", collist[2])
      data3$Pattern <- Pattern
      if(!is.null(result_list)){
      data3 <- dplyr::select(data3, Row.names, Pattern, FDR, FC_x, FC_y, everything())
      colnames(data3)[3] <- "padj"
      colnames(data3)[4] <- FC_xlab
      colnames(data3)[5] <- FC_ylab
      data3 <- data3[,- which(colnames(data3) == "MAP")]
      rownames(data3) <- data3$Row.names
      data3 <- data3[,-1]
      return(data3)
    }else{
    result <- dplyr::filter(data3, apply(data3[,2:(Cond_1 + Cond_2 + Cond_3)],1,mean) > basemean)
    sig <- rep(3, nrow(result))
    sig[which(result$FDR <= fdr & result$FC_x < log2(1/fc) & result$FC_y < log2(1/fc) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 2
    sig[which(result$FDR <= fdr & result$FC_x > log2(fc) & result$FC_y > log2(fc) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 1
    data3 <- data.frame(Row.names = result$Row.names, sig = sig, Pattern = result$Pattern, FC_x_axis = result$FC_x,
                        FC_y_axis = result$FC_y, padj = result$FDR)
    if((sum(sig == 1) >= 1) && (sum(sig == 2) >= 1)){
      new.levels <- c( paste0(specific,"_high"), paste0(specific,"_low"), "NS" )
      col = c("red","blue", "darkgray")}
    if((sum(sig == 1) >= 1) && (sum(sig == 2) == 0)){
      new.levels <- c(paste0(specific,"_high: "), "NS" )
      col = c("red", "darkgray")}
    if((sum(sig == 1) == 0) && (sum(sig == 2) >= 1)){
      new.levels <- c(paste0(specific,"_low: "), "NS" )
      col = c("blue", "darkgray")}
    if((sum(sig == 1) == 0) && (sum(sig == 2) == 0)){
      new.levels <- c("NS")
      col = "darkgray"}
    
    data3$sig <- factor(data3$sig, labels = new.levels)
    return(data3)
    }
  }
}

data_3degcount2 <- function(data3, Species, org){
  if(is.null(data3)){
    return(NULL)
  }else{
    if(length(unique(data3$sig)) == 1){
      data4 <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
      return(data4)
    } else {
      data4 <- dplyr::filter(data3, sig != "NS")
      if(str_detect(data4$Row.names[1], "ENS")){
        if(Species != "not selected"){
          my.symbols <- data4$Row.names
          gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data4 <- merge(data4, gene_IDs, by="Row.names")
        }
      }else{
        if(Species != "not selected"){
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
  }
}
cond3_scatter_plot <- function(data, data4, result_Condm, result_FDR, specific, 
                               fc, fdr, basemean, y_axis=NULL, x_axis=NULL,
                               GOI=NULL, heatmap = TRUE, Species){
  if(is.null(data)){
    return(NULL)
  }else{
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
    if(str_detect(rownames(data)[1], "ENS")){
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
    if(str_detect(data3$Row.names[1], "ENS")){
      if(Species != "not selected"){
        data3 <- merge(data3, data, by="Row.names")
      }
    }
    complete_data <- stats::na.omit(data3)
    labs_data <- subset(complete_data, padj <= fdr & Row.names !=
                          "" & (FC_x)*(FC_y) >= log2(fc))
    labs_data<-  labs_data[sort(labs_data$FC_xy, decreasing = T, index=T)$ix,]
    labs_data <- dplyr::filter(labs_data, sig != "NS")
    labs_data2 <- utils::head(labs_data, 20)
    font.label <- data.frame(size=5, color="black", face = "plain")
    set.seed(42)
    FC_x <- FC_y <- Row.names <- padj <- NULL
    p <- ggplot(data3, aes(x = FC_x, y = FC_y)) + geom_point(aes(color = sig),size = 0.1)
    p <- p  + geom_hline(yintercept = c(-log2(fc), log2(fc)), linetype = c(2, 2), color = c("black", "black"))+
      geom_vline(xintercept = c(-log2(fc), log2(fc)),linetype = c(2, 2), color = c("black", "black"))
    p <- p +
      theme_bw()+ scale_color_manual(values = col)+
      theme(legend.position = "top" , legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 12),
            axis.text.y= ggplot2::element_text(size = 12),
            text = ggplot2::element_text(size = 15),
            title = ggplot2::element_text(size = 15)) +
      xlab(FC_xlab) + ylab(FC_ylab)
    if((sum(sig == 1) >= 1) && (sum(sig == 2) >= 1)){
      p <- p + geom_point(data=dplyr::filter(data3, sig == new.levels[1]),color = "red", size= 0.25 )
      p <- p + geom_point(data=dplyr::filter(data3, sig == new.levels[2]),color = "blue", size= 0.25 )
    }
    if((sum(sig == 1) >= 1) && (sum(sig == 2) == 0)){
      p <- p + geom_point(data=dplyr::filter(data3, sig == new.levels[1]),color = "red", size= 0.25 )
    }
    if((sum(sig == 1) == 0) && (sum(sig == 2) >= 1)){
      p <- p + geom_point(data=dplyr::filter(data3, sig == new.levels[1]),color = "blue", size= 0.25 )
    }
    if (heatmap==TRUE) {
      if(!is.null(labs_data)) {
        if(str_detect(data3$Row.names[1], "ENS")){
          if(Species != "not selected"){
            p <- p + ggrepel::geom_text_repel(data = labs_data2, mapping = aes(label = Unique_ID),
                                              box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), 
                                              force = 1, fontface = font.label$face,
                                              size = font.label$size/2, color = font.label$color)
          }else{
            p <- p + ggrepel::geom_text_repel(data = labs_data2, mapping = aes(label = Row.names),
                                              box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), 
                                              force = 1, fontface = font.label$face,
                                              size = font.label$size/2, color = font.label$color)
          }
        }else{
        p <- p + ggrepel::geom_text_repel(data = labs_data2, mapping = aes(label = Row.names),
                                          box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), 
                                          force = 1, fontface = font.label$face,
                                          size = font.label$size/2, color = font.label$color)
        }
      } 
    }
    if(!is.null(y_axis) && !is.null(x_axis)){
      p <- p +  xlim(x_axis) + ylim(y_axis)
    }
    if(!is.null(GOI)) {
      for(name in GOI){
        if(str_detect(data3$Row.names[1], "ENS")){
          if(Species != "not selected"){
            data3$color[data3$Unique_ID == name] <- "GOI"
          }else{
            data3$color[data3$Row.names == name] <- "GOI"
          }
        }else{
          data3$color[data3$Row.names == name] <- "GOI"
        }
      }
      if(str_detect(data3$Row.names[1], "ENS")){
        if(Species != "not selected"){
          p <- p + geom_point(data=dplyr::filter(data3, color == "GOI"),color="green", size=1)
          p <- p + ggrepel::geom_text_repel(data = dplyr::filter(data3, color == "GOI"), mapping = aes(label = Unique_ID),
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1)
        }else{
          p <- p + geom_point(data=dplyr::filter(data3, color == "GOI"),color="green", size=1)
          p <- p + ggrepel::geom_text_repel(data = dplyr::filter(data3, color == "GOI"), mapping = aes(label = Row.names),
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1)
        }
      }else{
        p <- p + geom_point(data=dplyr::filter(data3, color == "GOI"),color="green", size=1)
        p <- p + ggrepel::geom_text_repel(data = dplyr::filter(data3, color == "GOI"), mapping = aes(label = Row.names),
                                          box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1)
      }
    }
    if(heatmap == TRUE){
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
                            show_row_names = F, show_row_dend = T,column_names_side = "top"))
    }
    
    p <- plot_grid(p, ht, rel_widths = c(2, 1))
    }
    return(p)
  }
}
enrichment3_1 <- function(data3, data4, Species, Gene_set, org, org_code, H_t2g){
  if(!is.null(Gene_set) && Species != "not selected"){
    cnet_list <- list()
    if(is.null(data4)){
      return(NULL)
    }else{
        withProgress(message = "dotplot",{
          for (name in unique(data3$sig)) {
            if (name != "NS"){
              if(is.null(H_t2g)){
                cnet_list <- NULL
              }else{
                H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
              em <- enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
              if (length(as.data.frame(em)$ID) == 0) {
                cnet1 <- NULL
              } else {
                cnet1 <- as.data.frame(setReadable(em, org, 'ENTREZID'))
                cnet1$Group <- name
                cnet1 <- cnet1[sort(cnet1$qvalue, decreasing = F, index=T)$ix,]
                cnet_list[[name]] = cnet1
              }
              }
            }
          }
          if (length(cnet_list) == 2){
            cnet1 <- cnet_list[[1]]
            cnet2 <- cnet_list[[2]]}
          if (length(cnet_list) == 1){
            cnet1 <- cnet_list[[1]]
            cnet2 <- NULL}
          if (length(cnet_list) == 0){
            cnet1 <- NULL
            cnet2 <- NULL}
          if((length(colnames(cnet1)) != 2 ) && (length(colnames(cnet2)) != 2 )) data <- rbind(cnet1, cnet2)
          if((length(colnames(cnet1)) == 2 ) && (length(colnames(cnet2)) != 2 )) data <- cnet2
          if((length(colnames(cnet1)) != 2 ) && (length(colnames(cnet2)) == 2 )) data <- cnet1
          if((length(colnames(cnet1)) == 2 ) && (length(colnames(cnet2)) == 2 )) data <- NULL
          if(!is.null(data)) data$GeneRatio <- parse_ratio(data$GeneRatio)
          return(data)
          incProgress()
        })
    }
  }
}
keggEnrichment2 <- function(data3, data4, Species, Gene_set, org, org_code, H_t2g){
  if(!is.null(Gene_set) && Species != "not selected"){
    if(is.null(data4)){
      return(NULL)
    }else{
      cnet_list <- list()
      cnet_list2 <- list()
        for (name in unique(data3$sig)) {
          if (name != "NS"){
            if(is.null(H_t2g)){
              cnet_list2 <- NULL
            }else{
              H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
            em <- enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) == 0) {
              cnet1 <- NULL
            } else {
              cnet1 <- as.data.frame(setReadable(em, org, 'ENTREZID'))
              cnet1$Group <- name
              cnet1 <- cnet1[sort(cnet1$qvalue, decreasing = F, index=T)$ix,]
              cnet1 <- cnet1[1:5,]
              cnet_list2[[name]] = cnet1
            }
            }
          }
        }
        if (length(cnet_list2) == 0){
          d <- NULL
        }else{
          if (length(cnet_list2) == 1) data <- cnet_list2[[1]]
          if (length(cnet_list2) == 2) data<- rbind(cnet_list2[[1]], cnet_list2[[2]])
          data <- dplyr::filter(data, !is.na(Group))
          data <- dplyr::filter(data, !is.na(Description))
          data$GeneRatio <- parse_ratio(data$GeneRatio)
          if ((length(data$Description) == 0) || length(which(!is.na(unique(data$qvalue))))==0) {
            d <- NULL
          } else{
            data$Description <- gsub("_", " ", data$Description)
            data <- dplyr::mutate(data, x = paste0(Group, 1/(-log10(eval(parse(text = "qvalue"))))))
            data$x <- gsub(":","", data$x)
            data <- dplyr::arrange(data, x)
            idx <- order(data[["x"]], decreasing = FALSE)
            data$Description <- factor(data$Description,
                                       levels=rev(unique(data$Description[idx])))
            d <- as.grob(ggplot(data, aes(x = Group,y= Description,color=qvalue,size=GeneRatio))+
                           geom_point() +
                           scale_color_continuous(low="red", high="blue",
                                                  guide=guide_colorbar(reverse=TRUE)) +
                           scale_size(range=c(1, 6))+ theme_dose(font.size=8)+ylab(NULL)+xlab(NULL) +
                           scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
          }}
          for (name in unique(data3$sig)) {
            if (name != "NS"){
              if(is.null(H_t2g)){
                cnet_list <- NULL
              }else{
                H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
              em <- enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
              if (length(as.data.frame(em)$ID) == 0) {
                cnet1 <- NULL
              } else cnet1 <- setReadable(em, org, 'ENTREZID')
              if ((length(as.data.frame(cnet1)$ID) == 0) || 
                  length(which(!is.na(unique(as.data.frame(cnet1)$qvalue))))==0) {
                c <- NULL
              } else{
                c <- cnetplot(cnet1, cex_label_gene = 0.7, cex_label_category = 0.75,
                              cex_category = 0.75, colorEdge = TRUE)
                c <- try(as.grob(c + guides(edge_color = "none")))
                if(length(class(c)) == 1){
                  if(class(c) == "try-error") c <- NULL
                }
                cnet_list[[name]] = c
              }
              }
            }
          }
      if (length(cnet_list) == 2){
        cnet1 <- cnet_list[[1]]
        cnet2 <- cnet_list[[2]]}
      if (length(cnet_list) == 1){
        cnet1 <- cnet_list[[1]]
        cnet2 <- NULL}
      if (length(cnet_list) == 0){
        cnet1 <- NULL
        cnet2 <- NULL}
      p <- plot_grid(d, cnet1, cnet2, nrow = 1)
    }
  }else return(NULL)
}
enrich_for_table <- function(data, H_t2g, Gene_set){
  if(length(as.data.frame(data)$Description) == 0 || is.null(H_t2g)){
    return(NULL)
  }else{
    colnames(data)[1] <- "gs_name"
    H_t2g <- H_t2g %>% distinct(gs_name, .keep_all = T)
    data2 <- left_join(data, H_t2g, by="gs_name")  %>% as.data.frame()
    if(Gene_set == "DoRothEA regulon (activator)" || Gene_set == "DoRothEA regulon (repressor)"){
      data3 <- data.frame(Group = data2$Group, Gene_set_name = data2$gs_name, Confidence = data2$confidence,
                          Count = data2$Count, GeneRatio = data2$GeneRatio, BgRatio = data2$BgRatio, pvalue = data2$pvalue, 
                          p.adjust = data2$p.adjust, qvalue = data2$qvalue, GeneSymbol = data2$geneID)
    }else{
      if(Gene_set == "Custom gene set"){
        data3 <- data.frame(Group = data2$Group, Gene_set_name = data2$gs_name,
                            Count = data2$Count, GeneRatio = data2$GeneRatio, BgRatio = data2$BgRatio, pvalue = data2$pvalue, 
                            p.adjust = data2$p.adjust, qvalue = data2$qvalue, GeneSymbol = data2$geneID)
      }else{
        data3 <- data.frame(Group = data2$Group, Gene_set_name = data2$gs_name, ID = data2$gs_id, Description = data2$gs_description,
                            Count = data2$Count, GeneRatio = data2$GeneRatio, BgRatio = data2$BgRatio, pvalue = data2$pvalue, 
                            p.adjust = data2$p.adjust, qvalue = data2$qvalue, GeneSymbol = data2$geneID)
        
      }
        return(data3) 
    }
  }
}
GeneList_for_enrichment <- function(Species, Gene_set, org, Custom_gene_list){
  if(Species != "not selected" || is.null(Gene_set) || is.null(org)){
    switch (Species,
            "Mus musculus" = species <- "Mus musculus",
            "Homo sapiens" = species <- "Homo sapiens",
            "Rattus norvegicus" = species <- "Rattus norvegicus",
            "Xenopus laevis" = species <- "Xenopus laevis",
            "Drosophila melanogaster" = species <- "Drosophila melanogaster",
            "Caenorhabditis elegans" = species <- "Caenorhabditis elegans")
    if(Gene_set == "MSigDB Hallmark"){
      H_t2g <- msigdbr(species = species, category = "H") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description) 
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="HALLMARK_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "KEGG"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:KEGG") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="KEGG_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "Transcription factor targets"){
      H_t2g <- msigdbr(species = species, category = "C3")
      H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "TFT:GTRD" | gs_subcat == "TFT:TFT_Legacy") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
    }
    if(Gene_set == "DoRothEA regulon (activator)"){
      H_t2g <- as_tibble(dorothea(species = Species,  type = "DoRothEA regulon (activator)", org = org)) %>%
        dplyr::select(gs_name, entrez_gene, confidence)
      H_t2g$entrez_gene <- as.integer(H_t2g$entrez_gene)
    }
    if(Gene_set == "DoRothEA regulon (repressor)"){
      H_t2g <- as_tibble(dorothea(species = Species,  type = "DoRothEA regulon (repressor)", org = org)) %>%
        dplyr::select(gs_name, entrez_gene, confidence)
      H_t2g$entrez_gene <- as.integer(H_t2g$entrez_gene)
    }
    if(Gene_set == "Reactome"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="REACTOME_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "miRNA target"){
      H_t2g <- msigdbr(species = species, category = "C3")
      H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "MIR:MIRDB" | gs_subcat == "MIR:MIR_Legacy") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
    }
    if(Gene_set == "GO biological process"){
      H_t2g <- msigdbr(species = species, category = "C5", subcategory = "GO:BP") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOBP_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "GO cellular component"){
      H_t2g <- msigdbr(species = species, category = "C5", subcategory = "GO:CC") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOCC_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "GO molecular function"){
      H_t2g <- msigdbr(species = species, category = "C5", subcategory = "GO:MF") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOMF_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "Human phenotype ontology"){
      H_t2g <- msigdbr(species = species, category = "C5", subcategory = "HPO") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="HP_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "WikiPathways"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:WIKIPATHWAYS") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="WP_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "PID (Pathway Interaction Database)"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:PID") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="PID_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "BioCarta"){
      H_t2g <- msigdbr(species = species, category = "C2", subcategory = "CP:BIOCARTA") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="BIOCARTA_", replacement = "")
      H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
    }
    if(Gene_set == "Custom gene set"){
      if(!is.null(Custom_gene_list)){
        H_t2g <- gene_list_convert_for_enrichment(data= Custom_gene_list, Species = Species)
        H_t2g <- data.frame(gs_name = H_t2g$Group, entrez_gene = H_t2g$ENTREZID)
        H_t2g$gs_name <- gsub(":", "_", H_t2g$gs_name)
      }else H_t2g <- NULL
    }
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Tnf", replacement = "TNF")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Tgf", replacement = "TGF")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Pi3k_akt_mtor", replacement = "PI3K_Akt_mTOR")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Il6_", replacement = "IL6_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Il2_", replacement = "IL2_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Kras", replacement = "KRas")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Uv_r", replacement = "UV_r")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Dna_", replacement = "DNA_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Rna_", replacement = "RNA_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Mrna_", replacement = "mRNA_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="E2f", replacement = "E2F")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="P53", replacement = "p53")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="G2m", replacement = "G2M")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Mtorc", replacement = "mTORC")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Ecm_", replacement = "ECM_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="Abc_", replacement = "ABC_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="No1_", replacement = "NO1_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="_mirna", replacement = "_miRNA")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="gtpase_", replacement = "GTPase_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="rho_", replacement = "Rho_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="jnk_", replacement = "JNK_")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="_jak", replacement = "_JAK")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="_stat", replacement = "_STAT")
    H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="_nfkb", replacement = "_NFkB")
    return(H_t2g)
  }else return(NULL)
}
GOIboxplot <- function(data){
  collist <- gsub("\\_.+$", "", colnames(data))
  collist <- unique(collist)
  data$Row.names <- rownames(data)
  data <- data %>% gather(key=sample, value=value,-Row.names)
  data$sample <- gsub("\\_.+$", "", data$sample)
  data$Row.names <- as.factor(data$Row.names)
  data$sample <- factor(data$sample,levels=collist,ordered=TRUE)
  data$value <- as.numeric(data$value)
  p <- ggpubr::ggboxplot(data, x = "sample", y = "value",
                         fill = "sample", scales = "free",
                         add = "jitter",
                         xlab = FALSE, ylab = "Normalized_count", ylim = c(0, NA))
  p <- (facet(p, facet.by = "Row.names",
              panel.labs.background = list(fill = "transparent", color = "transparent"),
              scales = "free", short.panel.labs = T, panel.labs.font = list(size=15))+ 
          theme(axis.text.x = element_blank(),
                panel.background = element_rect(fill = "transparent", size = 0.5),
                title = element_text(size = 10),text = element_text(size = 12),
                axis.title.y = element_text(size=15),legend.text = element_text(size=15),
                legend.title = element_blank()))
  return(p)
}
enrich_viewer_forMulti1 <- function(df, Species, org){
  if(is.null(df) || Species == "not selected"){
    return(NULL)
  }else{
    my.symbols <- df$GeneID
    if(str_detect(df$GeneID[1], "ENS")){
      gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                      keytype = "ENSEMBL",
                                      columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
      colnames(gene_IDs) <- c("GeneID","SYMBOL", "ENTREZID")
    }else{
      gene_IDs <- AnnotationDbi::select(org, keys = my.symbols,
                                        keytype = "SYMBOL",
                                        columns = c("ENTREZID", "SYMBOL"))
      colnames(gene_IDs) <- c("GeneID","ENTREZID")
    }
    gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
    data <- merge(df, gene_IDs, by="GeneID")
    return(data)
  }
}
enrich_viewer_forMulti2 <- function(df, Species, Gene_set, org, org_code, H_t2g){
  if(!is.null(Gene_set)){
    data3 <- enrich_viewer_forMulti1(df = df, Species = Species, org = org)
    if(is.null(data3)){
      return(NULL)
    }else{
        withProgress(message = "enrichment analysis",{
          if(is.null(H_t2g)){
            df <- NULL
          }else{
            H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
          df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
          for (name in unique(data3$Group)) {
            em <- enricher(data3$ENTREZID[data3$Group == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) != 0) {
              if(length(colnames(as.data.frame(em))) == 9){
                cnet1 <- as.data.frame(setReadable(em, org, 'ENTREZID'))
                cnet1$Group <- name
                df <- rbind(df, cnet1)
              }
            }
          }
          }
          if(length(df$ID) !=0){
            df$GeneRatio <- parse_ratio(df$GeneRatio)
            return(df)
          }else return(NULL)
        })
      }
  } 
}
enrich_keggGO_global <- function(formula_res, Gene_set){
    if(!is.null(Gene_set)){
      if(is.null(formula_res)){
        return(NULL)
      }else{
            formula_res <- clusterProfiler.dplyr::filter(formula_res, !is.na(qvalue))
            p1 <- as.grob(dotplot(formula_res, color ="qvalue", font.size = 10))
          p <- plot_grid(p1)
          return(p)
        }
      }
}
enrich_genelist <- function(data, Gene_set, H_t2g, org, showCategory=5){
    if(!is.null(Gene_set)){
      if(is.null(data)){
        return(NULL)
      }else{
          if(is.null(H_t2g)){
            df <- NULL
          }else{
            H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
          df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
          for (name in unique(data$Group)) {
            em <- enricher(data$ENTREZID[data$Group == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) != 0) {
              if(length(colnames(as.data.frame(em))) == 9){
                cnet1 <- as.data.frame(setReadable(em, org, 'ENTREZID'))
                cnet1$Group <- name
                cnet1 <- cnet1[sort(cnet1$pvalue, decreasing = F, index=T)$ix,]
                if (length(cnet1$pvalue) > showCategory){
                  cnet1 <- cnet1[1:showCategory,]
                }
                df <- rbind(df, cnet1)
              }
            }
          }
          }
          if ((length(df$Description) == 0) || length(which(!is.na(unique(df$qvalue)))) == 0) {
            p1 <- NULL
          } else{
            df$GeneRatio <- parse_ratio(df$GeneRatio)
            df <- dplyr::filter(df, !is.na(qvalue))
            df$Description <- gsub("_", " ", df$Description)
            df <- dplyr::mutate(df, x = paste0(Group, 1/(-log10(eval(parse(text = "qvalue"))))))
            df$x <- gsub(":","", df$x)
            df <- dplyr::arrange(df, x)
            idx <- order(df[["x"]], decreasing = FALSE)
            df$Description <- factor(df$Description,
                                     levels=rev(unique(df$Description[idx])))
            p1 <- as.grob(ggplot(df, aes(x = Group,y= Description,color=qvalue,size=GeneRatio))+
                            geom_point() +
                            scale_color_continuous(low="red", high="blue",
                                                   guide=guide_colorbar(reverse=TRUE)) +
                            scale_size(range=c(1, 6))+ theme_dose(font.size=8)+ylab(NULL)+xlab(NULL)+
                            scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
            p <- plot_grid(p1)
            return(p)
          }
        }
    }
}
cnet_global <- function(data, group, Gene_set, H_t2g, org, org_code,showCategory=5){
    if(!is.null(Gene_set)){
      if(is.null(data) || is.null(group)){
        return(NULL)
      }else{
        data2 <- dplyr::filter(data, Group == group)
          if(is.null(H_t2g)){
            kk1 <- NULL
          }else{
            H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
          kk1 <- try(enricher(data2$ENTREZID, TERM2GENE=H_t2g2, qvalueCutoff = 0.05))
          if (class(kk1) == "try-error") kk1 <- NA
          }
        if(length(as.data.frame(kk1)$ID) == 0){
          cnet1 <- NULL
        } else {
          cnet1 <- setReadable(kk1, org, 'ENTREZID')
        }
        if (length(as.data.frame(cnet1)$ID) == 0) {
          p2 <- NULL
        } else{
          p2 <- try(as.grob(cnetplot(cnet1,
                                     cex_label_gene = 0.7, cex_label_category = 0.75,showCategory = showCategory,
                                     cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none")))
          if(length(class(p2)) == 1){
            if(class(p2) == "try-error") p2 <- NULL
          }
        }
        p <- plot_grid(p2)
        return(p)
      }
    }
}
dotplot_for_output <- function(data, plot_genelist, Gene_set, Species){
  if(!is.null(Gene_set)){
    if(is.null(data)){
      return(NULL)
    }else{
      withProgress(message = "Plot results",{
        if(Species != "not selected"){
            print(plot_genelist)
        }
        incProgress(1)
      })
    }
  }
}
cnet_for_output <- function(data, plot_data, Gene_set, Species){
    if(!is.null(Gene_set)){
      if(is.null(data)){
        return(NULL)
      }else{
        if(Species != "not selected"){
          withProgress(message = "cnet plot",{
            p <- plot_data
            print(p)
            incProgress(1)
          })
        }else return(NULL)
      }
    }
}



