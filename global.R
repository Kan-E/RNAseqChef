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
library(shinyBS)
library(plotly,verbose=FALSE)
library('shinyjs', verbose = FALSE)
library(BiocManager)
library(clusterProfiler.dplyr)
library(dorothea)
options(repos = BiocManager::repositories())

msigdbr_species <- msigdbr_species()$species_name
read_df <- function(tmp){
  if(is.null(tmp)) {
    return(NULL)
  }else{
    if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
    if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
    if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
    return(df)
  }
}
dorothea <- function(species, type, org){
  if(species == "Mus musculus"){
    net <- dorothea::dorothea_mm
  }else{
    net <- dorothea::dorothea_hs
  }
  net2 <- net %>% filter(confidence != "D") %>% filter(confidence != "E")
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
  net3 <- data.frame(gs_name = net2$tf, entrez_gene = net2$ENTREZID, target = net2$target)
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
    net3 <- data.frame(gs_name = merge$gs_name, entrez_gene = merge$converted_entrez_gene)
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
data_3degcount1 <- function(data,result_Condm, result_FDR, specific, fc, fdr, basemean){
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
    if(specific == 1) {specific == collist[1]
      FC_xlab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[2], ")"))
      FC_ylab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[3], ")"))
      result_Condm$FC_x <- log2((result_Condm$C1 + 0.01)/(result_Condm$C2 + 0.01))
      result_Condm$FC_y <- log2((result_Condm$C1 + 0.01)/(result_Condm$C3 + 0.01))
      Pattern1 <- "Pattern4"
      Pattern2 <- "Pattern5"
    }
    if(specific == 2) {specific == collist[2]
      FC_xlab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[1], ")"))
      FC_ylab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[3], ")"))
      result_Condm$FC_x <- log2((result_Condm$C2 + 0.01)/(result_Condm$C1 + 0.01))
      result_Condm$FC_y <- log2((result_Condm$C2 + 0.01)/(result_Condm$C3 + 0.01))
      Pattern1 <- "Pattern3"
      Pattern2 <- "Pattern5"
    }
    if(specific == 3) {specific == collist[3]
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
      new.levels <- c( paste0(specific,"_up"), paste0(specific,"_down"), "NS" )
      col = c("red","blue", "darkgray")}
    if((sum(sig == 1) >= 1) && (sum(sig == 2) == 0)){
      new.levels <- c(paste0(specific,"_up: "), "NS" )
      col = c("red", "darkgray")}
    if((sum(sig == 1) == 0) && (sum(sig == 2) >= 1)){
      new.levels <- c(paste0(specific,"_down: "), "NS" )
      col = c("blue", "darkgray")}
    if((sum(sig == 1) == 0) && (sum(sig == 2) == 0)){
      new.levels <- c("NS")
      col = "darkgray"}
    
    data3$sig <- factor(data3$sig, labels = new.levels)
    return(data3)
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
cond3_scatter_plot <- function(data, data4, result_Condm, result_FDR, specific, fc, fdr, basemean){
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
    if(specific == 1) {specific == collist[1]
      FC_xlab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[2], ")"))
      FC_ylab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[3], ")"))
      result_Condm$FC_x <- log2((result_Condm$C1 + 0.01)/(result_Condm$C2 + 0.01))
      result_Condm$FC_y <- log2((result_Condm$C1 + 0.01)/(result_Condm$C3 + 0.01))
      Pattern1 <- "Pattern4"
      Pattern2 <- "Pattern5"
    }
    if(specific == 2) {specific == collist[2]
      FC_xlab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[1], ")"))
      FC_ylab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[3], ")"))
      result_Condm$FC_x <- log2((result_Condm$C2 + 0.01)/(result_Condm$C1 + 0.01))
      result_Condm$FC_y <- log2((result_Condm$C2 + 0.01)/(result_Condm$C3 + 0.01))
      Pattern1 <- "Pattern3"
      Pattern2 <- "Pattern5"
    }
    if(specific == 3) {specific == collist[3]
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
    complete_data <- stats::na.omit(data3)
    labs_data <- subset(complete_data, padj <= fdr & Row.names !=
                          "" & (FC_x)*(FC_y) >= log2(fc))
    labs_data<-  labs_data[sort(labs_data$FC_xy, decreasing = T, index=T)$ix,]
    labs_data <- dplyr::filter(labs_data, sig != "NS")
    labs_data2 <- utils::head(labs_data, 20)
    font.label <- data.frame(size=5, color="black", face = "plain")
    set.seed(42)
    FC_x <- FC_y <- sig <- Row.names <- padj <- NULL
    p <- ggplot(data3, aes(x = FC_x, y = FC_y)) + geom_point(aes(color = sig),size = 0.1)
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
    if(!is.null(labs_data)) {
      p <- p + ggrepel::geom_text_repel(data = labs_data2, mapping = aes(label = Row.names),
                                        box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), 
                                        force = 1, fontface = font.label$face,
                                        size = font.label$size/2, color = font.label$color)
    }
    p <- p  + geom_hline(yintercept = c(-log2(fc), log2(fc)), linetype = c(2, 2), color = c("black", "black"))+
      geom_vline(xintercept = c(-log2(fc), log2(fc)),linetype = c(2, 2), color = c("black", "black"))
    p <- p +
      theme_bw()+ scale_color_manual(values = col)+
      theme(legend.position = "top" , legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 10),
            axis.text.y= ggplot2::element_text(size = 10),
            text = ggplot2::element_text(size = 10),
            title = ggplot2::element_text(size = 10)) +
      xlab(FC_xlab) + ylab(FC_ylab)
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
    
    p <- plot_grid(p, ht, rel_widths = c(2, 1))
    return(p)
  }
}
enrichment3_1 <- function(data3, data4, Species, Gene_set, org, org_code, H_t2g){
  if(!is.null(Gene_set) && Species != "not selected"){
    cnet_list <- list()
    if(is.null(data4)){
      return(NULL)
    }else{
      if(Gene_set != "MSigDB Hallmark" && Gene_set != "Transcription factor targets" &&
         Gene_set != "DoRothEA regulon (activator)" && Gene_set != "DoRothEA regulon (repressor)"){
        if(Gene_set == "KEGG"){
          withProgress(message = "KEGG dotplot",{
            formula_res <- try(compareCluster(ENTREZID~sig, data=data4,fun="enrichKEGG", organism =org_code), silent = T)
            incProgress()
          })
        }
        if(Gene_set == "GO"){
          withProgress(message = "GO dotplot",{
            formula_res <- try(compareCluster(ENTREZID~sig, data=data4,fun="enrichGO", OrgDb=org), silent = T)
            incProgress()
          })
        }
        if (class(formula_res) == "try-error") formula_res <- NA
        return(formula_res)
      }else{
        withProgress(message = "dotplot",{
          for (name in unique(data3$sig)) {
            if (name != "NS"){
              em <- enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
              if (length(as.data.frame(em)$ID) == 0) {
                cnet1 <- NULL
              } else {
                cnet1 <- as.data.frame(setReadable(em, org, 'ENTREZID'))
                cnet1$Cluster <- name
                cnet1 <- cnet1[sort(cnet1$qvalue, decreasing = F, index=T)$ix,]
                cnet_list[[name]] = cnet1
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
          if((length(cnet1) = 0 ) && (length(cnet2) = 0)) {
            data <- NA
          }else if(is.null(cnet2)) {
            data <- cnet1
          }else data <- rbind(cnet1, cnet2)
          data["Description"] <- lapply(data["Description"], gsub, pattern="HALLMARK_", replacement = "")
          data$GeneRatio <- parse_ratio(data$GeneRatio)
          return(data)
          incProgress()
        })
      }
    }
  }
}
keggEnrichment2 <- function(data3, data4, formula_res, Species, Gene_set, org, org_code, H_t2g){
  if(!is.null(Gene_set) && Species != "not selected"){
    if(is.null(data4)){
      return(NULL)
    }else{
      cnet_list <- list()
      cnet_list2 <- list()
      if(Gene_set != "MSigDB Hallmark" && Gene_set != "Transcription factor targets" &&
         Gene_set != "DoRothEA regulon (activator)" && Gene_set != "DoRothEA regulon (repressor)"){
        if ((length(as.data.frame(formula_res)$Description) == 0) ||
            length(which(!is.na(unique(as.data.frame(formula_res)$qvalue)))) == 0) {
          d <- NULL
        } else{
          formula_res <- clusterProfiler.dplyr::filter(formula_res, !is.na(qvalue))
          d <- as.grob(dotplot(formula_res, showCategory=5, color ="qvalue" ,font.size=10))
        }
        for (name in unique(data3$sig)) {
          if (name != "NS"){
            if(Gene_set == "KEGG"){
              withProgress(message = "KEGG cnet plot",{
                kk1 <- enrichKEGG(data4$ENTREZID[data4$sig == name], organism =org_code)
                incProgress()
              })
            }
            if(Gene_set == "GO"){
              withProgress(mfessage = "GO cnet plot",{
                kk1 <- enrichGO(data4$ENTREZID[data4$sig == name], OrgDb=org)
                incProgress()
              })
            }
            if (is.null(kk1)) {
              cnet1 <- NULL
            } else cnet1 <- setReadable(kk1, org, 'ENTREZID')
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
      }else{
        for (name in unique(data3$sig)) {
          if (name != "NS"){
            em <- enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) == 0) {
              cnet1 <- NULL
            } else {
              cnet1 <- as.data.frame(setReadable(em, org, 'ENTREZID'))
              cnet1$Cluster <- name
              cnet1 <- cnet1[sort(cnet1$qvalue, decreasing = F, index=T)$ix,]
              cnet1 <- cnet1[1:5,]
              cnet_list2[[name]] = cnet1
            }
          }
        }
        if (length(cnet_list2) == 0){
          d <- NULL
        }else{
          if (length(cnet_list2) == 1) data <- cnet_list2[[1]]
          if (length(cnet_list2) == 2) data<- rbind(cnet_list2[[1]], cnet_list2[[2]])
          data <- dplyr::filter(data, !is.na(Cluster))
          data <- dplyr::filter(data, !is.na(Description))
          data["Description"] <- lapply(data["Description"], gsub, pattern="HALLMARK_", replacement = "")
          data$GeneRatio <- parse_ratio(data$GeneRatio)
          if ((length(data$Description) == 0) || length(which(!is.na(unique(data$qvalue))))==0) {
            d <- NULL
          } else{
            data <- dplyr::mutate(data, x = paste0(Cluster, eval(parse(text = "GeneRatio"))))
            data$x <- gsub(":","", data$x)
            data <- dplyr::arrange(data, x)
            idx <- order(data[["x"]], decreasing = FALSE)
            data$Description <- factor(data$Description,
                                       levels=rev(unique(data$Description[idx])))
            d <- as.grob(ggplot(data, aes(x = Cluster,y= Description,color=qvalue,size=GeneRatio))+
                           geom_point() +
                           scale_color_continuous(low="red", high="blue",
                                                  guide=guide_colorbar(reverse=TRUE)) +
                           scale_size(range=c(3, 8))+ theme_dose(font.size=8)+ylab(NULL))
          }}
        withProgress(message = "cnet plot",{
          for (name in unique(data3$sig)) {
            if (name != "NS"){
              em <- enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
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
          incProgress()
        })
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
GeneList_for_enrichment <- function(Species, Gene_set, org){
  if(Species != "not selected"){
    if(Gene_set != "MSigDB Hallmark" && Gene_set != "Transcription factor targets" &&
       Gene_set != "DoRothEA regulon (activator)" && Gene_set != "DoRothEA regulon (repressor)"){
      return(NULL)
    }else{
      switch (Species,
              "Mus musculus" = species <- "Mus musculus",
              "Homo sapiens" = species <- "Homo sapiens",
              "Rattus norvegicus" = species <- "Rattus norvegicus",
              "Xenopus laevis" = species <- "Xenopus laevis",
              "Drosophila melanogaster" = species <- "Drosophila melanogaster",
              "Caenorhabditis elegans" = species <- "Caenorhabditis elegans")
      if(Gene_set == "MSigDB Hallmark"){
        H_t2g <- msigdbr(species = species, category = "H") %>%
          dplyr::select(gs_name, entrez_gene) 
      }
      if(Gene_set == "Transcription factor targets"){
        H_t2g <- msigdbr(species = species, category = "C3")
        H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "TFT:GTRD" | gs_subcat == "TFT:TFT_Legacy") %>%
          dplyr::select(gs_name, entrez_gene)
      }
      if(Gene_set == "DoRothEA regulon (activator)"){
        H_t2g <- as_tibble(dorothea(species = Species,  type = "DoRothEA regulon (activator)", org = org)) %>%
          dplyr::select(gs_name, entrez_gene)
        H_t2g$entrez_gene <- as.integer(H_t2g$entrez_gene)
      }
      if(Gene_set == "DoRothEA regulon (repressor)"){
        H_t2g <- as_tibble(dorothea(species = Species,  type = "DoRothEA regulon (repressor)", org = org)) %>%
          dplyr::select(gs_name, entrez_gene)
        H_t2g$entrez_gene <- as.integer(H_t2g$entrez_gene)
      }
      return(H_t2g)
    }
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
              scales = "free", short.panel.labs = T)+
          theme(axis.text.x= element_text(size = 0),
                axis.text.y= element_text(size = 10),
                panel.background = element_rect(fill = "transparent", size = 0.5),
                title = element_text(size = 10),text = element_text(size = 20)))
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
      if(Gene_set != "MSigDB Hallmark" && Gene_set != "Transcription factor targets" &&
         Gene_set != "DoRothEA regulon (activator)" && Gene_set != "DoRothEA regulon (repressor)"){
        if(Gene_set == "KEGG"){
          withProgress(message = "KEGG enrichment analysis",{
            formula_res <- try(compareCluster(ENTREZID~Group, data=data3,
                                              fun="enrichKEGG", organism=org_code), silent = T)
            incProgress(1)
          })
        }
        if(Gene_set == "GO"){
          withProgress(message = "GO enrichment analysis",{
            formula_res <- try(compareCluster(ENTREZID~Group, data=data3,
                                              fun="enrichGO", OrgDb=org), silent =T)
            incProgress(1)
          })
        }
        if (length(as.data.frame(formula_res)$Description) == 0) {
          formula_res <- NULL
        }else{
          formula_res <-setReadable(formula_res, org, 'ENTREZID')
        }
        return(formula_res)
      }else{
        withProgress(message = "enrichment analysis",{
          df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
          for (name in unique(data3$Group)) {
            em <- enricher(data3$ENTREZID[data3$Group == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) != 0) {
              if(length(colnames(as.data.frame(em))) == 9){
                cnet1 <- as.data.frame(setReadable(em, org, 'ENTREZID'))
                cnet1$Group <- name
                df <- rbind(df, cnet1)
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
  } 
}