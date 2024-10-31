#unique idの\tを除く (write.table前)

library(shiny)
library(shinyBS, verbose = FALSE)
library('shinyjs', verbose = FALSE)
library(DT)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(org.Dm.eg.db)
library(org.Ce.eg.db)
library(org.Xl.eg.db)
library(AnnotationDbi)
library(org.Bt.eg.db)
library(org.Cf.eg.db)
library(org.Dr.eg.db)
library(org.Gg.eg.db)
library(org.Mmu.eg.db)
library(org.Pt.eg.db)
library(org.Sc.sgd.db)
library(org.At.tair.db)
library(tidyverse)
library(rstatix)
library(multcomp)
library(tools)
library(ggpubr)
library(venn)
library(ggrepel)
library(ggplotify)
library(gridExtra)
library(cowplot)
library(DESeq2)
library(EBSeq)
library(ggnewscale)
library(edgeR)
library(qvalue)
library(DEGreport)
library(msigdbr)
library(ComplexHeatmap)
library(plotly,verbose=FALSE)
library(BiocManager)
library(umap)
library(biomaRt)
library(colorspace)
library(pdftools)
library(magick)
library(clue)
library(ggrastr) ##devtools::install_github('VPetukhov/ggrastr')
options(repos = BiocManager::repositories())
file.copy("Rmd/pair_report.Rmd",file.path(tempdir(),"pair_report.Rmd"), overwrite = TRUE)
file.copy("Rmd/pair_batch_report.Rmd",file.path(tempdir(),"pair_batch_report.Rmd"), overwrite = TRUE)
file.copy("Rmd/3conditions_report.Rmd",file.path(tempdir(),"3conditions_report.Rmd"), overwrite = TRUE)
file.copy("Rmd/multi_report.Rmd",file.path(tempdir(),"multi_report.Rmd"), overwrite = TRUE)
file.copy("dds.rds",file.path(tempdir(),"dds.rds"), overwrite = TRUE)
msigdbr_species <- msigdbr::msigdbr_species()$species_name
gene_set_list <- c("MSigDB Hallmark", "KEGG", "Reactome", "PID (Pathway Interaction Database)",
                   "BioCarta","WikiPathways", "GO biological process", 
                   "GO cellular component","GO molecular function", "Human phenotype ontology", 
                   "DoRothEA regulon (activator)", "DoRothEA regulon (repressor)",
                   "Transcription factor targets", "miRNA target","Position")
biomart_data <- read.table("https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/non-model.txt",sep = "\t", row.names = 1,header = T,quote = "")
biomart_plants <- read.table("https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/non-model_plants.txt",sep = "\t",header = T,quote = "")
biomart_fungi <- read.table("https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/non-model_fungi.txt",sep = "\t",header = T,quote = "")
biomart_metazoa <- read.table("https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/non-model_metazoa.txt",sep = "\t",header = T)
no_orgDb_animals<-c(biomart_data$Scientific_common_name)
no_orgDb_plants<-c(biomart_plants$Scientific_common_name)
no_orgDb_fungi<-c(biomart_fungi$Scientific_common_name)
no_orgDb_metazoa<-c(biomart_metazoa$Scientific_common_name)
no_orgDb <- c(no_orgDb_animals,no_orgDb_metazoa,no_orgDb_plants,no_orgDb_fungi)
species_list <- c("not selected", "Homo sapiens", "Mus musculus", "Rattus norvegicus", 
                  "Drosophila melanogaster", "Caenorhabditis elegans","Bos taurus","Canis lupus familiaris",
                  "Danio rerio","Gallus gallus","Macaca mulatta",
                  "Pan troglodytes","Saccharomyces cerevisiae","Xenopus laevis","Arabidopsis thaliana",
                  no_orgDb)
species_list_nonmodel <- no_orgDb
read_df <- function(tmp, Species=NULL){
  if(is.null(tmp)) {
    return(NULL)
  }else{
    if(sum(is.element(c("csv","txt","tsv","xlsx"), tools::file_ext(tmp))) == 0) validate("Error: the file extension is in an unexpected format.")
    if(tools::file_ext(tmp) == "xlsx") {
      df2 <- readxl::read_xlsx(tmp) 
      df2 <- as.data.frame(df2)
      df <- try(data.frame(row.names = df2[,1]),silent = T)
      if(class(df) != "try-error") {
        if(dim(df2)[2] == 2){
          df <- data.frame(row.names = df2[,1],a = df2[,2])
          colnames(df)[1] <- colnames(df2)[2]
        }else{
          rownames(df2) <- df2[,1]
        df <- df2[,-1]
        colnames(df) <- gsub("-",".",colnames(df))
        }
      }
    }
    if(tools::file_ext(tmp) == "csv") df <- try(read.csv(tmp, header=TRUE, sep = ",", row.names = 1,quote = ""))
    if(tools::file_ext(tmp) == "txt" || tools::file_ext(tmp) == "tsv") df <- try(read.table(tmp, header=TRUE, sep = "\t", row.names = 1,quote = ""))
    if(class(df) == "try-error") {
      if(tools::file_ext(tmp) == "xlsx") df <- try(as.data.frame(readxl::read_xlsx(tmp)))
      if(tools::file_ext(tmp) == "csv") df <- try(read.csv(tmp, header=TRUE, sep = ",",quote = ""))
      if(tools::file_ext(tmp) == "txt" || tools::file_ext(tmp) == "tsv") df <- try(read.table(tmp, header=TRUE, sep = "\t",quote = ""))
      if(class(df) != "try-error") {
        validate("Error: There are duplicated genes or 'NA' in the uploaded data. Please fix them.")
      }else{
        validate(paste0("Error: the uploaded data is in an unexpected format. The original error message is as follows:\n",print(df)))
      }
    }else{
    rownames(df) = gsub("\"", "", rownames(df))
    rownames(df) = gsub(":", ".", rownames(df))
    rownames(df) = gsub("\\\\", ".", rownames(df))
    if(length(grep("SYMBOL", colnames(df))) != 0){
      df <- df[, - which(colnames(df) == "SYMBOL")]
    }
    if(length(colnames(df)) != 0){
    if(str_detect(colnames(df)[1], "^X\\.")){
    colnames(df) = str_sub(colnames(df), start = 3, end = -2) 
    }
    }
    if(length(grep("padj", colnames(df))) == 0 || length(grep("log2FoldChange", colnames(df))) == 0){
    df[is.na(df)] <- 0
    }
    if(dim(df)[2] != 0){
      if(length(grep("Protein.Ids", colnames(df))) != 0 || length(grep("First.Protein.Description", colnames(df)))){
      df <- df %>% distinct(Genes, .keep_all = T)
      df[df$Genes == "",]$Genes <- gsub("\\_.+$", "", df[df$Genes == "",]$Protein.Ids)
      rownames(df) <- df$Genes
      if(length(grep("Protein.Group", colnames(df))) != 0) df <- df[, - which(colnames(df) == "Protein.Group")]
      if(length(grep("Protein.Ids", colnames(df))) != 0) df <- df[, - which(colnames(df) == "Protein.Ids")]
      if(length(grep("Protein.Names", colnames(df))) != 0) df <- df[, - which(colnames(df) == "Protein.Names")]
      if(length(grep("First.Protein.Description", colnames(df))) != 0) df <- df[, - which(colnames(df) == "First.Protein.Description")]
      if(length(grep("Genes", colnames(df))) != 0) df <- df[, - which(colnames(df) == "Genes")]
      if(length(grep("Species", colnames(df))) != 0) df <- df[, - which(colnames(df) == "Species")]
      df[is.na(df)] <- 0
      df <- df[!str_detect(rownames(df), ";"),]
      df <- df[rownames(df) != "",]
    }
    }
    if(sum(!str_detect(rownames(df),"\\.")) == 0 & !str_detect(rownames(df)[1],"chr")) rownames(df) <- gsub("\\..+$", "",rownames(df))
    return(df)
    }
  }
}
read_gene_list <- function(tmp){
  if(is.null(tmp)) {
    return(NULL)
  }else{
    if(tools::file_ext(tmp) == "xlsx") {
      df <- try(readxl::read_xlsx(tmp))
      df <- as.data.frame(df)
    }
    if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",",quote = "")
    if(tools::file_ext(tmp) == "txt" || tools::file_ext(tmp) == "tsv") df <- read.table(tmp, header=TRUE, sep = "\t",quote = "")
    rownames(df) = gsub("\"", "", rownames(df))
    if(str_detect(colnames(df)[1], "^X\\.")){
      colnames(df) = str_sub(colnames(df), start = 3, end = -2) 
    }
    if(rownames(df)[1] == 1){
      df <- data.frame(Gene = df[,1], Group = df[,2])
    }else{
      df <- data.frame(Gene = rownames(df), Group = df[,1])
    }
    if(sum(!str_detect(df$Gene,"\\.")) == 0) df$Gene <- gsub("\\..+$", "",df$Gene)
    return(df)
  }
}
anno_rep <- function(row){
  if(!str_detect(colnames(row)[1], "_1") && !str_detect(colnames(row)[1], "_2") &&
     !str_detect(colnames(row)[1], "_3") && !str_detect(colnames(row)[1], "_rep1") && 
     !str_detect(colnames(row)[1], "_rep2") && !str_detect(colnames(row)[1], "_rep3")){
    colnames(row) <- gsub("\\.[0-9]+$", "", colnames(row))
    name_list <- colnames(row) %>% sort()
    row <- row %>% dplyr::select(all_of(name_list)) 
    unique_col <- unique(colnames(row))
    total <- 0
    for(i in 1:length(unique_col)){
      cond <- length(which(colnames(row) == unique_col[i]))
      for(k in 1:cond){
        colnames(row)[total + k] <- paste0(colnames(row)[total + k], "_", k)
      }
      total <- total + cond
    }
  }
  return(row)
}
anno_rep_meta <- function(meta){
  if(is.null(meta)) {
    return(NULL)
  }else{
  if(!str_detect(meta[1,1], "_1") && !str_detect(meta[1,1], "_2") && !str_detect(meta[1,1], "_3") && 
     !str_detect(meta[1,1], "_rep1") && !str_detect(meta[1,1], "_rep2") && !str_detect(meta[1,1], "_rep3")){
    meta[,1] <- gsub("\\.[0-9]+$", "", meta[,1])
    if(colnames(meta)[1] != "characteristics"){
    if(length(grep("characteristics", colnames(meta))) != 0){
      meta <- meta[, - which(colnames(meta) == "characteristics")]
    }
    }
    colnames(meta)[1] <- "characteristics"
    meta <- meta %>% dplyr::arrange(characteristics)
    unique_col <- unique(meta[,1])
    total <- 0
    for(i in 1:length(unique_col)){
      cond <- length(which(meta[,1] == unique_col[i]))
      for(k in 1:cond){
        meta[total + k,1] <- paste0(meta[total + k,1], "_", k)
      }
      total <- total + cond
    }
  }
  return(meta)
  }
}

gene_list_convert_for_enrichment <- function(gene_type,data, Ortholog,org,Species){
    if(is.null(data) || Species == "not selected"){
      return(NULL)
    }else{
      df <- data.frame(GeneID = data[,1], Group = data[,2])
      df$GeneID <- gsub("\\..*","", df$GeneID)
      my.symbols <- df$GeneID
      if(gene_type != "SYMBOL"){
        if(sum(is.element(no_orgDb, Species)) == 1){
          gene_IDs <- Ortholog
        }else{
        if(str_detect(my.symbols[1], "^AT.G")) key = "TAIR" else key = "ENSEMBL"
        gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                        keytype = key,
                                        columns = c(key,"SYMBOL", "ENTREZID"))
        }
        colnames(gene_IDs) <- c("GeneID","SYMBOL", "ENTREZID")
      }else{
        if(sum(is.element(no_orgDb, Species)) == 1){
          gene_IDs <- Ortholog
          gene_IDs <- gene_IDs[,-1]
        }else{
        gene_IDs <- AnnotationDbi::select(org, keys = my.symbols,
                                          keytype = "SYMBOL",
                                          columns = c("ENTREZID", "SYMBOL"))
        }
        colnames(gene_IDs) <- c("GeneID","ENTREZID")
      }
      gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
      data <- merge(df, gene_IDs, by="GeneID")
      return(data)
    }
}
dorothea <- function(species, confidence = "recommend",type){
  if(species == "Mus musculus"){
    net <- dorothea::dorothea_mm
    spe <- "Mus musculus"
  }else{
    net <- dorothea::dorothea_hs
    spe <- "Homo sapiens"
  }
  if(confidence == "recommend"){
  net2 <- net %>% dplyr::filter(confidence != "D") %>% dplyr::filter(confidence != "E")
  }else net2 <- net
  if(type == "DoRothEA regulon (activator)") net2 <- net2 %>% dplyr::filter(mor == 1)
  if(type == "DoRothEA regulon (repressor)") net2 <- net2 %>% dplyr::filter(mor == -1)
  my.symbols <- gsub("\\..*","", net2$target)
  gene_IDs<-AnnotationDbi::select(org(spe),keys = my.symbols,
                                  keytype = "SYMBOL",
                                  columns = c("SYMBOL", "ENTREZID"))
  colnames(gene_IDs) <- c("target", "ENTREZID")
  gene_IDs <- gene_IDs %>% distinct(target, .keep_all = T)
  gene_IDs <- na.omit(gene_IDs)
  net2 <- merge(net2, gene_IDs, by="target")
  net3 <- data.frame(gs_name = net2$tf, entrez_gene = net2$ENTREZID, target = net2$target, confidence = net2$confidence)
  net3 <- dplyr::arrange(net3, gs_name)
  if(species != "Mus musculus" && species != "Homo sapiens"){
    withProgress(message = paste0("Gene ID conversion from human to ", species, " for the regulon gene set. It takes a few minutes."),{
    genes <- net3$entrez_gene
    if(species == "Saccharomyces cerevisiae" || species == "Arabidopsis thaliana") validate("'Saccharomyces cerevisiae' cannot use DoRothEA regulon.")
    switch (species,
            "Rattus norvegicus" = set <- "rnorvegicus_gene_ensembl",
            "Xenopus tropicalis" = set <- "xtropicalis_gene_ensembl",
            "Drosophila melanogaster" = set <- "dmelanogaster_gene_ensembl",
            "Caenorhabditis elegans" = set <- "celegans_gene_ensembl",
            "Anolis carolinensis" = set <- "acarolinensis_gene_ensembl",
            "Bos taurus" = set <- "btaurus_gene_ensembl",
            "Canis lupus familiaris" = set <- "clfamiliaris_gene_ensembl",
            "Danio rerio" = set <- "drerio_gene_ensembl",
            "Equus caballus" = set <- "ecaballus_gene_ensembl",
            "Felis catus" = set <- "fcatus_gene_ensembl",
            "Gallus gallus" = set <- "ggallus_gene_ensembl",
            "Macaca mulatta" = set <- "mmulatta_gene_ensembl",
            "Monodelphis domestica" = set <- "mdomestica_gene_ensembl",
            "Ornithorhynchus anatinus" = set <- "oanatinus_gene_ensembl",
            "Pan troglodytes" = set <- "ptroglodytes_gene_ensembl")
    convert = useMart("ensembl", dataset = set, host="https://dec2021.archive.ensembl.org")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
    genes2 = getLDS(attributes = c("entrezgene_id"), filters = "entrezgene_id",
                   values = genes ,mart = human,
                 attributesL = c("entrezgene_id"),
                   martL = convert, uniqueRows=T)
    colnames(genes2) <- c("entrez_gene", "converted_entrez_gene")
    genes2 <- genes2 %>% distinct(converted_entrez_gene, .keep_all = T)
    merge <- merge(net3, genes2, by = "entrez_gene") 
    net3 <- data.frame(gs_name = merge$gs_name, entrez_gene = merge$converted_entrez_gene, confidence = merge$confidence)
    net3 <- dplyr::arrange(net3, gs_name)
    })
  }
  return(net3)
}
ensembl_archive <- c("https://dec2021.archive.ensembl.org",
                     "https://www.ensembl.org")
ensembl_archive_plants <- c("https://plants.ensembl.org",
                            "https://eg52-plants.ensembl.org")
ensembl_archive_fungi <- c("https://dec2021-fungi.ensembl.org",
                     "https://fungi.ensembl.org")
ensembl_archive_metazoa <- c("https://metazoa.ensembl.org",
                           "https://eg52-metazoa.ensembl.org")
orgDb_list <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus", 
                "Drosophila melanogaster", "Caenorhabditis elegans","Bos taurus","Canis lupus familiaris",
                "Danio rerio","Gallus gallus","Macaca mulatta","Pan troglodytes","Saccharomyces cerevisiae")
no_org_ID <- function(count=NULL,gene_list=NULL,Species,Ortholog,Biomart_archive){
  if(Species != "not selected"){
    if(!is.null(Ortholog)){
    if(sum(is.element(no_orgDb, Species)) == 1){
      withProgress(message = "preparing a gene annotation. It takes a few minutes.",{
        mart <- "ensembl"
        bio_data <- biomart_data
        if(Ortholog == "Arabidopsis thaliana") {
          mart <- "plants_mart"
          bio_data <- biomart_plants
        }
        if(Ortholog == "Saccharomyces cerevisiae") {
          mart <- "fungi_mart"
          bio_data <- biomart_fungi
        }
        if(Ortholog == "Drosophila melanogaster" || Ortholog == "Caenorhabditis elegans") {
          mart <- "metazoa_mart"
          bio_data <- biomart_metazoa
        }
        db <- useMart(mart, host=Biomart_archive)
        ensembl <- dplyr::filter(bio_data, Scientific_common_name == Species)$dataset
        use <- useDataset(ensembl, mart = db)
        genes_ensembl <- try(getBM(attributes = c("ensembl_gene_id","external_gene_name"), mart = use))
        if(Ortholog == "Mus musculus") ortho <- "mmusculus_gene_ensembl"
        if(Ortholog == "Homo sapiens") ortho <- "hsapiens_gene_ensembl"
        if(Ortholog == "Rattus norvegicus") ortho <- "rnorvegicus_gene_ensembl"
        if(Ortholog == "Drosophila melanogaster") ortho <- "dmelanogaster_eg_gene"
        if(Ortholog == "Caenorhabditis elegans") ortho <- "celegans_eg_gene"
        if(Ortholog == "Bos taurus") ortho <- "btaurus_gene_ensembl"
        if(Ortholog == "Canis lupus familiaris") ortho <- "clfamiliaris_gene_ensembl"
        if(Ortholog == "Danio rerio") ortho <- "drerio_gene_ensembl"
        if(Ortholog == "Gallus gallus") ortho <- "ggallus_gene_ensembl"
        if(Ortholog == "Macaca mulatta") ortho <- "mmulatta_gene_ensembl"
        if(Ortholog == "Pan troglodytes") ortho <- "ptroglodytes_gene_ensembl"
        if(Ortholog == "Saccharomyces cerevisiae") ortho <- "scerevisiae_eg_gene"
        if(Ortholog == "Arabidopsis thaliana") ortho <- "athaliana_eg_gene"
        if(class(genes_ensembl) == "try-error") {
          type <- "ENSEMBL"
        }else{
        if(is.null(gene_list)){
        count$ensembl_gene_id <- rownames(count)
        count$external_gene_name <- rownames(count)
        ENSEMBL <- merge(genes_ensembl,count,by="ensembl_gene_id")
        SYMBOL <- merge(genes_ensembl,count,by="external_gene_name")
        }else{
          gene_list$ensembl_gene_id <- gene_list[,1]
          gene_list$external_gene_name <- gene_list[,1]
          ENSEMBL <- merge(genes_ensembl,gene_list,by="ensembl_gene_id")
          SYMBOL <- merge(genes_ensembl,gene_list,by="external_gene_name")
        }
        if(dim(ENSEMBL)[1] > dim(SYMBOL)[1]) type <- "ENSEMBL" else type <- "SYMBOL"
        if(dim(ENSEMBL)[1] == 0 && dim(SYMBOL)[1] == 0) validate("Cannot identify gene IDs. Please check the 'Species' and use the 'Official gene symbol' or 'ENSEMBL ID' for gene names.")
        }
        if(type == "ENSEMBL"){
          attribute <- c("ensembl_gene_id")
          colname <- c("ENSEMBL","SYMBOL","ENTREZID")
          colname1<- c("ENSEMBL")
          genes <- getBM(attributes = c("ensembl_gene_id"), mart = use)
          filter <- "ensembl_gene_id"
        }else{
          attribute <- c("ensembl_gene_id","external_gene_name")
          colname <- c("Original_symbol","SYMBOL","ENTREZID")
          colname1<- c("Original_symbol")
          genes <- getBM(attributes = c("external_gene_name"), mart = use)
          filter <- "external_gene_name"
        }
        ortho_mart = useMart(mart, dataset = ortho, host=Biomart_archive)
        genes2 = try(getLDS(attributes = c("ensembl_gene_id"),
                        values = genes ,mart = use,filters = filter,
                        attributesL = c("external_gene_name","entrezgene_id"),
                        martL = ortho_mart, uniqueRows=T))
        if(class(genes2) == "try-error") {
          genes4 = try(getLDS(attributes = c("ensembl_gene_id"),
                              values = genes ,mart = use,filters = filter,
                              attributesL = c("external_gene_name","ensembl_gene_id"),
                              martL = ortho_mart, uniqueRows=T))
          if(class(genes4) == "try-error") {
            validate("biomart has encountered an unexpected server error.
                    Please try using a different 'biomart host' or try again later.")
          }else{
            if(Ortholog == "Drosophila melanogaster") org <- org.Dm.eg.db
            if(Ortholog == "Caenorhabditis elegans") org <- org.Ce.eg.db
            my.symbols <- genes4[,3]
            gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                            keytype = "ENSEMBL",
                                            columns = c("ENSEMBL","ENTREZID"))
            colnames(gene_IDs) <- c("Gene.stable.ID.1","ENTREZID")
            gene_IDs <- gene_IDs %>% distinct(Gene.stable.ID.1, .keep_all = T)
            genes3 <- merge(genes4, gene_IDs, by="Gene.stable.ID.1")
            genes2 <- genes3[,-1]
          }
        }
        colnames(genes) <- colname1
        colnames(genes2) <- colname
        gene3<-merge(genes,genes2,by=colname1,all=T)
        if(type == "ENSEMBL") org <- gene3 %>% distinct(ENSEMBL, .keep_all = T) else org <- gene3 %>% distinct(Original_symbol, .keep_all = T)
      })
    }
    return(org)
    }
  }
}
org <- function(Species, Ortholog=NULL){
  if(Species != "not selected"){
    if(sum(is.element(no_orgDb, Species)) == 0){
    switch (Species,
            "Mus musculus" = org <- org.Mm.eg.db,
            "Homo sapiens" = org <- org.Hs.eg.db,
            "Rattus norvegicus" = org <- org.Rn.eg.db,
            "Xenopus laevis" = org <- org.Xl.eg.db,
            "Drosophila melanogaster" = org <- org.Dm.eg.db,
            "Caenorhabditis elegans" = org <- org.Ce.eg.db,
            "Bos taurus" = org <- org.Bt.eg.db,
            "Canis lupus familiaris" = org <- org.Cf.eg.db,
            "Danio rerio" = org <- org.Dr.eg.db,
            "Gallus gallus" = org <- org.Gg.eg.db,
            "Macaca mulatta" = org <- org.Mmu.eg.db,
            "Pan troglodytes" = org <- org.Pt.eg.db,
            "Saccharomyces cerevisiae" = org <- org.Sc.sgd.db,
            "Arabidopsis thaliana" = org <- org.At.tair.db)
    }else{
    switch (Ortholog,
            "Mus musculus" = org <- org.Mm.eg.db,
            "Homo sapiens" = org <- org.Hs.eg.db,
            "Rattus norvegicus" = org <- org.Rn.eg.db,
            "Drosophila melanogaster" = org <- org.Dm.eg.db,
            "Caenorhabditis elegans" = org <- org.Ce.eg.db,
            "Bos taurus" = org <- org.Bt.eg.db,
            "Canis lupus familiaris" = org <- org.Cf.eg.db,
            "Danio rerio" = org <- org.Dr.eg.db,
            "Gallus gallus" = org <- org.Gg.eg.db,
            "Macaca mulatta" = org <- org.Mmu.eg.db,
            "Pan troglodytes" = org <- org.Pt.eg.db,
            "Saccharomyces cerevisiae" = org <- org.Sc.sgd.db,
            "Arabidopsis thaliana" = org <- org.At.tair.db)
    }
    return(org)
  }
}
org_code <- function(Species,Ortholog){
  if(Species != "not selected"){
    if(sum(is.element(no_orgDb, Species)) == 0){
    switch (Species,
            "Mus musculus" = org_code <- "mmu",
            "Homo sapiens" = org_code <- "hsa",
            "Rattus norvegicus" = org_code <- "rno",
            "Xenopus tropicalis" = org_code <- "xtr",
            "Xenopus laevis" = org_code <- "xla",
            "Drosophila melanogaster" = org_code <- "dme",
            "Caenorhabditis elegans" = org_code <- "cel",
            "Anolis carolinensis" = org_code <- "acs",
            "Bos taurus" = org_code <- "bta",
            "Canis lupus familiaris" = org_code <- "cfa",
            "Danio rerio" = org_code <- "dre",
            "Equus caballus" = org_code <- "ecb",
            "Felis catus" = org_code <- "fca",
            "Gallus gallus" = org_code <- "gga",
            "Macaca mulatta" = org_code <- "mcc",
            "Monodelphis domestica" = org_code <- "mdo",
            "Ornithorhynchus anatinus" = org_code <- "oaa",
            "Pan troglodytes" = org_code <- "ptr",
            "Heterocephalus glaber" = org_code <- "hgl",
            "Saccharomyces cerevisiae" = org_code <- "sce",
            "Arabidopsis thaliana" = org_code <- "ath")
    }else{
      switch (Ortholog,
              "Mus musculus" = org_code <- "mmu",
              "Homo sapiens" = org_code <- "hsa",
              "Rattus norvegicus" = org_code <- "rno",
              "Xenopus tropicalis" = org_code <- "xtr",
              "Xenopus laevis" = org_code <- "xla",
              "Drosophila melanogaster" = org_code <- "dme",
              "Caenorhabditis elegans" = org_code <- "cel",
              "Anolis carolinensis" = org_code <- "acs",
              "Bos taurus" = org_code <- "bta",
              "Canis lupus familiaris" = org_code <- "cfa",
              "Danio rerio" = org_code <- "dre",
              "Equus caballus" = org_code <- "ecb",
              "Felis catus" = org_code <- "fca",
              "Gallus gallus" = org_code <- "gga",
              "Macaca mulatta" = org_code <- "mcc",
              "Monodelphis domestica" = org_code <- "mdo",
              "Ornithorhynchus anatinus" = org_code <- "oaa",
              "Pan troglodytes" = org_code <- "ptr",
              "Heterocephalus glaber" = org_code <- "hgl",
              "Saccharomyces cerevisiae" = org_code <- "sce",
              "Arabidopsis thaliana" = org_code <- "ath")
    }
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
PCAplot <- function(data,legend=NULL){
  if(length(grep("SYMBOL", colnames(data))) != 0){
    data <- data[, - which(colnames(data) == "SYMBOL")]
  }
  if(length(grep("Unique_ID", colnames(data))) != 0){
    data <- data[, - which(colnames(data) == "Unique_ID")]
  }
  pca <- prcomp(data, scale. = T)
  label<- colnames(data)
  lab_x <- paste(summary(pca)$importance[2,1]*100,
                 "% of variance)", sep = "")
  lab_x <- paste("PC1 (", lab_x, sep = "")
  lab_y <- paste(summary(pca)$importance[2,2]*100,
                 "% of variance)", sep = "")
  lab_y <- paste("PC2 (", lab_y, sep = "")
  pca$rotation <- as.data.frame(pca$rotation)
  if(!is.null(legend)) {
    if(legend == "Legend"){
      legend_position <- "top" 
      label2<- NULL
    }else{
      legend_position <- "none"
      label2<- label
  }
  }
  g1 <- ggplot(pca$rotation,aes(x=pca$rotation[,1],
                                y=pca$rotation[,2],
                                col=gsub("\\_.+$", "", label), label = label2)) +
    geom_point()+
    theme(panel.background =element_rect(fill=NA,color=NA),
          panel.border = element_rect(fill = NA)) +
    xlab(lab_x) + ylab(lab_y) +
    theme(legend.position=legend_position, aspect.ratio=1)+ 
    guides(color=guide_legend(title=""))
  if(!is.null(legend)){
    if(legend == "Label") g1 <- g1 + geom_text_repel(show.legend = NULL)
  }
  rho <- cor(data,method="spearman")
  d <- dist(1-rho)
  mds <- try(as.data.frame(cmdscale(d)))
  if(class(mds) != "try-error"){
  g2 <- ggplot(mds, aes(x = mds[,1], y = mds[,2],
                        col = gsub("\\_.+$", "", label), label = label2)) +
    geom_point()+
    theme(panel.background =element_rect(fill=NA,color=NA),
          panel.border = element_rect(fill = NA)) +
    xlab("dim 1") + ylab("dim 2") +
    theme(legend.position=legend_position, aspect.ratio=1)+ 
    guides(color=guide_legend(title=""))
  }else g2 <- NULL
  if(!is.null(legend)){
    if(legend == "Label") g2 <- g2 + geom_text_repel(show.legend = NULL)
  }
  x <- NULL
  y <- NULL
  xend <- NULL
  yend <- NULL
  data.t <- t(data)
  hc <- hclust(dist(data.t), "ward.D2")
  dendr <- ggdendro::dendro_data(hc, type="rectangle")
  g3 <- ggplot() +
    geom_segment(data=ggdendro::segment(dendr),
                 aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data=ggdendro::label(dendr),
              aes(x, y, label=label, hjust=0), size=3) +
    theme(legend.position = legend_position,
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
umap_plot <- function(data, n_neighbors,lab=NULL){
  umap <- umap::umap(t(data),n_neighbors = n_neighbors, random_state = 123)
  data2 <- umap$layout %>% as.data.frame()
  label<- colnames(data)
  if(!is.null(lab)) {
    if(lab == "Legend"){
      label2<- NULL
    }else{
      label2<- label
    }
  }
  p<- ggplot(data2, mapping = aes(V1,V2, color = gsub("\\_.+$", "", label), label = label2))+
    geom_point()+xlab("UMAP_1") + ylab("UMAP_2")+
    theme(panel.background =element_rect(fill=NA,color=NA),panel.border = element_rect(fill = NA),
          aspect.ratio=1)+ 
    guides(color=guide_legend(title=""))
  if(!is.null(lab)){
    if(lab == "Label") p <- p + geom_text_repel(show.legend = NULL)
  }
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

data_3degcount1 <- function(gene_type,data,result_Condm, result_FDR, specific, fc, fdr, basemean,result_list=NULL){
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
      Pattern1 <- "pattern4"
      Pattern2 <- "pattern5"
    }
    if(specific == 2) {
      specific = collist[2]
      FC_xlab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[1], ")"))
      FC_ylab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[3], ")"))
      result_Condm$FC_x <- log2((result_Condm$C2 + 0.01)/(result_Condm$C1 + 0.01))
      result_Condm$FC_y <- log2((result_Condm$C2 + 0.01)/(result_Condm$C3 + 0.01))
      Pattern1 <- "pattern3"
      Pattern2 <- "pattern5"
    }
    if(specific == 3) {
      specific = collist[3]
      FC_xlab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[1], ")"))
      FC_ylab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[2], ")"))
      result_Condm$FC_x <- log2((result_Condm$C3 + 0.01)/(result_Condm$C1 + 0.01))
      result_Condm$FC_y <- log2((result_Condm$C3 + 0.01)/(result_Condm$C2 + 0.01))
      Pattern1 <- "pattern2"
      Pattern2 <- "pattern5"
    }
    print("FDR")
    print(head(result_FDR))
    result_FDR$FDR <- 1 - result_FDR$PPDE
    result <- merge(result_Condm, result_FDR, by=0)
    data$Row.names <- rownames(data)
    data2 <- merge(data, result, by="Row.names")
      data3 <- data2[,- which(colnames(data2) == "pattern1")]
      data3 <- data3[,- which(colnames(data3) == "pattern2")]
      data3 <- data3[,- which(colnames(data3) == "pattern3")]
      data3 <- data3[,- which(colnames(data3) == "pattern4")]
      data3 <- data3[,- which(colnames(data3) == "pattern5")]
      data3 <- data3[,- which(colnames(data3) == "C1")]
      data3 <- data3[,- which(colnames(data3) == "C2")]
      data3 <- data3[,- which(colnames(data3) == "C3")]
      data3 <- data3[,- which(colnames(data3) == "PPDE")]
      Pattern <- rep("cannot_be_classified", nrow(data3))
      Pattern[which(data3$MAP == "pattern1")] = paste(collist[1], "=", collist[2], "=", collist[3])
      Pattern[which(data3$MAP == "pattern2" & data3$FC_y > 0)] = paste(collist[1], "=", collist[2], ">", collist[3])
      Pattern[which(data3$MAP == "pattern2" & data3$FC_y < 0)] = paste(collist[1], "=", collist[2], "<", collist[3])
      Pattern[which(data3$MAP == "pattern3" & data3$FC_x > 0)] = paste(collist[1], "=", collist[3], ">", collist[2])
      Pattern[which(data3$MAP == "pattern3" & data3$FC_x < 0)] = paste(collist[1], "=", collist[3], "<", collist[2])
      Pattern[which(data3$MAP == "pattern4" & data3$FC_x > 0 & (data3$FC_y > 0))] = paste(collist[1], ">", collist[2], "=", collist[3])
      Pattern[which(data3$MAP == "pattern4" & data3$FC_x < 0 & (data3$FC_y < 0))] = paste(collist[1], "<", collist[2], "=", collist[3])
      Pattern[which(data3$MAP == "pattern5" & data3$FC_x > 0 & (data3$FC_y > 0) & ((data3$FC_x - data3$FC_y) > 0))] = paste(collist[1], ">", collist[2], ">", collist[3])
      Pattern[which(data3$MAP == "pattern5" & data3$FC_x > 0 & (data3$FC_y > 0) & ((data3$FC_x - data3$FC_y) < 0))] = paste(collist[1], ">", collist[3], ">", collist[2])
      Pattern[which(data3$MAP == "pattern5" & data3$FC_x < 0 & (data3$FC_y < 0) & ((data3$FC_x - data3$FC_y) > 0))] = paste(collist[1], "<", collist[3], "<", collist[2])
      Pattern[which(data3$MAP == "pattern5" & data3$FC_x < 0 & (data3$FC_y < 0) & ((data3$FC_x - data3$FC_y) < 0))] = paste(collist[1], "<", collist[2], "<", collist[3])
      Pattern[which(data3$MAP == "pattern5" & data3$FC_x > 0 & (data3$FC_y < 0))] = paste(collist[2], "<", collist[1], "<", collist[3])
      Pattern[which(data3$MAP == "pattern5" & data3$FC_x < 0 & (data3$FC_y > 0))] = paste(collist[3], "<", collist[1], "<", collist[2])
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

data_3degcount2 <- function(gene_type,data3, Species, Ortholog,org){
  if(is.null(data3)){
    return(NULL)
  }else{
    if(length(unique(data3$sig)) == 1){
      data4 <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
      return(data4)
    } else {
      data4 <- dplyr::filter(data3, sig != "NS")
      if(gene_type != "SYMBOL"){
        if(Species != "not selected"){
          if(sum(is.element(no_orgDb, Species)) == 1){
            gene_IDs <- Ortholog
          }else{
          my.symbols <- data4$Row.names
          if(str_detect(my.symbols[1], "^AT.G")) key = "TAIR" else key = "ENSEMBL"
          gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                          keytype = key,
                                          columns = c(key,"SYMBOL", "ENTREZID"))
          }
          colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data4 <- merge(data4, gene_IDs, by="Row.names")
        }
      }else{
        if(Species != "not selected"){
          if(sum(is.element(no_orgDb, Species)) == 1){
            gene_IDs <- Ortholog
            gene_IDs <- gene_IDs[,-1]
          }else{
          my.symbols <- data4$Row.names
          gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                          keytype = "SYMBOL",
                                          columns = c("SYMBOL", "ENTREZID"))
          }
          colnames(gene_IDs) <- c("Row.names", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data4 <- merge(data4, gene_IDs, by="Row.names")
        }
      }
      return(data4)
    }
  }
}
cond3_scatter_plot <- function(gene_type,data, data4, result_Condm, result_FDR, specific, 
                               fc, fdr, basemean, y_axis=NULL, x_axis=NULL,
                               GOI=NULL, heatmap = TRUE, Species,brush=NULL){
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
      Pattern1 <- "pattern4"
      Pattern2 <- "pattern5"
    }
    if(specific == 2) {
      specific = collist[2]
      FC_xlab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[1], ")"))
      FC_ylab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[3], ")"))
      result_Condm$FC_x <- log2((result_Condm$C2 + 0.01)/(result_Condm$C1 + 0.01))
      result_Condm$FC_y <- log2((result_Condm$C2 + 0.01)/(result_Condm$C3 + 0.01))
      Pattern1 <- "pattern3"
      Pattern2 <- "pattern5"
    }
    if(specific == 3) {
      specific = collist[3]
      FC_xlab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[1], ")"))
      FC_ylab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[2], ")"))
      result_Condm$FC_x <- log2((result_Condm$C3 + 0.01)/(result_Condm$C1 + 0.01))
      result_Condm$FC_y <- log2((result_Condm$C3 + 0.01)/(result_Condm$C2 + 0.01))
      Pattern1 <- "pattern2"
      Pattern2 <- "pattern5"
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
    complete_data <- stats::na.omit(data3)
    labs_data <- subset(complete_data, padj <= fdr & Row.names !=
                          "" & (FC_x)*(FC_y) >= log2(fc))
    labs_data<-  labs_data[sort(labs_data$FC_xy, decreasing = T, index=T)$ix,]
    labs_data <- dplyr::filter(labs_data, sig != "NS")
    labs_data2 <- utils::head(labs_data, 20)
    font.label <- data.frame(size=7.5, color="black", face = "bold.italic")
    set.seed(42)
    FC_x <- FC_y <- Row.names <- padj <- NULL
    p <- ggplot(data3, aes(x = FC_x, y = FC_y)) + geom_point(aes(color = sig),size = 0.1,alpha=.5)
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
        if(gene_type != "SYMBOL"){
          if(Species != "not selected"){
            p <- p + ggrepel::geom_text_repel(data = labs_data2, mapping = aes(label = Unique_ID),
                                              box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), 
                                              force = 1, fontface = font.label$face,show.legend = NULL,
                                              size = font.label$size/2, color = font.label$color)
          }else{
            p <- p + ggrepel::geom_text_repel(data = labs_data2, mapping = aes(label = Row.names),
                                              box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), 
                                              force = 1, fontface = font.label$face,show.legend = NULL,
                                              size = font.label$size/2, color = font.label$color)
          }
        }else{
        p <- p + ggrepel::geom_text_repel(data = labs_data2, mapping = aes(label = Row.names),
                                          box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), 
                                          force = 1, fontface = font.label$face,show.legend = NULL,
                                          size = font.label$size/2, color = font.label$color)
        }
      } 
    }
    if(!is.null(y_axis) && !is.null(x_axis)){
      p <- p +  xlim(x_axis) + ylim(y_axis)
    }
    if(!is.null(dim(brush)[1])){
    if(!is.null(GOI) || dim(brush)[1]!=0) {
      if(dim(brush)[1]!=0){
        if(gene_type != "SYMBOL"){
          if(Species != "not selected"){
            GOI <- brush$Unique_ID
          }else{
            GOI <- brush$Row.names
          }
        }else{
          GOI <- brush$Row.names
        }
      }
      for(name in GOI){
        if(gene_type != "SYMBOL"){
          if(Species != "not selected"){
            data3$color[data3$Unique_ID == name] <- "GOI"
          }else{
            data3$color[data3$Row.names == name] <- "GOI"
          }
        }else{
          data3$color[data3$Row.names == name] <- "GOI"
        }
      }
      if(gene_type != "SYMBOL"){
        if(Species != "not selected"){
          p <- p + geom_point(data=dplyr::filter(data3, color == "GOI"),color="green", size=1)
          p <- p + ggrepel::geom_label_repel(data = dplyr::filter(data3, color == "GOI"), mapping = aes(label = Unique_ID),alpha = 0.6,label.size = NA, 
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1, fontface = "bold.italic")
        }else{
          p <- p + geom_point(data=dplyr::filter(data3, color == "GOI"),color="green", size=1)
          p <- p + ggrepel::geom_label_repel(data = dplyr::filter(data3, color == "GOI"), mapping = aes(label = Row.names),alpha = 0.6,label.size = NA, 
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1, fontface = "bold.italic")
        }
      }else{
        p <- p + geom_point(data=dplyr::filter(data3, color == "GOI"),color="green", size=1)
        p <- p + ggrepel::geom_label_repel(data = dplyr::filter(data3, color == "GOI"), mapping = aes(label = Row.names),alpha = 0.6,label.size = NA, 
                                          box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1, fontface = "bold.italic")
      }
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
      data.z <- genefilter::genescale(data5, axis=1, method="Z")
      ht <- as.grob(Heatmap(data.z, name = "z-score", column_order = colnames(data.z),
                            clustering_method_columns = 'ward.D2',use_raster = TRUE,
                            show_row_names = F, show_row_dend = T,column_names_side = "top"))
    }
    
    p <- plot_grid(p, ht, rel_widths = c(2, 1))
    }
    return(p)
  }
}
cond3_scatter_range <- function(gene_type,data, data4, result_Condm, result_FDR, specific, 
                                fc, fdr, basemean){
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
      Pattern1 <- "pattern4"
      Pattern2 <- "pattern5"
    }
    if(specific == 2) {
      specific = collist[2]
      FC_xlab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[1], ")"))
      FC_ylab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[3], ")"))
      result_Condm$FC_x <- log2((result_Condm$C2 + 0.01)/(result_Condm$C1 + 0.01))
      result_Condm$FC_y <- log2((result_Condm$C2 + 0.01)/(result_Condm$C3 + 0.01))
      Pattern1 <- "pattern3"
      Pattern2 <- "pattern5"
    }
    if(specific == 3) {
      specific = collist[3]
      FC_xlab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[1], ")"))
      FC_ylab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[2], ")"))
      result_Condm$FC_x <- log2((result_Condm$C3 + 0.01)/(result_Condm$C1 + 0.01))
      result_Condm$FC_y <- log2((result_Condm$C3 + 0.01)/(result_Condm$C2 + 0.01))
      Pattern1 <- "pattern2"
      Pattern2 <- "pattern5"
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
    return(data3)
  }
}
v1_GetMultiFC<-function(EBMultiOut, collist,count,SmallNum = 0.01,EBSeq_mode=F){
  if(EBSeq_mode==F){
    norm_count <- as.data.frame(EBMultiOut$DataNorm)
    colnames(norm_count) <- colnames(count)
    for(name in unique(collist)){
      cond <- norm_count %>% dplyr::select(starts_with(name))
      cond_ave <- apply(cond,1,mean)
      if(name == unique(collist)[1]) df <- data.frame(cond_ave) else df <- data.frame(df, cond_ave)
    }
  colnames(df) <- c("C1","C2","C3")
  return(df)
  }
}
enrichment3_1 <- function(data3, data4, cnet_list2){
  if(!is.null(cnet_list2)){
        withProgress(message = "dotplot",{
          data <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          colnames(data) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
          for (name in names(cnet_list2)) {
              em <- cnet_list2[[name]]
              sum <- length(data4$ENTREZID[data4$sig == name])
              if (length(as.data.frame(em)$ID) == 0) {
                cnet1 <- NULL
              } else {
                cnet1 <- as.data.frame(em)
                cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
                cnet1 <- cnet1[sort(cnet1$qvalue, decreasing = F, index=T)$ix,]
                data <- rbind(data, cnet1)
              }
          }
          if(!is.null(data)) data$GeneRatio <- DOSE::parse_ratio(data$GeneRatio)
          return(data)
          incProgress()
        })
  }
}

keggEnrichment1 <- function(data3, data4, Species, Gene_set, org, H_t2g){
  if(!is.null(Gene_set) && Species != "not selected"){
    if(is.null(data4)){
      return(NULL)
    }else{
      cnet_list2 <- list()
      for (name in unique(data3$sig)) {
        if (name != "NS"){
          if(is.null(H_t2g)){
            cnet_list2 <- NULL
          }else{
            H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
            em <- clusterProfiler::enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) == 0) {
              cnet1 <- NULL
            } else {
              cnet1 <- clusterProfiler::setReadable(em, org, 'ENTREZID')
              cnet_list2[[name]] = cnet1
            }
          }
        }
      }
      return(cnet_list2)
    }
  }else return(NULL)
}
keggEnrichment1_xenopus <- function(data3, data4, Species, Gene_set, org, org_code){
  if(!is.null(Gene_set) && Species != "not selected"){
    if(is.null(data4)){
      return(NULL)
    }else{
      cnet_list2 <- list()
      for (name in unique(data3$sig)) {
        if (name != "NS"){
            if(Gene_set == "KEGG"){
              em <- enrichKEGG(data4$ENTREZID[data4$sig == name], organism = org_code, pvalueCutoff = 0.05,keyType = "ncbi-geneid")
            }
            if(Gene_set == "GO biological process"){
              em <- enrichGO(data4$ENTREZID[data4$sig == name], OrgDb = org, ont = "BP",pvalueCutoff = 0.05)
            }
            if(Gene_set == "GO cellular component"){
              em <- enrichGO(data4$ENTREZID[data4$sig == name], OrgDb= org, ont = "CC",pvalueCutoff = 0.05) 
            }
            if(Gene_set == "GO molecular function"){
              em <- enrichGO(data4$ENTREZID[data4$sig == name], OrgDb = org, ont = "MF",pvalueCutoff = 0.05) 
            }
            if (length(as.data.frame(em)$ID) == 0) {
              cnet1 <- NULL
            } else {
              cnet1 <- clusterProfiler::setReadable(em, org, 'ENTREZID')
              cnet_list2[[name]] = cnet1
            }
        }
      }
      return(cnet_list2)
    }
  }else return(NULL)
}

keggEnrichment2 <- function(data3, data4,cnet_list2){
  if(!is.null(cnet_list2)){
    data <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
    colnames(data) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
    cnet_list <- list()
    for (name in names(cnet_list2)) {
      sum <- length(data4$ENTREZID[data4$sig == name])
      em <- cnet_list2[[name]]
      if (length(as.data.frame(em)$ID) == 0) {
        cnet1 <- NULL
      } else {
        cnet1 <- em
        if ((length(as.data.frame(cnet1)$ID) == 0) || 
            length(which(!is.na(unique(as.data.frame(cnet1)$qvalue))))==0) {
          c <- NULL
        } else{
          c <- clusterProfiler::cnetplot(cnet1, cex_label_gene = 0.7, cex_label_category = 0.75,
                        cex_category = 0.75, colorEdge = TRUE)
          c <- try(as.grob(c + guides(edge_color = "none")))
          if(length(class(c)) == 1){
            if(class(c) == "try-error") c <- NULL
          }
          cnet_list[[name]] = c
        }
        cnet1 <- as.data.frame(em)
        cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
        cnet1 <- cnet1[sort(cnet1$qvalue, decreasing = F, index=T)$ix,]
        cnet1 <- cnet1[1:5,]
        data <- rbind(data, cnet1)
      }
    }
    data <- dplyr::filter(data, !is.na(Group))
    data <- dplyr::filter(data, !is.na(Description))
    data$GeneRatio <- DOSE::parse_ratio(data$GeneRatio)
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
                     scale_size(range=c(1, 6))+ DOSE::theme_dose(font.size=12)+ylab(NULL)+xlab(NULL) +
                     scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
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
  }else return(NULL)
}

enrich_for_table <- function(data, H_t2g, Gene_set){
  if(length(as.data.frame(data)$Description) == 0 || is.null(H_t2g)){
    return(NULL)
  }else{
    colnames(data)[1] <- "gs_name"
    H_t2g <- H_t2g %>% distinct(gs_name, .keep_all = T)
    data2 <- left_join(data, H_t2g, by="gs_name")  %>% as.data.frame()
    data2$Group <- gsub("\n"," ", data2$Group)
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
GeneList_for_enrichment <- function(Species, Ortholog,Gene_set, org, Custom_gene_list,Biomart_archive,gene_type=NULL){
  if(Species != "not selected" || is.null(Gene_set) || is.null(org)){
    if(Species %in% orgDb_list == TRUE) species <- Species else species <- Ortholog
    H_t2g <- NULL
      if(Gene_set == "MSigDB Hallmark"){
        H_t2g <- msigdbr::msigdbr(species = species, category = "H") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description) 
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="HALLMARK_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="P53", replacement = "p53")
      }
      if(Gene_set == "KEGG"){
        H_t2g <- msigdbr::msigdbr(species = species, category = "C2", subcategory = "CP:KEGG") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="KEGG_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
    if(Gene_set == "Position"){
      H_t2g <- msigdbr::msigdbr(species = species, category = "C1") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description) 
    }
      if(Gene_set == "Transcription factor targets"){
        H_t2g <- msigdbr::msigdbr(species = species, category = "C3")
        H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "TFT:GTRD" | gs_subcat == "TFT:TFT_Legacy") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      }
      if(Gene_set == "Reactome"){
        H_t2g <- msigdbr::msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="REACTOME_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
      if(Gene_set == "miRNA target"){
        H_t2g <- msigdbr::msigdbr(species = species, category = "C3")
        H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "MIR:MIRDB" | gs_subcat == "MIR:MIR_Legacy") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      }
      if(Gene_set == "GO biological process"){
        H_t2g <- msigdbr::msigdbr(species = species, category = "C5", subcategory = "GO:BP") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOBP_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
      if(Gene_set == "GO cellular component"){
        H_t2g <- msigdbr::msigdbr(species = species, category = "C5", subcategory = "GO:CC") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOCC_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
      if(Gene_set == "GO molecular function"){
        H_t2g <- msigdbr::msigdbr(species = species, category = "C5", subcategory = "GO:MF") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOMF_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
      if(Gene_set == "Human phenotype ontology"){
        H_t2g <- msigdbr::msigdbr(species = species, category = "C5", subcategory = "HPO") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="HP_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
      if(Gene_set == "WikiPathways"){
        H_t2g <- msigdbr::msigdbr(species = species, category = "C2", subcategory = "CP:WIKIPATHWAYS") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="WP_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
      if(Gene_set == "BioCarta"){
        H_t2g <- msigdbr::msigdbr(species = species, category = "C2", subcategory = "CP:BIOCARTA") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="BIOCARTA_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
      if(Gene_set == "Custom gene set"){
        if(!is.null(Custom_gene_list)){
          H_t2g <- gene_list_convert_for_enrichment(data= Custom_gene_list, Species = Species)
          H_t2g <- data.frame(gs_name = H_t2g$Group, entrez_gene = H_t2g$ENTREZID)
          H_t2g$gs_name <- gsub(":", "_", H_t2g$gs_name)
          H_t2g <- H_t2g %>% dplyr::filter(!is.na(entrez_gene))
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
      if(Gene_set == "DoRothEA regulon (activator)"){
        H_t2g <- as_tibble(dorothea(species = species, type = "DoRothEA regulon (activator)")) %>%
          dplyr::select(gs_name, entrez_gene, confidence)
        H_t2g$entrez_gene <- as.integer(H_t2g$entrez_gene)
      }
      if(Gene_set == "DoRothEA regulon (repressor)"){
        H_t2g <- as_tibble(dorothea(species = species, type = "DoRothEA regulon (repressor)")) %>%
          dplyr::select(gs_name, entrez_gene, confidence)
        H_t2g$entrez_gene <- as.integer(H_t2g$entrez_gene)
      }
      if(Gene_set == "PID (Pathway Interaction Database)"){
        H_t2g <- msigdbr::msigdbr(species = species, category = "C2", subcategory = "CP:PID") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="PID_", replacement = "")
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="PATHWAY", replacement = "pathway")
      }
      return(H_t2g)
  }else return(NULL)
}
GOI_color_palette<-c("default","Set1","Set2","Set3","Paired","Dark2","Accent","Spectral")
GOIboxplot <- function(data,statistical_test=NULL,plottype="Boxplot",ymin=0,
                       pair=NULL,ssGSEA=FALSE,color_design="new",color="default",rev="OFF"){
  print("GOIboxplot start")
  collist <- gsub("\\_.+$", "", colnames(data))
  collist <- unique(collist)
  rowlist <- rownames(data)
  data$Row.names <- rownames(data)
  data <- data %>% gather(key=sample, value=value,-Row.names)
  if(!is.null(pair)) data <- dplyr::left_join(data,pair)
  data$sample <- gsub("\\_.+$", "", data$sample)
  data$Row.names <- as.factor(data$Row.names)
  data$sample <- factor(data$sample,levels=collist,ordered=TRUE)
  data$value <- as.numeric(data$value)
  stat.test <-NULL
  if(!is.null(statistical_test) && statistical_test != "not_selected"){
    res <- data.frame(matrix(rep(NA, 11), nrow=1))[numeric(0), ]
    colnames(res) <- c("Gene", "group1", "group2", "term", "null.value","Std.Error","coefficients","t.value","p.adj","xmin", "xmax")
    if (statistical_test == "TukeyHSD"){
      stat.test <- data %>% 
        group_by(Row.names) %>% 
        tukey_hsd(value ~ sample) %>% 
        add_significance("p.adj") %>% 
        add_xy_position(scales = "free", step.increase = 0.2)
    }
    if(statistical_test == "Welch's t-test"){
      stat.test <- data %>% 
        group_by(Row.names) %>% 
        t_test(value ~ sample) %>% 
        add_significance() %>% 
        add_xy_position(scales = "free", step.increase = 0.2)
    }
    if(statistical_test == "Wilcoxon test"){
      stat.test <- data %>% 
        group_by(Row.names) %>% 
        wilcox_test(value ~ sample) %>% 
        add_significance() %>% 
        add_xy_position(scales = "free", step.increase = 0.2)
    }
    if(statistical_test == "Dunn's test"){
        stat.test <- data %>% 
          group_by(Row.names) %>% 
          dunn_test(value ~ sample) %>% 
          add_significance() %>% 
          add_xy_position(scales = "free", step.increase = 0.2)
    }
    if(statistical_test == "Dunnet's test"){
      for (name2 in rowlist){
        data2 <- dplyr::filter(data, Row.names == name2)
        dun <- stats::aov(value~sample, data2)
        dunnette <- multcomp::glht(model = dun, linfct=multcomp::mcp(sample="Dunnett"))
        dunnette2 <- try(summary(dunnette))
        p.adj <- c()
        coefficients <- c()
        Std.Error <- c()
        t.value <- c()
        group1 <- c()
        group2 <- c()
        term <- c()
        null.value <- c()
        xmin <- c()
        xmax <- c()
        if(length(class(dunnette2)) == 1){
        if(class(dunnette2) == "try-error"){
          for (i in 1:(length(collist)-1)){
            p.adj <- c(p.adj, NA)
            coefficients <- c(coefficients, NA)
            Std.Error <- c(Std.Error, NA)
            t.value <- c(t.value, NA)
            group1 <- c(group1, c(collist[1]))
            group2 <- c(group2, c(collist[i+1]))
            term <- c(term, c("sample"))
            null.value <- c(null.value, NA)
            xmin <- c(xmin, c(1))
            xmax <- c(xmax, c(i+1))
          }
        }}else{
        for (i in 1:(length(collist)-1)){
          p.adj <- c(p.adj, dunnette2[["test"]][["pvalues"]][i])
          coefficients <- c(coefficients, dunnette2[["test"]][["coefficients"]][i])
          Std.Error <- c(Std.Error, dunnette2[["test"]][["sigma"]][i])
          t.value <- c(t.value, dunnette2[["test"]][["tstat"]][i])
          group1 <- c(group1, c(collist[1]))
          group2 <- c(group2, c(collist[i+1]))
          term <- c(term, c("sample"))
          null.value <- c(null.value, 0)
          xmin <- c(xmin, c(1))
          xmax <- c(xmax, c(i+1))
        }
        }
        res2 <- data.frame(Row.names = name2, group1 = group1, group2 = group2, term = term,
                           null.value = null.value, Std.Error = Std.Error, coefficients = coefficients,
                           t.value = t.value, p.adj = p.adj, xmin = xmin, xmax = xmax)
        res <- rbind(res, res2)
      }
      res <- res %>% arrange(Row.names) %>% group_by(Row.names)
      stat.test2 <- data %>% 
        group_by(Row.names) %>% 
        get_y_position(value ~ sample, scales = "free", step.increase = 0.15, fun = "mean_se") %>% 
        dplyr::filter(group1 == collist[1])
      stat.test3 <- cbind(stat.test2,res[,-1:-3])
      stat.test3$Row.names <- as.factor(stat.test3$Row.names)
      stat.test <- stat.test3 %>% add_significance("p.adj")
    }
  }
  stat.test$group1 <- gsub("\n"," ",stat.test$group1)
  stat.test$group2 <- gsub("\n"," ",stat.test$group2)
  if(length(rowlist) > 200){
    p <- NULL
  }else{
    if(ssGSEA == FALSE) {
      if(is.na(ymin)) ylim <- NULL else ylim = c(ymin, NA)
      ylab = "Normalized_count"
    }else {
      ylim = NULL
      ylab = "ssGSEA score"
    }
    if (plottype == "Boxplot"){
      if(color_design=="new"){
      p <- ggplot(data, aes(x=sample,y=value))+
        geom_boxplot(aes(group=sample,colour=sample,fill=after_scale(alpha(colour,0.5))))+
        geom_jitter(alpha=1)+
        xlab(NULL)+ylab(ylab)+theme_classic()
      }else{
      p<- ggpubr::ggboxplot(data, x = "sample", y = "value",
                            fill = "sample", scales = "free",
                            add = "jitter",add.params = list(alpha=1),
                            xlab = FALSE, ylab = ylab)
      }
    }
    if (plottype == "Barplot"){
      if(color_design=="new"){
        data2 <- 
          data %>% 
          group_by(sample, Row.names) %>% 
          summarise(mean = mean(value), sd = sd(value), n = n(), se = sd / sqrt(n))
        p <- ggplot(data, aes(x=sample,y=value,colour=sample,fill=after_scale(alpha(colour,0.5))))+
          geom_bar(stat = "identity",data = data2,aes(x=sample,y=mean))+
          geom_errorbar(data = data2,aes(x=sample,y=mean,ymin = mean - se, ymax = mean + se,colour=sample), width = 0.3)+
          geom_jitter(color="black",alpha=1)+
          xlab(NULL)+ylab(ylab)+theme_classic()
      }else{
      p <- ggbarplot(data,x = "sample", y = "value", scales = "free",
                     facet.by = "Row.names", fill = "sample",add = c("mean_se", "jitter"),
                     add.params = list(size=0.5), xlab = FALSE)
      }
    }
    if (plottype == "Errorplot"){
      p <- ggerrorplot(data,x = "sample", y = "value",
                       scales = "free", add = "jitter", facet.by = "Row.names",
                       add.params = list(size=0.5), xlab = FALSE, error.plot = "errorbar") + 
        stat_summary(geom = "point", shape = 95,size = 5,col = "black", fun = "mean")
    }
    if (plottype == "Violin plot"){
      if(color_design=="new"){
        p <- ggplot(data, aes(x=sample,y=value))+
          geom_violin(aes(group=sample,colour=sample,fill=after_scale(alpha(colour,0.5))),trim=FALSE)+
          geom_boxplot(aes(group=sample), colour="black",fill="white",width = .2)+

          xlab(NULL)+ylab(ylab)+theme_classic()
      }else{
      p <- ggviolin(data,x = "sample", y = "value",
                    facet.by = "Row.names", fill = "sample",add = c("jitter","boxplot"),
                    add.params = list(size=0.5,fill = "white",alpha=1), xlab = FALSE,alpha = 0.5)
      }
    }
    if(!is.null(pair)){
      if(plottype == "Boxplot"){
      p <- ggplot(data, aes(x = sample, y = value)) + geom_boxplot(aes(fill=sample))+
        geom_line(aes(group = pair),alpha = .2) +
        geom_point() + theme_classic() + theme(legend.position = "top")+ 
        xlab(NULL)  + ylab(ylab)
      }
      if(plottype == "without boxplot"){
        p <- ggplot(data, aes(x = sample, y = value,group = pair)) + 
          geom_line() +
          geom_point(aes(color = sample)) + theme_classic() + theme(legend.position = "top")+ 
          xlab(NULL)  + ylab(ylab)+ 
          scale_color_manual(values=c("#00BFC4", "#F8766D"))
      }
    }
    if(color!="default"){
      if(is.null(rev)) validate("")
      if(rev == "ON") direction=-1 else direction=1
      if(color_design == "new") p <- p + scale_color_brewer(palette=color,direction = direction) else p <- p + scale_fill_brewer(palette=color,direction = direction)
    }else{
      if(length(unique(data$sample)) < 4){
        if(length(unique(data$sample)) ==2) color <- c("grey40", "brown2")
        if(length(unique(data$sample)) ==3) color <- c("grey40", "dodgerblue3", "brown2")
       if(color_design == "new")  p <- p + scale_color_manual(values=color) else  p <- p + scale_fill_manual(values=color)
      }
      
    }
  p <- (facet(p, facet.by = "Row.names",
              panel.labs.background = list(fill = "transparent", color = "transparent"),
              scales = "free", short.panel.labs = T, panel.labs.font = list(size=15, face = "italic"))+ 
          theme(axis.text.x = element_blank(),legend.position = "top",
                panel.background = element_rect(fill = "transparent", size = 0.5),
                title = element_text(size = 10),text = element_text(size = 12),
                axis.title.y = element_text(size=15),axis.text.y = element_text(size = 15),
                legend.text = element_text(size=15),legend.title = element_blank()))
  }
  if(!is.null(ylim)) p <- p + ylim(ylim)
  print("GOIboxplot end")
  if(!is.null(statistical_test) && statistical_test != "not_selected"){
    if(length(rowlist) <= 200) p <- p + stat_pvalue_manual(stat.test,hide.ns = T, size = 5)
    df <- list()
    df[["plot"]] <- p
    df[["statistical_test"]] <- stat.test
    return(df)
  }else return(p)
}
enrich_viewer_forMulti1 <- function(gene_type,df, Species, Ortholog, org){
  if(is.null(df) || Species == "not selected"){
    return(NULL)
  }else{
    my.symbols <- df$GeneID
    if(gene_type != "SYMBOL"){
      if(sum(is.element(no_orgDb, Species)) == 1){
        print(Ortholog)
        gene_IDs <- Ortholog
      }else{
      if(str_detect(my.symbols[1], "^AT.G")) key = "TAIR" else key = "ENSEMBL"
      gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                      keytype = key,
                                      columns = c(key,"SYMBOL", "ENTREZID"))
      }
      colnames(gene_IDs) <- c("GeneID","SYMBOL", "ENTREZID")
    }else{
      if(sum(is.element(no_orgDb, Species)) == 1){
        gene_IDs <- Ortholog
        gene_IDs <- gene_IDs[,-1]
      }else{
      gene_IDs <- AnnotationDbi::select(org, keys = my.symbols,
                                        keytype = "SYMBOL",
                                        columns = c("ENTREZID", "SYMBOL"))
      }
      colnames(gene_IDs) <- c("GeneID","ENTREZID")
    }
    gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
    data <- merge(df, gene_IDs, by="GeneID")
    return(data)
  }
}
enrich_viewer_forMulti2 <- function(gene_type,df, Species,Ortholog, Gene_set, org, org_code, H_t2g){
  if(!is.null(Gene_set)){
    data3 <- enrich_viewer_forMulti1(gene_type=gene_type,df = df, Species = Species, org = org,Ortholog=Ortholog)
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
            sum <- length(data3$ENTREZID[data3$Group == name])
            em <- clusterProfiler::enricher(data3$ENTREZID[data3$Group == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) != 0) {
              if(length(colnames(as.data.frame(em))) == 9){
                cnet1 <- as.data.frame(clusterProfiler::setReadable(em, org, 'ENTREZID'))
                cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
                df <- rbind(df, cnet1)
              }
            }
          }
          }
          if(length(df$ID) !=0){
            df$GeneRatio <- DOSE::parse_ratio(df$GeneRatio)
            return(df)
          }else return(NULL)
        })
      }
  } 
}
enrich_viewer_forMulti2_xenopus <- function(df, Species,Ortholog, Gene_set, org, org_code){
  if(!is.null(Gene_set)){
    data3 <- enrich_viewer_forMulti1(df = df, Species = Species,Ortholog = Ortholog, org = org)
    if(is.null(data3)){
      return(NULL)
    }else{
      withProgress(message = "enrichment analysis",{
          df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
          for (name in unique(data3$Group)) {
            sum <- length(data3$ENTREZID[data3$Group == name])
            if(Gene_set == "KEGG"){
              em <- enrichKEGG(data3$ENTREZID[data3$Group == name], organism = org_code, pvalueCutoff = 0.05,keyType = "ncbi-geneid")
            }
            if(Gene_set == "GO biological process"){
              em <- enrichGO(data3$ENTREZID[data3$Group == name], OrgDb = org, ont = "BP",pvalueCutoff = 0.05)
            }
            if(Gene_set == "GO cellular component"){
              em <- enrichGO(data3$ENTREZID[data3$Group == name], OrgDb= org, ont = "CC",pvalueCutoff = 0.05) 
            }
            if(Gene_set == "GO molecular function"){
              em <- enrichGO(data3$ENTREZID[data3$Group == name], OrgDb = org, ont = "MF",pvalueCutoff = 0.05) 
            }
            if (length(as.data.frame(em)$ID) != 0) {
              if(length(colnames(as.data.frame(em))) == 9){
                cnet1 <- as.data.frame(clusterProfiler::setReadable(em, org, 'ENTREZID'))
                cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
                df <- rbind(df, cnet1)
              }
            }
          }
        if(length(df$ID) !=0){
          df$GeneRatio <- DOSE::parse_ratio(df$GeneRatio)
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

enrich_gene_list <- function(data, Gene_set, H_t2g, org,org_code=NULL){
  if(!is.null(Gene_set)){
    if(is.null(data)){
      return(NULL)
    }else{
      if(is.null(H_t2g)){
        df <- NULL
      }else{
        H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
        df <- list()
        for (name in unique(data$Group)) {
          sum <- length(data$ENTREZID[data$Group == name])
          em <- clusterProfiler::enricher(data$ENTREZID[data$Group == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
          if (length(as.data.frame(em)$ID) != 0) {
            if(length(colnames(as.data.frame(em))) == 9){
              cnet1 <- clusterProfiler::setReadable(em, org, 'ENTREZID')
              df[[name]] <- cnet1
            }
          }
        }
      }
        return(df)
      }
    }
}
enrich_gene_list_xenopus <- function(data, Gene_set, org,org_code=NULL){
    if(is.null(data)){
      return(NULL)
    }else{
        df <- list()
        for (name in unique(data$Group)) {
          sum <- length(data$ENTREZID[data$Group == name])
          if(Gene_set == "KEGG"){
            em <- enrichKEGG(data$ENTREZID[data$Group == name], organism = org_code, pvalueCutoff = 0.05,keyType = "ncbi-geneid")
          }
          if(Gene_set == "GO biological process"){
            em <- enrichGO(data$ENTREZID[data$Group == name], OrgDb = org, ont = "BP",pvalueCutoff = 0.05)
          }
          if(Gene_set == "GO cellular component"){
            em <- enrichGO(data$ENTREZID[data$Group == name], OrgDb= org, ont = "CC",pvalueCutoff = 0.05) 
          }
          if(Gene_set == "GO molecular function"){
            em <- enrichGO(data$ENTREZID[data$Group == name], OrgDb = org, ont = "MF",pvalueCutoff = 0.05) 
          }
          if (length(as.data.frame(em)$ID) != 0) {
            if(length(colnames(as.data.frame(em))) == 9){
              cnet1 <- clusterProfiler::setReadable(em, org, 'ENTREZID')
              df[[name]] <- cnet1
            }
          }
        }
      return(df)
    }
}


enrich_genelist <- function(data, enrich_gene_list, showCategory=5,section=NULL,group_order=NULL){
      if(is.null(data) || is.null(enrich_gene_list)){
        return(NULL)
      }else{
          df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
          cluster_list <- c()
          for (name in names(enrich_gene_list)) {
            sum <- length(data$ENTREZID[data$Group == name])
            em <- enrich_gene_list[[name]]
            if (length(as.data.frame(em)$ID) != 0) {
              if(length(colnames(as.data.frame(em))) == 9){
                cnet1 <- as.data.frame(em)
                cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
                cluster_list <- c(cluster_list, paste(name, "\n","(",sum, ")",sep = ""))
                if(!is.null(group_order)) group_order[which(group_order == name)] <- paste(name, "\n","(",sum, ")",sep = "")
                cnet1 <- cnet1[sort(cnet1$pvalue, decreasing = F, index=T)$ix,]
                if (length(cnet1$pvalue) > showCategory){
                  cnet1 <- cnet1[1:showCategory,]
                }
                df <- rbind(df, cnet1)
              }
            }
          }
          if(!is.null(group_order)) group_order <- group_order[group_order %in% cluster_list]
          if ((length(df$Description) == 0) || length(which(!is.na(unique(df$qvalue)))) == 0) {
            p1 <- NULL
          } else{
            if(!is.null(section)){
            if(section == "enrichmentviewer"){
              df$Group <- gsub("_", " ", df$Group)
              if(!is.null(group_order)) group_order <- gsub("_", " ", group_order)
              for(i in 1:length(df$Group)){
                df$Group[i] <- paste(strwrap(df$Group[i], width = 15),collapse = "\n")
              }
              for(i in 1:length(unique(df$Group))){
                if(!is.null(group_order)) group_order[i] <- paste(strwrap(group_order[i], width = 15),collapse = "\n")
              }
              df$Group <- gsub(" \\(", "\n\\(", df$Group)
              if(!is.null(group_order)) group_order <- gsub(" \\(", "\n\\(", group_order)
            }
            if(section == "venn"){
              df$Group <- gsub(":", ": ", df$Group)
              if(!is.null(group_order)) group_order <- gsub(":", ": ", group_order)
              for(i in 1:length(df$Group)){
                df$Group[i] <- paste(strwrap(df$Group[i], width = 15),collapse = "\n")
              }
              for(i in 1:length(unique(df$Group))){
                if(!is.null(group_order)) group_order[i] <- paste(strwrap(group_order[i], width = 15),collapse = "\n")
              }
              df$Group <- gsub(" \\(", "\n\\(", df$Group)
              if(!is.null(group_order)) group_order <- gsub(" \\(", "\n\\(", group_order)
            }
            }
            if(!is.null(group_order)) {
              df$Group <- factor(df$Group, levels=group_order)
              df <- df %>% dplyr::arrange(Group) 
            }
            df$GeneRatio <- DOSE::parse_ratio(df$GeneRatio)
            df <- dplyr::filter(df, !is.na(qvalue))
            df$Description <- gsub("_", " ", df$Description)
            df <- dplyr::mutate(df, x = paste0(Group, 1/(-log10(eval(parse(text = "qvalue"))))))
            df$x <- gsub(":","", df$x)
            df <- dplyr::arrange(df, Group, x)
            idx <- order(df[["Group"]], df[["x"]], decreasing = FALSE)
            df$Description <- factor(df$Description,
                                     levels=rev(unique(df$Description[idx])))
            p1 <- as.grob(ggplot(df, aes(x = Group,y= Description,color=qvalue,size=GeneRatio))+
                            geom_point() +
                            scale_color_continuous(low="red", high="blue",
                                                   guide=guide_colorbar(reverse=TRUE)) +
                            scale_size(range=c(1, 6))+ DOSE::theme_dose(font.size=12)+ylab(NULL)+xlab(NULL)+
                            scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top"))
            p <- plot_grid(p1)
            return(p)
          }
        }
}
cnet_global <- function(data, group, enrich_gene_list, showCategory=5){
      if(is.null(data) || is.null(group) || is.null(enrich_gene_list)){
        return(NULL)
      }else{
        data2 <- dplyr::filter(data, Group == group)
          cnet1 <- enrich_gene_list[[group]]
        if(length(as.data.frame(cnet1)$ID) == 0) {
          p2 <- NULL
        }else{
          p2 <- try(as.grob(clusterProfiler::cnetplot(cnet1,
                                     cex_label_gene = 0.7, cex_label_category = 0.75, showCategory = showCategory,
                                     cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none")))
          if(length(class(p2)) == 1){
            if(class(p2) == "try-error") p2 <- NULL
          }
        }
        p <- plot_grid(p2)
        return(p)
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

getTargetSeq <- function(Species, upstream, downstream){
  library(TFBSTools)
  if(Species == "Mus musculus"){
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  }
  if(Species == "Homo sapiens"){
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  }
  x <- promoters(genes(txdb), upstream = upstream, downstream = downstream)
  return(x)
}

MotifAnalysis <- function(data, Species, org,x){
  withProgress(message = "Motif analysis takes about 2 min per group",{
    library(TFBSTools)
    library(monaLisa)
    library(GenomicRanges)
    library(BiocParallel)
    library(SummarizedExperiment)
    library(JASPAR2020)
    if(Species == "Mus musculus"){
      library(BSgenome.Mmusculus.UCSC.mm10)
      genome = BSgenome.Mmusculus.UCSC.mm10
      tax <- 10090
    }
    if(Species == "Homo sapiens"){
      library(BSgenome.Hsapiens.UCSC.hg19)
      genome = BSgenome.Hsapiens.UCSC.hg19
      tax <- 9606
    }
  pwms <- getMatrixSet(JASPAR2020,
                       opts = list(matrixtype = "PWM",
                                   tax_group = "vertebrates",
                                   species = tax
                                   ))
  df <- data.frame(GeneID = data[,1], Group = data[,2])
  df2 <- list()
  group_file <- length(unique(df$Group))
  perc <- 0
  for(name in unique(df$Group)){
    perc <- perc + 1
    data <- dplyr::filter(df, Group == name)
    my.symbols <- data$GeneID
    group.name <- paste(name, "\n(", length(my.symbols),")",sep = "")
    if(str_detect(my.symbols[1], "ENS") || str_detect(my.symbols[1], "FBgn") ||
       str_detect(my.symbols[1], "^AT.G")){
      if(sum(is.element(no_orgDb, Species)) == 1){
        gene_IDs <- org(Species)
        gene_IDs <- gene_IDs[,-2]
      }else{
      if(str_detect(my.symbols[1], "^AT.G")) key = "TAIR" else key = "ENSEMBL"
      gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                      keytype = key,
                                      columns = c("ENTREZID",key))
      }
      colnames(gene_IDs) <- c("ENSEMBL","gene_id")
    }else{
      if(sum(is.element(no_orgDb, Species)) == 1){
        gene_IDs <- org(Species)
        gene_IDs <- gene_IDs[,-1]
      }else{
      gene_IDs <- AnnotationDbi::select(org, keys = my.symbols,
                                        keytype = "SYMBOL",
                                        columns = c("SYMBOL","ENTREZID"))
      }
      colnames(gene_IDs) <- c("SYMBOL","gene_id")
    }
    y <- subset(x, gene_id %in% gene_IDs$gene_id)
    if(length(rownames(as.data.frame(y))) == 0) stop("Incorrect species")
    seq <- getSeq(genome, y)
    se <- calcBinnedMotifEnrR(seqs = seq,
                              pwmL = pwms,
                              background = "genome",
                              genome = genome,
                              genome.regions = subset(x, ! gene_id %in% gene_IDs$gene_id),
                              genome.oversample = 2,
                              BPPARAM = BiocParallel::SerialParam(RNGseed = 42),
                              verbose = TRUE)
      res <- data.frame(motif.id = elementMetadata(se)$motif.id, motif.name = elementMetadata(se)$motif.name,
                        motif.percentGC = elementMetadata(se)$motif.percentGC,
                        negLog10P = assay(se,"negLog10P"),negLog10Padj = assay(se,"negLog10Padj"), 
                        log2enr = assay(se,"log2enr"),pearsonResid = assay(se,"pearsonResid"),
                        expForegroundWgtWithHits = assay(se,"expForegroundWgtWithHits"),
                        sumForegroundWgtWithHits = assay(se,"sumForegroundWgtWithHits"),
                        sumBackgroundWgtWithHits = assay(se,"sumBackgroundWgtWithHits"),
                        Group = group.name)
      df2[[name]] <- res
    incProgress(1/group_file, message = paste("Finish motif analysis of Group '", name, "', ", perc, "/", group_file,sep = ""))
  }
  return(df2)
})
}

MotifRegion <- function(data, target_motif, Species, x){
  if(Species == "Mus musculus"){
    genome = BSgenome.Mmusculus.UCSC.mm10
  }
  if(Species == "Homo sapiens"){
    genome = BSgenome.Hsapiens.UCSC.hg19
  }
  df <- data.frame(GeneID = data[,1], Group = data[,2])
  target_motif$Group <- gsub(" ", "\n", target_motif$Group)
  name <- gsub("\\\n.+$", "", target_motif$Group)
  data <- dplyr::filter(df, Group %in% name)
  my.symbols <- data$GeneID
  if(str_detect(my.symbols[1], "ENS") || str_detect(my.symbols[1], "FBgn") ||
     str_detect(my.symbols[1], "^AT.G")){
    if(sum(is.element(no_orgDb, Species)) == 1){
      gene_IDs <- org(Species)
      gene_IDs <- data.frame(gene_id = gene_IDs$ENTREZID, ENSEMBL = gene_IDs$ENSEMBL)
    }else{
    if(str_detect(my.symbols[1], "^AT.G")) key = "TAIR" else key = "ENSEMBL"
    gene_IDs<-AnnotationDbi::select(org(Species),keys = my.symbols,
                                    keytype = key,
                                    columns = c("ENTREZID",key))
    }
    colnames(gene_IDs) <- c("gene_id","ENSEMBL")
  }else{
    if(sum(is.element(no_orgDb, Species)) == 1){
      gene_IDs <- org(Species)
      gene_IDs <- gene_IDs[,-1]
    }else{
    gene_IDs <- AnnotationDbi::select(org(Species), keys = my.symbols,
                                      keytype = "SYMBOL",
                                      columns = c("SYMBOL","ENTREZID"))
    }
    colnames(gene_IDs) <- c("SYMBOL","gene_id")
  }
  y <- subset(x, gene_id %in% gene_IDs$gene_id)
  if(length(rownames(as.data.frame(y))) == 0) stop("Incorrect species")
  seq <- getSeq(genome, y)
  pfm <- getMatrixByID(JASPAR2020,target_motif$motif.id)
  pwm <- toPWM(pfm)
  res <- findMotifHits(query = pwm,
                       subject = seq,
                       min.score = 6.0,
                       method = "matchPWM",
                       BPPARAM = BiocParallel::SerialParam()) %>% as.data.frame()
  my.symbols <- as.character(res$seqnames)
  gene_IDs <- AnnotationDbi::select(org(Species), keys = my.symbols,
                                    keytype = "ENTREZID",
                                    columns = c("SYMBOL","ENTREZID"))
  colnames(gene_IDs) <- c("seqnames","SYMBOL")
  res2<-merge(gene_IDs,res,by="seqnames")
  res2 <- res2[,-1]
  return(res2)
}

Motifplot <- function(df2, showCategory=5,padj,data,group_order){
  df <- data.frame(matrix(rep(NA, 11), nrow=1))[numeric(0), ]
  data <- data.frame(GeneID = data[,1], Group = data[,2])
  for(name in names(df2)){
    res <- df2[[name]]
    print(name)
    data2 <- dplyr::filter(data, Group == name)
    my.symbols <- data2$GeneID
    if(!is.null(group_order)) group_order[which(group_order == name)] <- paste(name, "\n(", length(my.symbols),")",sep = "")
    res <- dplyr::filter(res, X1 > -log10(padj))
    res <- res %>% dplyr::arrange(-X1.1)
    if(length(rownames(res)) > showCategory){
      res <- res[1:showCategory,]
    }
    df <- rbind(df, res)
  }
  colnames(df) <- c("motif.id", "motif.name","motif.percentGC", "negLog10P", "negLog10Padj", "log2enr",
                    "pearsonResid", "expForegroundWgtWithHits", "sumForegroundWgtWithHits", "sumBackgroundWgtWithHits",
                    "Group")
  if(length(df$motif.id) == 0){
    return(NULL)
  }else{
    df$Group <- gsub("_", " ", df$Group)
    print(unique(df$Group))
    
    if(!is.null(group_order)) group_order <- gsub("_", " ", group_order)
    print(group_order)
    if(!is.null(group_order)) df$Group <- factor(df$Group,levels=group_order)
  df$padj <- 10^(-df$negLog10Padj)
  df <- dplyr::mutate(df, x = paste0(Group, 1/-log10(eval(parse(text = "padj")))))
  df$x <- gsub(":","", df$x)
  df <- dplyr::arrange(df, x)
  idx <- order(df[["Group"]], df[["x"]], decreasing = FALSE)
  df$motif.name <- factor(df$motif.name,
                          levels=rev(unique(df$motif.name[idx])))
  d <- ggplot(df, aes(x = Group,y= motif.name,color=padj,size=log2enr))+
    geom_point() +
    scale_color_continuous(low="red", high="blue",
                           guide=guide_colorbar(reverse=TRUE)) +
    scale_size(range=c(1, 6))+ DOSE::theme_dose(font.size=15)+ylab(NULL)+xlab(NULL) +
    scale_y_discrete(labels = label_wrap_gen(30)) + scale_x_discrete(position = "top")+
    theme(plot.margin=margin(l=-0.75,unit="cm"))
  
  df <- df %>% distinct(motif.id, .keep_all = T)
  width.seqlogo = 2
  highlight <- NULL
  clres <- FALSE
  optsL <- list(ID = df$motif.id)
  pfm1 <- TFBSTools::getMatrixSet(JASPAR2020, opts = optsL)
  maxwidth <- max(vapply(TFBSTools::Matrix(pfm1), ncol, 0L))
  grobL <- lapply(pfm1, seqLogoGrob, xmax = maxwidth, xjust = "center")
  hmSeqlogo <- HeatmapAnnotation(logo = annoSeqlogo(grobL = grobL, 
                                                    which = "row", space = unit(1, "mm"),
                                                    width = unit(width.seqlogo, "inch")), 
                                 show_legend = FALSE, show_annotation_name = FALSE, 
                                 which = "row")
  tmp <- matrix(rep(NA, length(df$motif.id)),ncol = 1, 
                dimnames = list(df$motif.name, NULL))
  
  hmMotifs <- Heatmap(matrix = tmp, name = "names", width = unit(0, "inch"), 
                      na_col = NA, col = c(`TRUE` = "green3",`FALSE` = "white"), 
                      cluster_rows = clres, show_row_dend = show_dendrogram, 
                      cluster_columns = FALSE, show_row_names = TRUE, row_names_side = "left", 
                      show_column_names = FALSE, show_heatmap_legend = FALSE,
                      left_annotation = hmSeqlogo)
  h <- grid.grabExpr(print(hmMotifs),wrap.grobs=TRUE)
  p <- plot_grid(plot_grid(NULL, h, ncol = 1, rel_heights = c(0.05:10)),as.grob(d))
  
  return(p)
  }
}

GOIheatmap <- function(data.z, show_row_names = TRUE, type = NULL, GOI = NULL, all=FALSE){
  if(length(rownames(data.z)) <= 50) {
    if(!is.null(type)) {if(type == "ALL") show_row_names = FALSE else show_row_names = TRUE}
  }else{
    show_row_names = FALSE
  }
  if(!is.null(all)){
  if(all == TRUE) show_row_names = TRUE
  }
  ht <- Heatmap(data.z, name = "z-score",column_order = colnames(data.z),
                clustering_method_columns = 'ward.D2',use_raster = TRUE,
                show_row_names = show_row_names, show_row_dend = F,column_names_side = "top",
                row_names_gp = gpar(fontface = "italic"))
  if(!is.null(type)) {
  if(type == "ALL" && !is.null(GOI) && all != TRUE) {
    indexes <- which(rownames(data.z) %in% GOI)
    labels <- rownames(data.z)[indexes]
    ht <- ht + rowAnnotation(
      link = anno_mark(at = indexes, labels = labels,which="row",link_width = unit(1, "cm"),
                       labels_gp = gpar(fontface = "italic")),
      width = unit(1, "cm") + max_text_width(labels)
    )
  }
  }
  return(ht)
}

ensembl2symbol <- function(data,Species,Ortholog,gene_type,org, merge=TRUE, rowname = TRUE){
  if(Species != "not selected"){
  if(gene_type != "SYMBOL"){
  data <- as.data.frame(data)
  if(sum(is.element(no_orgDb, Species)) == 1){
    gene_IDs <- try(Ortholog[,-3])
    if(length(class(gene_IDs)) == 1){
      if(class(gene_IDs) == "try-error") validate("biomart has encountered an unexpected server error.
                                                \nPlease try using a different 'biomart host' or try again later.")
    }
  }else{
    if(gene_type == "ENSEMBL"){
      my.symbols <- gsub("\\..*","", rownames(data))
      if(str_detect(my.symbols[1], "^AT.G")) key = "TAIR" else key = "ENSEMBL"
      gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                      keytype = key,
                                      columns = c(key,"SYMBOL"))
    }
  }
    colnames(gene_IDs) <- c("Row.names","SYMBOL")
    data$Row.names <- gsub("\\..*","", rownames(data))
    gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
    if(merge == TRUE){
      data2 <- merge(data, gene_IDs, by="Row.names")
      rownames(data2) <- data2$Row.names
      data <- data2[,-1]
    }else {
      if(rowname == TRUE) rownames(gene_IDs) <- gene_IDs$Row.names
      data <- gene_IDs
    }
    return(data) 
  }
  }
  return(data)
  }

gene_type <- function(my.symbols,org,Species){
  if(Species != "not selected"){
    if(sum(is.element(no_orgDb, Species)) != 1){
    ENSEMBL<-try(AnnotationDbi::select(org,keys = my.symbols,
                                       keytype = "ENSEMBL",
                                       columns = c("ENSEMBL", "ENTREZID")))
    SYMBOL <-try(AnnotationDbi::select(org,keys = my.symbols,
                                       keytype = "SYMBOL",
                                       columns = c("SYMBOL", "ENTREZID")))
    if(class(ENSEMBL) == "try-error" && class(SYMBOL) != "try-error") {type <- "SYMBOL"
    }else if(class(ENSEMBL) != "try-error" && class(SYMBOL) == "try-error") {type <- "ENSEMBL"
    }else if(class(ENSEMBL) == "try-error" && class(SYMBOL) == "try-error") {validate("Cannot identify gene IDs. Please check the 'Species' and use the 'Official gene symbol' or 'ENSEMBL ID' for gene names.")
    }else{
      if(dim(ENSEMBL)[1] > dim(SYMBOL)[1]) type <- "ENSEMBL" else type <- "SYMBOL"
    }
    }else type <- "non-model organism"
  }else type <- "not selected"
  return(type)
}

consensus_kmeans = function(mat, centers, km_repeats) {
  partition_list = lapply(seq_len(km_repeats), function(i) {
    as.cl_hard_partition(kmeans(mat, centers))
  })
  partition_list = cl_ensemble(list = partition_list)
  partition_consensus = cl_consensus(partition_list)
  as.vector(cl_class_ids(partition_consensus)) 
}

corr_plot_pair <- function(data,corr_color,GOI_x,GOI_y){
  p1 <- NULL
  p2 <- NULL
  if(corr_color == ""){
    p1 <- ggplot(data, aes(x=log10(.data[[GOI_x]]+1),y=log10(.data[[GOI_y]]+1))) +
      geom_smooth(method=lm, se=FALSE, color='#2C3E50',linetype="dashed",size=0.5)
  }else if(corr_color == "sample_name"){
    label <- gsub("\\_.+$", "", rownames(data))
    p1 <- ggplot(data, aes(x=log10(.data[[GOI_x]]+1),y=log10(.data[[GOI_y]]+1), col=label)) +
      geom_smooth(method=lm, se=FALSE, color='#2C3E50',linetype="dashed",size=0.5)
  }else {
    p2 <- ggplot(data, aes(x=log10(.data[[GOI_x]]+1),y=log10(.data[[GOI_y]]+1), col=log10(.data[[corr_color]]))) +
      geom_smooth(method=lm, se=FALSE, color='#2C3E50',linetype="dashed",size=0.5)
  }
  if(!is.null(p1)){
    p <- p1 +
      geom_point()+ 
      theme_bw()+ 
      xlab(paste(strwrap(paste0("log10(", GOI_x," + 1)"), width = 30),collapse = "\n"))+
      ylab(paste(strwrap(paste0("log10(", GOI_y," + 1)"), width = 30),collapse = "\n"))+
      theme(legend.position = "top" , legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 12),
            axis.text.y= ggplot2::element_text(size = 12),
            text = ggplot2::element_text(size = 15),
            title = ggplot2::element_text(size = 15),
            plot.title = element_text(size = 15))
  }
  if(!is.null(p2)){
    p <- p2 +
      geom_point()+
      scale_color_continuous(low="blue", high="red")+ 
      theme_bw()+
      ggtitle(paste(strwrap(paste0("color = log10(", corr_color," + 1)"), width = 30),collapse = "\n"))+ 
      xlab(paste(strwrap(paste0("log10(", GOI_x," + 1)"), width = 30),collapse = "\n"))+
      ylab(paste(strwrap(paste0("log10(", GOI_y," + 1)"), width = 30),collapse = "\n"))+
      theme(legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 12),
            axis.text.y= ggplot2::element_text(size = 12),
            text = ggplot2::element_text(size = 15),
            title = ggplot2::element_text(size = 15),
            plot.title = element_text(size = 15))
  }
  return(p)
}


