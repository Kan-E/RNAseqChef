## limma, 数字が頭文字だとエラー

## limma, 数字から始まるサンプル名によるエラー
## proteinIDへの対応
## Receptor-Ligand interaction


##isoform経過：gene_type, no_org_ID, ensembl2geneID完了
library(shiny)
library(shinyBS, verbose = FALSE)
library('shinyjs', verbose = FALSE)
library(DT)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)
library(tidyr)
library(rstatix)
library(multcomp)
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
library(DEGreport)
library(ComplexHeatmap)
library(biomaRt)
library(limma)
library(colorspace)
library(clue)
library(ggrastr) ##devtools::install_github('VPetukhov/ggrastr')
library(statmod)
library(eulerr)
orgdb_package_names <- c(
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "org.Rn.eg.db",
  "org.Dm.eg.db",
  "org.Ce.eg.db",
  "org.Xl.eg.db",
  "org.Bt.eg.db",
  "org.Cf.eg.db",
  "org.Dr.eg.db",
  "org.Gg.eg.db",
  "org.Mmu.eg.db",
  "org.Pt.eg.db",
  "org.Sc.sgd.db",
  "org.At.tair.db"
)

orgdb_package_by_species <- c(
  "Homo sapiens" = "org.Hs.eg.db",
  "Mus musculus" = "org.Mm.eg.db",
  "Rattus norvegicus" = "org.Rn.eg.db",
  "Xenopus laevis" = "org.Xl.eg.db",
  "Drosophila melanogaster" = "org.Dm.eg.db",
  "Caenorhabditis elegans" = "org.Ce.eg.db",
  "Bos taurus" = "org.Bt.eg.db",
  "Canis lupus familiaris" = "org.Cf.eg.db",
  "Danio rerio" = "org.Dr.eg.db",
  "Gallus gallus" = "org.Gg.eg.db",
  "Macaca mulatta" = "org.Mmu.eg.db",
  "Pan troglodytes" = "org.Pt.eg.db",
  "Saccharomyces cerevisiae" = "org.Sc.sgd.db",
  "Arabidopsis thaliana" = "org.At.tair.db"
)

rsconnect_orgdb_dependency_hints <- function() {
  if (FALSE) {
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    library(org.Rn.eg.db)
    library(org.Xl.eg.db)
    library(org.Dm.eg.db)
    library(org.Ce.eg.db)
    library(org.Bt.eg.db)
    library(org.Cf.eg.db)
    library(org.Dr.eg.db)
    library(org.Gg.eg.db)
    library(org.Mmu.eg.db)
    library(org.Pt.eg.db)
    library(org.Sc.sgd.db)
    library(org.At.tair.db)
  }
  invisible(NULL)
}

rsconnect_orgdb_dependency_hints()

rsconnect_dynamic_dependency_hints <- function() {
  if (FALSE) {
    library(enrichplot)
    library(clusterProfiler)
    library(KEGGREST)
    library(GO.db)
    library(GSVA)
    library(digest)
    library(pdftools)
    library(png)
    library(rtracklayer)
    library(rmarkdown)
    library(dorothea)
    library(DOSE)
    library(ggdendro)
    library(umap)
    library(genefilter)
    library(IHW)
    library(qvalue)
    library(BiocParallel)
    library(SummarizedExperiment)
    library(magick)
  }
  invisible(NULL)
}

rsconnect_dynamic_dependency_hints()

load_orgdb_package <- function(package_name) {
  package_name <- as.character(package_name)[1]
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop("Package not installed: ", package_name)
  }
  getExportedValue(package_name, package_name)
}

register_lazy_orgdb_packages <- function(env) {
  for (package_name in orgdb_package_names) {
    local({
      pkg <- package_name
      delayedAssign(
        pkg,
        load_orgdb_package(pkg),
        assign.env = env
      )
    })
  }
}

register_lazy_orgdb_packages(environment())
bioc_repos <- tryCatch(
  suppressWarnings(BiocManager::repositories()),
  error = function(e) getOption("repos")
)
options(repos = bioc_repos)
app_root <- normalizePath(getwd(), mustWork = TRUE)

app_file <- function(...) {
  file.path(app_root, ...)
}

copy_to_temp_if_exists <- function(source_path, target_name = basename(source_path)) {
  source_file <- app_file(source_path)
  if (file.exists(source_file)) {
    file.copy(source_file, file.path(tempdir(), target_name), overwrite = TRUE)
  }
}

read_local_or_remote_table <- function(local_path, remote_path = NULL, ...) {
  local_file <- app_file(local_path)
  if (file.exists(local_file)) {
    return(read.table(local_file, ...))
  }
  if (!is.null(remote_path)) {
    return(read.table(remote_path, ...))
  }
  stop("Missing file: ", local_path)
}

example_data_path <- function(local_path, remote_path = NULL) {
  local_file <- app_file(local_path)
  if (file.exists(local_file)) {
    return(local_file)
  }
  remote_path
}

copy_to_temp_if_exists("Rmd/pair_report.Rmd")
copy_to_temp_if_exists("Rmd/pair_batch_report.Rmd")
copy_to_temp_if_exists("Rmd/3conditions_report.Rmd")
copy_to_temp_if_exists("Rmd/multi_report.Rmd")
copy_to_temp_if_exists("dds.rds")
msigdbr_species <- msigdbr::msigdbr_species()$species_name
gene_set_list <- c("MSigDB Hallmark", "KEGG", "Reactome", "PID (Pathway Interaction Database)",
                   "BioCarta","WikiPathways", "GO biological process", 
                   "GO cellular component","GO molecular function", "Human phenotype ontology", 
                   "DoRothEA regulon (activator)", "DoRothEA regulon (repressor)",
                   "Transcription factor targets", "miRNA target","Position",
                   "CGP (chemical and genetic pertubations)","ImmuneSigDB","Macrophage (444 gene sets from ImmuneSigDB)","VAX (vaccine response)",
                   "Cell type signature","Custom gene set")
biomart_data <- read_local_or_remote_table("data/non-model.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/non-model.txt", sep = "\t", row.names = 1, header = T, quote = "")
biomart_plants <- read_local_or_remote_table("data/non-model_plants.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/non-model_plants.txt", sep = "\t", header = T, quote = "")
biomart_fungi <- read_local_or_remote_table("data/non-model_fungi.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/non-model_fungi.txt", sep = "\t", header = T, quote = "")
biomart_metazoa <- read_local_or_remote_table("data/non-model_metazoa.txt", "https://raw.githubusercontent.com/Kan-E/RNAseqChef/main/data/non-model_metazoa.txt", sep = "\t", header = T)
no_orgDb_animals<-c(biomart_data$Scientific_common_name)
no_orgDb_plants<-c(biomart_plants$Scientific_common_name)
no_orgDb_fungi<-c(biomart_fungi$Scientific_common_name)
no_orgDb_metazoa<-c(biomart_metazoa$Scientific_common_name)
force_orgdb_species <- c("Arabidopsis thaliana")
build_non_model_species_aliases <- function(...) {
  tables <- list(...)
  alias_tables <- lapply(tables, function(df) {
    canonical <- as.character(df$Scientific_common_name)
    cols <- intersect(c("Scientific_common_name", "Scientific_name", "Common_name", "description"), colnames(df))
    do.call(rbind, lapply(cols, function(col) {
      data.frame(alias = as.character(df[[col]]), canonical = canonical, stringsAsFactors = FALSE)
    }))
  })
  alias_df <- unique(do.call(rbind, alias_tables))
  alias_df <- alias_df[!is.na(alias_df$alias) & nzchar(alias_df$alias), , drop = FALSE]
  alias_df$alias_lower <- tolower(alias_df$alias)
  alias_df
}
species_alias_df <- build_non_model_species_aliases(biomart_data, biomart_metazoa, biomart_plants, biomart_fungi)
species_alias_df <- species_alias_df[!species_alias_df$canonical %in% force_orgdb_species, , drop = FALSE]
species_list_nonmodel <- setdiff(unique(c(no_orgDb_animals,no_orgDb_metazoa,no_orgDb_plants,no_orgDb_fungi)), force_orgdb_species)
no_orgDb <- unique(as.character(species_alias_df$alias))
species_list <- c("not selected", "Homo sapiens", "Mus musculus", "Rattus norvegicus", 
                  "Drosophila melanogaster", "Caenorhabditis elegans","Bos taurus","Canis lupus familiaris",
                  "Danio rerio","Gallus gallus","Macaca mulatta",
                  "Pan troglodytes","Saccharomyces cerevisiae","Xenopus laevis","Arabidopsis thaliana",
                  species_list_nonmodel)
normalize_species_name <- function(species){
  if(length(species) != 1 || is.na(species) || !nzchar(species)) return(species)
  if(species %in% species_list_nonmodel) return(species)
  idx <- match(species, species_alias_df$alias)
  if(!is.na(idx)) return(species_alias_df$canonical[[idx]])
  idx <- match(tolower(species), species_alias_df$alias_lower)
  if(!is.na(idx)) return(species_alias_df$canonical[[idx]])
  species
}
normalize_species_input <- function(species){
  if(is.null(species) || !length(species) || all(is.na(species))) return("not selected")
  species <- as.character(species[[1]])
  if(!nzchar(species)) return("not selected")
  normalize_species_name(species)
}

normalize_ortholog_input <- function(ortholog){
  if(is.null(ortholog) || !length(ortholog) || all(is.na(ortholog))) return(NULL)
  ortholog <- as.character(ortholog[[1]])
  if(!nzchar(ortholog) || ortholog == "not selected") return(NULL)
  ortholog
}
sgd_symbol_column <- function(org_obj) {
  if (!is.null(org_obj) &&
      !is.null(org_obj$packageName) &&
      identical(org_obj$packageName, "org.Sc.sgd.db")) {
    "GENENAME"
  } else {
    "SYMBOL"
  }
}
gene_primary_key <- function(my.symbols) {
  if (length(my.symbols) == 0 || is.na(my.symbols[1]) || !nzchar(my.symbols[1])) {
    return("ENSEMBL")
  }
  first_id <- gsub("\\..*$", "", as.character(my.symbols[1]))
  if (str_detect(first_id, "^AT[1-5MC]G[0-9]+$")) {
    return("TAIR")
  }
  if (str_detect(first_id, "^WBGene[0-9]+$")) {
    return("WORMBASE")
  }
  if (str_detect(first_id, "^S[0-9]{9}$")) {
    return("SGD")
  }
  "ENSEMBL"
}
is_ensembl_like_id <- function(my.symbols, org_obj = NULL) {
  if (length(my.symbols) == 0 || is.na(my.symbols[1]) || !nzchar(my.symbols[1])) {
    return(FALSE)
  }
  first_id <- gsub("\\..*$", "", as.character(my.symbols[1]))
  if (str_detect(first_id, "^ENS") ||
      str_detect(first_id, "^FBgn") ||
      str_detect(first_id, "^AT[1-5MC]G[0-9]+$") ||
      str_detect(first_id, "^WBGene[0-9]+$")) {
    return(TRUE)
  }
  if (!is.null(org_obj) &&
      !is.null(org_obj$packageName) &&
      identical(org_obj$packageName, "org.Sc.sgd.db")) {
    return(str_detect(first_id, "^Y[A-Z]{2}[0-9]{3}[CW](?:-[A-Z])?$") ||
             str_detect(first_id, "^S[0-9]{9}$"))
  }
  FALSE
}
is_transcript_like_id <- function(my.symbols) {
  if (length(my.symbols) == 0) {
    return(FALSE)
  }
  first_id <- as.character(my.symbols[1])
  if (is.na(first_id) || !nzchar(first_id)) {
    return(FALSE)
  }
  first_id <- gsub("\\..*$", "", first_id)
  str_detect(first_id, "^ENS[A-Z0-9]*T[0-9]+$") ||
    str_detect(first_id, "^FBtr") ||
    str_detect(first_id, "^(NM|NR|XM|XR)_[0-9]+$")
}
prepare_ranked_gene_list <- function(data, value_col = "log2FoldChange", id_col = "ENTREZID") {
  if (is.null(data) || !is.data.frame(data) || !nrow(data) ||
      !all(c(value_col, id_col) %in% colnames(data))) {
    return(NULL)
  }
  stats <- suppressWarnings(as.numeric(data[[value_col]]))
  ids <- as.character(data[[id_col]])
  keep <- !is.na(stats) & is.finite(stats) & !is.na(ids) & nzchar(ids)
  if (!any(keep)) {
    return(NULL)
  }
  ranked <- data.frame(
    gene_id = ids[keep],
    stat = stats[keep],
    stringsAsFactors = FALSE
  )
  ranked <- ranked[order(abs(ranked$stat), decreasing = TRUE), , drop = FALSE]
  ranked <- ranked[!duplicated(ranked$gene_id), , drop = FALSE]
  geneList <- ranked$stat
  names(geneList) <- ranked$gene_id
  sort(geneList, decreasing = TRUE)
}
gene_set_cache <- new.env(parent = emptyenv())
read_by_extension <- function(tmp, use_row_names = FALSE, na_strings = NULL){
  ext <- tools::file_ext(tmp)
  if(ext == "xlsx"){
    df <- try(as.data.frame(readxl::read_xlsx(tmp)), silent = TRUE)
    if(use_row_names && !inherits(df, "try-error")){
      original_colnames <- colnames(df)
      rowname_try <- try(data.frame(row.names = df[,1]), silent = TRUE)
      if(!inherits(rowname_try, "try-error")){
        if(dim(df)[2] == 2){
          df <- data.frame(row.names = df[,1], a = df[,2])
          colnames(df)[1] <- original_colnames[2]
        }else{
          rownames(df) <- df[,1]
          df <- df[,-1]
          colnames(df) <- gsub("-",".",colnames(df))
        }
      }
    }
    return(df)
  }
  if(ext == "csv"){
    return(try(read.csv(tmp, header = TRUE, sep = ",", row.names = if(use_row_names) 1 else NULL, quote = "", na.strings = na_strings), silent = TRUE))
  }
  if(ext == "txt" || ext == "tsv"){
    return(try(read.table(tmp, header = TRUE, sep = "\t", row.names = if(use_row_names) 1 else NULL, quote = "", na.strings = na_strings), silent = TRUE))
  }
  try(stop("unsupported extension"), silent = TRUE)
}

clean_x_prefix <- function(df){
  if(length(colnames(df)) != 0 && str_detect(colnames(df)[1], "^X\\.")){
    colnames(df) <- sub("^X\\.", "", colnames(df))
  }
  df
}

read_df <- function(tmp, Species=NULL){
  if(is.null(tmp)) {
    return(NULL)
  }else{
    if(sum(is.element(c("csv","txt","tsv","xlsx"), tools::file_ext(tmp))) == 0) validate("Error: the file extension is in an unexpected format.")
    df <- read_by_extension(tmp, use_row_names = TRUE)
    if(inherits(df, "try-error") || length(grep("Protein.Ids", colnames(df))) != 0 || length(grep("First.Protein.Description", colnames(df)))) {
      df <- read_by_extension(tmp, use_row_names = FALSE)
      if(!inherits(df, "try-error")) {
        if(dim(df)[2] != 0){
          if(length(grep("Protein.Ids", colnames(df))) != 0 || length(grep("First.Protein.Description", colnames(df)))){
            df <- df %>% distinct(Genes, .keep_all = T)
            df[is.na(df)] <- 0
            if(length(grep("Protein.Ids", colnames(df))) != 0) df[df$Genes == "",]$Genes <- gsub("\\_.+$", "", df[df$Genes == "",]$Protein.Ids)
            rownames(df) <- df$Genes
            if(length(grep("Protein.Group", colnames(df))) != 0) df <- df[, - which(colnames(df) == "Protein.Group")]
            if(length(grep("Protein.Ids", colnames(df))) != 0) df <- df[, - which(colnames(df) == "Protein.Ids")]
            if(length(grep("Protein.Names", colnames(df))) != 0) df <- df[, - which(colnames(df) == "Protein.Names")]
            if(length(grep("First.Protein.Description", colnames(df))) != 0) df <- df[, - which(colnames(df) == "First.Protein.Description")]
            if(length(grep("Genes", colnames(df))) != 0) df <- df[, - which(colnames(df) == "Genes")]
            if(length(grep("Species", colnames(df))) != 0) df <- df[, - which(colnames(df) == "Species")]
            df <- df[!str_detect(rownames(df), ";"),]
            df <- df[rownames(df) != "",]
          }else validate("Error: There are duplicated genes or 'NA' in the uploaded data. Please fix them.")
        }
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
    df <- clean_x_prefix(df)
    if(length(grep("padj", colnames(df))) == 0 || length(grep("log2FoldChange", colnames(df))) == 0){
    df[is.na(df)] <- 0
    }
    if(sum(!str_detect(rownames(df),"\\.")) == 0  & !str_detect(rownames(df)[1],"chr")) rownames(df) <- gsub("\\..+$", "",rownames(df))
    return(df)
    }
  }
}
read_gene_list <- function(tmp){
  if(is.null(tmp)) {
    return(NULL)
  }else{
    df <- read_by_extension(tmp, use_row_names = FALSE)
    if(inherits(df, "try-error")) validate("Error: failed to load the uploaded gene list.")
    rownames(df) = gsub("\"", "", rownames(df))
    df <- clean_x_prefix(df)
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
    for(i in seq_along(unique_col)){
      cond <- length(which(colnames(row) == unique_col[i]))
      for(k in seq_len(cond)){
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
    for(i in seq_along(unique_col)){
      cond <- length(which(meta[,1] == unique_col[i]))
      for(k in seq_len(cond)){
        meta[total + k,1] <- paste0(meta[total + k,1], "_", k)
      }
      total <- total + cond
    }
  }
  return(meta)
  }
}

gene_list_convert_for_enrichment <- function(gene_type,data, Ortholog,Isoform,org,Species){
    if(is.null(data) || Species == "not selected"){
      return(NULL)
    }else{
      df <- data.frame(GeneID = data[,1], Group = data[,2])
      df$GeneID <- gsub("\\..*","", df$GeneID)
      my.symbols <- df$GeneID
      if(gene_type != "SYMBOL"){
        if(sum(is.element(no_orgDb, Species)) == 1){
          gene_IDs <- Ortholog
        }else if(gene_type == "isoform"){
          gene_IDs <- Isoform
        }else{
        key <- gene_primary_key(my.symbols)
        SYMBOL <- sgd_symbol_column(org)
        gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                        keytype = key,
                                        columns = c(key,SYMBOL, "ENTREZID"))
        }
        colnames(gene_IDs) <- c("GeneID","SYMBOL", "ENTREZID")
      }else{
        if(sum(is.element(no_orgDb, Species)) == 1){
          gene_IDs <- Ortholog[,-1]
        }else if(gene_type == "isoform"){
          gene_IDs <- Isoform[,-1]
        }else{
          SYMBOL <- sgd_symbol_column(org)
        gene_IDs <- AnnotationDbi::select(org, keys = my.symbols,
                                          keytype = SYMBOL,
                                          columns = c("ENTREZID", SYMBOL))
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
    convert <- with_biomart_timeout(useMart("ensembl", dataset = set, host = "https://dec2021.archive.ensembl.org"))
    human <- with_biomart_timeout(useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org"))
    genes2 <- with_biomart_timeout(getLDS(attributes = c("entrezgene_id"), filters = "entrezgene_id",
                                          values = genes, mart = human,
                                          attributesL = c("entrezgene_id"),
                                          martL = convert, uniqueRows = TRUE))
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
biomart_timeout_seconds <- 15

with_biomart_timeout <- function(expr, seconds = biomart_timeout_seconds) {
  old_timeout <- getOption("timeout")
  options(timeout = max(1L, as.integer(seconds)))
  setTimeLimit(elapsed = seconds, transient = TRUE)
  on.exit({
    options(timeout = old_timeout)
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  }, add = TRUE)
  force(expr)
}

biomart_host_candidates <- function(mart, host) {
  if (identical(mart, "plants_mart")) {
    return(unique(c(host, ensembl_archive_plants)))
  }
  if (identical(mart, "fungi_mart")) {
    return(unique(c(host, ensembl_archive_fungi)))
  }
  if (identical(mart, "metazoa_mart")) {
    return(unique(c(host, ensembl_archive_metazoa)))
  }
  unique(c(host))
}

biomart_prefers_martservice <- function(mart) {
  identical(mart, "plants_mart")
}

pick_first_match <- function(x) {
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

ortholog_dataset_name <- function(species) {
  switch(species,
         "Mus musculus" = "mmusculus_gene_ensembl",
         "Homo sapiens" = "hsapiens_gene_ensembl",
         "Rattus norvegicus" = "rnorvegicus_gene_ensembl",
         "Drosophila melanogaster" = "dmelanogaster_eg_gene",
         "Caenorhabditis elegans" = "celegans_eg_gene",
         "Bos taurus" = "btaurus_gene_ensembl",
         "Canis lupus familiaris" = "clfamiliaris_gene_ensembl",
         "Danio rerio" = "drerio_gene_ensembl",
         "Gallus gallus" = "ggallus_gene_ensembl",
         "Macaca mulatta" = "mmulatta_gene_ensembl",
         "Pan troglodytes" = "ptroglodytes_gene_ensembl",
         "Saccharomyces cerevisiae" = "scerevisiae_eg_gene",
         "Arabidopsis thaliana" = "athaliana_eg_gene",
         NULL)
}

ortholog_orgdb <- function(species) {
  switch(species,
         "Mus musculus" = org.Mm.eg.db,
         "Homo sapiens" = org.Hs.eg.db,
         "Rattus norvegicus" = org.Rn.eg.db,
         "Drosophila melanogaster" = org.Dm.eg.db,
         "Caenorhabditis elegans" = org.Ce.eg.db,
         "Bos taurus" = org.Bt.eg.db,
         "Canis lupus familiaris" = org.Cf.eg.db,
         "Danio rerio" = org.Dr.eg.db,
         "Gallus gallus" = org.Gg.eg.db,
         "Macaca mulatta" = org.Mmu.eg.db,
         "Pan troglodytes" = org.Pt.eg.db,
         "Saccharomyces cerevisiae" = org.Sc.sgd.db,
         "Arabidopsis thaliana" = org.At.tair.db,
         NULL)
}

ortholog_keytype <- function(species) {
  switch(species,
         "Arabidopsis thaliana" = "TAIR",
         "Saccharomyces cerevisiae" = "ORF",
         "ENSEMBL")
}

ortholog_dataset_prefix <- function(dataset) {
  if (is.null(dataset) || !nzchar(dataset)) {
    return(dataset)
  }
  if (grepl("_eg_gene$", dataset)) {
    return(sub("_gene$", "", dataset))
  }
  sub("_gene_ensembl$", "", dataset)
}

direct_martservice_query <- function(host, mart, dataset, attributes, filter_name = NULL, values = NULL, unique_rows = TRUE) {
  host_candidates <- biomart_host_candidates(mart, host)
  if (length(attributes) == 0) {
    return(data.frame())
  }

  attr_xml <- paste0("<Attribute name = '", attributes, "'/>", collapse = "")
  filter_xml <- ""
  if (!is.null(filter_name) && length(filter_name) == 1 && !is.null(values)) {
    values <- unique(as.character(values))
    values <- values[!is.na(values) & nzchar(values)]
    filter_xml <- paste0("<Filter name = '", filter_name, "' value = '", paste(values, collapse = ","), "' />")
  }
  xml <- paste0(
    "<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query>",
    "<Query virtualSchemaName = '", mart, "' uniqueRows = '", if (unique_rows) "1" else "0",
    "' count='0' datasetConfigVersion='0.6' header='1' formatter='TSV'>",
    "<Dataset name = '", dataset, "'>", attr_xml, filter_xml, "</Dataset></Query>"
  )

  last_error <- NULL
  for (candidate in host_candidates) {
    service_url <- paste0(sub("/+$", "", candidate), "/biomart/martservice?query=", utils::URLencode(xml, reserved = TRUE))
    result <- try(with_biomart_timeout(utils::read.delim(service_url,
                                                         header = TRUE,
                                                         sep = "\t",
                                                         quote = "",
                                                         stringsAsFactors = FALSE,
                                                         check.names = FALSE,
                                                         fill = TRUE)), silent = TRUE)
    if (inherits(result, "try-error")) {
      last_error <- result
      next
    }
    if (ncol(result) == length(attributes)) {
      colnames(result) <- attributes
    }
    if (ncol(result) == 0) {
      result <- as.data.frame(matrix(nrow = 0, ncol = length(attributes)))
      colnames(result) <- attributes
    }
    return(result)
  }
  if (!is.null(last_error)) {
    stop(last_error)
  }
  stop("biomart query failed")
}

open_biomart_dataset <- function(mart, dataset, host) {
  host_candidates <- biomart_host_candidates(mart, host)

  primary <- NULL
  if (identical(mart, "ensembl")) {
    for (candidate in host_candidates) {
      primary <- try(with_biomart_timeout(useMart(mart, dataset = dataset, host = candidate)), silent = TRUE)
      if (!inherits(primary, "try-error")) {
        return(primary)
      }
    }
  } else {
    for (candidate in host_candidates) {
      db <- try(with_biomart_timeout(useMart(mart, host = candidate)), silent = TRUE)
      if (inherits(db, "try-error")) {
        primary <- db
        next
      }
      primary <- try(with_biomart_timeout(useDataset(dataset, mart = db)), silent = TRUE)
      if (!inherits(primary, "try-error")) {
        return(primary)
      }
    }
  }

  if (identical(mart, "ensembl")) {
    fallback <- try(with_biomart_timeout(useEnsembl(biomart = "genes", dataset = dataset)), silent = TRUE)
    if (!inherits(fallback, "try-error")) {
      return(fallback)
    }
  }
  primary
}

direct_homolog_map <- function(source_mart, source_ids, source_filter, source_output_col, ortholog_species,
                               mart_name = "ensembl", dataset = NULL, host = NULL) {
  ortho_dataset <- ortholog_dataset_name(ortholog_species)
  ortho_org <- ortholog_orgdb(ortholog_species)
  ortho_key <- ortholog_keytype(ortholog_species)
  if (is.null(ortho_dataset) || is.null(ortho_org)) {
    return(NULL)
  }

  attr_names <- try(with_biomart_timeout(listAttributes(source_mart)$name), silent = TRUE)
  if (inherits(attr_names, "try-error")) {
    return(NULL)
  }

  prefix <- ortholog_dataset_prefix(ortho_dataset)
  homolog_gene_attr <- pick_first_match(grep(paste0("^", prefix, "_homolog_ensembl_gene$"), attr_names, value = TRUE))
  homolog_symbol_attr <- pick_first_match(grep(paste0("^", prefix, "_homolog_(associated_gene_name|gene_name|external_gene_name)$"), attr_names, value = TRUE))

  query_attrs <- unique(c(source_filter, homolog_symbol_attr, homolog_gene_attr))
  query_attrs <- query_attrs[!is.na(query_attrs) & nzchar(query_attrs)]
  if (length(query_attrs) <= 1) {
    return(NULL)
  }

  if (!is.null(dataset) && !is.null(host) && biomart_prefers_martservice(mart_name)) {
    homolog_map <- try(direct_martservice_query(host = host,
                                                mart = mart_name,
                                                dataset = dataset,
                                                attributes = query_attrs,
                                                filter_name = source_filter,
                                                values = source_ids), silent = TRUE)
  } else {
    homolog_map <- try(with_biomart_timeout(getBM(attributes = query_attrs,
                                                  filters = source_filter,
                                                  values = source_ids,
                                                  mart = source_mart)), silent = TRUE)
    if ((inherits(homolog_map, "try-error") || nrow(homolog_map) == 0) &&
        !is.null(dataset) && !is.null(host) && mart_name != "ensembl") {
      homolog_map <- try(direct_martservice_query(host = host,
                                                  mart = mart_name,
                                                  dataset = dataset,
                                                  attributes = query_attrs,
                                                  filter_name = source_filter,
                                                  values = source_ids), silent = TRUE)
    }
  }
  if (inherits(homolog_map, "try-error") || nrow(homolog_map) == 0) {
    return(NULL)
  }

  homolog_map <- homolog_map %>%
    distinct(.data[[source_filter]], .keep_all = TRUE)
  colnames(homolog_map)[colnames(homolog_map) == source_filter] <- source_output_col

  if (!is.na(homolog_symbol_attr) && homolog_symbol_attr %in% colnames(homolog_map)) {
    colnames(homolog_map)[colnames(homolog_map) == homolog_symbol_attr] <- "SYMBOL"
  } else {
    homolog_map$SYMBOL <- NA_character_
  }

  if (!is.na(homolog_gene_attr) && homolog_gene_attr %in% colnames(homolog_map)) {
    colnames(homolog_map)[colnames(homolog_map) == homolog_gene_attr] <- "ORTHOLOG_ENSEMBL"
    ortho_keys <- unique(homolog_map$ORTHOLOG_ENSEMBL)
    ortho_keys <- ortho_keys[!is.na(ortho_keys) & nzchar(ortho_keys)]
    if (length(ortho_keys) > 0) {
      entrez_map <- AnnotationDbi::select(ortho_org,
                                          keys = ortho_keys,
                                          keytype = ortho_key,
                                          columns = c(ortho_key, "ENTREZID"))
      if (!is.null(entrez_map) && nrow(entrez_map) > 0) {
        colnames(entrez_map)[1] <- "ORTHOLOG_ENSEMBL"
        entrez_map <- entrez_map %>%
          distinct(ORTHOLOG_ENSEMBL, .keep_all = TRUE)
        homolog_map <- merge(homolog_map, entrez_map,
                             by = "ORTHOLOG_ENSEMBL", all.x = TRUE)
      }
    }
  }

  if (!"ENTREZID" %in% colnames(homolog_map)) {
    homolog_map$ENTREZID <- NA_character_
  }

  homolog_map %>%
    dplyr::select(all_of(source_output_col), SYMBOL, ENTREZID) %>%
    distinct(.data[[source_output_col]], .keep_all = TRUE)
}

resolve_non_model_species <- function(species, bio_data) {
  species <- normalize_species_name(species)
  candidates <- as.character(bio_data$Scientific_common_name)
  if (species %in% candidates) species else species
}

no_org_ID <- function(count=NULL,gene_list=NULL,Species,Ortholog,Biomart_archive,RNA_type="gene_level"){
  Species <- normalize_species_input(Species)
  Ortholog <- normalize_ortholog_input(Ortholog)
  if(Species != "not selected"){
    Species <- normalize_species_name(Species)
    if(!is.null(Ortholog)){
      if(sum(is.element(no_orgDb, Species)) == 1){
        org_df <- withProgress(message = "preparing a gene annotation. It takes a few minutes.",{
          if(RNA_type == "gene_level"){
            ENS <- "ensembl_gene_id"
            ENT1 <- "external_gene_name"
            ENT2 <- "external_gene_name"
            ENT3 <- "external_gene_name"
          }else{
            ENS <- "ensembl_transcript_id"
            ENT1 <- c("refseq_mrna","refseq_ncrna")
            ENT2 <- "refseq"
            ENT3 <- c("refseq_mrna","refseq_ncrna","external_gene_name")
          }

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

          species_name <- resolve_non_model_species(Species, bio_data)
          ensembl <- dplyr::filter(bio_data, Scientific_common_name == species_name)$dataset
          use <- open_biomart_dataset(mart, ensembl, Biomart_archive)
          if(inherits(use, "try-error")) {
            validate("biomart has encountered an unexpected server error.
                    Please try using a different 'biomart host' or try again later.")
          }
          source_ids <- if (is.null(gene_list)) rownames(count) else gene_list[[1]]
          source_ids <- unique(as.character(source_ids))
          source_ids <- source_ids[!is.na(source_ids) & nzchar(source_ids)]
          validate(need(length(source_ids) != 0, "No valid gene IDs were found in the input."))

          if(RNA_type == "gene_level"){
            if (mart == "ensembl" || !biomart_prefers_martservice(mart)) {
              genes_ensembl_match <- try(with_biomart_timeout(getBM(attributes = c(ENS,ENT1),
                                                                    filters = ENS,
                                                                    values = source_ids,
                                                                    mart = use)), silent = TRUE)
              genes_symbol_match <- try(with_biomart_timeout(getBM(attributes = c(ENS,ENT1),
                                                                   filters = ENT1,
                                                                   values = source_ids,
                                                                   mart = use)), silent = TRUE)
              if (mart != "ensembl") {
                if (inherits(genes_ensembl_match, "try-error") || nrow(genes_ensembl_match) == 0) {
                  genes_ensembl_match <- try(direct_martservice_query(host = Biomart_archive,
                                                                      mart = mart,
                                                                      dataset = ensembl,
                                                                      attributes = c(ENS, ENT1),
                                                                      filter_name = ENS,
                                                                      values = source_ids), silent = TRUE)
                }
                if (inherits(genes_symbol_match, "try-error") || nrow(genes_symbol_match) == 0) {
                  genes_symbol_match <- try(direct_martservice_query(host = Biomart_archive,
                                                                     mart = mart,
                                                                     dataset = ensembl,
                                                                     attributes = c(ENS, ENT1),
                                                                     filter_name = ENT1,
                                                                     values = source_ids), silent = TRUE)
                }
              }
            } else {
              genes_ensembl_match <- try(direct_martservice_query(host = Biomart_archive,
                                                                  mart = mart,
                                                                  dataset = ensembl,
                                                                  attributes = c(ENS, ENT1),
                                                                  filter_name = ENS,
                                                                  values = source_ids), silent = TRUE)
              genes_symbol_match <- try(direct_martservice_query(host = Biomart_archive,
                                                                 mart = mart,
                                                                 dataset = ensembl,
                                                                 attributes = c(ENS, ENT1),
                                                                 filter_name = ENT1,
                                                                 values = source_ids), silent = TRUE)
              if (inherits(genes_ensembl_match, "try-error") || nrow(genes_ensembl_match) == 0) {
                genes_ensembl_match <- try(with_biomart_timeout(getBM(attributes = c(ENS,ENT1),
                                                                      filters = ENS,
                                                                      values = source_ids,
                                                                      mart = use)), silent = TRUE)
              }
              if (inherits(genes_symbol_match, "try-error") || nrow(genes_symbol_match) == 0) {
                genes_symbol_match <- try(with_biomart_timeout(getBM(attributes = c(ENS,ENT1),
                                                                     filters = ENT1,
                                                                     values = source_ids,
                                                                     mart = use)), silent = TRUE)
              }
            }
          }else{
            genes_ensembl <- try(with_biomart_timeout(getBM(attributes = c(ENS,ENT1), mart = use)), silent = TRUE)
          }

          ortho <- ortholog_dataset_name(Ortholog)
          if(RNA_type == "gene_level"){
            ENSEMBL_n <- if(inherits(genes_ensembl_match, "try-error")) 0 else nrow(genes_ensembl_match)
            SYMBOL_n <- if(inherits(genes_symbol_match, "try-error")) 0 else nrow(genes_symbol_match)
            if(ENSEMBL_n == 0 && SYMBOL_n == 0) validate("Cannot identify gene IDs. Please check the 'Species' and use the 'Official gene symbol' or 'ENSEMBL ID' for gene names.")
            if(ENSEMBL_n >= SYMBOL_n) type <- "ENSEMBL" else type <- "SYMBOL"
          }else if(inherits(genes_ensembl, "try-error")) {
            type <- "ENSEMBL"
          }else{
            if(RNA_type != "gene_level"){
              genes_ensembl <- dplyr::mutate(genes_ensembl, !!ENT2 := gsub("NA","",paste0(genes_ensembl$refseq_mrna,genes_ensembl$refseq_ncrna)))
            }
            if(is.null(gene_list)){
              count <- count %>% dplyr::mutate(!!ENS := rownames(count))
              count <- count %>% dplyr::mutate(!!ENT2 := rownames(count))
              ENSEMBL <- merge(genes_ensembl,count,by=ENS)
              SYMBOL <- merge(genes_ensembl,count,by=ENT2)
            }else{
              gene_list <- gene_list %>% dplyr::mutate(!!ENS := gene_list[,1])
              gene_list <- gene_list %>% dplyr::mutate(!!ENT2 := gene_list[,1])
              ENSEMBL <- merge(genes_ensembl,gene_list,by=ENS)
              SYMBOL <- merge(genes_ensembl,gene_list,by=ENT2)
            }
            if(dim(ENSEMBL)[1] > dim(SYMBOL)[1]) type <- "ENSEMBL" else type <- "SYMBOL"
            if(dim(ENSEMBL)[1] == 0 && dim(SYMBOL)[1] == 0) validate("Cannot identify gene IDs. Please check the 'Species' and use the 'Official gene symbol' or 'ENSEMBL ID' for gene names.")
          }
          if(type == "ENSEMBL"){
            if(RNA_type == "gene_level"){
          colname <- c("ENSEMBL","SYMBOL","ENTREZID")
          colname1<- c("ENSEMBL")
          }else{
            colname <- c("ENSEMBL","refseq_mrna","refseq_ncrna","SYMBOL","ENTREZID")
            colname1<- c("ENSEMBL")
          }
          genes <- data.frame(source_id = source_ids, stringsAsFactors = FALSE)
          colnames(genes)[1] <- ENS
          filter <- ENS
        }else{
            colname <- c("Original_symbol","SYMBOL","ENTREZID")
            colname1<- c("Original_symbol")
            genes <- data.frame(source_id = source_ids, stringsAsFactors = FALSE)
            colnames(genes)[1] <- ENT2
            filter <- ENT1
        }

          genes2_direct <- NULL
          if (RNA_type == "gene_level" && (mart != "ensembl" && Ortholog == "Arabidopsis thaliana" || mart == "ensembl")) {
            genes2_direct <- direct_homolog_map(source_mart = use,
                                                source_ids = source_ids,
                                                source_filter = filter,
                                                source_output_col = colname1[[1]],
                                                ortholog_species = Ortholog,
                                                mart_name = mart,
                                                dataset = ensembl,
                                                host = Biomart_archive)
          }
          if (!is.null(genes2_direct) && nrow(genes2_direct) != 0) {
            genes2 <- genes2_direct
          } else if (mart != "ensembl" && Ortholog == "Arabidopsis thaliana") {
            validate("biomart has encountered an unexpected server error.
                    Please try using a different 'biomart host' or try again later.")
          } else {
            ortho_mart <- open_biomart_dataset(mart, ortho, Biomart_archive)
            if(inherits(ortho_mart, "try-error")) {
              validate("biomart has encountered an unexpected server error.
                    Please try using a different 'biomart host' or try again later.")
            }
            genes2 <- try(with_biomart_timeout(getLDS(attributes = c(ENS),
                                                      values = source_ids, mart = use, filters = filter,
                                                      attributesL = c(ENT3,"entrezgene_id"),
                                                      martL = ortho_mart, uniqueRows = TRUE)), silent = TRUE)
            if(!inherits(genes2, "try-error")) {
              if(RNA_type != "gene_level"){
                colnames(genes2) <- colname
                genes2 <- dplyr::mutate(genes2, !!ENT2 := gsub("NA","",paste0(genes2$refseq_mrna,genes2$refseq_ncrna)))
                genes2 <- genes2 %>% dplyr::select(all_of(colname[1]), refseq, "SYMBOL","ENTREZID")
                colname <- colnames(genes2)
              }
            }else{
              genes4 <- try(with_biomart_timeout(getLDS(attributes = c(ENS),
                                                        values = source_ids, mart = use, filters = filter,
                                                        attributesL = c(ENT3,ENS),
                                                        martL = ortho_mart, uniqueRows = TRUE)), silent = TRUE)
              if(inherits(genes4, "try-error")) {
                genes2 <- direct_homolog_map(source_mart = use,
                                             source_ids = source_ids,
                                             source_filter = filter,
                                             source_output_col = colname1[[1]],
                                             ortholog_species = Ortholog,
                                             mart_name = mart,
                                             dataset = ensembl,
                                             host = Biomart_archive)
                if (is.null(genes2) || nrow(genes2) == 0) {
                  validate("biomart has encountered an unexpected server error.
                    Please try using a different 'biomart host' or try again later.")
                }
              }else{
              if(Ortholog == "Drosophila melanogaster") org <- org.Dm.eg.db
              if(Ortholog == "Caenorhabditis elegans") org <- org.Ce.eg.db
              my.symbols <- genes4[,3]
              if(RNA_type == "gene_level"){
                ENS_dc <- "ENSEMBL"
                name_dc <- "Gene.stable.ID.1"
              }else{
                ENS_dc <- "ENSEMBLTRANS"
                name_dc <- "Transcript.stable.ID.1"
              }
              gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                              keytype = ENS_dc,
                                              columns = c(ENS_dc,"ENTREZID"))
              colnames(gene_IDs) <- c(name_dc,"ENTREZID")
              gene_IDs <- gene_IDs %>% distinct(!!as.name(name_dc), .keep_all = TRUE)
              genes3 <- merge(genes4, gene_IDs, by=name_dc)
              genes2 <- genes3[,-1]
              genes2 <- dplyr::mutate(genes2, !!ENT2 := gsub("NA","",paste0(genes2$RefSeq.mRNA.ID,genes2$RefSeq.ncRNA.ID)))
              genes2 <- genes2 %>% dplyr::select(all_of(colname[1]), refseq, "Gene.name","ENTREZID")
              colname <- c(colname[1],"refseq","SYMBOL","ENTREZID")
              colname <- colnames(genes2)
            }
          }
          }
          colnames(genes) <- colname1
          colnames(genes2) <- colname
          gene3<-merge(genes,genes2,by=colname1,all=TRUE)
          if(type == "ENSEMBL") gene3 %>% distinct(ENSEMBL, .keep_all = TRUE) else gene3 %>% distinct(Original_symbol, .keep_all = TRUE)
        })
      }
      return(org_df)
    }
  }
}
#isoform----
isoform_ID <- function(count=NULL,gene_list=NULL,Species,Ortholog,Biomart_archive,RNA_type="gene_level"){
  Species <- normalize_species_input(Species)
  Ortholog <- normalize_ortholog_input(Ortholog)
  if(Species != "not selected"){
    if(RNA_type == "transcript_level"){
      withProgress(message = "preparing a gene annotation. It takes a few minutes.",{
        if(sum(is.element(no_orgDb, Species)) == 1){
          gene_id <- no_org_ID(count=count,gene_list=gene_list,Species=Species,Ortholog=Ortholog,
                               Biomart_archive=Biomart_archive,RNA_type=RNA_type)
        }else{
          mart <- "ensembl"
          bio_data <- biomart_data
          if(Species == "Arabidopsis thaliana") {
            mart <- "plants_mart"
            bio_data <- biomart_plants
          }
          if(Species == "Saccharomyces cerevisiae") {
            mart <- "fungi_mart"
            bio_data <- biomart_fungi
          }
          if(Species == "Drosophila melanogaster" || Species == "Caenorhabditis elegans") {
            mart <- "metazoa_mart"
            bio_data <- biomart_metazoa
          }
          switch (Species,
                  "Homo sapiens" = ensembl <- "hsapiens_gene_ensembl",
                  "Mus musculus" = ensembl <- "mmusculus_gene_ensembl",
                  "Rattus norvegicus" = ensembl <- "rnorvegicus_gene_ensembl",
                  "Xenopus tropicalis" = ensembl <- "xtropicalis_gene_ensembl",
                  "Drosophila melanogaster" = ensembl <- "dmelanogaster_gene_ensembl",
                  "Caenorhabditis elegans" = ensembl <- "celegans_gene_ensembl",
                  "Anolis carolinensis" = ensembl <- "acarolinensis_gene_ensembl",
                  "Bos taurus" = ensembl <- "btaurus_gene_ensembl",
                  "Canis lupus familiaris" = ensembl <- "clfamiliaris_gene_ensembl",
                  "Danio rerio" = ensembl <- "drerio_gene_ensembl",
                  "Equus caballus" = ensembl <- "ecaballus_gene_ensembl",
                  "Felis catus" = ensembl <- "fcatus_gene_ensembl",
                  "Gallus gallus" = ensembl <- "ggallus_gene_ensembl",
                  "Macaca mulatta" = ensembl <- "mmulatta_gene_ensembl",
                  "Monodelphis domestica" = ensembl <- "mdomestica_gene_ensembl",
                  "Ornithorhynchus anatinus" = ensembl <- "oanatinus_gene_ensembl",
                  "Pan troglodytes" = ensembl <- "ptroglodytes_gene_ensembl")
          use <- open_biomart_dataset(mart, ensembl, Biomart_archive)
          if(inherits(use, "try-error")) {
            validate("biomart has encountered an unexpected server error. Please try using a different 'biomart host' or try again later.")
          }
          genes_ensembl <- try(with_biomart_timeout(getBM(attributes = c("ensembl_transcript_id","external_gene_name","entrezgene_id"), mart = use)), silent = TRUE)
          if (!inherits(genes_ensembl, "try-error")) {
            genes_ensembl <- genes_ensembl %>%
              distinct(ensembl_transcript_id, .keep_all = TRUE) %>%
              dplyr::filter(!is.na(ensembl_transcript_id))
          } else {
            genes_ensembl <- NULL
          }
          genes_refseq <- try(with_biomart_timeout(getBM(attributes = c("refseq_mrna","refseq_ncrna","external_gene_name","entrezgene_id"), mart = use)), silent = TRUE)
          if(inherits(genes_refseq, "try-error")) {
            if (is.null(genes_ensembl)) {
              validate("biomart has encountered an unexpected server error. Please try using a different 'biomart host' or try again later.")
            }
            type <- "ENSEMBL"
          }else{
            genes_refseq <- dplyr::mutate(genes_refseq, refseq = gsub("NA","",paste0(genes_refseq$refseq_mrna,genes_refseq$refseq_ncrna))) %>%
              dplyr::filter(!is.na(refseq)) %>%
              distinct(refseq, .keep_all = TRUE) %>%
              dplyr::select(refseq,external_gene_name,entrezgene_id)
            if(is.null(gene_list)){
              count <- count %>% dplyr::mutate(ensembl_transcript_id = rownames(count))
              count <- count %>% dplyr::mutate(refseq = rownames(count))
              table <- count
            }else{
              gene_list <- gene_list %>% dplyr::mutate(ensembl_transcript_id = gene_list[,1])
              gene_list <- gene_list %>% dplyr::mutate(refseq = gene_list[,1])
              table <- gene_list
            }
            ENSEMBL <- if (is.null(genes_ensembl)) table[0, , drop = FALSE] else merge(genes_ensembl,table,by="ensembl_transcript_id")
            REFSEQ <- merge(genes_refseq,table,by="refseq")
            if(is.null(genes_ensembl) || dim(ENSEMBL)[1] <= dim(REFSEQ)[1]) type <- "REFSEQ" else type <- "ENSEMBL"
            if(dim(ENSEMBL)[1] == 0 && dim(REFSEQ)[1] == 0) validate("Cannot identify gene IDs. Please check the 'Species' and use the 'Official gene symbol' or 'ENSEMBL ID' for gene names.")
          }
          table$check <- TRUE
          if(type == "ENSEMBL"){
            validate(need(!is.null(genes_ensembl), "biomart failed to retrieve transcript annotation. Please try a different 'biomart host' or try again later."))
            gene_id <- merge(genes_ensembl,table,by="ensembl_transcript_id",all=TRUE) %>%
              distinct(ensembl_transcript_id, .keep_all = TRUE) %>%
              dplyr::filter(check==TRUE) %>%
              dplyr::select(ensembl_transcript_id, external_gene_name,entrezgene_id)
          }else{
            gene_id <- merge(genes_refseq,table,by="refseq",all=TRUE) %>%
              distinct(refseq, .keep_all = TRUE) %>%
              dplyr::filter(check==TRUE) %>%
              dplyr::select(refseq, external_gene_name,entrezgene_id)
          }
          colnames(gene_id) <- c("Transcript_ID","SYMBOL","ENTREZID")
        }
      })
      return(gene_id)
    }
  }
}

isoform_balance <- function(data,norm_count){
  filtered_count <- norm_count %>% 
    dplyr::filter(!is.na(SYMBOL)) %>% 
    dplyr::filter(SYMBOL != "") %>% 
    dplyr::filter(SYMBOL %in% data$SYMBOL)
  collist <- gsub("\\_.+$", "", colnames(data))
  collist <- collist[-which(collist == "SYMBOL")]
  for(i in 1:length(unique(collist))){
    cond1 <- unique(collist)[i]
    cond1_ave <- filtered_count %>% dplyr::select(starts_with(cond1)) %>% apply(1,mean)
    table2 <- data.frame(cond=cond1_ave)
    if(i != 1) table <- cbind(table,table2) else table <- table2
  } 
  table <- data.frame(table,SYMBOL=filtered_count$SYMBOL,row.names = rownames(filtered_count))
  res <- data.frame(matrix(rep(NA, 13), nrow=1))[numeric(0), ]
  for(name in unique(table$SYMBOL)){
    table2 <- table %>% dplyr::filter(SYMBOL==name) %>% dplyr::select(-SYMBOL)
    if(dim(table2)[1] != 1){
      cond_sumFC <- log2((colSums(table2)[1]+0.01)/(colSums(table2)[2]+0.01))
      chisq <- chisq.test(table2)
      res2 <- data.frame(SYMBOL = name, 
                         statistic=chisq$statistic,
                         method="Chi-square test",
                         pvalue=chisq$p.value,
                         transcript_n = dim(table2)[1],
                         GeneLevel_FC = abs(cond_sumFC),
                         row.names = NULL)
      res <- rbind(res,res2)
    }
  }
  res$padj<- p.adjust(res$pvalue,method = "BH")
  res <- res %>% dplyr::arrange(padj)
  df <- list()
  df[["res"]] <- res
  df[["table"]] <- table
  return(df)
}
barplot_forTranscript <- function(table,name){
  table2 <- table %>% dplyr::filter(SYMBOL %in% name) %>% 
    rownames_to_column(var = "Row.names") %>% 
    dplyr::group_by(SYMBOL) %>% dplyr::mutate(n=paste0("No.",row_number())) %>%
    dplyr::ungroup(SYMBOL)
  table3 <- table2  %>% dplyr::select(-SYMBOL,-n) %>% 
    gather(key=sample, value=value,-Row.names)
  SYMBOL <- data.frame(SYMBOL = table2$SYMBOL,Row.names = table2$Row.names,n=table2$n)
  table3 <- merge(SYMBOL,table3,by="Row.names")
  g <- ggplot(table3, aes(x = sample, y = value, fill = n)) + 
    geom_bar(position="stack",stat = "identity") + theme_bw()
  p <- facet(g,facet.by = "SYMBOL",scales = "free",panel.labs.background = list(fill = "transparent", color = "transparent"),
             short.panel.labs = T, panel.labs.font = list(size=15, face = "italic"))+ 
    theme(axis.text.x = element_blank(),legend.position = "none",
          panel.background = element_rect(fill = "transparent", size = 0.5),
          title = element_text(size = 10),text = element_text(size = 12),
          axis.title.y = element_text(size=15),legend.text = element_text(size=15),
          legend.title = element_blank())+ 
    xlab(NULL) + ylab("Normalized count")
  return(p)
}
#-----------
org <- function(Species, Ortholog=NULL){
  Species <- normalize_species_input(Species)
  Ortholog <- normalize_ortholog_input(Ortholog)
  if(Species == "not selected") return(NULL)
  target_species <- if(sum(is.element(no_orgDb, Species)) == 0) Species else Ortholog
  package_name <- unname(orgdb_package_by_species[target_species])
  if(is.na(package_name) || !nzchar(package_name)) return(NULL)
  load_orgdb_package(package_name)
}
org_code <- function(Species,Ortholog){
  Species <- normalize_species_input(Species)
  Ortholog <- normalize_ortholog_input(Ortholog)
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
    X <- as.matrix(t(data))
    X <- X[, apply(X, 2, sd, na.rm = TRUE) > 0, drop = FALSE]  
    pca <- prcomp(X, scale. = T)
    label<- colnames(data)
    label<- gsub("\\_.+$", "", label)
    lab_x <- paste(summary(pca)$importance[2,1]*100,
                   "% of variance)", sep = "")
    lab_x <- paste("PC1 (", lab_x, sep = "")
    lab_y <- paste(summary(pca)$importance[2,2]*100,
                   "% of variance)", sep = "")
    lab_y <- paste("PC2 (", lab_y, sep = "")
    pca$x <- as.data.frame(pca$x)
    return(pca$x)
  } 
}
PCAplot <- function(data,legend=NULL){
  if(length(grep("SYMBOL", colnames(data))) != 0){
    data <- data[, - which(colnames(data) == "SYMBOL")]
  }
  if(length(grep("Unique_ID", colnames(data))) != 0){
    data <- data[, - which(colnames(data) == "Unique_ID")]
  }
  X <- as.matrix(t(data))
  X <- X[, apply(X, 2, sd, na.rm = TRUE) > 0, drop = FALSE]  
  pca <- prcomp(X, scale. = T)
  label<- colnames(data)
  lab_x <- paste(summary(pca)$importance[2,1]*100,
                 "% of variance)", sep = "")
  lab_x <- paste("PC1 (", lab_x, sep = "")
  lab_y <- paste(summary(pca)$importance[2,2]*100,
                 "% of variance)", sep = "")
  lab_y <- paste("PC2 (", lab_y, sep = "")
  pca$x <- as.data.frame(pca$x)
  legend <- if(is.null(legend) || !length(legend)){
    NULL
  }else{
    legend_value <- as.character(legend[[1]])
    if(is.na(legend_value) || !nzchar(legend_value)) NULL else legend_value
  }
  legend_position <- "none"
  label2 <- NULL
  if(!is.null(legend)) {
    if(identical(legend, "Legend")){
      legend_position <- "top" 
      label2 <- NULL
    }else{
      legend_position <- "none"
      label2 <- label
    }
  }
  g1 <- ggplot(pca$x,aes(x=pca$x[,1],
                                y=pca$x[,2],
                                col=gsub("\\_.+$", "", label), label = label2)) +
    geom_point()+
    theme(panel.background =element_rect(fill=NA,color=NA),
          panel.border = element_rect(fill = NA)) +
    xlab(lab_x) + ylab(lab_y) +
    theme(legend.position=legend_position, aspect.ratio=1)+ 
    guides(color=guide_legend(title=""))
  if(!is.null(legend)){
    if(identical(legend, "Label")) g1 <- g1 + geom_text_repel(show.legend = NULL)
  }
  rho <- cor(data,method="spearman")
  d <- as.dist(1-rho)
  mds <- try(as.data.frame(cmdscale(d)))
  if(!inherits(mds, "try-error")){
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
    if(identical(legend, "Label")) g2 <- g2 + geom_text_repel(show.legend = NULL)
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
umap_neighbor_info <- function(data, default = 15L){
  item_count <- if(is.null(data)) 0L else ncol(data)
  if(is.na(item_count) || item_count < 3L){
    return(list(valid = FALSE, value = NA_integer_, max = NA_integer_,
                message = "umap: at least 3 samples are required"))
  }
  max_neighbors <- item_count - 1L
  value <- min(as.integer(default), max_neighbors)
  value <- max(2L, value)
  list(valid = TRUE, value = value, max = max_neighbors, message = NULL)
}

umap_error_message <- function(data, n_neighbors = NULL){
  default_neighbors <- if(is.null(n_neighbors)) 15L else as.integer(n_neighbors)
  info <- umap_neighbor_info(data, default = default_neighbors)
  if(!info$valid){
    return(info$message)
  }
  if(is.null(n_neighbors)){
    return(NULL)
  }
  if(is.na(as.integer(n_neighbors)) || as.integer(n_neighbors) > info$max){
    return("umap: number of neighbors must be smaller than number of items")
  }
  NULL
}

umap_plot <- function(data, n_neighbors,lab=NULL){
  err <- umap_error_message(data, n_neighbors)
  if(!is.null(err)){
    stop(err, call. = FALSE)
  }
  info <- umap_neighbor_info(data, default = n_neighbors)
  n_neighbors <- info$value
  umap <- umap::umap(t(data),n_neighbors = n_neighbors, random_state = 123)
  data2 <- umap$layout %>% as.data.frame()
  label <- colnames(data)
  label2 <- NULL
  if(!is.null(lab) && identical(lab, "Label")){
    label2 <- label
  }
  p <- ggplot(data2, mapping = aes(V1,V2, color = gsub("\\_.+$", "", label), label = label2))+
    geom_point()+xlab("UMAP_1") + ylab("UMAP_2")+
    theme(panel.background =element_rect(fill=NA,color=NA),panel.border = element_rect(fill = NA),
          aspect.ratio=1)+ 
    guides(color=guide_legend(title=""))
  if(!is.null(lab) && identical(lab, "Label")){
    p <- p + geom_text_repel(show.legend = NULL)
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

data_3degcount2 <- function(gene_type,data3, Species, Ortholog,Isoform,org){
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
          }else if(gene_type == "isoform"){
            gene_IDs <- Isoform
          }else{
          my.symbols <- data4$Row.names
          key <- gene_primary_key(my.symbols)
          SYMBOL <- sgd_symbol_column(org)
          gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                          keytype = key,
                                          columns = c(key,SYMBOL, "ENTREZID"))
          }
          colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data4 <- merge(data4, gene_IDs, by="Row.names")
        }
      }else{
        if(Species != "not selected"){
          if(sum(is.element(no_orgDb, Species)) == 1){
            gene_IDs <- Ortholog[,-1]
          }else if(gene_type == "isoform"){
            gene_IDs <- Isoform[,-1]
          }else{
          my.symbols <- data4$Row.names
          SYMBOL <- sgd_symbol_column(org)
          gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                          keytype = SYMBOL,
                                          columns = c(SYMBOL, "ENTREZID"))
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
                               fc, fdr, basemean, y_axis=NULL, x_axis=NULL,id_cut,
                               GOI=NULL, heatmap = TRUE, Species,brush=NULL,
                               GOI_color_type="default",cond3_pathway_color_gene=NULL){
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
    if(GOI_color_type == "default"){
      sig <- rep(3, nrow(result))
      sig[which(result$FDR <= fdr & result$FC_x < log2(1/fc) & result$FC_y < log2(1/fc) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 2
      sig[which(result$FDR <= fdr & result$FC_x > log2(fc) & result$FC_y > log2(fc) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 1
      df <- NULL
    }else{
      sig <- rep(3, nrow(result))
      df <- cond3_pathway_color_gene
      for(name in rownames(df)){
        sig[which(result$Row.names == name)] = 4
        sig[which(result$Row.names == name & result$FDR <= fdr & result$FC_x < log2(1/fc) & result$FC_y < log2(1/fc) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 2
        sig[which(result$Row.names == name & result$FDR <= fdr & result$FC_x > log2(fc) & result$FC_y > log2(fc) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 1
      }
    }
    data3 <- data.frame(Row.names = result$Row.names, FC_x = result$FC_x,
                        FC_y = result$FC_y, padj = result$FDR, sig = sig, FC_xy = result$FC_x * result$FC_y)
    high_label <- NULL
    low_label <- NULL
    if(sum(sig == 1) >= 1) high_label <- paste0(specific, "_high: ", sum(sig == 1))
    if(sum(sig == 2) >= 1) low_label <- paste0(specific, "_low: ", sum(sig == 2))
    label_map <- c()
    if(!is.null(high_label)) label_map["1"] <- high_label
    if(!is.null(low_label)) label_map["2"] <- low_label
    label_map["3"] <- "NS"
    if(GOI_color_type != "default") label_map["4"] <- "others"
    color_map <- c()
    if(!is.null(high_label)) color_map[high_label] <- "red"
    if(!is.null(low_label)) color_map[low_label] <- "blue"
    color_map["NS"] <- if (GOI_color_type == "default") "darkgray" else "gray1"
    if(GOI_color_type != "default") color_map["others"] <- "lightgray"
    legend_breaks <- unname(label_map)
    data3$sig <- factor(data3$sig, levels = as.integer(names(label_map)), labels = legend_breaks)
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
      theme_bw()+ scale_color_manual(values = color_map, breaks = legend_breaks)+
      theme(legend.position = "top" , legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 12),
            axis.text.y= ggplot2::element_text(size = 12),
            text = ggplot2::element_text(size = 15),
            title = ggplot2::element_text(size = 15)) +
      xlab(FC_xlab) + ylab(FC_ylab)
    if((sum(sig == 1) >= 1) && (sum(sig == 2) >= 1)){
      if(GOI_color_type != "default") p <- p + geom_point(data=dplyr::filter(data3, sig == label_map[["4"]]),color = "gray1", size= 0.4 )
      p <- p + geom_point(data=dplyr::filter(data3, sig == label_map[["1"]]),color = "red", size= 0.4 )
      p <- p + geom_point(data=dplyr::filter(data3, sig == label_map[["2"]]),color = "blue", size= 0.4 )
    }
    if((sum(sig == 1) >= 1) && (sum(sig == 2) == 0)){
      if(GOI_color_type != "default") p <- p + geom_point(data=dplyr::filter(data3, sig == label_map[["4"]]),color = "gray1", size= 0.4 )
      p <- p + geom_point(data=dplyr::filter(data3, sig == label_map[["1"]]),color = "red", size= 0.4 )
    }
    if((sum(sig == 1) == 0) && (sum(sig == 2) >= 1)){
      if(GOI_color_type != "default") p <- p + geom_point(data=dplyr::filter(data3, sig == label_map[["4"]]),color = "gray1", size= 0.4 )
      p <- p + geom_point(data=dplyr::filter(data3, sig == label_map[["2"]]),color = "blue", size= 0.4 )
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
      if(GOI_color_type == "default"){
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
      }else{
        for(name2 in GOI){
          if(gsub(".+\\s", "", name2) %in% rownames(df)){
            if(gene_type != "SYMBOL"){
              if(Species != "not selected"){
                data3$color[data3$Unique_ID == name2] <- "GOI"
              }else{
                data3$color[data3$Row.names == name2] <- "GOI"
              }
            }else{
              data3$color[data3$Row.names == name2] <- "GOI"
            }
          }
        } 
      }
      if(length(grep("color", colnames(data3))) != 0) {
      if(gene_type != "SYMBOL"){
        if(Species != "not selected"){
          if(id_cut) {
            id_list <- gsub("\\\n.+$", "", data3$Unique_ID)
            dup_list <- unique(id_list[duplicated(id_list)])
            for(i in 1:length(data3$Unique_ID)){
              if(gsub("\\\n.+$", "", data3$Unique_ID[i]) == "NA") {
                data3$Unique_ID[i] <- gsub(".+\\s", "", data3$Unique_ID[i])
              }else if(! gsub("\\\n.+$", "", data3$Unique_ID[i]) %in% dup_list) {
                data3$Unique_ID[i] <- gsub("\\\n.+$", "", data3$Unique_ID[i])
              }
            }
          }
          p <- p + geom_point(data=dplyr::filter(data3, color == "GOI"),color="green", size=1)
          p <- p + ggrepel::geom_text_repel(data = dplyr::filter(data3, color == "GOI"), mapping = aes(label = Unique_ID),
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"),
                                            force = 1, fontface = "bold.italic",
                                            bg.color = grDevices::adjustcolor("white", alpha.f = 0.6), bg.r = 0.15)
        }else{
          p <- p + geom_point(data=dplyr::filter(data3, color == "GOI"),color="green", size=1)
          p <- p + ggrepel::geom_text_repel(data = dplyr::filter(data3, color == "GOI"), mapping = aes(label = Row.names),
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"),
                                            force = 1, fontface = "bold.italic",
                                            bg.color = grDevices::adjustcolor("white", alpha.f = 0.6), bg.r = 0.15)
        }
      }else{
        p <- p + geom_point(data=dplyr::filter(data3, color == "GOI"),color="green", size=1)
        p <- p + ggrepel::geom_text_repel(data = dplyr::filter(data3, color == "GOI"), mapping = aes(label = Row.names),
                                          box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"),
                                          force = 1, fontface = "bold.italic",
                                          bg.color = grDevices::adjustcolor("white", alpha.f = 0.6), bg.r = 0.15)
      }
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
          data <- NULL
          for (name in names(cnet_list2)) {
              em <- cnet_list2[[name]]
              sum <- length(data4$ENTREZID[data4$sig == name])
              if (length(as.data.frame(em)$ID) == 0) {
                cnet1 <- NULL
              } else {
                cnet1 <- as.data.frame(em)
                cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
                cnet1 <- cnet1[sort(cnet1$qvalue, decreasing = F, index=T)$ix,]
                data <- dplyr::bind_rows(data, cnet1)
              }
          }
          if(!is.null(data) && nrow(data) > 0) data$GeneRatio <- DOSE::parse_ratio(data$GeneRatio)
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
              em <- clusterProfiler::enrichKEGG(data4$ENTREZID[data4$sig == name], organism = org_code, pvalueCutoff = 0.05,keyType = "ncbi-geneid")
            }
            if(Gene_set == "GO biological process"){
              em <- clusterProfiler::enrichGO(data4$ENTREZID[data4$sig == name], OrgDb = org, ont = "BP",pvalueCutoff = 0.05)
            }
            if(Gene_set == "GO cellular component"){
              em <- clusterProfiler::enrichGO(data4$ENTREZID[data4$sig == name], OrgDb= org, ont = "CC",pvalueCutoff = 0.05) 
            }
            if(Gene_set == "GO molecular function"){
              em <- clusterProfiler::enrichGO(data4$ENTREZID[data4$sig == name], OrgDb = org, ont = "MF",pvalueCutoff = 0.05) 
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

build_cnetplot <- function(cnet, foldChange = NULL, showCategory = NULL,
                           cex_label_gene = 0.7, cex_label_category = 0.75,
                           cex_category = 0.75, color_edge = "category") {
  args_base <- list(x = cnet)
  if (!is.null(showCategory)) args_base$showCategory <- showCategory
  if (!is.null(foldChange)) args_base$foldChange <- foldChange
  args_core <- args_base
  if (!is.null(color_edge)) args_core$color_edge <- color_edge
  args_base_label <- c(
    args_base,
    list(
      cex_label_gene = cex_label_gene,
      cex_label_category = cex_label_category,
      cex_category = cex_category
    )
  )
  args_label <- c(
    args_core,
    list(
      cex_label_gene = cex_label_gene,
      cex_label_category = cex_label_category,
      cex_category = cex_category
    )
  )

  candidates <- list(
    list(pkg = "enrichplot", args = c(args_base_label, list(color.params = list(edge = color_edge)))),
    list(pkg = "clusterProfiler", args = c(args_base_label, list(color.params = list(edge = color_edge)))),
    list(pkg = "enrichplot", args = c(args_base_label, list(colorEdge = TRUE))),
    list(pkg = "clusterProfiler", args = c(args_base_label, list(colorEdge = TRUE))),
    list(pkg = "enrichplot", args = args_base_label),
    list(pkg = "clusterProfiler", args = args_base_label),
    list(pkg = "enrichplot", args = c(args_base, list(color.params = list(edge = color_edge)))),
    list(pkg = "clusterProfiler", args = c(args_base, list(color.params = list(edge = color_edge)))),
    list(pkg = "enrichplot", args = c(args_base, list(colorEdge = TRUE))),
    list(pkg = "clusterProfiler", args = c(args_base, list(colorEdge = TRUE))),
    list(pkg = "enrichplot", args = args_core),
    list(pkg = "clusterProfiler", args = args_core),
    list(pkg = "enrichplot", args = args_base),
    list(pkg = "clusterProfiler", args = args_base),
    list(pkg = "enrichplot", args = args_label),
    list(pkg = "clusterProfiler", args = args_label),
    list(pkg = "enrichplot", args = args_base),
    list(pkg = "clusterProfiler", args = args_base)
  )

  for (candidate in candidates) {
    if (!requireNamespace(candidate$pkg, quietly = TRUE)) next
    plot_obj <- try(do.call(get("cnetplot", envir = asNamespace(candidate$pkg)), candidate$args), silent = TRUE)
    if (!inherits(plot_obj, "try-error")) {
      return(plot_obj)
    }
  }

  NULL
}

as_cnet_grob <- function(plot_obj, remove_all_legend = FALSE) {
  if (is.null(plot_obj)) return(NULL)

  plot_no_legend <- if (remove_all_legend) {
    try(plot_obj + theme(legend.position = "none"), silent = TRUE)
  } else {
    plot_edge_hidden <- try(plot_obj + guides(edge_colour = "none"), silent = TRUE)
    if (inherits(plot_edge_hidden, "try-error")) {
      plot_edge_hidden <- try(plot_obj + guides(edge_color = "none"), silent = TRUE)
    }
    plot_edge_hidden
  }

  if (!inherits(plot_no_legend, "try-error")) {
    grob_obj <- try(as.grob(plot_no_legend), silent = TRUE)
    if (!inherits(grob_obj, "try-error")) return(grob_obj)
  }

  grob_obj <- try(as.grob(plot_obj), silent = TRUE)
  if (inherits(grob_obj, "try-error")) return(NULL)
  grob_obj
}

keggEnrichment2 <- function(data3, data4,cnet_list2){
  if(!is.null(cnet_list2)){
    data <- NULL
    cnet_list <- list()
    for (name in names(cnet_list2)) {
      sum <- length(data4$ENTREZID[data4$sig == name])
      em <- cnet_list2[[name]]
      if (length(as.data.frame(em)$ID) == 0) {
        cnet1 <- NULL
      } else {
        cnet1 <- em
        cnet_df <- as.data.frame(cnet1)
        if(!"qvalue" %in% colnames(cnet_df) && "p.adjust" %in% colnames(cnet_df)){
          cnet_df$qvalue <- cnet_df$p.adjust
        }
        if ((length(as.data.frame(cnet1)$ID) == 0) || 
            !("qvalue" %in% colnames(cnet_df)) ||
            length(which(!is.na(unique(cnet_df$qvalue))))==0) {
          c <- NULL
        } else{
          c_plot <- build_cnetplot(cnet1,
                                   cex_label_gene = 0.7,
                                   cex_label_category = 0.75,
                                   cex_category = 0.75)
          c <- as_cnet_grob(c_plot, remove_all_legend = FALSE)
          cnet_list[[name]] = c
        }
        cnet1 <- cnet_df
        cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
        cnet1 <- cnet1[sort(cnet1$qvalue, decreasing = F, index=T)$ix,]
        cnet1 <- cnet1[1:5,]
        data <- dplyr::bind_rows(data, cnet1)
      }
    }
    if(is.null(data) || nrow(data) == 0) return(NULL)
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
empty_enrichment_gene_set <- function() {
  data.frame(
    gs_name = character(),
    entrez_gene = integer(),
    gs_id = character(),
    gs_description = character(),
    stringsAsFactors = FALSE
  )
}

msigdbr_subcollection_column <- function(df) {
  if ("gs_subcollection" %in% colnames(df)) {
    return("gs_subcollection")
  }
  if ("gs_subcat" %in% colnames(df)) {
    return("gs_subcat")
  }
  NULL
}

standardize_msigdbr_table <- function(df) {
  if (is.null(df) || !nrow(df)) {
    return(empty_enrichment_gene_set())
  }
  if (!"entrez_gene" %in% colnames(df) && "ncbi_gene" %in% colnames(df)) {
    df$entrez_gene <- df$ncbi_gene
  }
  if (!"gs_id" %in% colnames(df)) {
    df$gs_id <- NA_character_
  }
  if (!"gs_description" %in% colnames(df)) {
    df$gs_description <- NA_character_
  }
  as.data.frame(df)
}

msigdbr_fetch_compat <- function(species, collection, subcollection = NULL) {
  msigdbr_args <- names(formals(msigdbr::msigdbr))
  supports_new_args <- all(c("collection", "subcollection") %in% msigdbr_args)
  if (supports_new_args) {
    if (is.null(subcollection)) {
      res <- try(msigdbr::msigdbr(species = species, collection = collection), silent = TRUE)
    } else {
      res <- try(msigdbr::msigdbr(species = species, collection = collection, subcollection = subcollection), silent = TRUE)
    }
  } else {
    if (is.null(subcollection)) {
      res <- try(msigdbr::msigdbr(species = species, category = collection), silent = TRUE)
    } else {
      res <- try(msigdbr::msigdbr(species = species, category = collection, subcategory = subcollection), silent = TRUE)
    }
  }
  if (inherits(res, "try-error") || is.null(res) || !nrow(res)) {
    return(empty_enrichment_gene_set())
  }
  standardize_msigdbr_table(res)
}

msigdbr_fetch_with_fallbacks <- function(species, collection, subcollections) {
  for (subcollection in subcollections) {
    res <- msigdbr_fetch_compat(species = species, collection = collection, subcollection = subcollection)
    if (nrow(res)) {
      return(res)
    }
  }
  empty_enrichment_gene_set()
}

msigdbr_filter_compat <- function(df, subcollections) {
  if (is.null(df) || !nrow(df)) {
    return(empty_enrichment_gene_set())
  }
  sub_col <- msigdbr_subcollection_column(df)
  if (is.null(sub_col)) {
    return(empty_enrichment_gene_set())
  }
  filtered <- df[df[[sub_col]] %in% subcollections, , drop = FALSE]
  if (!nrow(filtered)) {
    return(empty_enrichment_gene_set())
  }
  filtered
}

GeneList_for_enrichment <- function(Species, Ortholog,Gene_set, org, Custom_gene_list,Biomart_archive,gene_type=NULL){
  if(Species != "not selected" && !is.null(Gene_set) && !is.null(org)){
    if(Species != "Xenopus laevis" && Ortholog != "Arabidopsis thaliana" && Species != "Arabidopsis thaliana"){
    if(Species %in% orgDb_list == TRUE) species <- Species else species <- Ortholog
    ortholog_key <- if (is.null(Ortholog)) "" else Ortholog
    cache_key <- paste(species, ortholog_key, Gene_set, sep = "|")
    if(Gene_set != "Custom gene set" && exists(cache_key, envir = gene_set_cache, inherits = FALSE)) {
      return(get(cache_key, envir = gene_set_cache, inherits = FALSE))
    }
    H_t2g <- NULL
      if(Gene_set == "MSigDB Hallmark"){
        H_t2g <- msigdbr_fetch_compat(species = species, collection = "H") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description) 
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="HALLMARK_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="P53", replacement = "p53")
      }
      if(Gene_set == "KEGG"){
        H_t2g <- msigdbr_fetch_with_fallbacks(
          species = species,
          collection = "C2",
          subcollections = c("CP:KEGG", "CP:KEGG_LEGACY", "CP:KEGG_MEDICUS")
        ) %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="KEGG_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
    if(Gene_set == "Position"){
      H_t2g <- msigdbr_fetch_compat(species = species, collection = "C1") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description) 
    }
      if(Gene_set == "Transcription factor targets"){
        H_t2g <- msigdbr_fetch_compat(species = species, collection = "C3")
        H_t2g <- msigdbr_filter_compat(H_t2g, c("TFT:GTRD", "TFT:TFT_LEGACY")) %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      }
    if(Gene_set == "CGP (chemical and genetic pertubations)"){
      H_t2g <- msigdbr_fetch_compat(species = species, collection = "C2", subcollection = "CGP") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
    }
    if(Gene_set == "ImmuneSigDB"){
      H_t2g <- msigdbr_fetch_compat(species = species, collection = "C7")
      H_t2g <- msigdbr_filter_compat(H_t2g, c("C7", "IMMUNESIGDB")) %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
    }
    if(Gene_set == "Macrophage (444 gene sets from ImmuneSigDB)"){
      H_t2g <- msigdbr_fetch_compat(species = species, collection = "C7")
      H_t2g <- msigdbr_filter_compat(H_t2g, c("C7", "IMMUNESIGDB")) %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description) %>% as.data.frame()
      macrophages <- H_t2g %>% dplyr::filter(grepl(x = gs_name,pattern="MACROPHAGE"))
      BMDMs <- H_t2g %>% dplyr::filter(grepl(x = gs_name,pattern="BMDM"))
      H_t2g <- rbind(macrophages, BMDMs)
    }
    if(Gene_set == "VAX (vaccine response)"){
      H_t2g <- msigdbr_fetch_compat(species = species, collection = "C7")
      H_t2g <- msigdbr_filter_compat(H_t2g, c("C7", "VAX")) %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
    }
    if(Gene_set == "Cell type signature"){
      H_t2g <- msigdbr_fetch_compat(species = species, collection = "C8") %>%
        dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
    }
      if(Gene_set == "Reactome"){
        H_t2g <- msigdbr_fetch_compat(species = species, collection = "C2", subcollection = "CP:REACTOME") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="REACTOME_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
      if(Gene_set == "miRNA target"){
        H_t2g <- msigdbr_fetch_compat(species = species, collection = "C3")
        H_t2g <- msigdbr_filter_compat(H_t2g, c("MIR:MIRDB", "MIR:MIR_LEGACY")) %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
      }
      if(Gene_set == "GO biological process"){
        H_t2g <- msigdbr_fetch_compat(species = species, collection = "C5", subcollection = "GO:BP") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOBP_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
      if(Gene_set == "GO cellular component"){
        H_t2g <- msigdbr_fetch_compat(species = species, collection = "C5", subcollection = "GO:CC") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOCC_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
      if(Gene_set == "GO molecular function"){
        H_t2g <- msigdbr_fetch_compat(species = species, collection = "C5", subcollection = "GO:MF") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="GOMF_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
      if(Gene_set == "Human phenotype ontology"){
        H_t2g <- msigdbr_fetch_compat(species = species, collection = "C5", subcollection = "HPO") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="HP_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
      if(Gene_set == "WikiPathways"){
        H_t2g <- msigdbr_fetch_compat(species = species, collection = "C2", subcollection = "CP:WIKIPATHWAYS") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="WP_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
      if(Gene_set == "BioCarta"){
        H_t2g <- msigdbr_fetch_compat(species = species, collection = "C2", subcollection = "CP:BIOCARTA") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="BIOCARTA_", replacement = "")
        H_t2g$gs_name <- H_t2g$gs_name %>% str_to_lower() %>% str_to_title()
      }
      if(Gene_set == "Custom gene set"){
        if(!is.null(Custom_gene_list)){
          H_t2g <- gene_list_convert_for_enrichment(data= Custom_gene_list, Species = Species,gene_type=gene_type,org=org)
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
        H_t2g <- msigdbr_fetch_compat(species = species, collection = "C2", subcollection = "CP:PID") %>%
          dplyr::select(gs_name, entrez_gene, gs_id, gs_description)
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="PID_", replacement = "")
        H_t2g["gs_name"] <- lapply(H_t2g["gs_name"], gsub, pattern="PATHWAY", replacement = "pathway")
      }
      H_t2g <- H_t2g %>% as.data.frame()
      H_t2g$gs_name <- gsub("_"," ",H_t2g$gs_name)
      if(Gene_set != "Custom gene set") {
        assign(cache_key, H_t2g, envir = gene_set_cache)
      }
    }else{
      print("plant")
      H_t2g <- At_Xl_path(Species,Gene_set)
    }
    print(head(H_t2g))
    return(H_t2g)
  }else return(NULL)
}
GOI_color_palette<-c("default","Set1","Set2","Set3","Paired","Dark2","Accent","Spectral")
GOIboxplot <- function(data,statistical_test=NULL,plottype="Boxplot",ymin=0, ylabel=NULL,
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
        if(inherits(dunnette2, "try-error")){
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
        }else{
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
      if(is.null(ylabel)) ylab = "Normalized_count" else ylab = ylabel
    }else {
      ylim = NULL
      if(is.null(ylabel)) ylab = "ssGSEA score" else ylab = ylabel
    }
    if (plottype == "Boxplot"){
      if(color_design=="new"){
      p <- ggplot(data, aes(x=sample,y=value))+
        geom_boxplot(aes(group=sample,colour=sample,fill=after_scale(alpha(colour,0.5))),outlier.shape = NA)+
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
      p <- ggplot(data, aes(x = sample, y = value)) + geom_boxplot(aes(fill=sample),outlier.shape = NA)+
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
enrich_viewer_forMulti1 <- function(gene_type,df, Species, Ortholog,Isoform, org){
  if(is.null(df) || Species == "not selected"){
    return(NULL)
  }else{
    my.symbols <- df$GeneID
    if(gene_type != "SYMBOL"){
      if(sum(is.element(no_orgDb, Species)) == 1){
        gene_IDs <- Ortholog
      }else if(gene_type == "isoform"){
        gene_IDs <- Isoform
      }else{
      key <- gene_primary_key(my.symbols)
      SYMBOL <- sgd_symbol_column(org)
      gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                      keytype = key,
                                      columns = c(key,SYMBOL, "ENTREZID"))
      }
      colnames(gene_IDs) <- c("GeneID","SYMBOL", "ENTREZID")
    }else{
      if(sum(is.element(no_orgDb, Species)) == 1){
        gene_IDs <- Ortholog[,-1]
      }else if(gene_type == "isoform"){
        gene_IDs <- Isoform[,-1]
      }else{
      SYMBOL <- sgd_symbol_column(org)
      gene_IDs <- AnnotationDbi::select(org, keys = my.symbols,
                                        keytype = SYMBOL,
                                        columns = c("ENTREZID", SYMBOL))
      }
      colnames(gene_IDs) <- c("GeneID","ENTREZID")
    }
    gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
    data <- merge(df, gene_IDs, by="GeneID")
    return(data)
  }
}
enrich_viewer_forMulti2 <- function(gene_type,df, Species,Ortholog,Isoform, Gene_set, org, org_code, H_t2g){
  if(!is.null(Gene_set)){
    data3 <- enrich_viewer_forMulti1(gene_type=gene_type,df = df, Species = Species, org = org,Ortholog=Ortholog,Isoform=Isoform)
    if(is.null(data3)){
      return(NULL)
    }else{
        withProgress(message = "enrichment analysis",{
          if(is.null(H_t2g)){
            df <- NULL
          }else{
            H_t2g2 <- H_t2g %>% dplyr::select(gs_name, entrez_gene)
          df <- NULL
          for (name in unique(data3$Group)) {
            sum <- length(data3$ENTREZID[data3$Group == name])
            em <- clusterProfiler::enricher(data3$ENTREZID[data3$Group == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) != 0) {
              cnet1 <- as.data.frame(clusterProfiler::setReadable(em, org, 'ENTREZID'))
              if(!"qvalue" %in% colnames(cnet1) && "p.adjust" %in% colnames(cnet1)) cnet1$qvalue <- cnet1$p.adjust
              cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
              df <- dplyr::bind_rows(df, cnet1)
            }
          }
          }
          if(!is.null(df) && nrow(df) != 0){
            df$GeneRatio <- DOSE::parse_ratio(df$GeneRatio)
            return(df)
          }else return(NULL)
        })
      }
  } 
}
enrich_viewer_forMulti2_xenopus <- function(df, Species,Ortholog,Isoform, Gene_set, org, org_code){
  if(!is.null(Gene_set)){
    data3 <- enrich_viewer_forMulti1(df = df, Species = Species,Ortholog = Ortholog,Isoform=Isoform, org = org)
    if(is.null(data3)){
      return(NULL)
    }else{
      withProgress(message = "enrichment analysis",{
          df <- NULL
          for (name in unique(data3$Group)) {
            sum <- length(data3$ENTREZID[data3$Group == name])
            if(Gene_set == "KEGG"){
              em <- clusterProfiler::enrichKEGG(data3$ENTREZID[data3$Group == name], organism = org_code, pvalueCutoff = 0.05,keyType = "ncbi-geneid")
            }
            if(Gene_set == "GO biological process"){
              em <- clusterProfiler::enrichGO(data3$ENTREZID[data3$Group == name], OrgDb = org, ont = "BP",pvalueCutoff = 0.05)
            }
            if(Gene_set == "GO cellular component"){
              em <- clusterProfiler::enrichGO(data3$ENTREZID[data3$Group == name], OrgDb= org, ont = "CC",pvalueCutoff = 0.05) 
            }
            if(Gene_set == "GO molecular function"){
              em <- clusterProfiler::enrichGO(data3$ENTREZID[data3$Group == name], OrgDb = org, ont = "MF",pvalueCutoff = 0.05) 
            }
            if (length(as.data.frame(em)$ID) != 0) {
              cnet1 <- as.data.frame(clusterProfiler::setReadable(em, org, 'ENTREZID'))
              if(!"qvalue" %in% colnames(cnet1) && "p.adjust" %in% colnames(cnet1)) cnet1$qvalue <- cnet1$p.adjust
              cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
              df <- dplyr::bind_rows(df, cnet1)
            }
          }
        if(!is.null(df) && nrow(df) != 0){
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
            formula_res <- dplyr::filter(formula_res, !is.na(qvalue))
            p1 <- as.grob(enrichplot::dotplot(formula_res, color ="qvalue", font.size = 10))
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
          em <- clusterProfiler::enricher(data$ENTREZID[data$Group == name], TERM2GENE=H_t2g2, pvalueCutoff = 0.5)
          if (length(as.data.frame(em)$ID) != 0) {
            cnet1 <- clusterProfiler::setReadable(em, org, 'ENTREZID')
            df[[name]] <- cnet1
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
            em <- clusterProfiler::enrichKEGG(data$ENTREZID[data$Group == name], organism = org_code, pvalueCutoff = 0.05,keyType = "ncbi-geneid")
          }
          if(Gene_set == "GO biological process"){
            em <- clusterProfiler::enrichGO(data$ENTREZID[data$Group == name], OrgDb = org, ont = "BP",pvalueCutoff = 0.05)
          }
          if(Gene_set == "GO cellular component"){
            em <- clusterProfiler::enrichGO(data$ENTREZID[data$Group == name], OrgDb= org, ont = "CC",pvalueCutoff = 0.05) 
          }
          if(Gene_set == "GO molecular function"){
            em <- clusterProfiler::enrichGO(data$ENTREZID[data$Group == name], OrgDb = org, ont = "MF",pvalueCutoff = 0.05) 
          }
          if (length(as.data.frame(em)$ID) != 0) {
            cnet1 <- clusterProfiler::setReadable(em, org, 'ENTREZID')
            df[[name]] <- cnet1
          }
        }
      return(df)
    }
}


enrich_genelist <- function(data, enrich_gene_list, showCategory=5,section=NULL,group_order=NULL){
      if(is.null(data) || is.null(enrich_gene_list)){
        return(NULL)
      }else{
          df <- NULL
          cluster_list <- c()
          for (name in names(enrich_gene_list)) {
            sum <- length(data$ENTREZID[data$Group == name])
            em <- enrich_gene_list[[name]]
            if (length(as.data.frame(em)$ID) != 0) {
              cnet1 <- as.data.frame(em)
              if(!"qvalue" %in% colnames(cnet1) && "p.adjust" %in% colnames(cnet1)) cnet1$qvalue <- cnet1$p.adjust
              cnet1$Group <- paste(name, "\n","(",sum, ")",sep = "")
              cluster_list <- c(cluster_list, paste(name, "\n","(",sum, ")",sep = ""))
              if(!is.null(group_order)) group_order[which(group_order == name)] <- paste(name, "\n","(",sum, ")",sep = "")
              cnet1 <- cnet1[sort(cnet1$pvalue, decreasing = F, index=T)$ix,]
              if (length(cnet1$pvalue) > showCategory){
                cnet1 <- cnet1[1:showCategory,]
              }
              df <- dplyr::bind_rows(df, cnet1)
            }
          }

          if(is.null(df) || nrow(df) == 0) return(NULL)
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
          p2_plot <- build_cnetplot(cnet1,
                                    showCategory = showCategory,
                                    cex_label_gene = 0.7,
                                    cex_label_category = 0.75,
                                    cex_category = 0.75)
          p2 <- as_cnet_grob(p2_plot, remove_all_legend = FALSE)
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

ensembl2symbol <- function(data,Species,Ortholog,Isoform,gene_type,org, merge=TRUE, rowname = TRUE){
  if(Species != "not selected"){
  if(gene_type != "SYMBOL"){
  data <- as.data.frame(data)
  if(sum(is.element(no_orgDb, Species)) == 1){
    gene_IDs <- try(Ortholog[,-3])
    if(inherits(gene_IDs, "try-error")) validate("biomart has encountered an unexpected server error.
                                                \nPlease try using a different 'biomart host' or try again later.")
  }else{
    if(gene_type == "ENSEMBL"){
      my.symbols <- gsub("\\..*","", rownames(data))
      key <- gene_primary_key(my.symbols)
      SYMBOL <- sgd_symbol_column(org)
      gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                      keytype = key,
                                      columns = c(key,SYMBOL))
    }
    if(gene_type == "isoform"){
      gene_IDs <- try(Isoform[,-3])
      if(inherits(gene_IDs, "try-error")) validate("biomart has encountered an unexpected server error.
                                                \nPlease try using a different 'biomart host' or try again later.")
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

gene_type <- function(my.symbols,org,Species,RNA_type="gene_level"){
    Species <- normalize_species_input(Species)
    if(Species != "not selected"){
  if(RNA_type != "gene_level" && is_transcript_like_id(my.symbols)) {
    type <- "isoform"
  } else {
    ENS <- "ENSEMBL"
    ENT <- "SYMBOL"
    if(sum(is.element(no_orgDb, Species)) != 1){
      key <- gene_primary_key(my.symbols)
      SYMBOLa <- sgd_symbol_column(org)
    ENSEMBL<-try(AnnotationDbi::select(org,keys = my.symbols,
                                       keytype = key,
                                       columns = c(key, "ENTREZID")))
    SYMBOL <-try(AnnotationDbi::select(org,keys = my.symbols,
                                       keytype = SYMBOLa,
                                       columns = c(SYMBOLa, "ENTREZID")))
    if(inherits(ENSEMBL, "try-error") && !inherits(SYMBOL, "try-error")) {type <- "SYMBOL"
    }else if(!inherits(ENSEMBL, "try-error") && inherits(SYMBOL, "try-error")) {type <- "ENSEMBL"
    }else if(inherits(ENSEMBL, "try-error") && inherits(SYMBOL, "try-error")) {validate("Cannot identify gene IDs. Please check the 'Species' and use the 'Official gene symbol' or 'ENSEMBL ID' for gene names.")
    }else{
      if(dim(ENSEMBL)[1] > dim(SYMBOL)[1]) type <- "ENSEMBL" else type <- "SYMBOL"
    }
    }else type <- "non-model organism"
  }
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

corr_plot_pair <- function(data,corr_color,GOI_x,GOI_y,logscale=TRUE){
  p1 <- NULL
  p2 <- NULL
  if(corr_color == ""){
    if(logscale) p1 <- ggplot(data, aes(x=log10(.data[[GOI_x]]+1),y=log10(.data[[GOI_y]]+1))) else p1 <- ggplot(data, aes(x=.data[[GOI_x]],y=.data[[GOI_y]]))
    p1 <- p1 + geom_smooth(method=lm, se=FALSE, color='#2C3E50',linetype="dashed",size=0.5)
  }else if(corr_color == "sample_name"){
    label <- gsub("\\_.+$", "", rownames(data))
    if(logscale) p1 <- ggplot(data, aes(x=log10(.data[[GOI_x]]+1),y=log10(.data[[GOI_y]]+1), col=label)) else p1 <- ggplot(data, aes(x=.data[[GOI_x]],y=.data[[GOI_y]], col=label))
    p1 <- p1 + geom_smooth(method=lm, se=FALSE, color='#2C3E50',linetype="dashed",size=0.5)
  }else {
    if(logscale) p1 <- ggplot(data, aes(x=log10(.data[[GOI_x]]+1),y=log10(.data[[GOI_y]]+1), col=log10(.data[[corr_color]]))) else p1 <- ggplot(data, aes(x=.data[[GOI_x]],y=.data[[GOI_y]], col=.data[[corr_color]]))
    p1 <- p1 + geom_smooth(method=lm, se=FALSE, color='#2C3E50',linetype="dashed",size=0.5)
  }
  if(!is.null(p1)){
    p <- p1 +
      geom_point()+ 
      theme_bw()+ 
      theme(legend.position = "top" , legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 12),
            axis.text.y= ggplot2::element_text(size = 12),
            text = ggplot2::element_text(size = 15),
            title = ggplot2::element_text(size = 15),
            plot.title = element_text(size = 15))
    if(logscale) p<- p + 
        xlab(paste(strwrap(paste0("log10(", GOI_x," + 1)"), width = 30),collapse = "\n"))+
        ylab(paste(strwrap(paste0("log10(", GOI_y," + 1)"), width = 30),collapse = "\n"))
  }
  if(!is.null(p2)){
    p <- p2 +
      geom_point()+
      scale_color_continuous(low="blue", high="red")+ 
      theme_bw()+
      theme(legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 12),
            axis.text.y= ggplot2::element_text(size = 12),
            text = ggplot2::element_text(size = 15),
            title = ggplot2::element_text(size = 15),
            plot.title = element_text(size = 15))
    if(logscale) p <- p + ggtitle(paste(strwrap(paste0("color = log10(", corr_color," + 1)"), width = 30),collapse = "\n"))+ 
        xlab(paste(strwrap(paste0("log10(", GOI_x," + 1)"), width = 30),collapse = "\n"))+
        ylab(paste(strwrap(paste0("log10(", GOI_y," + 1)"), width = 30),collapse = "\n"))
  }
  return(p)
}


ssGSEA <- function(norm_count, gene_set,org,gene_type,Species,Ortholog){
  set.seed(12345)
  genesbyGeneSet <- split(gene_set$GeneID,gene_set$gs_name)
  if(length(grep("SYMBOL", colnames(norm_count))) != 0) norm_count <- norm_count[, - which(colnames(norm_count) == "SYMBOL")]
  
  ssgseaPar <- GSVA::ssgseaParam(as.matrix(norm_count),genesbyGeneSet)
  ssgsea.score <- GSVA::gsva(ssgseaPar)
  return(ssgsea.score)
}

geneid_convert_ssGSEA <- function(norm_count,gene_type,Species,Ortholog,Isoform,org){
  my.symbols <- rownames(norm_count)
  if(gene_type != "SYMBOL"){
    if(sum(is.element(no_orgDb, Species)) == 1){
      gene_IDs <- Ortholog
    }else if(gene_type == "isoform"){
      gene_IDs <- Isoform
    }else{
      key <- gene_primary_key(my.symbols)
      SYMBOL <- sgd_symbol_column(org)
      gene_IDs<-AnnotationDbi::select(org,keys = my.symbols,
                                      keytype = SYMBOL,
                                      columns = c(SYMBOL,key, "ENTREZID"))
    }
    colnames(gene_IDs) <- c("SYMBOL","GeneID", "ENTREZID")
  }else{
    if(sum(is.element(no_orgDb, Species)) == 1){
      gene_IDs <- Ortholog[,-1]
    }else if(gene_type == "isoform"){
      gene_IDs <- Isoform
    }else{
      SYMBOL <- sgd_symbol_column(org)
      gene_IDs <- AnnotationDbi::select(org, keys = my.symbols,
                                        keytype = SYMBOL,
                                        columns = c(SYMBOL,"ENTREZID"))
    }
    colnames(gene_IDs) <- c("GeneID","ENTREZID")
  }
  gene_IDs <- na.omit(gene_IDs)
  gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
  rownames(gene_IDs) <- gene_IDs$GeneID

  return(gene_IDs)
}
At_Xl_path <- function(Species,Gene_set){
  if(Gene_set == "KEGG"){
    ##KEGG
    if(Species == "Arabidopsis thaliana") {
      species <- "ath" 
      gsub_name <- ' - Arabidopsis thaliana \\(thale cress\\)'
    }else if(Species == "Xenopus laevis") {
      species <- "xla"
      gsub_name <- ' - Xenopus laevis \\(African clawed frog\\)'
    }else validate("")
    library(KEGGREST)
    print("keggrest")
    pathwayList <- keggList("pathway", species)
    df <- data.frame(pathway = names(pathwayList), gs_name = pathwayList)
    df$gs_name <- gsub(gsub_name,"",df$gs_name)
    kegg <- keggLink(species, "pathway")
    df2 <- data.frame(pathway = names(kegg), entrez_gene = kegg)
    df2$pathway <- gsub("path:", "", df2$pathway)
    df2$entrez_gene <- as.character(gsub(paste0(species,":"), "", df2$entrez_gene))
    pathwayList <- merge(df,df2,by="pathway")
    print(pathwayList)
  }
  
  if(Species == "Arabidopsis thaliana") {
    if(Gene_set == "ARACYC"){
      ##ARACYC
      pathwayList <- toTable(org.At.tairARACYC)
    }else if(str_detect(Gene_set, "GO")){
      pathwayList <- toTable(org.At.tairGO2TAIR)[,1:2]
      tb = AnnotationDbi::toTable(GO.db::GOTERM)[,-1]
      pathwayList <- merge(pathwayList,tb,by="go_id")
      colnames(pathwayList)[2:3] <- c("entrez_gene","gs_name")
    }
  }else if(Species == "Xenopus laevis") {
    ##GO
    if(str_detect(Gene_set, "GO")){
      pathwayList <- toTable(org.Xl.egGO2EG)[,1:2]
      tb = AnnotationDbi::toTable(GO.db::GOTERM)[,-1]
      pathwayList <- merge(pathwayList,tb,by="go_id")
      colnames(pathwayList)[2:3] <- c("entrez_gene","gs_name")
    }
  }else validate("")
  return(pathwayList)
}
