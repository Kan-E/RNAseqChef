library(dorothea)
library(AnnotationDbi)
library(dplyr)

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
return(net3)
}
