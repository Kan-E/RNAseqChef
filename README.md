[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7095217.svg)](https://doi.org/10.5281/zenodo.7095217)

# RNAseqChef
RNAseqChef, an RNA-seq data controller highlighting gene expression features, is a web-based application ([https://imeg-ku.shinyapps.io/RNAseqChef/](https://imeg-ku.shinyapps.io/RNAseqChef/)) for automated, systematic, and integrated RNA-seq differential expression analysis. RNAseqChef is designed for wet-bench scientists with little computational programming skill to dissect multiple RNA-seq datasets quickly. <br>

# Manual
Manual：           [https://github.com/Kan-E/RNAseqChef/wiki](https://github.com/Kan-E/RNAseqChef/wiki) <br>
Manual (Japanese)：[https://kan-e.github.io/RNAseqChef_manual_japanese/](https://kan-e.github.io/RNAseqChef_manual_japanese/) <br>

# Local installation
## Method 1 (Docker is required)
- Download Docker
- Run the following commands once to get the docker image of RNAseqChef<br>
```
docker pull omicschef/rnaseqchef:v1.0.9
```
You may now run RNAseqChef with just one command in the command line:
```
docker run --rm -p 3838:3838 omicschef/rnaseqchef:v1.0.9
```
Please access http://localhost:3838 in your browser.

## Method 2 (R environment setup is required)
- Download R and RStudio (In the case of macOS, additionally install XQuartz and Xcode)
- Run the following commands once<br>

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
install.packages("https://cran.r-project.org/src/contrib/Archive/lasso2/lasso2_1.2-22.tar.gz",repos = NULL, type = "source")

pkgs <- c("shiny","DT","readxl","rstatix","multcomp","tidyverse","ggpubr","venn","ggrepel",
"ggdendro","ggplotify","gridExtra","cowplot","DESeq2","EBSeq","ggnewscale","edgeR","IHW",
"qvalue","org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","org.Xl.eg.db","org.Dm.eg.db",
"org.Ce.eg.db","AnnotationDbi","clusterProfiler","enrichplot","DOSE","msigdbr","genefilter",
"ComplexHeatmap","shinyBS","plotly","shinyjs","DEGreport","devtools","dorothea","umap", "biomaRt",
"monaLisa","GenomicRanges","BiocParallel","SummarizedExperiment","JASPAR2020",
"TxDb.Mmusculus.UCSC.mm10.knownGene","TxDb.Hsapiens.UCSC.hg19.knownGene",
"BSgenome.Mmusculus.UCSC.mm10","BSgenome.Hsapiens.UCSC.hg19","TFBSTools", "org.Bt.eg.db",
"org.Dr.eg.db","org.Cf.eg.db","org.Gg.eg.db","org.Mmu.eg.db","org.Pt.eg.db","org.Sc.sgd.db",
"org.At.tair.db","colorspace","magick","pdftools","clue")
options(repos = BiocManager::repositories())
for(pkg in pkgs) if (!require(pkg, character.only = T)){
    BiocManager::install(pkg, update = F)
}
devtools::install_github("YuLab-SMU/clusterProfiler.dplyr")
devtools::install_github('VPetukhov/ggrastr')
```

You may now run RNAseqChef with just one command in R:<br>
```
shiny::runGitHub("RNAseqChef", "Kan-E")
```

# Citation
When the user publishes the results from RNAseqChef analysis, please cite our original paper.
- Etoh K. & Nakao M. A web-based integrative transcriptome analysis, RNAseqChef, uncovers cell/tissue type-dependent action of sulforaphane. _JBC_, 299(6), 104810 (2023).
https://doi.org/10.1016/j.jbc.2023.104810

# Reference
Shiny framework
- Winston Chang, Joe Cheng, JJ Allaire, Carson Sievert, Barret Schloerke, Yihui
  Xie, Jeff Allen, Jonathan McPherson, Alan Dipert and Barbara Borges (2021).
  shiny: Web Application Framework for R. R package version 1.7.1.
  https://CRAN.R-project.org/package=shiny
- Eric Bailey (2022). shinyBS: Twitter Bootstrap Components for Shiny. R package
  version 0.61.1. https://CRAN.R-project.org/package=shinyBS
- Yihui Xie, Joe Cheng and Xianying Tan (2022). DT: A Wrapper of the JavaScript
  Library 'DataTables'. R package version 0.23.
  https://CRAN.R-project.org/package=DT

EBSeq (for ebseq)
- Ning Leng and Christina Kendziorski (2020). EBSeq: An R package for gene and isoform
  differential expression analysis of RNA-seq data. R package version 1.30.0.
  
DESeq2 (for deseq2)
- Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for
  RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)

edgeR (for edger)
- Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential
  expression analysis of digital gene expression data. Bioinformatics 26, 139-140

IHW, Independent hypothesis weighting, and qvalue (for fdr control method of deseq2 and edger)
- Nikolaos Ignatiadis, Bernd Klaus, Judith Zaugg and Wolfgang Huber (2016): Data-driven hypothesis
  weighting increases detection power in genome-scale multiple testing. Nature Methods 13:577,
  doi: 10.1038/nmeth.3885
- John D. Storey, Andrew J. Bass, Alan Dabney and David Robinson (2021). qvalue: Q-value
  estimation for false discovery rate control. R package version 2.26.0.
  http://github.com/jdstorey/qvalue

ggdendro (for dendrograms)
- Andrie de Vries and Brian D. Ripley (2020). ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'. R package version 0.1.22. https://CRAN.R-project.org/package=ggdendro

umap (for UMAP)
- Konopka T (2022). _umap: Uniform Manifold Approximation and Projection_. R package
  version 0.2.8.0, <https://CRAN.R-project.org/package=umap>.

clusterProfiler, DOSE, msigdbr, dorothea (for enrichment analysis)
- T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141
- Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015 31(4):608-609
- Dolgalev I (2022). _msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format_. R
  package version 7.5.1, <https://CRAN.R-project.org/package=msigdbr>.
- Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. 'Benchmark and integration of resources for the estimation of human transcription factor activities.' Genome Research. 2019. DOI: 10.1101/gr.240663.118.

AnnotationDbi,AnnotationHub, org.Hs.eg.db, org.Mm.eg.db, org.Rn.eg.db, org.Xl.eg.db, org.Dm.eg.db, org.Ce.eg.db, org.Bt.eg.db, org.Cf.eg.db, org.Dr.eg.db, org.Gg.eg.db, org.Mmu.eg.db, org.Pt.eg.db, org.Sc.sgd.db, and org.Ss.eg.db (for genome wide annotation)
- Hervé Pagès, Marc Carlson, Seth Falcon and Nianhua Li (2020). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.52.0. https://bioconductor.org/packages/AnnotationDbi
- Marc Carlson (2020). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.12.0.
- Marc Carlson (2020). org.Mm.eg.db: Genome wide annotation for Mouse. R package version 3.12.0.
- Marc Carlson (2022). org.Rn.eg.db: Genome wide annotation for Rat. R package version 3.15.0.
- Marc Carlson (2022). org.Xl.eg.db: Genome wide annotation for Xenopus. R package version 3.15.0.
- Marc Carlson (2022). org.Dm.eg.db: Genome wide annotation for Fly. R package version 3.15.0.
- Marc Carlson (2022). org.Ce.eg.db: Genome wide annotation for Worm. R package version 3.15.0.
- Marc Carlson (2022). org.Bt.eg.db: Genome wide annotation for Bovine. R package version 3.15.0.
- Marc Carlson (2022). org.Cf.eg.db: Genome wide annotation for Canine. R package version 3.15.0.
- Marc Carlson (2022). org.Dr.eg.db: Genome wide annotation for Zebrafish. R package version 3.15.0.
- Marc Carlson (2022). org.Gg.eg.db: Genome wide annotation for Chicken. R package version 3.15.0.
- Marc Carlson (2022). org.Mmu.eg.db: Genome wide annotation for Rhesus. R package version 3.15.0.
- Marc Carlson (2022). org.Pt.eg.db: Genome wide annotation for Chimp. R package version 3.15.0.
- Marc Carlson (2022). org.Sc.sgd.db: Genome wide annotation for Yeast. R package version 3.15.0.

biomaRt (for ortholog)
- Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package
  biomaRt. Steffen Durinck, Paul T. Spellman, Ewan Birney and Wolfgang Huber, Nature Protocols 4,
  1184-1191 (2009).

genefilter (for z-score normalization)
- R. Gentleman, V. Carey, W. Huber and F. Hahne (2021). genefilter: methods for filtering genes from high-throughput experiments. R package version 1.72.1.

ComplexHeatmap (for heatmap and k-means clustering)
- Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.

ggplot2 and ggpubr (for boxplot and scater plot)
- H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
- Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.4.0. https://CRAN.R-project.org/package=ggpubr

DEGreport (for divisive clustering analysis)
- Pantano L (2022). DEGreport: Report of DEG analysis. R package version 1.32.0, http://lpantano.github.io/DEGreport

venn (for venn diagram analysis)
- Adrian Dusa (2021). venn: Draw Venn Diagrams. R package version 1.10. https://CRAN.R-project.org/package=venn

GenomicRanges, TxDb.Mmusculus.UCSC.mm10.knownGene, TxDb.Hsapiens.UCSC.hg19.knownGene, BSgenome.Mmusculus.UCSC.mm10, and BSgenome.Hsapiens.UCSC.hg19 (for promoter sequence)
- Lawrence M, Huber W, Pag\`es H, Aboyoun P, Carlson M, et al. (2013) Software for Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118. doi:10.1371/journal.pcbi.1003118"
- Team BC, Maintainer BP (2019). _TxDb.Mmusculus.UCSC.mm10.knownGene: Annotation package for TxDb object(s)_. R package version 3.10.0.
- Team TBD (2021). _BSgenome.Mmusculus.UCSC.mm10: Full genome sequences for Mus musculus (UCSC version mm10, based on GRCm38.p6)_. R package version 1.4.3.
- Carlson M, Maintainer BP (2015). _TxDb.Hsapiens.UCSC.hg19.knownGene: Annotation package for TxDb object(s)_. R package version 3.2.2.
- Team TBD (2020). _BSgenome.Hsapiens.UCSC.hg19: Full genome sequences for Homo sapiens (UCSC version hg19, based on GRCh37.p13)_. R package version 1.4.3.

monaLisa, TFBSTools, BiocParallel, SummarizedExperiment, and JASPAR2020 (for promoter motif analysis)
- Machlab D, Burger L, Soneson C, Rijli FM, Schübeler D, Stadler MB. monaLisa: an R/Bioconductor package for identifying regulatory motifs. Bioinformatics (2022).
- Tan, G., and Lenhard, B. (2016). TFBSTools: an R/bioconductor package for transcription factor binding site analysis. Bioinformatics 32, 1555-1556.
- Morgan M, Wang J, Obenchain V, Lang M, Thompson R, Turaga N (2022). _BiocParallel:Bioconductor facilities for parallel evaluation_. R package version 1.30.3, <https://github.com/Bioconductor/BiocParallel>.
- Morgan M, Obenchain V, Hester J, Pagès H (2022). _SummarizedExperiment:SummarizedExperiment container_. R package version 1.26.1, <https://bioconductor.org/packages/SummarizedExperiment>.
- Baranasic D (2020). _JASPAR2020: Data package for JASPAR database (version 2020)_. R package version 0.99.10, <http://jaspar.genereg.net/>.

dplyr and tidyr (for data manipulation)
- Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.7. https://CRAN.R-project.org/package=dplyr
- Hadley Wickham (2021). tidyr: Tidy Messy Data. R package version 1.1.3. https://CRAN.R-project.org/package=tidyr

# License
This shiny code is licensed under the GPLv3. Please see the file [LICENSE.md](https://github.com/Kan-E/RNAseqChef/blob/main/LICENSE.md) for information.<br>
```
RNAseqChef, an RNA-seq data controller highlighting gene expression features
Shiny App for automated, systematic, and integrated RNA-seq differential expression analysis
Copyright (C) 2022  Kan Etoh

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

You may contact the author of this code, Kan Etoh, at <kaneto@kumamoto-u.ac.jp>
```
# Author

Kan Etoh
<[kaneto@kumamoto-u.ac.jp](mailto:kaneto@kumamoto-u.ac.jp)>
