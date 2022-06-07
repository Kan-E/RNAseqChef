# RNAseqChef
RNAseqChef, an RNA-seq data controller highlighting gene expression features, is a web-based application (https://kan-e.shinyapps.io/RNAseqChef/) for automated, systematic, and integrated RNA-seq differential expression analysis. RNAseqChef is designed for wet-bench scientists with little computational programming skill to dissect multiple RNA-seq datasets quickly. <br>

RNAseqChef consists of mainly five panels including, `Pair-wise DEG`, `3 conditions DEG`, `Venn diagram`, `Normalized count data`, and `Enrichment viewer`. `Pair-wise DEG` and `3 conditions DEG` have functions for single dataset DEG analysis such as DEG detection, clustering analysis, gene expression profiling, and enrichment analysis. `Venn diagram`, `normalized count data`, and `Enrichment viewer` have functions for multiple dataset analysis such as Venn diagram analysis, k-means clustering, and gene-set extraction and profiling.

# Local installation
To run this app locally on your machine, R environment setup is required.
- Download R and RStudio (In the case of macOS, additionally install XQuartz and Xcode)
- Run the following commands once
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

pkgs <- c("shiny","DT","gdata","rstatix","multcomp","tidyverse","ggpubr","venn","ggrepel","ggdendro","ggplotify","gridExtra","cowplot","DESeq2","EBSeq","ggnewscale","edgeR","IHW","qvalue","org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","org.Xl.eg.db","org.Dm.eg.db","org.Ce.eg.db","AnnotationDbi","clusterProfiler","enrichplot","DOSE","msigdbr","genefilter","ComplexHeatmap","shinyBS","plotly","shinyjs","devtools")

for(pkg in pkgs) if (!require(pkg, character.only = T)){
    BiocManager::install(pkg, update = F)
}
```

You may now run RNAseqChef with just one command in R:
```
shiny::runGitHub("RNAseqChef", "Kan-E")
```

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

clusterProfiler, DOSE, msigdbr (for enrichment analysis)
- T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141
- Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015 31(4):608-609
- Dolgalev I (2022). _msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format_. R
  package version 7.5.1, <https://CRAN.R-project.org/package=msigdbr>.

AnnotationDbi, org.Hs.eg.db, org.Mm.eg.db, org.Rn.eg.db, org.Xl.eg.db, org.Dm.eg.db, and org.Ce.eg.db (for genome wide annotation)
- Hervé Pagès, Marc Carlson, Seth Falcon and Nianhua Li (2020). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.52.0. https://bioconductor.org/packages/AnnotationDbi
- Marc Carlson (2020). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.12.0.
- Marc Carlson (2020). org.Mm.eg.db: Genome wide annotation for Mouse. R package version 3.12.0.
- Marc Carlson (2022). org.Rn.eg.db: Genome wide annotation for Rat. R package version 3.15.0.
- Marc Carlson (2022). org.Xl.eg.db: Genome wide annotation for Worm. R package version 3.15.0.
- Marc Carlson (2022). org.Dm.eg.db: Genome wide annotation for Rat. R package version 3.15.0.
- Marc Carlson (2022). org.Ce.eg.db: Genome wide annotation for Worm. R package version 3.15.0.

genefilter (for z-score normalization)
- R. Gentleman, V. Carey, W. Huber and F. Hahne (2021). genefilter: methods for filtering genes from high-throughput experiments. R package version 1.72.1.

ComplexHeatmap (for heatmap and k-means clustering)
- Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.

ggplot2 and ggpubr (for boxplot and scater plot)
- H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
- Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.4.0. https://CRAN.R-project.org/package=ggpubr

venn (for venn diagram analysis)
- Adrian Dusa (2021). venn: Draw Venn Diagrams. R package version 1.10. https://CRAN.R-project.org/package=venn

dplyr and tidyr (for data manipulation)
- Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.7. https://CRAN.R-project.org/package=dplyr
- Hadley Wickham (2021). tidyr: Tidy Messy Data. R package version 1.1.3. https://CRAN.R-project.org/package=tidyr

# License
This shiny code is licensed under the GPLv3. Please see the file LICENSE.md for information.
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
<kaneto@kumamoto-u.ac.jp>
