FROM rocker/shiny-verse:latest
MAINTAINER Kan Etoh
RUN apt-get update && apt-get install -y \
    build-essential \
    libglpk40 \
    libbz2-dev \
    liblzma-dev \
    libgsl-dev \
    r-cran-gsl \
    wget \
    git \
    libmagick++-dev \
    libpoppler-dev \
    libpoppler-cpp-dev
RUN R -e "install.packages('BiocManager')" && \
    R -e "BiocManager::install('shiny', update = F)" && \
    R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/lasso2/lasso2_1.2-22.tar.gz',repos = NULL, type = 'source')" && \
    R -e "BiocManager::install('DT', update = F)" && \
    R -e "BiocManager::install('readxl', update = F)" && \
    R -e "BiocManager::install('rstatix', update = F)" && \
    R -e "BiocManager::install('multcomp', update = F)" && \
    R -e "BiocManager::install('venn', update = F)" && \
    R -e "BiocManager::install('ggrepel', update = F)" && \
    R -e "BiocManager::install('ggdendro', update = F)" && \
    R -e "BiocManager::install('ggplotify', update = F)" && \
    R -e "BiocManager::install('gridExtra', update = F)" && \
    R -e "BiocManager::install('cowplot', update = F)" && \
    R -e "BiocManager::install('DESeq2', update = F)" && \
    R -e "BiocManager::install('EBSeq', update = F)" && \
    R -e "BiocManager::install('ggnewscale', update = F)" && \
    R -e "BiocManager::install('edgeR', update = F)" && \
    R -e "BiocManager::install('IHW', update = F)" && \
    R -e "BiocManager::install('qvalue', update = F)" && \
    R -e "BiocManager::install('org.Hs.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Mm.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Rn.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Xl.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Dm.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Ce.eg.db', update = F)" && \
    R -e "BiocManager::install('DOSE', update = F)" && \
    R -e "BiocManager::install('msigdbr', update = F)" && \
    R -e "BiocManager::install('genefilter', update = F)" && \
    R -e "BiocManager::install('ComplexHeatmap', update = F)" && \
    R -e "BiocManager::install('shinyBS', update = F)" && \
    R -e "BiocManager::install('plotly', update = F)" && \
    R -e "BiocManager::install('shinyjs', update = F)" && \
    R -e "BiocManager::install('DEGreport', update = F)" && \
    R -e "BiocManager::install('dorothea', update = F)" && \
    R -e "BiocManager::install('umap', update = F)" && \
    R -e "BiocManager::install('ggpubr', update = F)" && \
    R -e "BiocManager::install('biomaRt', update = F)" && \
    R -e  "devtools::install_github('YuLab-SMU/clusterProfiler.dplyr')" && \
    R -e "BiocManager::install('enrichplot', update = F)" && \
    R -e "BiocManager::install('clusterProfiler', update = F)" && \
    R -e "BiocManager::install('monaLisa', update = F)" && \
    R -e "BiocManager::install('GenomicRanges', update = F)" && \
    R -e "BiocManager::install('BiocParallel', update = F)" && \
    R -e "BiocManager::install('SummarizedExperiment', update = F)" && \
    R -e "BiocManager::install('JASPAR2020', update = F)" && \
    R -e "BiocManager::install('TFBSTools', update = F)" && \
    R -e "BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene', update = F)" && \
    R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene', update = F)" && \
    R -e "BiocManager::install('BSgenome.Mmusculus.UCSC.mm10', update = F)" && \
    R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19', update = F)" && \
    R -e "BiocManager::install('org.Bt.eg.db')" && \
    R -e "BiocManager::install('org.Cf.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Dr.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Gg.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Mmu.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Pt.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Sc.sgd.db', update = F)" && \
    R -e "BiocManager::install('org.At.tair.db', update = F)" && \
    R -e "BiocManager::install('colorspace', update = F)" && \
    R -e "BiocManager::install('magick', update = F)" && \
    R -e "BiocManager::install('pdftools', update = F)"
RUN sudo rm -rf /srv/shiny-server/sample-apps /srv/shiny-server/index.html /srv/shiny-server/01_hello /srv/shiny-server/02_text /srv/shiny-server/03_reactivity /srv/shiny-server/04_mpg /srv/shiny-server/05_sliders /srv/shiny-server/06_tabsets /srv/shiny-server/07_widgets /srv/shiny-server/08_html /srv/shiny-server/09_upload /srv/shiny-server/10_download /srv/shiny-server/11_timer
RUN mkdir -p /srv/shiny-server/RNAseqChef
COPY ui.R /srv/shiny-server/RNAseqChef/
COPY server.R /srv/shiny-server/RNAseqChef/
COPY global.R /srv/shiny-server/RNAseqChef/
COPY google-analytics.html /srv/shiny-server/RNAseqChef/
COPY www /srv/shiny-server/RNAseqChef/www/
COPY dds.rds /srv/shiny-server/RNAseqChef/
COPY Rmd /srv/shiny-server/RNAseqChef/Rmd/
COPY navAppend.js /srv/shiny-server/RNAseqChef/
COPY shiny-server.conf /etc/shiny-server/
RUN chown -R shiny:shiny /srv/shiny-server
EXPOSE 3838
CMD exec shiny-server >> /var/log/shiny-server.log 2>&1

