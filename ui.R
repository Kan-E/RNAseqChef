library(shiny)
library(DT)
library(gdata)
library(rstatix)
library(multcomp)
library(ggcorrplot)
library(tidyverse)
library(tools)
library(ggpubr)
library(venn)
library(tidyverse)
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
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(org.Dm.eg.db)
library(org.Ce.eg.db)
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
options(repos = BiocManager::repositories())

shinyUI(fluidPage(
  tags$head(includeHTML(("google-analytics.html"))),
  tags$style(
    type = 'text/css',
    # add the name of the tab you want to use as title in data-value
    HTML(
      ".container-fluid > .nav > li >
                        a[data-value='Title'] {font-size: 20px}"
    )
  ),
  navbarPage(
    footer=p(hr(),p("ShinyApp created by Kan Etoh",align="center",width=4),
             p(("Copyright (C) 2022, code licensed under MIT"),align="center",width=4),
             p(("Code available on Github:"),a("https://github.com/Kan-E/RNAseqChef",href="https://github.com/Kan-E/RNAseqChef"),align="center",width=4)),
    "",
    id='navBar',
    tabPanel("RNAseqChef" ,value='Title', icon = icon("utensils"),
             fluidRow(    
               column(12,
                      h1("RNAseqChef",align="center"),br(),
                      p("RNAseqChef, an RNA-seq data controller highlighting gene expression features, is a web-based application for automated, systematic, and integrated RNA-seq differential expression analysis.",
                        align="center"),
                      br(),
                      column(6, 
                             h4("Pair-wise DEG"),
                             "Detects and visualizes differentially expressed genes",br(),br(),
                             img(src="Pair-wise_DEG.png", width = 600,height = 300), br(),br(),
                             h4("3 conditions DEG"),
                             "Detects and visualizes differentially expressed genes by EBSeq multi-comparison analysis",br(),br(),
                             img(src="3cond_DEG.png", width = 600,height = 475)),
                      column(6, 
                             h4("Venn diagram"),
                             "Detects and visualizes the overlap between DEGs from multiple datasets",br(),br(),
                             img(src="Venn.png", width = 600,height = 230),br(),br(),br(),br(),
                             h4("Normalized count analysis"),
                             "identifies similar samples and gene expression pattern by clustering methods",br(),br(),
                             img(src="Normalized.png", width = 550,height = 500))
               )
             )
    ),
    # pair-wise -------------------------------------
    tabPanel("Pair-wise DEG",
             sidebarLayout(
               # sidebar---------------------------------
               sidebarPanel(
                 radioButtons('data_file_type','Input:',
                              c('Row_count_matrix'="Row1",
                                'Option: Row_count_matrix + Metadata'="Row2"
                              ),selected = "Row1"),
                 # Conditional panels appear based on input.data_file_type selection
                 conditionalPanel(condition="input.data_file_type=='Row1'",
                                  strong("Count matrix format: "),br(),
                                  "The replication number is represented by the underbar.",br(),
                                  "Do not use it for anything else.", br(),
                                  img(src="input_format1.png", height = 125, width = 204), br(),
                                  fileInput("file3",
                                            label = "Select a row count matrix file (txt, csv)",
                                            accept = c("txt", "csv"),
                                            multiple = FALSE,
                                            width = "80%")
                 ),
                 conditionalPanel(condition="input.data_file_type=='Row2'",
                                  strong("Count matrix format: "),br(),
                                  "You can use the matrix file whose column name is accession number, and extract the colums you want to analyze by using",
                                  "the metadata.", br(),
                                  "The replication number is represented by the underbar in the characteristics of metadata.",br(),br(),
                                  img(src="input_format2.png", height = 302, width = 253), br(),
                                  fileInput("file1",
                                            label = "Select a row count matrix file (txt, csv)",
                                            accept = c("txt", "csv"),
                                            multiple = FALSE,
                                            width = "80%"),
                                  fileInput("file2",
                                            label = "Select a metadata file to define samples for the following analysis",
                                            accept = c("txt", "csv"),
                                            multiple = FALSE,
                                            width = "80%")
                 ),
                 radioButtons('DEG_method','DEG analysis method:',
                              c('DESeq2'="DESeq2",
                                'EBSeq'="EBSeq",
                                'edgeR'="edgeR"
                              ),selected = "DESeq2"),
                 conditionalPanel(condition="input.DEG_method=='DESeq2'",
                                  selectInput("FDR_method", "FDR method", c("BH", "Qvalue", "IHW"), selected = "BH")
                 ),
                 conditionalPanel(condition="input.DEG_method=='edgeR'",
                                  selectInput("FDR_method", "FDR method", c("BH", "Qvalue", "IHW"), selected = "BH")
                 ),
                 conditionalPanel(condition="input.DEG_method=='edgeR'"),
                 fluidRow(
                   column(6, selectInput("Species", "Species", c("not selected", "human", "mouse", "rat", "fly", "worm"), selected = "not selected"))),
                 h4("Cut-off conditions:"),
                 fluidRow(
                   column(4, numericInput("fc", "Fold Change", min   = 1, max   = NA, value = 2)),
                   column(4, numericInput("fdr", "FDR", min   = 0, max   = 1, value = 0.05)),
                   column(4, numericInput("basemean", "Basemean", min   = 0, max   = NA, value = 0))
                 ),
                 "Option: Normalized count file:",br(),
                 "You can use normalized count (e.g. TPM count) for basemean cutoff and boxplot.",
                 fileInput("norm_file1",
                           label = "Select a normalized count file",
                           accept = c("txt", "csv"),
                           multiple = FALSE,
                           width = "80%"),
                 actionButton("goButton", "example data"),
                 tags$head(tags$style("#goButton{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                           tags$style("
          body {
            padding: 0 !important;
          }"
                           ))
               ), #sidebarPanel
               
               # Main Panel -------------------------------------
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Input Data",
                            bsCollapse(id="input_collapse_panel",open="Row_count_panel",multiple = TRUE,
                                       bsCollapsePanel(title="Row_count_matrix:",
                                                       value="Row_count_panel",
                                                       dataTableOutput('Row_count_matrix')
                                       ),
                                       bsCollapsePanel(title="Metadata:",
                                                       value="Metadata_panel",
                                                       dataTableOutput('Metadata')
                                       ),
                                       bsCollapsePanel(title="Defined_row_count_matrix:",
                                                       value="D_row_count_matrix_panel",
                                                       fluidRow(
                                                         column(4, downloadButton("download_pair_d_row_count", "Download difined row count"))
                                                       ),
                                                       dataTableOutput('D_Row_count_matrix')
                                       )
                            )
                   ),
                   tabPanel("Result overview",
                            fluidRow(
                              column(4, downloadButton("download_pair_PCA", "Download clustering analysis"))
                            ),
                            plotOutput("PCA"),
                            fluidRow(
                              column(4, downloadButton("download_pair_MA", "Download MA-plot"))
                            ),
                            plotOutput("MA"),
                            bsCollapse(id="input_collapse_pair_DEG",open="DEG_panel",multiple = TRUE,
                                       bsCollapsePanel(title="DEG_result:",
                                                       value="DEG_panel",
                                                       fluidRow(
                                                         column(4, downloadButton("download_pair_DEG_result", "Download DEG result"))
                                                       ),
                                                       dataTableOutput("DEG_result")
                                       ),
                                       bsCollapsePanel(title="DEG_count_up:",
                                                       value="deg_count_up_panel",
                                                       fluidRow(
                                                         column(4, downloadButton("download_pair_deg_count_up", "Download DEG count up"))
                                                       ),
                                                       dataTableOutput("pair_deg_up")
                                       ),
                                       bsCollapsePanel(title="DEG_count_down:",
                                                       value="deg_count_down_panel",
                                                       fluidRow(
                                                         column(4, downloadButton("download_pair_deg_count_down", "Download DEG count down"))
                                                       ),
                                                       dataTableOutput("pair_deg_down")
                                       ),
                                       bsCollapsePanel(title="Normalized_Count_matrix:",
                                                       value="norm_panel",
                                                       fluidRow(
                                                         column(4, downloadButton("download_pair_norm_count", "Download normalized count matrix"))
                                                       ),
                                                       dataTableOutput("Normalized_Count_matrix")
                                       ),
                                       bsCollapsePanel(title="PCA:",
                                                       value="PCA_panel",
                                                       fluidRow(
                                                         column(4, downloadButton("download_pair_PCA_table", "Download PCA table"))
                                                       ),
                                                       dataTableOutput("pair_PCA_data")
                                       )
                            )),
                   tabPanel("GOI profiling",
                            fluidRow(
                              column(4, downloadButton("download_pair_volcano", "Download volcano plot")),
                              column(4, downloadButton("download_pair_GOIheatmap", "Download heatmap"))
                            ),
                            fluidRow(
                              column(4, htmlOutput("GOI")),
                              column(4, htmlOutput("volcano_x")),
                              column(4, htmlOutput("volcano_y"))
                            ),
                            fluidRow(
                              column(8, plotOutput("volcano1")),
                              column(4, plotOutput("GOIheatmap"))
                            ),
                            div(
                              plotOutput("GOIbox", height = "100%"),
                              style = "height: calc(100vh  - 100px)"
                            ),
                            fluidRow(
                              column(4, downloadButton("download_pair_GOIbox", "Download boxplot"))
                            )),
                   tabPanel("Enrichment analysis",
                            fluidRow(
                              column(4, textOutput("Spe1"),
                                     tags$head(tags$style("#Spe1{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                              column(4, htmlOutput("Gene_set")),
                              column(4, downloadButton("download_pair_enrichment", "Download"))
                            ),
                            plotOutput("enrichment1"),
                            plotOutput("enrichment2"),
                            bsCollapse(id="input_collapse_pair_enrich",open="ORA_panel",multiple = TRUE,
                                       bsCollapsePanel(title="Enrichment result:",
                                                       value="ORA_panel",
                                                       fluidRow(
                                                         column(4, downloadButton("download_pair_enrichment_table", "Download enrichment result"))
                                                       ),
                                                       dataTableOutput('pair_enrichment_result')
                                       ),
                                       bsCollapsePanel(title="GSEA result:",
                                                       value="GSEA_panel",
                                                       fluidRow(
                                                         column(4, downloadButton("download_pair_GSEA_table", "Download GSEA result"))
                                                       ),
                                                       dataTableOutput('pair_GSEA_result')
                                       )
                            )
                   )
                 )
               ) # main panel
             ) #sidebarLayout
    ), #tabPanel
    # 3conditions -------------------------------------
    tabPanel("3 conditions DEG",
             sidebarLayout(
               
               # sidebar_3conditions---------------------------------
               sidebarPanel(
                 radioButtons('data_file_type2','Input:',
                              c('Row_count_matrix'="Row3",
                                'Option: Row_count_matrix + Metadata'="Row4"
                              ),selected = "Row3"),
                 # Conditional panels appear based on input.data_file_type selection
                 conditionalPanel(condition="input.data_file_type2=='Row3'",
                                  strong("Count matrix format: "),br(),
                                  "The replication number is represented by the underbar.",br(),
                                  "Do not use it for anything else.", br(),
                                  img(src="input_format1.png", height = 125, width = 204), br(),
                                  fileInput("file4",
                                            label = "Select a row count matrix file (txt, csv)",
                                            accept = c("txt", "csv"),
                                            multiple = FALSE,
                                            width = "80%")
                 ),
                 conditionalPanel(condition="input.data_file_type2=='Row4'",
                                  strong("Count matrix format: "),br(),
                                  "You can use the matrix file whose column name is accession number, and extract the colums you want to analyze by using",
                                  "the metadata.", br(),
                                  "The replication number is represented by the underbar in the characteristics of metadata.",br(),br(),
                                  img(src="input_format2.png", height = 302, width = 253), br(),
                                  fileInput("file5",
                                            label = "Select a row count matrix file (txt, csv)",
                                            accept = c("txt", "csv"),
                                            multiple = FALSE,
                                            width = "80%"),
                                  fileInput("file6",
                                            label = "Select a metadata file to define samples for the following analysis",
                                            accept = c("txt", "csv"),
                                            multiple = FALSE,
                                            width = "80%")
                 ),
                 fluidRow(
                   column(6, selectInput("Species2", "Species", c("not selected", "human", "mouse", "rat", "fly", "worm"), selected = "not selected"))),
                 h4("Cut-off conditions:"),
                 fluidRow(
                   column(4, numericInput("fc2", "Fold Change", min   = 0, max   = NA, value = 2)),
                   column(4, numericInput("fdr2", "FDR", min   = 0, max   = NA, value = 0.05)),
                   column(4, numericInput("basemean2", "Basemean", min   = 0, max   = NA, value = 0))
                 ),
                 "Option: Normalized count file:",br(),
                 "You can use normalized count (e.g. TPM count) for basemean cutoff and boxplot.",
                 fileInput("norm_file2",
                           label = "Select a normalized count file",
                           accept = c("txt", "csv"),
                           multiple = FALSE,
                           width = "80%"),
                 actionButton("goButton2", "example data"),
                 tags$head(tags$style("#goButton{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                           tags$style("
          body {
            padding: 0 !important;
          }"
                           ))
               ), #sidebarPanel
               
               # Main Panel -------------------------------------
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Input 3 conditions Data",
                            bsCollapse(id="input_collapse_panel2",open="Row_count_panel2",multiple = FALSE,
                                       bsCollapsePanel(title="Row_count_matrix:",
                                                       value="Row_count_panel2",
                                                       dataTableOutput('Row_count_matrix2')
                                       ),
                                       bsCollapsePanel(title="Metadata:",
                                                       value="Metadata_panel2",
                                                       dataTableOutput('Metadata2')
                                       ),
                                       bsCollapsePanel(title="Defined_row_count_matrix:",
                                                       value="D_row_count_matrix_panel2",
                                                       fluidRow(
                                                         column(4, downloadButton("download_cond3_d_row_count", "Download difined row count"))
                                                       ),
                                                       dataTableOutput('D_Row_count_matrix2')
                                       )
                            )
                   ),
                   tabPanel("Result overview",
                            fluidRow(
                              column(4, downloadButton("download_3cond_PCA", "Download clustering analysis"))
                            ),
                            plotOutput("PCA2"),
                            fluidRow(
                              column(4, downloadButton("download_3cond_scatter1", "Download scatter plot1")),
                              column(4, downloadButton("download_3cond_DEG_table1", "Download DEG_result1"))
                            ),
                            plotOutput("scatter_1"),
                            dataTableOutput("DEG_result2_1"),
                            fluidRow(
                              column(4, downloadButton("download_3cond_scatter2", "Download scatter plot2")),
                              column(4, downloadButton("download_3cond_DEG_table2", "Download DEG_result2"))
                            ),
                            plotOutput("scatter_2"),
                            dataTableOutput("DEG_result2_2"),
                            fluidRow(
                              column(4, downloadButton("download_3cond_scatter3", "Download scatter plot3")),
                              column(4, downloadButton("download_3cond_DEG_table3", "Download DEG_result3"))
                            ),
                            plotOutput("scatter_3"),
                            dataTableOutput("DEG_result2_3"),
                            bsCollapse(id="input_collapse_3_DEG",open="norm_panel2",multiple = TRUE,
                                       bsCollapsePanel(title="Normalized_Count_matrix:",
                                                       value="norm_panel2",
                                                       fluidRow(
                                                         column(4, downloadButton("download_cond3_norm_count", "Download normalized count"))
                                                       ),
                                                       dataTableOutput("Normalized_Count_matrix2")
                                       ),
                                       bsCollapsePanel(title="PCA:",
                                                       value="PCA3_panel",
                                                       fluidRow(
                                                         column(4, downloadButton("download_cond3_pca_table", "Download PCA table"))
                                                       ),
                                                       dataTableOutput("PCA3_data")
                                       )
                            )),
                   tabPanel("GOI profiling",
                            fluidRow(
                              column(4, downloadButton("download_3cond_scatter", "Download scatter plot")),
                              column(4, downloadButton("download_3cond_GOIheat", "Download heatmap"))
                            ),
                            fluidRow(
                              column(4, htmlOutput("GOI2"))
                            ),
                            plotOutput("cond3_GOIheatmap"),
                            div(
                              plotOutput("cond3_GOIboxplot", height = "100%"),
                              style = "height: calc(100vh  - 100px)"
                            ),
                            fluidRow(
                              column(4, downloadButton("download_3cond_GOIbox", "Download boxplot"))
                            )
                   ),
                   tabPanel("Enrichment analysis",
                            fluidRow(
                              column(4, textOutput("Spe2"),
                                     tags$head(tags$style("#Spe2{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                              column(4, htmlOutput("Gene_set2")),
                              column(4, downloadButton("download_3enrich", "Download"))
                            ),
                            plotOutput("keggenrichment2_1"),
                            plotOutput("keggenrichment2_2"),
                            plotOutput("keggenrichment2_3"),
                            bsCollapse(id="input_collapse_3_enrich",open="ORA3_1_panel",multiple = TRUE,
                                       bsCollapsePanel(title="Enrichment result1:",
                                                       value="ORA3_1_panel",
                                                       dataTableOutput('enrichment3_result_1')
                                       ),
                                       bsCollapsePanel(title="Enrichment result2:",
                                                       value="ORA3_2_panel",
                                                       dataTableOutput('enrichment3_result_2')
                                       ),
                                       bsCollapsePanel(title="Enrichment result3:",
                                                       value="ORA3_3_panel",
                                                       dataTableOutput('enrichment3_result_3')
                                       )
                            ))
                 )
               ) # main panel
             ) #sidebarLayout
    ), #tabPanel
    # venn diagram analysis---------------------------------
    tabPanel("Venn diagram",
             # titlePanel(h5("Upload Files")),
             sidebarLayout(
               # venn diagram analysis---------------------------------
               sidebarPanel(
                 fileInput(
                   inputId = "files", 
                   label = "Select gene list files (txt)", 
                   multiple = TRUE,
                   accept = c("text/txt",
                              "text/tab-separated-values,text/plain",
                              ".txt")
                 ),
                 "Count extraction of intersection",
                 fileInput(
                   inputId = "file_for_venn", 
                   label = "Choose a normalized count file (txt)", 
                   multiple = TRUE,
                   accept = c("text/txt",
                              "text/tab-separated-values,text/plain",
                              ".txt")
                 ),
                 br(),br(),
                 h4("Input for integrated heatmap"),
                 fileInput(
                   inputId = "countfiles", 
                   label = "Select normalized count files (txt)", 
                   multiple = TRUE,
                   accept = c("text/txt",
                              "text/tab-separated-values,text/plain",
                              ".txt")
                 ),
                 fluidRow(
                   column(6, selectInput(
                     inputId = "pre_zscoring", 
                     label = "Option: Pre-zscoring", 
                     multiple = FALSE,choices = c("TRUE", "FALSE"), selected = "TRUE"))
                 ),
                 actionButton("goButton_venn", "example data"),
                 tags$head(tags$style("#goButton{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                           tags$style("
          body {
            padding: 0 !important;
          }"
                           )
                 ) #sidebarPanel
               ),
               
               # Main Panel -------------------------------------
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Input gene lists",
                            plotOutput("venn"),
                            bsCollapse(id="venn_collapse_panel",open="venn_result_panel",multiple = TRUE,
                                       bsCollapsePanel(title="venn result:",
                                                       value="venn_result_panel",
                                                       downloadButton("download_venn_result", "Download venn result"),
                                                       dataTableOutput("venn_result")
                                       ),
                                       bsCollapsePanel(title="intersection count table:",
                                                       value="intersection_count_panel",
                                                       fluidRow(
                                                         column(4, htmlOutput("select_file1")),
                                                         column(4, downloadButton("download_intersection_count_table", "Download intersection count table"))
                                                       ),
                                                       dataTableOutput("intersection_count")
                                       )
                            )
                   ),
                   tabPanel("integrated_heatmap",
                            fluidRow(
                              column(4, htmlOutput("select_file2")),
                              column(4, downloadButton("download_integrated_heatmap", "download integrated heatmap"))
                            ),
                            plotOutput("intheatmap"),
                            bsCollapse(id="int_result_collapse_panel",open="integrated_count_panel",multiple = TRUE,
                                       bsCollapsePanel(title="integrated normalized count:",
                                                       value="integrated_count_panel",
                                                       downloadButton("download_integrated_count_table", "Download integrated count table"),
                                                       dataTableOutput("integrated_count_table")
                                       ),
                                       bsCollapsePanel(title="integrated zscored normalized count:",
                                                       value="integrated_z_count_panel",
                                                       downloadButton("download_integrated_z_count_table", "Download integrated zscored count table"),
                                                       dataTableOutput("integrated_count_z_table")
                                       )
                            )
                   )
                 )
               ) #sidebarLayout
             ) #tabPanel
    ),
    # Normalized count data analysis -------------------------------------
    tabPanel("Normalized count analysis",
             # titlePanel(h5("Upload Files")),
             sidebarLayout(
               # Normalized count data analysis---------------------------------
               sidebarPanel(
                 radioButtons('data_file_type3','Input:',
                              c('Normalized_count_matrix'="Row5",
                                'Option: Normalized_count_matrix + Metadata'="Row6"
                              ),selected = "Row5"),
                 # Conditional panels appear based on input.data_file_type selection
                 conditionalPanel(condition="input.data_file_type3=='Row5'",
                                  strong("Count matrix format: "),br(),
                                  "The replication number is represented by the underbar.",br(),
                                  "Do not use it for anything else.", br(),
                                  fileInput("file7",
                                            label = "Select a row count matrix file (txt, csv)",
                                            accept = c("txt", "csv"),
                                            multiple = FALSE,
                                            width = "80%")
                 ),
                 conditionalPanel(condition="input.data_file_type3=='Row6'",
                                  strong("Count matrix format: "),br(),
                                  "You can use the matrix file whose column name is accession number, and extract the colums you want to analyze by using",
                                  "the metadata.", br(),
                                  "The replication number is represented by the underbar in the characteristics of metadata.",br(),br(),
                                  img(src="input_format2.png", height = 302, width = 253), br(),
                                  fileInput("file8",
                                            label = "Select a normalized count matrix file (txt, csv)",
                                            accept = c("txt", "csv"),
                                            multiple = FALSE,
                                            width = "80%"),
                                  fileInput("file9",
                                            label = "Select a metadata file to define samples for the following analysis",
                                            accept = c("txt", "csv"),
                                            multiple = FALSE,
                                            width = "80%")
                 ),
                 fileInput("file10",
                           label = "Option: Select a gene list file for GOI extraction",
                           accept = c("txt", "csv", "xlsx"),
                           multiple = FALSE,
                           width = "80%"),
                 h4("Cut-off conditions:"),
                 fluidRow(
                   column(4, numericInput("basemean3", "Basemean", min   = 0, max   = NA, value = 0))
                 ),
                 fluidRow(
                   column(6, selectInput("Species3", "Species", c("not selected", "human", "mouse", "rat", "fly", "worm"), selected = "not selected"))),
                 actionButton("goButton3", "example data"),
                 tags$head(tags$style("#goButton{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                           tags$style("
          body {
            padding: 0 !important;
          }"
                           )
                 ) #sidebarPanel
               ),
               
               # Main Panel -------------------------------------
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Input Normalized Count Data",
                            bsCollapse(id="norm_input_collapse_panel",open="Normalized_count_panel",multiple = FALSE,
                                       bsCollapsePanel(title="Norm_count_matrix:",
                                                       value="Norm_count_panel",
                                                       dataTableOutput('norm_count_input1')
                                       ),
                                       bsCollapsePanel(title="Metadata:",
                                                       value="norm_Metadata_panel",
                                                       dataTableOutput('norm_Metadata')
                                       ),
                                       bsCollapsePanel(title="Gene list:",
                                                       value="gene_list_panel",
                                                       dataTableOutput('Gene_list')
                                       ),
                                       bsCollapsePanel(title="Defined_normalized_count_matrix:",
                                                       value="D_norm_count_matrix_panel",
                                                       fluidRow(
                                                         column(4, downloadButton("download_d_norm_count", "Download difined normalized count"))
                                                       ),
                                                       dataTableOutput('d_norm_count')
                                       )
                            )
                   ),
                   tabPanel("Clustering",
                            plotOutput("norm_PCA"),
                            fluidRow(
                              column(3, downloadButton("download_norm_PCA", "Download Clustering")),
                              column(3, downloadButton("download_norm_heatmap", "Download heatmap")),
                              column(6, plotOutput("norm_heatmap"))
                            ),
                            bsCollapse(id="norm_result_collapse_panel",open="cutoff_count_panel",multiple = TRUE,
                                       bsCollapsePanel(title="basemean-cutoff_Norm_count_matrix:",
                                                       value="cutoff_count_panel",
                                                       dataTableOutput("d_norm_count_cutoff")
                                       ),
                                       bsCollapsePanel(title="PCA_table:",
                                                       value="norm_PCA_panel",
                                                       dataTableOutput("norm_PCA_table")
                                       )
                            ),
                            fluidRow(
                              column(4, downloadButton("download_norm_pca_table", "Download PCA table"))
                            ),
                   ),
                   tabPanel("GOI profiling",
                            fluidRow(
                              column(4, downloadButton("download_norm_GOIheat", "Download heatmap"))
                            ),
                            fluidRow(
                              column(4, htmlOutput("GOI3")),
                              column(8, plotOutput("norm_GOIheatmap"))
                            ),
                            div(
                              plotOutput("norm_GOIboxplot", height = "100%"),
                              style = "height: calc(100vh  - 100px)"
                            ),
                            fluidRow(
                              column(4, downloadButton("download_norm_GOIbox", "Download boxplot"))
                            )
                   ),
                   tabPanel("k-means clustering",
                            fluidRow(
                              column(4, htmlOutput("norm_kmeans_num"),
                                     downloadButton("download_norm_kmeans_heatmap", "Download heatmap")),
                              column(8, plotOutput("norm_kmeans_heatmap"))
                            ),
                            downloadButton("download_norm_kmeans_cluster", "Download k-means clustering result"),
                            dataTableOutput("norm_kmeans_count_table")
                   )
                 )
               )
             ) #sidebarLayout
    ), #tabPanel
    #Instruction--------------------------
    tabPanel("Reference",
             fluidRow(    
               column(10,
                      h2("To cite RNAseqChef:"),
                      p("Web tool:",br(), "Kan Etoh: RNAseqChef (2022)", a("https://kan-e.shinyapps.io/RNAseqChef/",href="https://kan-e.shinyapps.io/RNAseqChef/")),
                      p("Source code:",br(), "Kan Etoh: RNAseqChef: web application for automated, systematic, and integrated RNA-seq differential expression analysis. (2022)", a("https://github.com/Kan-E/RNAseqChef/",href="https://github.com/Kan-E/RNAseqChef/")),br(),
                      br(),
                      h2("Reference:"),
                      "- Winston Chang, Joe Cheng, JJ Allaire, Carson Sievert, Barret Schloerke, Yihui Xie, Jeff Allen, Jonathan McPherson, Alan Dipert and Barbara Borges (2021). shiny: Web Application Framework for R. R package version 1.7.1. https://CRAN.R-project.org/package=shiny",br(),
                      "- Ning Leng and Christina Kendziorski (2020). EBSeq: An R package for gene and isoform
  differential expression analysis of RNA-seq data. R package version 1.30.0.",br(),
                      "- Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for
  RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)",br(),
                      "- Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential
  expression analysis of digital gene expression data. Bioinformatics 26, 139-140",br(),
                      "- Nikolaos Ignatiadis, Bernd Klaus, Judith Zaugg and Wolfgang Huber (2016): Data-driven hypothesis
  weighting increases detection power in genome-scale multiple testing. Nature Methods 13:577,
  doi: 10.1038/nmeth.3885",br(),
                      "- John D. Storey, Andrew J. Bass, Alan Dabney and David Robinson (2021). qvalue: Q-value
  estimation for false discovery rate control. R package version 2.26.0.
  http://github.com/jdstorey/qvalue",br(),
                      "- Andrie de Vries and Brian D. Ripley (2020). ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'. R package version 0.1.22. https://CRAN.R-project.org/package=ggdendro",br(),
                      "- T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141",br(),
                      "- Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015 31(4):608-609",br(),
                      "- Hervé Pagès, Marc Carlson, Seth Falcon and Nianhua Li (2020). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.52.0. https://bioconductor.org/packages/AnnotationDbi",br(),
                      "- Marc Carlson (2020). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.12.0.",br(),
                      "- Marc Carlson (2020). org.Mm.eg.db: Genome wide annotation for Mouse. R package version 3.12.0.",br(),
                      "- R. Gentleman, V. Carey, W. Huber and F. Hahne (2021). genefilter: methods for filtering genes from high-throughput experiments. R package version 1.72.1.",br(),
                      "- Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.",br(),
                      "- H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.", br(),
                      "- Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.4.0. https://CRAN.R-project.org/package=ggpubr",br(),
                      "- Adrian Dusa (2021). venn: Draw Venn Diagrams. R package version 1.10. https://CRAN.R-project.org/package=venn",br(),
                      "- Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.7. https://CRAN.R-project.org/package=dplyr",br(),
                      "- Hadley Wickham (2021). tidyr: Tidy Messy Data. R package version 1.1.3. https://CRAN.R-project.org/package=tidyr"
               )
             )
    )
  ) 
)
)
