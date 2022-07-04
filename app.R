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
options(repos = BiocManager::repositories())
msigdbr_species <- msigdbr_species()$species_name
popoverTempate <- 
  '<div class="popover popover-lg" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content"></div></div>'


ui<- fluidPage(
  tags$head(includeHTML(("google-analytics.html"))),
  tags$style(
    type = 'text/css',
    HTML(
      ".container-fluid > .nav > li >
                        a[data-value='Title'] {font-size: 20px}",
      ".navbar{margin:3px;}"
    )
  ),
  tags$style(".popover{
            max-width: 500px;
          }"),
  navbarPage(
    footer=p(hr(),p("ShinyApp created by Kan Etoh",align="center",width=4),
             p(("Copyright (C) 2022, code licensed under GPLv3"),align="center",width=4),
             p(("Code available on Github:"),a("https://github.com/Kan-E/RNAseqChef",href="https://github.com/Kan-E/RNAseqChef"),align="center",width=4)),
  "",
  id='navBar',
  tabPanel("RNAseqChef" ,value='Title', icon = icon("utensils"),
           fluidRow(
             column(12,
                    br(),br(),
                    h1(strong("RNAseqChef"),align="center"),br(),
                    p("RNAseqChef, an RNA-seq data controller highlighting gene expression features, is a web-based application for automated, systematic, and integrated RNA-seq differential expression analysis.",
                      align="center"),br(),br(),style={'background-color:beige;font-size: 16px;'},
                    ),
             column(12,
                    column(6,hr(),
                           h4(strong("Pair-wise DEG")),
                           "Detects and visualizes differentially expressed genes",
                           img(src="Pair-wise_DEG.png", width = 600,height = 300),br(),hr(),
                           h4(strong("3 conditions DEG")),
                           "Detects and visualizes differentially expressed genes by EBSeq multi-comparison analysis",
                           img(src="3cond_DEG.png", width = 600,height = 475),br(),hr(),
                           h4(strong("Multi DEG")),
                           "Detects and visualizes differentially expressed genes by DESeq2 LRT following clustering analysis",
                           img(src="Multi DEG.png", width = 600,height = 682)),
                    column(6,hr(),
                           h4(strong("Venn diagram")),
                           "Detects and visualizes the overlap between DEGs from multiple datasets",
                           img(src="Venn.png", width = 600,height = 230),br(),hr(),
                           h4(strong("Normalized count analysis")),
                           "identifies similar samples and gene expression patterns by clustering methods",
                           img(src="Normalized.png", width = 550,height = 500),hr(),
                           h4(strong("Enrichment viewer")),
                           "determines and visualizes biological functions of gene set of interest",
                           img(src="enrichment_viewer.png", width = 550,height = 250)
                           )
             )
           )
  ),
  # pair-wise -------------------------------------
  tabPanel("Pair-wise DEG",
           sidebarLayout(
             # sidebar---------------------------------
             sidebarPanel(
               radioButtons('data_file_type','Input:',
                            c('Raw_count_matrix'="Row1",
                              'Option: Raw_count_matrix + Metadata'="Row2",
                              'Option: Batch mode (not displayed in the output panel)'="Row11"
                            ),selected = "Row1"),
               conditionalPanel(condition="input.data_file_type=='Row1'",
                                fileInput("file3",
                                          strong(
                                            span("Select a raw count matrix file (txt, csv)"),
                                            span(icon("info-circle"), id = "icon1", 
                                                 options = list(template = popoverTempate))
                                          ),
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%"),
                                bsPopover("icon1", "Count matrix format:", 
                                          content=paste(strong("The replication number"), "is represented by", strong("the underline"),".<br>", strong("Do not use it for anything else"),".<br><br>",
                                                        img(src="input_format1.png", width = 400,height = 250)), 
                                          placement = "right",options = list(container = "body")),
               ),
               conditionalPanel(condition="input.data_file_type=='Row2'",
                                fileInput("file1",
                                          strong(
                                            span("Select a raw count matrix file (txt, csv)"),
                                            span(icon("info-circle"), id = "icon2", 
                                                 options = list(template = popoverTempate))
                                          ),
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%"),
                                fileInput("file2",
                                          "Select a metadata file to define samples for the following analysis",
                                          accept = c("txt", "csv","xlsx"),
                                          multiple = FALSE,
                                          width = "80%"),
                                bsPopover("icon2", "Metadata format:", 
                                          content=paste("The first column is", strong("the sample name"), "used in the raw count data.<br>", 
                                                        "The second column is", strong("the corresponding sample name"), "that matches the sample name in the first column.<br><br>",
                                                        img(src="input_format2.png", width = 400,height = 400)),
                                          placement = "right",options = list(container = "body")),
               ),
               conditionalPanel(condition="input.data_file_type=='Row11'",
                                fileInput("file11",
                                          strong(
                                            span("Select raw count matrix files (txt, csv)"),
                                            span(icon("info-circle"), id = "icon3", 
                                                 options = list(template = popoverTempate))
                                          ),
                                          accept = c("txt", "csv"),
                                          multiple = TRUE,
                                          width = "80%"),
                                bsPopover("icon3", "Count matrix format:", 
                                          content=paste(strong("The replication number"), "is represented by", strong("the underline"),".<br>", strong("Do not use it for anything else"),".<br>",
                                                        "There is no limitation to the number of uploaded files.<br><br>",
                                                        img(src="input_format1.png", width = 400,height = 250)), 
                                          placement = "right",options = list(container = "body")),
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
                 column(6, selectInput("Species", "Species", c("not selected", "Homo sapiens", "Mus musculus", "Rattus norvegicus", "Xenopus laevis", 
                                                               "Drosophila melanogaster", "Caenorhabditis elegans"), selected = "not selected"))),
               h4("Cut-off conditions:"),
               fluidRow(
                 column(4, numericInput("fc", "Fold Change", min   = 1, max   = NA, value = 2)),
                 column(4, numericInput("fdr", "FDR", min   = 0, max   = 1, value = 0.05)),
                 column(4, numericInput("basemean", "Basemean", min   = 0, max   = NA, value = 0))
               ),
               fileInput("norm_file1",
                         strong(
                           span("Option: Select a normalized count file"),
                           span(icon("info-circle"), id = "icon4", 
                                options = list(template = popoverTempate))
                         ),
                         accept = c("txt", "csv"),
                         multiple = TRUE,
                         width = "80%"),
               bsPopover("icon4", "Option: Normalized count file:", 
                         content="You can use a normalized count data (e.g. TPM count) for basemean cutoff and boxplot.", 
                         placement = "right",options = list(container = "body")),
               actionButton("goButton", "example data (mouse)"),
               tags$head(tags$style("#goButton{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                         tags$style("
          body {
            padding: 0 !important;
          }"
                         )),
               conditionalPanel(condition="input.data_file_type=='Row11'",
                                fluidRow(column(7),
                                  column(5, downloadButton("downloadData", "Download zip"))
                                )
               ),
             ), #sidebarPanel

             # Main Panel -------------------------------------
             mainPanel(
               tabsetPanel(
                 type = "tabs",
                 tabPanel("Input Data",
                          bsCollapse(id="input_collapse_panel",open="Row_count_panel",multiple = TRUE,
                                     bsCollapsePanel(title="Raw_count_matrix:",
                                                     value="Row_count_panel",
                                                     dataTableOutput('Row_count_matrix')
                                     ),
                                     bsCollapsePanel(title="Metadata:",
                                                     value="Metadata_panel",
                                                     dataTableOutput('Metadata')
                                     ),
                                     bsCollapsePanel(title="Defined_raw_count_matrix:",
                                                     value="D_row_count_matrix_panel",
                                                     fluidRow(
                                                       column(4, downloadButton("download_pair_d_row_count", "Download defined raw count"))
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
                            c('Raw_count_matrix'="Row3",
                              'Option: Raw_count_matrix + Metadata'="Row4"
                            ),selected = "Row3"),
               # Conditional panels appear based on input.data_file_type selection
               conditionalPanel(condition="input.data_file_type2=='Row3'",
                                fileInput("file4",
                                          strong(
                                            span("Select a raw count matrix file (txt, csv)"),
                                            span(icon("info-circle"), id = "icon5", 
                                                 options = list(template = popoverTempate))
                                          ),
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%"),
                                bsPopover("icon5", "Count matrix format:", 
                                          content=paste(strong("The replication number"), "is represented by", strong("the underline"),".<br>", strong("Do not use it for anything else"),".<br><br>",
                                                        img(src="input_format1.png", width = 400,height = 250)), 
                                          placement = "right",options = list(container = "body")),
               ),
               conditionalPanel(condition="input.data_file_type2=='Row4'",
                                fileInput("file5",
                                          strong(
                                            span("Select a raw count matrix file (txt, csv)"),
                                            span(icon("info-circle"), id = "icon6", 
                                                 options = list(template = popoverTempate))
                                          ),
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%"),
                                fileInput("file6",
                                          label = "Select a metadata file to define samples for the following analysis",
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%"),
                                bsPopover("icon6", "Metadata format:", 
                                          content=paste("The first column is", strong("the sample name"), "used in the raw count data.<br>", 
                                                        "The second column is", strong("the corresponding sample name"), "that matches the sample name in the first column.<br><br>",
                                                        img(src="input_format2.png", width = 400,height = 400)),
                                          placement = "right",options = list(container = "body")),
               ),
               fluidRow(
                 column(6, selectInput("Species2", "Species", c("not selected", "Homo sapiens", "Mus musculus", "Rattus norvegicus", "Xenopus laevis", 
                                                                "Drosophila melanogaster", "Caenorhabditis elegans"), selected = "not selected"))),
               h4("Cut-off conditions:"),
               fluidRow(
                 column(4, numericInput("fc2", "Fold Change", min   = 1, max   = NA, value = 2)),
                 column(4, numericInput("fdr2", "FDR", min   = 0, max   = 0.05, value = 0.05)),
                 column(4, numericInput("basemean2", "Basemean", min   = 0, max   = NA, value = 1))
               ),
               fileInput("norm_file2",
                         strong(
                           span("Option: Select a normalized count file"),
                           span(icon("info-circle"), id = "icon7", 
                                options = list(template = popoverTempate))
                         ),
                         accept = c("txt", "csv"),
                         multiple = FALSE,
                         width = "80%"),
               bsPopover("icon7", "Option: Normalized count file:", 
                         content="You can use a normalized count data (e.g. TPM count) for basemean cutoff and boxplot.", 
                         placement = "right",options = list(container = "body")),
               actionButton("goButton2", "example data (mouse)"),
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
                                     bsCollapsePanel(title="Raw_count_matrix:",
                                                     value="Row_count_panel2",
                                                     dataTableOutput('Row_count_matrix2')
                                     ),
                                     bsCollapsePanel(title="Metadata:",
                                                     value="Metadata_panel2",
                                                     dataTableOutput('Metadata2')
                                     ),
                                     bsCollapsePanel(title="Defined_raw_count_matrix:",
                                                     value="D_row_count_matrix_panel2",
                                                     fluidRow(
                                                       column(4, downloadButton("download_cond3_d_row_count", "Download defined row count"))
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
                            column(4, downloadButton("download_cond3_enrichment", "Download"))
                          ),
                          plotOutput("keggenrichment2_1"),
                          plotOutput("keggenrichment2_2"),
                          plotOutput("keggenrichment2_3"),
                          bsCollapse(id="input_collapse_3_enrich",open="ORA3_1_panel",multiple = TRUE,
                                     bsCollapsePanel(title="Enrichment result1:",
                                                     value="ORA3_1_panel",
                                                     fluidRow(
                                                       column(4, downloadButton("download_cond3_enrichment_table1", "Download enrichment result table 1"))
                                                     ),
                                                     dataTableOutput('enrichment3_result_1')
                                     ),
                                     bsCollapsePanel(title="Enrichment result2:",
                                                     value="ORA3_2_panel",
                                                     fluidRow(
                                                       column(4, downloadButton("download_cond3_enrichment_table2", "Download enrichment result table 2"))
                                                     ),
                                                     dataTableOutput('enrichment3_result_2')
                                     ),
                                     bsCollapsePanel(title="Enrichment result3:",
                                                     value="ORA3_3_panel",
                                                     fluidRow(
                                                       column(4, downloadButton("download_cond3_enrichment_table3", "Download enrichment result table 3"))
                                                     ),
                                                     dataTableOutput('enrichment3_result_3')
                                     )
                          ))
               )
             ) # main panel
           ) #sidebarLayout
  ), #tabPanel
  # Multi DEG -------------------------------------
  tabPanel("Multi DEG",
           sidebarLayout(
             # sidebar---------------------------------
             sidebarPanel(
               radioButtons('multi_data_file_type','Input:',
                            c('Raw_count_matrix (One-factor multi-condition)'="Row1",
                              'Raw_count_matrix + metadata (Two-factor multi-condition)'="Row2"
                            ),selected = "Row1"),
               conditionalPanel(condition="input.multi_data_file_type=='Row1'",
                                fileInput("multi_file1",
                                          strong(
                                            span("Select a raw count matrix file (txt, csv)"),
                                            span(icon("info-circle"), id = "icon8", 
                                                 options = list(template = popoverTempate))
                                          ),
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%"),
                                bsPopover("icon8", "Count matrix format:", 
                                          content=paste(strong("The replication number"), "is represented by", strong("the underline"),".<br>", strong("Do not use it for anything else"),".<br><br>",
                                                        img(src="input_format1.png", width = 400,height = 250)), 
                                          placement = "right",options = list(container = "body")),
               ),
               conditionalPanel(condition="input.multi_data_file_type=='Row2'",
                                fileInput("multi_file2",
                                          strong(
                                            span("Select a raw count matrix file (txt, csv)"),
                                            span(icon("info-circle"), id = "icon9", 
                                                 options = list(template = popoverTempate))
                                          ),
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%"),
                                fileInput("multi_file3",
                                          label = "Select a metadata file to define samples for the following analysis",
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%"),
                                bsPopover("icon9", "Metadata format:", 
                                          content=paste("The first column is", strong("the sample name"), "used in the raw count data.<br>", 
                                                        "The second column is", strong("the corresponding sample name"), "that matches the sample name in the first column.<br><br>",
                                                        img(src="input_format2.png", width = 400,height = 400)),
                                          placement = "right",options = list(container = "body")),
               ),
               fluidRow(column(6,  selectInput("FDR_method6", "FDR method", c("BH", "Qvalue", "IHW"), selected = "BH")),
                 column(6, selectInput("Species6", "Species", c("not selected", "Homo sapiens", "Mus musculus", "Rattus norvegicus", "Xenopus laevis", 
                                                                "Drosophila melanogaster", "Caenorhabditis elegans"), selected = "not selected"))),
               h4("Cut-off conditions:"),
               fluidRow(
                 column(4, numericInput("fc6", "Fold Change", min   = 0, max   = NA, value = 1.5)),
                 column(4, numericInput("fdr6", "FDR", min   = 0, max   = 1, value = 0.05)),
                 column(4, numericInput("basemean6", "Basemean", min   = 0, max   = NA, value = 0))
               ),
               fileInput("multi_norm_file1",
                         strong(
                           span("Option: Select a normalized count file"),
                           span(icon("info-circle"), id = "icon10", 
                                options = list(template = popoverTempate))
                         ),
                         accept = c("txt", "csv"),
                         multiple = TRUE,
                         width = "80%"),
               bsPopover("icon10", "Option: Normalized count file:", 
                         content="You can use a normalized count data (e.g. TPM count) for basemean cutoff and boxplot.", 
                         placement = "right",options = list(container = "body")),
               actionButton("goButton6", "example data (mouse)"),
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
                          bsCollapse(id="multi_input_collapse_panel",open="multi_Row_count_panel",multiple = TRUE,
                                     bsCollapsePanel(title="Raw_count_matrix:",
                                                     value="multi_Row_count_panel",
                                                     dataTableOutput('multi_Row_count_matrix')
                                     ),
                                     bsCollapsePanel(title="Metadata:",
                                                     value="multi_Metadata_panel",
                                                     dataTableOutput('multi_Metadata')
                                     ),
                                     bsCollapsePanel(title="Defined_raw_count_matrix:",
                                                     value="multi_d_Row_count_panel",
                                                     dataTableOutput('multi_d_Row_count_matrix')
                                     )
                          )
                 ),
                 tabPanel("Result overview",
                          fluidRow(
                            column(4, downloadButton("download_multi_PCA", "Download clustering analysis"))
                          ),
                          div(
                            plotOutput("multi_PCA", height = "100%"),
                            style = "height: calc(100vh  - 100px)"
                          ),
                          bsCollapse(id="input_collapse_multi_DEG",open="DEG_panel",multiple = TRUE,
                                     bsCollapsePanel(title="DEG_result:",
                                                     value="DEG_panel",
                                                     fluidRow(
                                                       column(4, downloadButton("download_multi_DEG_result", "Download DEG result"))
                                                     ),
                                                     dataTableOutput("multi_DEG_result")
                                     ),
                                     bsCollapsePanel(title="Normalized_Count_matrix:",
                                                     value="multi_norm_panel",
                                                     fluidRow(
                                                       column(4, downloadButton("download_multi_norm_count", "Download normalized count matrix"))
                                                     ),
                                                     dataTableOutput("multi_Normalized_Count_matrix")
                                     ),
                                     bsCollapsePanel(title="PCA:",
                                                     value="multi_PCA_panel",
                                                     fluidRow(
                                                       column(4, downloadButton("download_multi_PCA_table", "Download PCA table"))
                                                     ),
                                                     dataTableOutput("multi_PCA_data")
                                     )
                          )),
                 tabPanel("Divisive clustering",
                          fluidRow(
                            column(4, downloadButton("download_multi_boxplot", "Download boxplots"))
                          ),
                          fluidRow(
                            column(3, htmlOutput("selectFC")),
                            column(3, textOutput("multi_DEG_total1"),
                                   tags$head(tags$style("#multi_DEG_total1{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"))),
                            column(6, htmlOutput("topP")),
                          ),
                          div(
                            plotOutput("multi_boxplot", height = "100%"),
                            style = "height: calc(100vh  - 100px)"
                          ),
                          bsCollapsePanel(title="DEG_pattern:",
                                          value="multi_deg_pattern_panel",
                                          fluidRow(
                                            column(4, downloadButton("download_deg_pattern_list", "Download DEG pattern list"))
                                          ),
                                          dataTableOutput("multi_pattern1_list")
                          ),
                          bsCollapsePanel(title="DEG_pattern_normalized_count_data:",
                                          value="multi_deg_pattern_count_panel",
                                          fluidRow(
                                            column(4, htmlOutput("multi_select_file1")),
                                            column(4, downloadButton("download_deg_pattern_count", "Download DEG pattern normalized count data"))
                                          ),
                                          dataTableOutput("multi_pattern1_count")
                          ),
                          fluidRow(
                            column(4, htmlOutput("multi_whichGroup1_1")),
                            column(4, htmlOutput("Gene_set7")),
                            column(4, downloadButton("download_multi_cluster_enrichment", "Download dot plot"))
                          ),
                          fluidRow(
                            column(4, textOutput("multi_Spe"),
                                   tags$head(tags$style("#multi_Spe{color: red;
                                         font-size: 20px;
                                         font-style: bold;
                                      }")))
                          ),
                          plotOutput("multi_enrichment3"),
                          fluidRow(
                            column(4, htmlOutput("multi_whichGroup1_2")),
                            column(4, downloadButton("download_multi_enrichment_cnet", "Download cnet plot"))
                          ),
                          plotOutput("multi_enrichment4"),
                          fluidRow(
                            column(4, downloadButton("download_multi_enrichment_table", "Download enrichment result"))
                          ),
                          dataTableOutput('multi_enrichment_result')
                 ),
                 tabPanel("k-means clustering",
                          fluidRow(
                            column(3, htmlOutput("selectFC2")),
                            column(3, textOutput("multi_DEG_total2"),
                                   tags$head(tags$style("#multi_DEG_total2{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"))),
                            column(6, htmlOutput("topP2"))
                          ),
                          fluidRow(
                            column(4, htmlOutput("multi_kmeans_num"),
                                   downloadButton("download_multi_kmeans_heatmap", "Download heatmap"),
                                   downloadButton("download_multi_kmeans_boxplot", "Download boxplots")),
                            column(8, plotOutput("multi_kmeans_heatmap"))
                          ),
                          div(
                            plotOutput("multi_kmeans_boxplot", height = "100%"),
                            style = "height: calc(100vh  - 100px)"
                          ),
                          bsCollapsePanel(title="kmeans_result:",
                                          value="multi_deg_kmeans_pattern_count_panel",
                                          fluidRow(
                                            column(4, downloadButton("download_multi_kmeans_cluster", "Download k-means clustering result"))
                                          ),
                                          dataTableOutput("multi_kmeans_count_table")
                          ),
                          bsCollapsePanel(title="cluster_normalized_count_data:",
                                          value="multi_deg_kmeans_pattern_count_panel",
                                          fluidRow(
                                            column(4, htmlOutput("multi_select_file2")),
                                            column(4, downloadButton("download_deg_kmeans_pattern_count", "Download DEG pattern normalized count data"))
                                          ),
                                          dataTableOutput("multi_pattern2_count")
                          ),
                          
                          fluidRow(
                            column(4, htmlOutput("multi_whichGroup2_1")),
                            column(4, htmlOutput("Gene_set8")),
                            column(4, downloadButton("download_multi_cluster_enrichment2", "Download dot plot"))
                          ),
                          fluidRow(
                            column(4, textOutput("multi_Spe2"),
                                   tags$head(tags$style("#multi_Spe2{color: red;
                                         font-size: 20px;
                                         font-style: bold;
                                      }")))
                          ),
                          plotOutput("multi_enrichment5"),
                          fluidRow(
                            column(4, htmlOutput("multi_whichGroup2_2")),
                            column(4, downloadButton("download_multi_enrichment_cnet2", "Download cnet plot"))
                          ),
                          plotOutput("multi_enrichment6"),
                          fluidRow(
                            column(4, downloadButton("download_multi_enrichment_table2", "Download enrichment result"))
                          ),
                          dataTableOutput('multi_enrichment_result2')
                 ),
                 tabPanel("GSEA",
                          fluidRow(
                            column(3, htmlOutput("selectEnrich_pair")),
                            column(4, htmlOutput("Gene_set6")),
                            column(4, downloadButton("download_multi_enrichment", "Download"))
                          ),
                          fluidRow(
                            column(4, textOutput("multi_Spe1"),
                                   tags$head(tags$style("#multi_Spe1{color: red;
                                         font-size: 20px;
                                         font-style: bold;
                                      }")))
                          ),
                          plotOutput("multi_enrichment1"),
                          bsCollapse(id="input_collapse_multi_enrich",open="ORA_panel",multiple = TRUE,
                                     bsCollapsePanel(title="GSEA result:",
                                                     value="multi_GSEA_panel",
                                                     fluidRow(
                                                       column(4, downloadButton("download_multi_GSEA_table", "Download GSEA result"))
                                                     ),
                                                     dataTableOutput('multi_GSEA_result')
                                     )
                          )
                 )
               )
             ) # main panel
           ) #sidebarLayout
  ), #tabPanel
  # venn diagram analysis---------------------------------
  tabPanel("Venn diagram",
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
                          fluidRow(
                            column(4, downloadButton("download_vennplot", "Download venn diagram"))
                          ),
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
                                fileInput("file7",
                                          strong(
                                            span("Select a normalized count matrix file (txt, csv)"),
                                            span(icon("info-circle"), id = "icon11", 
                                                 options = list(template = popoverTempate))
                                          ),
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%"),
                                bsPopover("icon11", "Count matrix format:", 
                                          content=paste(strong("The replication number"), "is represented by", strong("the underline"),".<br>", strong("Do not use it for anything else"),".<br><br>",
                                                        img(src="input_format1.png", width = 400,height = 250)), 
                                          placement = "right",options = list(container = "body")),
               ),
               conditionalPanel(condition="input.data_file_type3=='Row6'",
                                fileInput("file8",
                                          strong(
                                            span("Select a normalized count matrix file (txt, csv)"),
                                            span(icon("info-circle"), id = "icon12", 
                                                 options = list(template = popoverTempate))
                                          ),
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%"),
                                fileInput("file9",
                                          label = "Select a metadata file to define samples for the following analysis",
                                          accept = c("txt", "csv"),
                                          multiple = FALSE,
                                          width = "80%"),
                                bsPopover("icon12", "Metadata format:", 
                                          content=paste("The first column is", strong("the sample name"), "used in the raw count data.<br>", 
                                                        "The second column is", strong("the corresponding sample name"), "that matches the sample name in the first column.<br><br>",
                                                        img(src="input_format2.png", width = 400,height = 400)), 
                                          placement = "right",options = list(container = "body")),
               ),
               fluidRow(
                 column(6, selectInput("Species3", "Species", c("not selected", "Homo sapiens", "Mus musculus", "Rattus norvegicus", "Xenopus laevis",
                                                                "Drosophila melanogaster", "Caenorhabditis elegans"), selected = "not selected")),
                 column(6, selectInput("normHeat", "Heatmap in clustering panel", c("OFF", "Draw"), selected = "OFF"))),
               h4("Filter option 1:"),
               fileInput("file10",
                         label = "Select a gene list file for gene extraction",
                         accept = c("txt", "csv", "xlsx"),
                         multiple = FALSE,
                         width = "80%"),
               h4("Filter option 2:"),
               fluidRow(
                 column(4, numericInput("basemean3", "Basemean", min   = 0, max   = NA, value = 0),
                       ),
                 column(8,  "If the gene number is >10,000 after filtering, heatmap and k-means clustering are not recommended.", br(),
                        "This is due to a limitation of the server memory.")
               ),
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
                                                     column(4, downloadButton("download_d_norm_count", "Download defined normalized count"))
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
                                 htmlOutput("kmeans_cv"),
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
  # enrichment viewer -------------------------------------
  tabPanel("Enrichment viewer",
           sidebarLayout(
             # enrichment viewer---------------------------------
             sidebarPanel(
               fileInput("enrich_data_file",
                         strong(
                           span("Select a gene list file (txt, csv, xlsx)"),
                           span(icon("info-circle"), id = "icon13", 
                                options = list(template = popoverTempate))
                         ),
                         accept = c("txt", "csv", "xlsx"),
                         multiple = FALSE,
                         width = "80%"),
               bsPopover("icon13", "Gene list format: ", 
                         content=paste("First column must be", strong("gene name"), "(Gene symbol or ENSEMBL ID).<br>", "Second column must be", strong("group or cluster name"),".<br>",
                                       "You can use result files of venn diagram analysis and k-means clustering as input.<br><br>", 
                                       img(src="input_format_enrich.png", width = 250,height = 400)), 
                         placement = "right",options = list(container = "body")),
               fluidRow(
                 column(6, selectInput("Species4", "Species", c("not selected", "Homo sapiens", "Mus musculus", "Rattus norvegicus", "Xenopus laevis",
                                                                "Drosophila melanogaster", "Caenorhabditis elegans"), selected = "not selected"))),
               actionButton("goButton4", "example data (mouse)"),
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
                 tabPanel("Input gene list for enrichment analysis",
                          dataTableOutput('enrichment_input')
                 ),
                 tabPanel("Enrichment analysis",
                          fluidRow(
                            column(4, textOutput("Spe3"),
                                   tags$head(tags$style("#Spe3{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                            column(4, htmlOutput("Gene_set3")),
                            column(4, downloadButton("download_enrichment", "Download dot plot"))
                          ),
                          plotOutput("enrichment3"),
                          fluidRow(
                            column(4, htmlOutput("whichGroup")),
                            column(4, downloadButton("download_enrichment_cnet", "Download cnet plot"))
                          ),
                          plotOutput("enrichment4"),
                          fluidRow(
                            column(4, downloadButton("download_enrichment_table", "Download enrichment result"))
                          ),
                          dataTableOutput('enrichment_result')
                 )
               )
             )
           ) #sidebarLayout
  ),
  #Instruction--------------------------
  navbarMenu("More",
             tabPanel("Volcano navi",
                      sidebarLayout(
                        # volcano navi---------------------------------
                        sidebarPanel(
                          fileInput("deg_file1",
                                    strong(
                                      span("Select a pair-wise DEG result file"),
                                      span(icon("info-circle"), id = "icon14", 
                                           options = list(template = popoverTempate))
                                    ),
                                    accept = c("txt", "csv"),
                                    multiple = FALSE,
                                    width = "80%"),
                          bsPopover("icon14", "File format: ", 
                                    content=paste("First column must be gene name (Gene symbol or ENSEMBL ID).<br>", 
                                                  "The file must contain", strong("log2FoldChange"), "and", strong("padj"), "columns.<br>",
                                                  "You can use a pair-wise DEG result file as input.<br><br>", 
                                                  img(src="input_format_volcano.png", width = 480,height = 230)), 
                                    placement = "right",options = list(container = "body")),
                          radioButtons('volcano_inputType','Reverse number signs of log2FoldChange:',
                                       c('ON'="reverseON",
                                         'OFF'="reverseOFF"
                                       ),selected = "reverseON"),
                          fileInput("deg_file2",
                                    label = "Option: select a normalized count file",
                                    accept = c("txt", "csv"),
                                    multiple = FALSE,
                                    width = "80%"),
                          fluidRow(
                            column(6, selectInput("Species5", "Species", c("not selected", "human", "mouse", "rat", "fly", "worm"), selected = "not selected"))
                          ),
                          fluidRow(
                            column(4, numericInput("fc4", "Fold Change", min   = 1, max   = NA, value = 2)),
                            column(4, numericInput("fdr4", "FDR", min   = 0, max   = 1, value = 0.05))
                          ),
                          actionButton("goButton5", "example data (mouse)"),
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
                            tabPanel("Input DEG result",
                                     bsCollapse(id="DEG_input_collapse_panel",open="DEG_panel",multiple = FALSE,
                                                bsCollapsePanel(title="DEG_result:",
                                                                value="DEG_panel",
                                                                dataTableOutput('DEG_input1')
                                                ),
                                                bsCollapsePanel(title="Option: normalized_count_matrix:",
                                                                value="normalized_count_matrix_panel",
                                                                dataTableOutput('norm_count_forDEG')
                                                )
                                     )
                            ),
                            tabPanel("GOI profiling",
                                     fluidRow(
                                       column(4, downloadButton("download_volcano_navi", "Download volcano plot")),
                                       column(4, downloadButton("download_deg_heatmap", "Download heatmap"))
                                     ),
                                     fluidRow(
                                       column(4, htmlOutput("degGOI")),
                                       column(4, htmlOutput("deg_volcano_x")),
                                       column(4, htmlOutput("deg_volcano_y"))
                                     ),
                                     fluidRow(
                                       column(8, plotOutput("deg_volcano1")),
                                       column(4, plotOutput("deg_GOIheatmap"))
                                     ),
                                     div(
                                       plotOutput("deg_GOIbox", height = "100%"),
                                       style = "height: calc(100vh  - 100px)"
                                     ),
                                     fluidRow(
                                       column(4, downloadButton("download_deg_GOIbox", "Download boxplot"))
                                     ))
                          )
                        )
                      ) #sidebarLayout
             ),
             tabPanel("MSigDB gene list",
                      sidebarLayout(
                        # MSigDB gene list---------------------------------
                        sidebarPanel(
                          selectInput("msigdbr_Species", "Species", c("", msigdbr_species), selected = "not selected"),
                          selectizeInput("msigdbr_gene_set", label="Name of gene set", choices = '')
                          #sidebarPanel
                        ),
                        
                        # Main Panel -------------------------------------
                        mainPanel(
                          fluidRow(
                            column(4, downloadButton("download_msigdbr_list", "Download"))
                          ),
                          dataTableOutput('msigdbr_geneset')
                        )
                      ) #sidebarLayout
             ),
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
                               "- Pantano L (2022). DEGreport: Report of DEG analysis. R package version 1.32.0, http://lpantano.github.io/DEGreport", br(),
                               "- Andrie de Vries and Brian D. Ripley (2020). ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'. R package version 0.1.22. https://CRAN.R-project.org/package=ggdendro",br(),
                               "- T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141",br(),
                               "- Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015 31(4):608-609",br(),
                               "- Dolgalev I (2022). _msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format_. R
  package version 7.5.1, <https://CRAN.R-project.org/package=msigdbr>.", br(),
                               "- Hervé Pagès, Marc Carlson, Seth Falcon and Nianhua Li (2020). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.52.0. https://bioconductor.org/packages/AnnotationDbi",br(),
                               "- Marc Carlson (2020). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.12.0.",br(),
                               "- Marc Carlson (2020). org.Mm.eg.db: Genome wide annotation for Mouse. R package version 3.12.0.",br(),
                               "- Marc Carlson (2020). org.Rn.eg.db: Genome wide annotation for Human. R package version 3.12.0.",br(),
                               "- Marc Carlson (2020). org.Xl.eg.db: Genome wide annotation for Human. R package version 3.12.0.",br(),
                               "- Marc Carlson (2020). org.Dm.eg.db: Genome wide annotation for Mouse. R package version 3.12.0.",br(),
                               "- Marc Carlson (2020). org.Ce.eg.db: Genome wide annotation for Human. R package version 3.12.0.",br(),
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
# server ---------------------------------
server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  # pair-wise ------------------------------------------------------------------------------
org1 <- reactive({
  if(input$Species != "not selected"){
    switch (input$Species,
            "Mus musculus" = org <- org.Mm.eg.db,
            "Homo sapiens" = org <- org.Hs.eg.db,
            "Rattus norvegicus" = org <- org.Rn.eg.db,
            "Xenopus laevis" = org <- org.Xl.eg.db,
            "Drosophila melanogaster" = org <- org.Dm.eg.db,
            "Caenorhabditis elegans" = org <- org.Ce.eg.db)
    return(org)
  }
})
org_code1 <- reactive({
  if(input$Species != "not selected"){
    switch (input$Species,
            "Mus musculus" = org_code <- "mmu",
            "Homo sapiens" = org_code <- "hsa",
            "Rattus norvegicus" = org_code <- "rno",
            "Xenopus laevis" = org_code <- "xla",
            "Drosophila melanogaster" = org_code <- "dme",
            "Caenorhabditis elegans" = org_code <- "cel")
    return(org_code)
  }
})

row_count_matrix <- reactive({
  withProgress(message = "Importing row count matrix, please wait",{
    if (input$data_file_type == "Row1"){
      tmp <- input$file3$datapath
      if(is.null(input$file3) && input$goButton > 0 )  tmp = "data/example1.txt"
      if(is.null(tmp)) {
        return(NULL)
      }else{
        if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
        if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
        if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
        return(df)
      }
    }
    if (input$data_file_type == "Row2"){
      tmp <- input$file1$datapath
      if(is.null(input$file1) && input$goButton == 0) return(NULL)
      if(is.null(input$file1) && input$goButton > 0 )  tmp = "data/example2.csv"
      if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
      if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
      if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
      return(df)
    }
  })
})
  metadata <- reactive({
    if (input$data_file_type != "Row2"){
      return(NULL)
    }else{
    tmp <- input$file2$datapath
    if(is.null(input$file2) && input$goButton == 0) return(NULL)
    if(is.null(input$file2) && input$goButton > 0 )  tmp = "data/example3.csv"
    if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
    if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
    if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
    rownames(df) <- gsub("-",".",rownames(df))
    return(df)
    }
  })
  norm_count_matrix <- reactive({
    if(is.null(input$norm_file1)){
      return(NULL)
    }else{
      if(length(input$norm_file1[, 1]) == 1){
        upload <- read.table(file = input$norm_file1[[1, 'datapath']], header = T ,sep="\t", row.names = 1)
      }else{
        upload = list()
        for(nr in 1:length(input$norm_file1[, 1])){
          df <- read.table(file = input$norm_file1[[nr, 'datapath']], header = T ,sep="\t", row.names = 1)
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
        return(row)
      }
    }else{
      meta <- metadata()
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
                           font.label = c("bold", 6),font.legend = "bold",
                           font.main = "bold",xlab = xlab,
                           ggtheme = ggplot2::theme_minimal(),
                           select.top.method = "fc"))
    data2 <- data_degcount2()
    if(is.null(data2)){
      ht <- NULL
    }else{
      data.z <- genescale(data2[,8:(7 + Cond_1 + Cond_2)], axis=1, method="Z")
      ht <- as.grob(Heatmap(data.z, name = "z-score",column_order = colnames(data.z),
                            clustering_method_columns = 'ward.D2',
                            show_row_names = F, show_row_dend = F))
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
        pdf(file, height = 4, width = 7)
        print(ma_heatmap_plot())
        dev.off()
        incProgress(1)
      })
    }
  )



  # pair-wise volcano--------------
  GOI_list <- reactive({
    withProgress(message = "Preparing GOI list",{
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
    incProgress(1)
    })
  })

  output$GOI <- renderUI({
    if(is.null(d_row_count_matrix())){
      return(NULL)
    }else{
    selectizeInput("GOI", "genes of interest (GOI)", c(GOI_list()),multiple = TRUE, options = list(delimiter = " ", create = T))
    }
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
            axis.text.x= ggplot2::element_text(size = 10),
            axis.text.y= ggplot2::element_text(size = 10),
            text = ggplot2::element_text(size = 10),
            title = ggplot2::element_text(size = 10)) +
      xlab("log2 fold change") + ylab("-log10(padj)") +
      xlim(input$xrange)+
      ylim(c(0, input$yrange))
    if(!is.null(label_data)) {
      if(str_detect(rownames(count)[1], "ENS")){
        if(input$Species != "not selected"){
          v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
          v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Unique_ID),
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1)
        }else{
          v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
          v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Row.names),
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1)
        }
      }else{
        v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
        v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Row.names),
                                          box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1)
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
        pdf(file, height = 5, width = 5)
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
      ht <- Heatmap(data.z, name = "z-score",column_order = colnames(data.z),
                    clustering_method_columns = 'ward.D2',
                    show_row_names = T, show_row_dend = F)
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
        pdf(file, height = 10, width = 7)
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
      collist <- gsub("\\_.+$", "", colnames(data3))
      collist <- unique(collist)
      data3$Row.names <- rownames(data3)
      data3 <- data3 %>% gather(key=sample, value=value,-Row.names)
      data3$sample <- gsub("\\_.+$", "", data3$sample)
      data3$Row.names <- as.factor(data3$Row.names)
      data3$sample <- factor(data3$sample,levels=collist,ordered=TRUE)
      data3$value <- as.numeric(data3$value)
      p <- ggpubr::ggboxplot(data3, x = "sample", y = "value",
                             fill = "sample", scales = "free",
                             add = "jitter",
                             xlab = FALSE, ylab = "Normalized_count", ylim = c(0, NA))
      p <- (facet(p, facet.by = "Row.names",
                  panel.labs.background = list(fill = "transparent", color = "transparent"),
                  scales = "free", short.panel.labs = T)+
              theme(axis.text.x= element_text(size = 0),
                    axis.text.y= element_text(size = 10),
                    panel.background = element_rect(fill = "transparent", size = 0.5),
                    title = element_text(size = 10),text = element_text(size = 20))
            + scale_fill_manual(values=c("gray", "#ff8082")))
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
        if ((length(rowlist) > 81) && (length(rowlist) <= 200))
        {pdf_hsize <- 22.5
        pdf_wsize <- 22.5}
        if ((length(rowlist) > 64) && (length(rowlist) <= 81))
        {pdf_hsize <- 20.25
        pdf_wsize <- 20.25}
        if ((length(rowlist) > 49) && (length(rowlist) <= 64))
        {pdf_hsize <- 18
        pdf_wsize <- 18}
        if ((length(rowlist) > 36) && (length(rowlist) <= 49))
        {pdf_hsize <- 15.75
        pdf_wsize <- 15.75}
        if ((length(rowlist) > 25) && (length(rowlist) <= 36))
        {pdf_hsize <- 13.5
        pdf_wsize <- 13.5}
        if ((length(rowlist) > 16) && (length(rowlist) <= 25))
        {pdf_hsize <- 11.5
        pdf_wsize <- 11.5}
        if ((length(rowlist) > 12) && (length(rowlist) <= 16))
        {pdf_hsize <- 9
        pdf_wsize <- 9}
        if ((length(rowlist) > 9) && (length(rowlist) <= 12))
        {pdf_hsize <- 7.5
        pdf_wsize <- 9}
        if ((length(rowlist) > 6) && (length(rowlist) <= 9))
        {pdf_hsize <- 7.5
        pdf_wsize <- 6.75}
        if ((length(rowlist) > 4) && (length(rowlist) <= 6))
        {pdf_hsize <- 6
        pdf_wsize <- 9}
        if (length(rowlist) == 4)
        {pdf_hsize <- 6
        pdf_wsize <- 6}
        if (length(rowlist) == 3)
        {pdf_hsize <- 3
        pdf_wsize <- 9}
        if (length(rowlist) == 2)
        {pdf_hsize <- 3
        pdf_wsize <- 6}
        if (length(rowlist) == 1)
        {pdf_hsize <- 3
        pdf_wsize <- 3}
        if (length(rowlist) > 200)
        {pdf_hsize <- 30
        pdf_wsize <- 30}
        pdf(file, height = pdf_hsize, width = pdf_wsize)
        print(pair_GOIbox())
        dev.off()
        incProgress(1)
      })
    }
  )

  # pair-wise PCA ------------------------------------------------------------------------------
  pairPCAdata <- reactive({
    if(is.null(d_row_count_matrix())){
      return(NULL)
    }else{
    data <- deg_norm_count()
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
  })

  pair_pca_plot <- reactive({
    data <- deg_norm_count()
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
  })

  output$download_pair_PCA = downloadHandler(
    filename = function(){
      paste0(download_pair_overview_dir(), "_PCA-MDS-dendrogram.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        pdf(file, height = 3.5, width = 9)
        print(pair_pca_plot())
        dev.off()
        incProgress(1)
      })
    }
  )

  output$PCA <- renderPlot({
    if(is.null(d_row_count_matrix())){
      return(NULL)
    }else{
    print(pair_pca_plot())
    }
  })

  output$pair_PCA_data <- DT::renderDataTable({
    pairPCAdata()
  })

  output$download_pair_PCA_table = downloadHandler(
    filename = function() {
      paste0(download_pair_overview_dir(), "_PCA_table.txt")
      },
    content = function(file){write.table(pairPCAdata(), file, row.names = T, sep = "\t", quote = F)}
  )
  # pair-wise enrichment ------------------------------------------------------------------------------
  output$Spe1 <- renderText({
    if(input$Species == "not selected") print("Please select 'Species'")
  })
  Hallmark_set <- reactive({
    if(input$Species != "not selected"){
    if(input$Gene_set != "MSigDB Hallmark" && input$Gene_set != "Transcription factor targets"){
      return(NULL)
      }else{
    switch (input$Species,
            "Mus musculus" = species <- "Mus musculus",
            "Homo sapiens" = species <- "Homo sapiens",
            "Rattus norvegicus" = species <- "Rattus norvegicus",
            "Xenopus laevis" = species <- "Xenopus laevis",
            "Drosophila melanogaster" = species <- "Drosophila melanogaster",
            "Caenorhabditis elegans" = species <- "Caenorhabditis elegans")
        if(input$Gene_set == "MSigDB Hallmark"){
          H_t2g <- msigdbr(species = species, category = "H") %>%
            dplyr::select(gs_name, entrez_gene) 
        }
        if(input$Gene_set == "Transcription factor targets"){
    H_t2g <- msigdbr(species = species, category = "C3")
    H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "TFT:GTRD" | gs_subcat == "TFT:TFT_Legacy") %>%
      dplyr::select(gs_name, entrez_gene)
        }
    return(H_t2g)
      }
    }else return(NULL)
  })

  enrichment_1_1 <- reactive({
    if(!is.null(input$Gene_set) && input$Species != "not selected"){
    data3 <- data_degcount2()
    if(input$Gene_set != "MSigDB Hallmark" && input$Gene_set != "Transcription factor targets"){
        if(input$Gene_set == "KEGG"){
          withProgress(message = "KEGG enrichment analysis",{
          formula_res <- try(compareCluster(ENTREZID~group, data=data3,
                                            fun="enrichKEGG", organism=org_code1()), silent = T)
          incProgress(1)
          })
        }
        if(input$Gene_set == "GO"){
          withProgress(message = "GO enrichment analysis",{
          formula_res <- try(compareCluster(ENTREZID~group, data=data3,
                                            fun="enrichGO", OrgDb=org1()), silent =T)
          incProgress(1)
          })
        }
        if (class(formula_res) == "try-error") {
          formula_res <- NULL
          return(formula_res)
        }else{
          if (length(as.data.frame(formula_res)$ID) == 0) {
            formula_res <- NULL
          }else{
          formula_res <- setReadable(formula_res, org1(), 'ENTREZID')
          }
          return(formula_res)
        }
      }else{
        withProgress(message = "enrichment analysis",{
        H_t2g <- Hallmark_set()
        em_up <- try(enricher(dplyr::filter(data3, group == "Up")$ENTREZID, TERM2GENE=H_t2g, pvalueCutoff = 0.05))
        em_down <- try(enricher(dplyr::filter(data3, group == "Down")$ENTREZID, TERM2GENE=H_t2g, pvalueCutoff = 0.05))
        if (class(em_up) == "try-error") em_up <- NA
        if (class(em_down) == "try-error") em_down <- NA
        if (length(as.data.frame(em_up)$ID) == 0) {
          em_up <- NA
        } else{
          em_up <- setReadable(em_up, org1(), 'ENTREZID')
        }
        if (length(as.data.frame(em_down)$ID) == 0) {
          em_down <- NA
        } else{
          em_down <- setReadable(em_down, org1(), 'ENTREZID')
        }
        if ((length(as.data.frame(em_up)$ID) == 0) && (length(as.data.frame(em_down)$ID) == 0))  {
          return(NULL)
        } else{
          group1 <- as.data.frame(em_up)
          group2 <- as.data.frame(em_down)
          if ((length(colnames(group1)) != 1) && (length(colnames(group2)) != 1))  {
            group1$Cluster <- "Up"
            group2$Cluster <- "Down"
            data<- rbind(group1, group2)
          }
          if((length(colnames(group1)) == 1) && (length(colnames(group2)) != 1)){
            group2$Cluster <- "Down"
            data <- group2
          }
          if((length(colnames(group1)) != 1) && (length(colnames(group2)) == 1)){
            group1$Cluster <- "Up"
            data <- group1
          }
          if((length(colnames(group1)) == 1) && (length(colnames(group2)) == 1)){
            data <- NA
          }
          if(length(colnames(data)) != 0) {
          data["Description"] <- lapply(data["Description"], gsub, pattern="HALLMARK_", replacement = "")
          data$GeneRatio <- parse_ratio(data$GeneRatio)
          return(data)
          }else return(NULL)
        }
        incProgress(1)
        })
        }
      }else{return(NULL)}
  })

  enrichment_1_gsea <- reactive({
    if(!is.null(input$Gene_set) && input$Species != "not selected"){
    data <- data_degcount()
    data3 <- data_degcount2()
    count <- deg_norm_count()
    data <- na.omit(data)
    geneList <- data$log2FoldChange
    names(geneList) = as.character(data$ENTREZID)
    geneList <- sort(geneList, decreasing = TRUE)
    if(input$Gene_set != "MSigDB Hallmark" && input$Gene_set != "Transcription factor targets"){
        if(input$Gene_set == "KEGG"){
          withProgress(message = "KEGG GSEA",{
          kk3 <- gseKEGG(geneList = geneList, pvalueCutoff = 0.05,
                         organism = org_code1(), keyType = "kegg",
                         exponent = 1, eps = 0, pAdjustMethod = "BH",
                         minGSSize = 50, maxGSSize = 500, by = "fgsea",
                         use_internal_data = FALSE, verbose = F)
          })
        }
        if(input$Gene_set == "GO"){
          withProgress(message = "GO GSEA takes a few minutes",{
          kk3 <- gseGO(geneList = geneList, pvalueCutoff = 0.05,
                       OrgDb = org1(), exponent = 1, eps = 0,
                       pAdjustMethod = "BH", minGSSize = 50,
                       maxGSSize = 500, by = "fgsea", verbose = F)
          incProgress(1)
          })
        }
        if (length(as.data.frame(kk3)$ID) == 0) {
          kk3 <- NULL
        } else{
          kk3 <- setReadable(kk3, org1(), 'ENTREZID')
        }
        return(kk3)
      }else{
        withProgress(message = "GSEA",{
        H_t2g <- Hallmark_set()
        em3 <- GSEA(geneList, TERM2GENE = H_t2g,pvalueCutoff = 0.05,exponent = 1, eps = 0, pAdjustMethod = "BH",
                    minGSSize = 50, maxGSSize = 500,by = "fgsea",verbose = F)
        if (length(as.data.frame(em3)$ID) == 0) {
          em4 <- NA
        } else{
          em4 <- setReadable(em3, org1(), 'ENTREZID')
        }
        return(em4)
        incProgress(1)
        })
      }
    }else return(NULL)
  })


  # pair-wise enrichment plot ------------------------------------------------------------------------------
  pair_enrich1_keggGO <- reactive({
    if(!is.null(input$Gene_set) && input$Species != "not selected"){
      if(input$Gene_set != "MSigDB Hallmark" && input$Gene_set != "Transcription factor targets"){
      count <- deg_norm_count()
      formula_res <- enrichment_1_1()
      if ((length(as.data.frame(formula_res)$Description) == 0) ||
          length(which(!is.na(unique(as.data.frame(formula_res)$qvalue)))) == 0) {
        p1 <- NULL
      } else{
        formula_res <- clusterProfiler.dplyr::filter(formula_res, !is.na(qvalue))
        p1 <- as.grob(dotplot(formula_res, color ="qvalue", font.size = 10))
      }
      kk4 <- enrichment_1_gsea()
      if (length(as.data.frame(kk4)$ID) == 0) {
        p4 <- NULL
      } else{
        if (length(as.data.frame(kk4)$ID) >= 5){
          p4 <- as.grob(gseaplot2(kk4, 1:5, pvalue_table = F))
        }else{
          p4 <- as.grob(gseaplot2(kk4, 1:length(kk4$ID), pvalue_table = F))
        }
      }
      p <- plot_grid(p1, p4, nrow = 1)
      return(p)
    }
    }else return(NULL)
  })

  pair_enrich1_H <- reactive({
    if(!is.null(input$Gene_set) && input$Species != "not selected"){
    if(input$Gene_set == "MSigDB Hallmark" || input$Gene_set == "Transcription factor targets"){
      count <- deg_norm_count()
      formula_res <- enrichment_1_1()
      data3 <- data_degcount2()
      H_t2g <- Hallmark_set()
      em_up <- enricher(dplyr::filter(data3, group == "Up")$ENTREZID, TERM2GENE=H_t2g, pvalueCutoff = 0.05)
      em_down <- enricher(dplyr::filter(data3, group == "Down")$ENTREZID, TERM2GENE=H_t2g, pvalueCutoff = 0.05)
      if (length(as.data.frame(em_up)$ID) == 0) {
        em_up <- NA
      } else{
        em_up <- setReadable(em_up, org1(), 'ENTREZID')
      }
      if (length(as.data.frame(em_down)$ID) == 0) {
        em_down <- NA
      } else{
        em_down <- setReadable(em_down, org1(), 'ENTREZID')
      }
      if ((length(em_up) == 0) && (length(em_down) == 0))  {
        p1 <- NULL
      } else{
        group1 <- as.data.frame(em_up)
        group2 <- as.data.frame(em_down)
        group1$Cluster <- "Up"
        group2$Cluster <- "Down"
        if(length(as.data.frame(em_up)$ID) != 0){
          if (length(group1$pvalue) > 5){
            group1 <- group1[sort(group1$pvalue, decreasing = F, index=T)$ix,]
            group1 <- group1[1:5,]
          }
        }
        if(length(as.data.frame(em_down)$ID) != 0){
          if (length(group2$pvalue) > 5){
            group2 <- group2[sort(group2$pvalue, decreasing = F, index=T)$ix,]
            group2 <- group2[1:5,]
          }
        }
        if ((length(colnames(group1)) != 2) && (length(colnames(group2)) != 2))  {
          data<- rbind(group1, group2)
        }
        if((length(colnames(group1)) == 2) && (length(colnames(group2)) != 2)){
          data <- group2
        }
        if((length(colnames(group1)) != 2) && (length(colnames(group2)) == 2)){
          data <- group1
        }
        if((length(colnames(group1)) == 2) && (length(colnames(group2)) == 2)){
          data <- NA
        }
        if(length(colnames(data)) != 0){
          data["Description"] <- lapply(data["Description"], gsub, pattern="HALLMARK_", replacement = "")
          data$GeneRatio <- parse_ratio(data$GeneRatio)
        if ((length(data$Description) == 0) || length(which(!is.na(unique(data$qvalue))))==0) {
          p1 <- NULL
        } else{
          p1 <- as.grob(ggplot(data, aes(x = Cluster,y=reorder(Description, GeneRatio)))+
                          geom_point(aes(color=qvalue,size=GeneRatio)) +
                          scale_color_continuous(low="red", high="blue",
                                                 guide=guide_colorbar(reverse=TRUE)) +
                          scale_size(range=c(3, 8))+ theme_dose(font.size=8)+ylab(NULL))
        }}else p1 <- NULL
      }
      em3 <- enrichment_1_gsea()
      if (length(as.data.frame(em3)$ID) == 0) {
        p4 <- NULL
      } else{
        if (length(as.data.frame(em3)$ID) >= 5){
          p4 <- as.grob(gseaplot2(em3, 1:5, pvalue_table = F))
        }else{
          p4 <- as.grob(gseaplot2(em3, 1:length(em3$ID), pvalue_table = F))
        }
      }
      p <- plot_grid(p1, p4, nrow = 1)
      return(p)
    }
    }else return(NULL)
  })

  output$enrichment1 <- renderPlot({
    if(!is.null(input$Gene_set)){
      if(is.null(d_row_count_matrix())){
        return(NULL)
      }else{
        withProgress(message = "Plot results",{
    if(input$Species != "not selected"){
      if(input$Gene_set != "MSigDB Hallmark" && input$Gene_set != "Transcription factor targets"){
      print(pair_enrich1_keggGO())
      }else{
        print(pair_enrich1_H())
      }
    }
          incProgress(1)
})
      }
    }
  })

  pair_enrich2 <- reactive({
    if(!is.null(input$Gene_set) && input$Species != "not selected"){
    data <- data_degcount()
    data3 <- data_degcount2()
    count <- deg_norm_count()
    upgene <- data3[data3$log2FoldChange > log(input$fc, 2),]
    downgene <- data3[data3$log2FoldChange < log(1/input$fc, 2),]
    if(input$Gene_set != "MSigDB Hallmark" && input$Gene_set != "Transcription factor targets"){
        if(input$Gene_set == "KEGG"){
          kk1 <- enrichKEGG(upgene$ENTREZID, organism =org_code1(),
                            pvalueCutoff = 0.05)
          kk2 <- enrichKEGG(downgene$ENTREZID, organism =org_code1(),
                            pvalueCutoff = 0.05)
        }
        if(input$Gene_set == "GO"){
          kk1 <- enrichGO(upgene$ENTREZID, OrgDb = org1(),
                          pvalueCutoff = 0.05)
          kk2 <- enrichGO(downgene$ENTREZID, OrgDb = org1(),
                          pvalueCutoff = 0.05)
        }
    }else{
      H_t2g <- Hallmark_set()
      kk1 <- try(enricher(dplyr::filter(data3, group == "Up")$ENTREZID, TERM2GENE=H_t2g, pvalueCutoff = 0.05))
      kk2 <- try(enricher(dplyr::filter(data3, group == "Down")$ENTREZID, TERM2GENE=H_t2g, pvalueCutoff = 0.05))
      if (class(kk1) == "try-error") kk1 <- NA
      if (class(kk2) == "try-error") kk2 <- NA
    }
      if(length(as.data.frame(kk1)$ID) == 0){
        cnet1 <- NULL
      } else {
        cnet1 <- setReadable(kk1, org1(), 'ENTREZID')
      }
      if(length(as.data.frame(kk2)$ID) == 0){
        cnet2 <- NULL
      } else {
        cnet2 <- setReadable(kk2, org1(), 'ENTREZID')
      }
      if (length(as.data.frame(cnet1)$ID) == 0) {
        p2 <- NULL
      } else{
        geneList_up <- upgene$log2FoldChange
        names(geneList_up) = as.character(upgene$ENTREZID)
        p2 <- try(as.grob(cnetplot(cnet1, foldChange=geneList_up,
                                   cex_label_gene = 0.7, cex_label_category = 0.75,
                                   cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none")))
        if(length(class(p2)) == 1){
          if(class(p2) == "try-error") p2 <- NULL
        }else{p2 <- as.grob(cnetplot(cnet1, foldChange=geneList_up,
                                         cex_label_gene = 0.7, cex_label_category = 0.75,
                                         cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none"))}
      }
      if (length(as.data.frame(cnet2)$ID) == 0) {
        p3 <- NULL
      } else{
        geneList_down <- downgene$log2FoldChange
        names(geneList_down) = as.character(downgene$ENTREZID)
        p3 <- try(as.grob(cnetplot(cnet2, foldChange=geneList_down,
                                   cex_label_gene = 0.7, cex_label_category = 0.75,
                                   cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none")))
        if(length(class(p3)) == 1){
          if(class(p3) == "try-error") p3 <- NULL
        }else{p3 <- as.grob(cnetplot(cnet2, foldChange=geneList_down,
                                         cex_label_gene = 0.7, cex_label_category = 0.75,
                                         cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none"))}
      }
      p <- plot_grid(p2, p3, nrow = 1)
      return(p)
    }else return(NULL)
  })

  output$enrichment2 <- renderPlot({
    if(!is.null(input$Gene_set)){
    if(is.null(d_row_count_matrix())){
      return(NULL)
    }else{
    if(input$Species != "not selected"){
      withProgress(message = "cnet plot",{
          p <- pair_enrich2()
        print(p)
        incProgress(1)
      })
    }else return(NULL)
    }
    }
  })

  output$download_pair_enrichment = downloadHandler(
    filename = function(){
      paste(download_pair_overview_dir(), paste0(input$Gene_set,"-enrichment.pdf"), sep="_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$Gene_set != "MSigDB Hallmark" && input$Gene_set != "Transcription factor targets"){
          p1 <- pair_enrich1_keggGO()
        }else{
          p1 <- pair_enrich1_H()
        }
        p2 <- pair_enrich2()
        pdf(file, height = 10, width = 12)
        print(plot_grid(p1, p2, nrow =2))
        dev.off()
        incProgress(1)
      })
    }
  )

  output$Gene_set <- renderUI({
    selectInput('Gene_set', 'Gene Set', c("KEGG", "GO", "MSigDB Hallmark", "Transcription factor targets"))
  })

  output$pair_enrichment_result <- DT::renderDataTable({
    as.data.frame(enrichment_1_1())
  })

  output$pair_GSEA_result <- DT::renderDataTable({
    as.data.frame(enrichment_1_gsea())
  })

  output$download_pair_enrichment_table = downloadHandler(
    filename = function() {
      paste(download_pair_overview_dir(), paste0(input$Gene_set,"-enrichment.txt"), sep="_")
    },
    content = function(file){write.table(as.data.frame(enrichment_1_1()), file, row.names = F, sep = "\t", quote = F)}
  )
  output$download_pair_GSEA_table = downloadHandler(
    filename = function() {
      paste(download_pair_overview_dir(), paste0(input$Gene_set,"-GSEA.txt"), sep="_")
    },
    content = function(file){write.table(as.data.frame(enrichment_1_gsea()), file, row.names = F, sep = "\t", quote = F)}
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
        df["day0"] <- list(read.table("data/day0.txt",header = T, row.names = 1))
        df["day1"] <- list(read.table("data/day1.txt",header = T, row.names = 1))
        df["day5"] <- list(read.table("data/day5.txt",header = T, row.names = 1))
        return(df)
      }
      return(NULL)
    }else{
      for(nr in 1:length(input$file11[, 1])){
        df <- read.table(file = input$file11[[nr, 'datapath']], header = T ,sep="\t", row.names = 1)
        upload[gsub("\\..+$", "", input$file11[nr,]$name)] <- list(df)
      }
      return(upload)
    }
  })
  
  
  dds_batch <- reactive({
    count_files <- batch_files()
    count_list <- list()
    for (name in names(count_files)) {
      count <- count_files[[name]]
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
      count_list[name] <- list(dds)
    }
    return(count_list)
  })
  
  
  deg_result_batch <- reactive({
    count_files <- batch_files()
    count_list <- list()
    for (name in names(count_files)) {
      count <- count_files[[name]]
      collist <- gsub("\\_.+$", "", colnames(count))
      if (input$DEG_method == "DESeq2") {
        withProgress(message = "DESeq2",{
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
          incProgress(1)
        })
      }
      if (input$DEG_method == "edgeR") {
        withProgress(message = "edgeR",{
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
          incProgress(1)
        })
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
      count_list[name] <- list(res)
    }
    return(count_list)
  })
  
  
  
  deg_norm_count_batch <- reactive({
    if(!is.null(norm_count_matrix())){
      return(norm_count_matrix())
    }else{
      count_files <- batch_files()
      count_list <- list()
      for (name in names(count_files)) {
        count <- count_files[[name]]
        collist <- gsub("\\_.+$", "", colnames(count))
        group <- data.frame(con = factor(collist))
        if (input$DEG_method == "DESeq2") {
          dds <- dds_batch()[[name]]
          contrast <- c("con", unique(collist))
          normalized_counts <- counts(dds, normalized=TRUE)
        }
        if (input$DEG_method == "edgeR") {
          dds <- dds_batch()[[name]]
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
        count_list[name] <- list(normalized_counts)
      }
      return(count_list)
    }
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
                                show_row_names = F, show_row_dend = F))
        }
        p <- plot_grid(m1, ht, rel_widths = c(2, 1))
        malist[[name]] <- p
      }
      return(malist)
    }
  })
  
  pair_pca_plot_batch <- reactive({
    data1 <- deg_norm_count_batch()
    pca_list <- list()
    for (name in names(data1)) {
      data <- data1[[name]]
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
        for(name in names(files)){
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
          write.table(deg, DEG, quote = F, row.names = T, sep = "\t")
          write.table(norm, norm_count, quote = F, row.names = T, sep = "\t")
          write.table(up, UP, quote = F, row.names = T, sep = "\t")
          write.table(down, DOWN, quote = F, row.names = T, sep = "\t")
          pdf(PCA, height = 3.5, width = 9)
          print(pca_r[[name]])
          dev.off()
          pdf(MA, height = 4, width = 7)
          print(ma_r[[name]])
          dev.off()
        }
        zip(zipfile=fname, files=fs)
      })
    },
    contentType = "application/zip"
  )
  
  #multi DEG--------------------------
  org6 <- reactive({
    if(input$Species6 != "not selected"){
      switch (input$Species6,
              "Mus musculus" = org <- org.Mm.eg.db,
              "Homo sapiens" = org <- org.Hs.eg.db,
              "Rattus norvegicus" = org <- org.Rn.eg.db,
              "Xenopus laevis" = org <- org.Xl.eg.db,
              "Drosophila melanogaster" = org <- org.Dm.eg.db,
              "Caenorhabditis elegans" = org <- org.Ce.eg.db)
      return(org)
    }
  })
  org_code6 <- reactive({
    if(input$Species6 != "not selected"){
      switch (input$Species6,
              "Mus musculus" = org_code <- "mmu",
              "Homo sapiens" = org_code <- "hsa",
              "Rattus norvegicus" = org_code <- "rno",
              "Xenopus laevis" = org_code <- "xla",
              "Drosophila melanogaster" = org_code <- "dme",
              "Caenorhabditis elegans" = org_code <- "cel")
      return(org_code)
    }
  })
  
  multi_row_count_matrix <- reactive({
    withProgress(message = "Importing row count matrix, please wait",{
      if (input$multi_data_file_type == "Row1"){
        tmp <- input$multi_file1$datapath
        if(is.null(input$multi_file1) && input$goButton6 > 0 )  tmp = "data/example4.txt"
        if(is.null(tmp)) {
          return(NULL)
        }else{
          if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
          if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
          if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
          return(df)
        }
      }
      if (input$multi_data_file_type == "Row2"){
        tmp <- input$multi_file2$datapath
        if(is.null(input$multi_file2) && input$goButton6 == 0) return(NULL)
        if(is.null(input$multi_file2) && input$goButton6 > 0 )  tmp = "data/example8.txt"
        if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
        if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
        if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
        return(df)
      }
    })
  })
  multi_metadata <- reactive({
    if (input$multi_data_file_type != "Row2"){
      return(NULL)
    }else{
      tmp <- input$multi_file3$datapath
      if(is.null(input$multi_file3) && input$goButton6 == 0) return(NULL)
      if(is.null(input$multi_file3) && input$goButton6 > 0 )  tmp = "data/example5.csv"
      if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
      if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
      if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
      rownames(df) <- gsub("-",".",rownames(df))
      df[,1] <- gsub("\\_.+$", "", df[,1])
      df[,2] <- gsub("\\_.+$", "", df[,2])
      return(df)
    }
  })
  
  multi_d_row_count_matrix <- reactive({
    row <- multi_row_count_matrix()
    if (input$multi_data_file_type == "Row1"){
      if(is.null(row)) {
        return(NULL)
      }else{
        return(row)
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
        upload <- read.table(file = input$multi_norm_file1[[1, 'datapath']], header = T ,sep="\t", row.names = 1)
      }else{
        upload = list()
        for(nr in 1:length(input$multi_norm_file1[, 1])){
          df <- read.table(file = input$multi_norm_file1[[nr, 'datapath']], header = T ,sep="\t", row.names = 1)
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
      print(paste0("The DEG number after the filteration: ", length(multi_pattern1()$gene)))
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
                     selected = c(multi_select()[1], multi_select()[2]), multiple = TRUE, 
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
          clusters <- degPatterns(cluster_rlog, metadata = meta, time = "condition", plot = FALSE)
        }else{
          meta <- data.frame(condition=factor(meta[,1]), type=factor(meta[,2]))
          meta[,2] <- factor(meta[,2], levels = unique(meta[,2]),ordered = TRUE)
          rownames(meta) <- colnames(count)
          clusters <- degPatterns(cluster_rlog, metadata = meta, time = "condition", col=colnames(meta)[2], plot = FALSE)
        }
        
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
            theme_bw()+ theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
        }
        else{
          p <- degPlotCluster(table, time = "condition", color=colnames(meta)[2], process = TRUE)+ 
            scale_color_brewer(palette = "Set1", direction=-1)+
            theme_bw()+ theme(legend.position = "top")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
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
      selectInput("multi_selectfile1", "cluster_list", choices = c("not selected", unique(clusters$cluster)), selected = "not selected",multiple = F)
    }
  })
  
  multi_pattern_extract <- reactive({
    data <- as.data.frame(multi_deg_norm_count())
    count <- multi_d_row_count_matrix()
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
      clusters
    }
  })
  output$multi_pattern1_count <- DT::renderDataTable({
    multi_pattern_extract()
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
          if ((clusterNumber > 81) && (clusterNumber <= 200))
          {pdf_hsize <- 15
          pdf_wsize <- 15}
          if ((clusterNumber > 64) && (clusterNumber <= 81))
          {pdf_hsize <- 14
          pdf_wsize <- 13.5}
          if ((clusterNumber > 49) && (clusterNumber <= 64))
          {pdf_hsize <- 13
          pdf_wsize <- 12}
          if ((clusterNumber > 36) && (clusterNumber <= 49))
          {pdf_hsize <- 11
          pdf_wsize <- 10.5}
          if ((clusterNumber > 25) && (clusterNumber <= 36))
          {pdf_hsize <- 10
          pdf_wsize <- 9}
          if ((clusterNumber > 16) && (clusterNumber <= 25))
          {pdf_hsize <- 8
          pdf_wsize <- 7.5}
          if ((clusterNumber > 12) && (clusterNumber <= 16))
          {pdf_hsize <- 7
          pdf_wsize <- 6}
          if ((clusterNumber > 9) && (clusterNumber <= 12))
          {pdf_hsize <- 7
          pdf_wsize <- 6}
          if ((clusterNumber > 6) && (clusterNumber <= 9))
          {pdf_hsize <- 6
          pdf_wsize <- 4.5}
          if ((clusterNumber > 4) && (clusterNumber <= 6))
          {pdf_hsize <- 6
          pdf_wsize <- 6}
          if (clusterNumber == 4)
          {pdf_hsize <- 6
          pdf_wsize <- 4}
          if (clusterNumber == 3)
          {pdf_hsize <- 4
          pdf_wsize <- 6}
          if (clusterNumber == 2)
          {pdf_hsize <- 4
          pdf_wsize <- 4}
          if (clusterNumber == 1)
          {pdf_hsize <- 4
          pdf_wsize <- 2}
          if (clusterNumber > 200)
          {pdf_hsize <- 30
          pdf_wsize <- 30}
          pdf(file, height = pdf_hsize, width = pdf_wsize)
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
                     selected = c(multi_select()[1], multi_select()[2]), multiple = TRUE, 
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
      print(paste0("The DEG number after the filteration: ", length(rownames(multi_deg_count1()))))
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
                      show_row_names = F)
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
            theme_bw()+ theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
        }
        else{
          p <- degPlotCluster(table, time = "condition", color=colnames(meta)[2], process = TRUE)+ 
            scale_color_brewer(palette = "Set1", direction=-1)+
            theme_bw()+ theme(legend.position = "top")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
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
          if ((clusterNumber > 81) && (clusterNumber <= 200))
          {pdf_hsize <- 15
          pdf_wsize <- 15}
          if ((clusterNumber > 64) && (clusterNumber <= 81))
          {pdf_hsize <- 14
          pdf_wsize <- 13.5}
          if ((clusterNumber > 49) && (clusterNumber <= 64))
          {pdf_hsize <- 13
          pdf_wsize <- 12}
          if ((clusterNumber > 36) && (clusterNumber <= 49))
          {pdf_hsize <- 11
          pdf_wsize <- 10.5}
          if ((clusterNumber > 25) && (clusterNumber <= 36))
          {pdf_hsize <- 10
          pdf_wsize <- 9}
          if ((clusterNumber > 16) && (clusterNumber <= 25))
          {pdf_hsize <- 8
          pdf_wsize <- 7.5}
          if ((clusterNumber > 12) && (clusterNumber <= 16))
          {pdf_hsize <- 7
          pdf_wsize <- 6}
          if ((clusterNumber > 9) && (clusterNumber <= 12))
          {pdf_hsize <- 7
          pdf_wsize <- 6}
          if ((clusterNumber > 6) && (clusterNumber <= 9))
          {pdf_hsize <- 6
          pdf_wsize <- 4.5}
          if ((clusterNumber > 4) && (clusterNumber <= 6))
          {pdf_hsize <- 6
          pdf_wsize <- 6}
          if (clusterNumber == 4)
          {pdf_hsize <- 6
          pdf_wsize <- 4}
          if (clusterNumber == 3)
          {pdf_hsize <- 4
          pdf_wsize <- 6}
          if (clusterNumber == 2)
          {pdf_hsize <- 4
          pdf_wsize <- 4}
          if (clusterNumber == 1)
          {pdf_hsize <- 4
          pdf_wsize <- 2}
          if (clusterNumber > 200)
          {pdf_hsize <- 30
          pdf_wsize <- 30}
          pdf(file, height = pdf_hsize, width = pdf_wsize)
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
    multi_kmeans_cluster()
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
        pdf(file, height = 10, width = 7)
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
      selectInput("multi_selectfile2", "cluster_list", choices = c("not selected", unique(clusters$Cluster)), selected = "not selected",multiple = F)
    }
  })
  
  multi_kmeans_pattern_extract <- reactive({
    count <- multi_d_row_count_matrix()
    clusters <- multi_kmeans_cluster()
    if(is.null(clusters)){
      return(NULL)
    }else{
      if(input$multi_selectfile2 == "not selected" || is.null(input$multi_selectfile2)){
        return(NULL)
      }else{
        cluster_name <- input$multi_selectfile2
        clusterCount <- dplyr::filter(clusters, Cluster == cluster_name)
        clusterCount <- clusterCount[,-1]
        return(clusterCount)
      }
    }
  })
  
  output$multi_pattern2_count <- DT::renderDataTable({
    multi_kmeans_pattern_extract()
  })
  
  output$download_deg_kmeans_pattern_count = downloadHandler(
    filename = function() {
      paste(paste(download_multi_overview_dir(),input$multi_selectfile2, sep = "_"), 
            "kmeans_count_table.txt", sep = "_")
    },
    content = function(file){write.table(multi_kmeans_pattern_extract(), file, row.names = T, sep = "\t", quote = F)}
  )
  
  #Multi DEG enrichment------------
  output$multi_Spe1 <- renderText({
    if(input$Species6 == "not selected") print("Please select 'Species'")
  })
  multi_Hallmark_set <- reactive({
    if(input$Species6 != "not selected"){
      if(input$Gene_set6 != "MSigDB Hallmark" && input$Gene_set6 != "Transcription factor targets"){
        return(NULL)
      }else{
        switch (input$Species6,
                "Mus musculus" = species <- "Mus musculus",
                "Homo sapiens" = species <- "Homo sapiens",
                "Rattus norvegicus" = species <- "Rattus norvegicus",
                "Xenopus laevis" = species <- "Xenopus laevis",
                "Drosophila melanogaster" = species <- "Drosophila melanogaster",
                "Caenorhabditis elegans" = species <- "Caenorhabditis elegans")
        if(input$Gene_set6 == "MSigDB Hallmark"){
          H_t2g <- msigdbr(species = species, category = "H") %>%
            dplyr::select(gs_name, entrez_gene) 
        }
        if(input$Gene_set6 == "Transcription factor targets"){
          H_t2g <- msigdbr(species = species, category = "C3")
          H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "TFT:GTRD" | gs_subcat == "TFT:TFT_Legacy") %>%
            dplyr::select(gs_name, entrez_gene)
        }
        return(H_t2g)
      }
    }else return(NULL)
  })
  
  output$selectEnrich_pair <- renderUI({
    if(is.null(multi_deg_result())){
      return(NULL)
    }else{
      selectizeInput("selectEnrich_pair", "Select a pair for GSEA", c(multi_select()),
                     selected = c(multi_select()[1], multi_select()[2]), multiple = TRUE, 
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
          my.symbols <- rownames(sig_res_LRT)
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
      if(input$Gene_set6 != "MSigDB Hallmark" && input$Gene_set6 != "Transcription factor targets"){
        if(input$Gene_set6 == "KEGG"){
          withProgress(message = "KEGG GSEA",{
            kk3 <- gseKEGG(geneList = geneList, pvalueCutoff = 0.05,
                           organism = org_code6(), keyType = "kegg",
                           exponent = 1, eps = 0, pAdjustMethod = "BH",
                           minGSSize = 50, maxGSSize = 500, by = "fgsea",
                           use_internal_data = FALSE, verbose = F)
          })
        }
        if(input$Gene_set6 == "GO"){
          withProgress(message = "GO GSEA takes a few minutes",{
            kk3 <- gseGO(geneList = geneList, pvalueCutoff = 0.05,
                         OrgDb = org6(), exponent = 1, eps = 0,
                         pAdjustMethod = "BH", minGSSize = 50,
                         maxGSSize = 500, by = "fgsea", verbose = F)
            incProgress(1)
          })
        }
        if (length(as.data.frame(kk3)$ID) == 0) {
          kk3 <- NULL
        } else{
          kk3 <- setReadable(kk3, org6(), 'ENTREZID')
        }
        return(kk3)
      }else{
        withProgress(message = "GSEA",{
          H_t2g <- multi_Hallmark_set()
          em3 <- GSEA(geneList, TERM2GENE = H_t2g,pvalueCutoff = 0.05,exponent = 1, eps = 0, pAdjustMethod = "BH",
                      minGSSize = 50, maxGSSize = 500,by = "fgsea",verbose = F)
          if (length(as.data.frame(em3)$ID) == 0) {
            em4 <- NA
          } else{
            em4 <- setReadable(em3, org6(), 'ENTREZID')
          }
          return(em4)
          incProgress(1)
        })
      }
    }else return(NULL)
  })
  
  
  # Multi DEG enrichment plot ------------------------------------------------------------------------------
  multi_enrich1_keggGO <- reactive({
    if(!is.null(input$Gene_set6) && input$Species6 != "not selected"){
      if(input$Gene_set6 != "MSigDB Hallmark" && input$Gene_set6 != "Transcription factor targets"){
        kk4 <- multi_enrichment_1_gsea()
        if (length(as.data.frame(kk4)$ID) == 0) {
          p4 <- NULL
        } else{
          if (length(as.data.frame(kk4)$ID) >= 5){
            p4 <- as.grob(gseaplot2(kk4, 1:5, pvalue_table = F))
          }else{
            p4 <- as.grob(gseaplot2(kk4, 1:length(kk4$ID), pvalue_table = F))
          }
        }
        p <- plot_grid(p4, nrow = 1)
        return(p)
      }
    }else return(NULL)
  })
  
  multi_enrich1_H <- reactive({
    if(!is.null(input$Gene_set6) && input$Species6 != "not selected"){
      if(input$Gene_set6 == "MSigDB Hallmark" || input$Gene_set6 == "Transcription factor targets"){
        H_t2g <- multi_Hallmark_set()
        em3 <- multi_enrichment_1_gsea()
        if (length(as.data.frame(em3)$ID) == 0) {
          p4 <- NULL
        } else{
          if (length(as.data.frame(em3)$ID) >= 5){
            p4 <- as.grob(gseaplot2(em3, 1:5, pvalue_table = F))
          }else{
            p4 <- as.grob(gseaplot2(em3, 1:length(em3$ID), pvalue_table = F))
          }
        }
        p <- plot_grid(p4, nrow = 1)
        return(p)
      }
    }else return(NULL)
  })
  
  output$multi_enrichment1 <- renderPlot({
    if(!is.null(input$Gene_set6)){
      if(is.null(multi_enrich_pairFC())){
        return(NULL)
      }else{
        withProgress(message = "Plot results",{
          if(input$Species6 != "not selected"){
            if(input$Gene_set6 != "MSigDB Hallmark" && input$Gene_set6 != "Transcription factor targets"){
              print(multi_enrich1_keggGO())
            }else{
              print(multi_enrich1_H())
            }
          }
          incProgress(1)
        })
      }
    }
  })
  
  output$Gene_set6 <- renderUI({
    selectInput('Gene_set6', 'Gene Set', c("KEGG", "GO", "MSigDB Hallmark", "Transcription factor targets"))
  })
  
  output$multi_GSEA_result <- DT::renderDataTable({
    if(is.null(dds) || length(input$selectEnrich_pair) != 2){
      return(NULL)
    }else{
      as.data.frame(multi_enrichment_1_gsea())
    }
  })
  
  output$download_multi_enrichment = downloadHandler(
    filename = function(){
      paste(download_multi_overview_dir(), paste0(input$Gene_set6,"-GSEA.pdf"), sep="_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$Gene_set6 != "MSigDB Hallmark" && input$Gene_set6 != "Transcription factor targets"){
          p1 <- multi_enrich1_keggGO()
        }else{
          p1 <- multi_enrich1_H()
        }
        pdf(file, height = 5, width = 7)
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
    content = function(file){write.table(as.data.frame(multi_enrichment_1_gsea()), file, row.names = F, sep = "\t", quote = F)}
  )
  
  #multi DEG enrichment 2--------
  multi_Hallmark_set2 <- reactive({
    if(input$Species6 != "not selected"){
      if(input$Gene_set7 != "MSigDB Hallmark" && input$Gene_set7 != "Transcription factor targets"){
        return(NULL)
      }else{
        switch (input$Species6,
                "Mus musculus" = species <- "Mus musculus",
                "Homo sapiens" = species <- "Homo sapiens",
                "Rattus norvegicus" = species <- "Rattus norvegicus",
                "Xenopus laevis" = species <- "Xenopus laevis",
                "Drosophila melanogaster" = species <- "Drosophila melanogaster",
                "Caenorhabditis elegans" = species <- "Caenorhabditis elegans")
        if(input$Gene_set7 == "MSigDB Hallmark"){
          H_t2g <- msigdbr(species = species, category = "H") %>%
            dplyr::select(gs_name, entrez_gene) 
        }
        if(input$Gene_set7 == "Transcription factor targets"){
          H_t2g <- msigdbr(species = species, category = "C3")
          H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "TFT:GTRD" | gs_subcat == "TFT:TFT_Legacy") %>%
            dplyr::select(gs_name, entrez_gene)
        }
        return(H_t2g)
      }
    }else return(NULL)
  })
  
  multi_Hallmark_set3 <- reactive({
    if(input$Species6 != "not selected"){
      if(input$Gene_set8 != "MSigDB Hallmark" && input$Gene_set8 != "Transcription factor targets"){
        return(NULL)
      }else{
        switch (input$Species6,
                "Mus musculus" = species <- "Mus musculus",
                "Homo sapiens" = species <- "Homo sapiens",
                "Rattus norvegicus" = species <- "Rattus norvegicus",
                "Xenopus laevis" = species <- "Xenopus laevis",
                "Drosophila melanogaster" = species <- "Drosophila melanogaster",
                "Caenorhabditis elegans" = species <- "Caenorhabditis elegans")
        if(input$Gene_set8 == "MSigDB Hallmark"){
          H_t2g <- msigdbr(species = species, category = "H") %>%
            dplyr::select(gs_name, entrez_gene) 
        }
        if(input$Gene_set8 == "Transcription factor targets"){
          H_t2g <- msigdbr(species = species, category = "C3")
          H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "TFT:GTRD" | gs_subcat == "TFT:TFT_Legacy") %>%
            dplyr::select(gs_name, entrez_gene)
        }
        return(H_t2g)
      }
    }else return(NULL)
  })
  
  
  output$Gene_set7 <- renderUI({
    selectInput('Gene_set7', 'Gene Set', c("KEGG", "GO", "MSigDB Hallmark", "Transcription factor targets"))
  })
  
  output$Gene_set8 <- renderUI({
    selectInput('Gene_set8', 'Gene Set', c("KEGG", "GO", "MSigDB Hallmark", "Transcription factor targets"))
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
  
  multi_enrich_viewer1 <- reactive({
    df <- multi_enrich_input1()
    if(is.null(df) || input$Species6 == "not selected"){
      return(NULL)
    }else{
      my.symbols <- df$GeneID
      if(str_detect(df$GeneID[1], "ENS")){
        gene_IDs<-AnnotationDbi::select(org6(),keys = my.symbols,
                                        keytype = "ENSEMBL",
                                        columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
        colnames(gene_IDs) <- c("GeneID","SYMBOL", "ENTREZID")
      }else{
        gene_IDs <- AnnotationDbi::select(org6(), keys = my.symbols,
                                          keytype = "SYMBOL",
                                          columns = c("ENTREZID", "SYMBOL"))
        colnames(gene_IDs) <- c("GeneID","ENTREZID")
      }
      gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
      data <- merge(df, gene_IDs, by="GeneID")
      return(data)
    }
  })
  
  multi_enrich_viewer11 <- reactive({
    df <- multi_enrich_input2()
    if(is.null(df) || input$Species6 == "not selected"){
      return(NULL)
    }else{
      my.symbols <- df$GeneID
      if(str_detect(df$GeneID[1], "ENS")){
        gene_IDs<-AnnotationDbi::select(org6(),keys = my.symbols,
                                        keytype = "ENSEMBL",
                                        columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
        colnames(gene_IDs) <- c("GeneID","SYMBOL", "ENTREZID")
      }else{
        gene_IDs <- AnnotationDbi::select(org6(), keys = my.symbols,
                                          keytype = "SYMBOL",
                                          columns = c("ENTREZID", "SYMBOL"))
        colnames(gene_IDs) <- c("GeneID","ENTREZID")
      }
      gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
      data <- merge(df, gene_IDs, by="GeneID")
      return(data)
    }
  })
  
  multi_enrich_viewer2 <- reactive({
    if(!is.null(input$Gene_set7)){
      data3 <- multi_enrich_viewer1()
      if(is.null(data3)){
        return(NULL)
      }else{
        if(input$Gene_set7 != "MSigDB Hallmark" && input$Gene_set7 != "Transcription factor targets"){
          if(input$Gene_set7 == "KEGG"){
            withProgress(message = "KEGG enrichment analysis",{
              formula_res <- try(compareCluster(ENTREZID~Group, data=data3,
                                                fun="enrichKEGG", organism=org_code6()), silent = T)
              incProgress(1)
            })
          }
          if(input$Gene_set7 == "GO"){
            withProgress(message = "GO enrichment analysis",{
              formula_res <- try(compareCluster(ENTREZID~Group, data=data3,
                                                fun="enrichGO", OrgDb=org6()), silent =T)
              incProgress(1)
            })
          }
          if (length(as.data.frame(formula_res)$Description) == 0) {
            formula_res <- NULL
          }else{
            formula_res <-setReadable(formula_res, org6(), 'ENTREZID')
          }
          return(formula_res)
        }else{
          withProgress(message = "enrichment analysis",{
            H_t2g <- multi_Hallmark_set2()
            df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
            colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
            for (name in unique(data3$Group)) {
              em <- enricher(data3$ENTREZID[data3$Group == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
              if (length(as.data.frame(em)$ID) != 0) {
                if(length(colnames(as.data.frame(em))) == 9){
                  cnet1 <- as.data.frame(setReadable(em, org6(), 'ENTREZID'))
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
  })
  
  multi_enrich_viewer12 <- reactive({
    if(!is.null(input$Gene_set8)){
      data3 <- multi_enrich_viewer11()
      if(is.null(data3)){
        return(NULL)
      }else{
        if(input$Gene_set8 != "MSigDB Hallmark" && input$Gene_set8 != "Transcription factor targets"){
          if(input$Gene_set8 == "KEGG"){
            withProgress(message = "KEGG enrichment analysis",{
              formula_res <- try(compareCluster(ENTREZID~Group, data=data3,
                                                fun="enrichKEGG", organism=org_code6()), silent = T)
              incProgress(1)
            })
          }
          if(input$Gene_set8 == "GO"){
            withProgress(message = "GO enrichment analysis",{
              formula_res <- try(compareCluster(ENTREZID~Group, data=data3,
                                                fun="enrichGO", OrgDb=org6()), silent =T)
              incProgress(1)
            })
          }
          if (length(as.data.frame(formula_res)$Description) == 0) {
            formula_res <- NULL
          }else{
            formula_res <-setReadable(formula_res, org6(), 'ENTREZID')
          }
          return(formula_res)
        }else{
          withProgress(message = "enrichment analysis",{
            H_t2g <- multi_Hallmark_set3()
            df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
            colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
            for (name in unique(data3$Group)) {
              em <- enricher(data3$ENTREZID[data3$Group == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
              if (length(as.data.frame(em)$ID) != 0) {
                if(length(colnames(as.data.frame(em))) == 9){
                  cnet1 <- as.data.frame(setReadable(em, org6(), 'ENTREZID'))
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
  })
  
  multi_enrich_keggGO <- reactive({
    if(!is.null(input$Gene_set7)){
      formula_res <- multi_enrich_viewer2()
      if(is.null(formula_res)){
        return(NULL)
      }else{
        if(input$Gene_set7 != "MSigDB Hallmark" && input$Gene_set7 != "Transcription factor targets"){
          if ((length(as.data.frame(formula_res)$Description) == 0) ||
              length(which(!is.na(unique(as.data.frame(formula_res)$qvalue)))) == 0) {
            p1 <- NULL
          } else{
            formula_res <- clusterProfiler.dplyr::filter(formula_res, !is.na(qvalue))
            p1 <- as.grob(dotplot(formula_res, color ="qvalue", font.size = 10))
          }
          p <- plot_grid(p1)
          return(p)
        }
      }
    }
  })
  
  multi_enrich_keggGO2 <- reactive({
    if(!is.null(input$Gene_set8)){
      formula_res <- multi_enrich_viewer12()
      if(is.null(formula_res)){
        return(NULL)
      }else{
        if(input$Gene_set8 != "MSigDB Hallmark" && input$Gene_set8 != "Transcription factor targets"){
          if ((length(as.data.frame(formula_res)$Description) == 0) ||
              length(which(!is.na(unique(as.data.frame(formula_res)$qvalue)))) == 0) {
            p1 <- NULL
          } else{
            formula_res <- clusterProfiler.dplyr::filter(formula_res, !is.na(qvalue))
            p1 <- as.grob(dotplot(formula_res, color ="qvalue", font.size = 10))
          }
          p <- plot_grid(p1)
          return(p)
        }
      }
    }
  })
  
  multi_enrich_H <- reactive({
    if(!is.null(input$Gene_set7)){
      data3 <- multi_enrich_viewer1()
      if(is.null(data3)){
        return(NULL)
      }else{
        if(input$Gene_set7 != "MSigDB Hallmark" && input$Gene_set7 != "Transcription factor targets"){
          return(NULL)
        }else{
          H_t2g <- multi_Hallmark_set2()
          df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
          for (name in unique(data3$Group)) {
            em <- enricher(data3$ENTREZID[data3$Group == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) != 0) {
              if(length(colnames(as.data.frame(em))) == 9){
                cnet1 <- as.data.frame(setReadable(em, org6(), 'ENTREZID'))
                cnet1$Group <- name
                cnet1 <- cnet1[sort(cnet1$pvalue, decreasing = F, index=T)$ix,]
                if (length(cnet1$pvalue) > 5){
                  cnet1 <- cnet1[1:5,]
                }
                df <- rbind(df, cnet1)
              }
            }
          }
          if ((length(df$Description) == 0) || length(which(!is.na(unique(df$qvalue)))) == 0) {
            p1 <- NULL
          } else{
            df["Description"] <- lapply(df["Description"], gsub, pattern="HALLMARK_", replacement = "") 
            df$GeneRatio <- parse_ratio(df$GeneRatio)
            df <- dplyr::filter(df, !is.na(qvalue))
            df <- dplyr::mutate(df, x = paste0(Group, eval(parse(text = "GeneRatio"))))
            df$x <- gsub(":","", df$x)
            df <- dplyr::arrange(df, x)
            idx <- order(df[["x"]], decreasing = FALSE)
            df$Description <- factor(df$Description,
                                     levels=rev(unique(df$Description[idx])))
            p1 <- as.grob(ggplot(df, aes(x = Group,y= Description,color=qvalue,size=GeneRatio))+
                            geom_point() +
                            scale_color_continuous(low="red", high="blue",
                                                   guide=guide_colorbar(reverse=TRUE)) +
                            scale_size(range=c(3, 8))+ theme_dose(font.size=8)+ylab(NULL))
            p <- plot_grid(p1)
            return(p)
          }
        }}
    }
  })
  
  multi_enrich_H2 <- reactive({
    if(!is.null(input$Gene_set8)){
      data3 <- multi_enrich_viewer11()
      if(is.null(data3)){
        return(NULL)
      }else{
        if(input$Gene_set8 != "MSigDB Hallmark" && input$Gene_set8 != "Transcription factor targets"){
          return(NULL)
        }else{
          H_t2g <- multi_Hallmark_set3()
          df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
          for (name in unique(data3$Group)) {
            em <- enricher(data3$ENTREZID[data3$Group == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) != 0) {
              if(length(colnames(as.data.frame(em))) == 9){
                cnet1 <- as.data.frame(setReadable(em, org6(), 'ENTREZID'))
                cnet1$Group <- name
                cnet1 <- cnet1[sort(cnet1$pvalue, decreasing = F, index=T)$ix,]
                if (length(cnet1$pvalue) > 5){
                  cnet1 <- cnet1[1:5,]
                }
                df <- rbind(df, cnet1)
              }
            }
          }
          if ((length(df$Description) == 0) || length(which(!is.na(unique(df$qvalue)))) == 0) {
            p1 <- NULL
          } else{
            df["Description"] <- lapply(df["Description"], gsub, pattern="HALLMARK_", replacement = "") 
            df$GeneRatio <- parse_ratio(df$GeneRatio)
            df <- dplyr::filter(df, !is.na(qvalue))
            df <- dplyr::mutate(df, x = paste0(Group, eval(parse(text = "GeneRatio"))))
            df$x <- gsub(":","", df$x)
            df <- dplyr::arrange(df, x)
            idx <- order(df[["x"]], decreasing = FALSE)
            df$Description <- factor(df$Description,
                                     levels=rev(unique(df$Description[idx])))
            p1 <- as.grob(ggplot(df, aes(x = Group,y= Description,color=qvalue,size=GeneRatio))+
                            geom_point() +
                            scale_color_continuous(low="red", high="blue",
                                                   guide=guide_colorbar(reverse=TRUE)) +
                            scale_size(range=c(3, 8))+ theme_dose(font.size=8)+ylab(NULL))
            p <- plot_grid(p1)
            return(p)
          }
        }}
    }
  })
  
  output$multi_enrichment3 <- renderPlot({
    if(!is.null(input$Gene_set7)){
      if(is.null(multi_enrich_viewer2())){
        return(NULL)
      }else{
        withProgress(message = "Plot results",{
          if(input$Species6 != "not selected"){
            if(input$Gene_set7 != "MSigDB Hallmark" && input$Gene_set7 != "Transcription factor targets"){
              print(multi_enrich_keggGO())
            }else{
              print(multi_enrich_H())
            }
          }
          incProgress(1)
        })
      }
    }
  })
  
  output$multi_enrichment5 <- renderPlot({
    if(!is.null(input$Gene_set8)){
      if(is.null(multi_enrich_viewer12())){
        return(NULL)
      }else{
        withProgress(message = "Plot results",{
          if(input$Species6 != "not selected"){
            if(input$Gene_set8 != "MSigDB Hallmark" && input$Gene_set8 != "Transcription factor targets"){
              print(multi_enrich_keggGO2())
            }else{
              print(multi_enrich_H2())
            }
          }
          incProgress(1)
        })
      }
    }
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
  
  multi_enrich_viewer3 <- reactive({
    df <- multi_enrich_input3()
    if(is.null(df) || input$Species6 == "not selected"){
      return(NULL)
    }else{
      my.symbols <- df$GeneID
      if(str_detect(df$GeneID[1], "ENS")){
        gene_IDs<-AnnotationDbi::select(org6(),keys = my.symbols,
                                        keytype = "ENSEMBL",
                                        columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
        colnames(gene_IDs) <- c("GeneID","SYMBOL", "ENTREZID")
      }else{
        gene_IDs <- AnnotationDbi::select(org6(), keys = my.symbols,
                                          keytype = "SYMBOL",
                                          columns = c("ENTREZID", "SYMBOL"))
        colnames(gene_IDs) <- c("GeneID","ENTREZID")
      }
      gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
      data <- merge(df, gene_IDs, by="GeneID")
      return(data)
    }
  })
  
  multi_enrich_viewer13 <- reactive({
    df <- multi_enrich_input4()
    if(is.null(df) || input$Species6 == "not selected"){
      return(NULL)
    }else{
      my.symbols <- df$GeneID
      if(str_detect(df$GeneID[1], "ENS")){
        gene_IDs<-AnnotationDbi::select(org6(),keys = my.symbols,
                                        keytype = "ENSEMBL",
                                        columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
        colnames(gene_IDs) <- c("GeneID","SYMBOL", "ENTREZID")
      }else{
        gene_IDs <- AnnotationDbi::select(org6(), keys = my.symbols,
                                          keytype = "SYMBOL",
                                          columns = c("ENTREZID", "SYMBOL"))
        colnames(gene_IDs) <- c("GeneID","ENTREZID")
      }
      gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
      data <- merge(df, gene_IDs, by="GeneID")
      return(data)
    }
  })
  
  multi_enrich2 <- reactive({
    if(!is.null(input$Gene_set7)){
      data <- multi_enrich_viewer3()
      group <- input$multi_whichGroup1_2
      if(is.null(data) || is.null(group)){
        return(NULL)
      }else{
        data2 <- dplyr::filter(data, Group == group)
        if(input$Gene_set7 != "MSigDB Hallmark" && input$Gene_set7 != "Transcription factor targets"){
          if(input$Gene_set7 == "KEGG"){
            kk1 <- enrichKEGG(data2$ENTREZID, organism =org_code6(),
                              pvalueCutoff = 0.05)
          }
          if(input$Gene_set7 == "GO"){
            kk1 <- enrichGO(data2$ENTREZID, OrgDb = org6(),
                            pvalueCutoff = 0.05)
          }
        }else{
          H_t2g <- multi_Hallmark_set2()
          kk1 <- try(enricher(data2$ENTREZID, TERM2GENE=H_t2g, pvalueCutoff = 0.05))
          if (class(kk1) == "try-error") kk1 <- NA
        }
        if(length(as.data.frame(kk1)$ID) == 0){
          cnet1 <- NULL
        } else {
          cnet1 <- setReadable(kk1, org6(), 'ENTREZID')
        }
        if (length(as.data.frame(cnet1)$ID) == 0) {
          p2 <- NULL
        } else{
          p2 <- try(as.grob(cnetplot(cnet1,
                                     cex_label_gene = 0.7, cex_label_category = 0.75,
                                     cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none")))
          if(length(class(p2)) == 1){
            if(class(p2) == "try-error") p2 <- NULL
          }else{p2 <- as.grob(cnetplot(cnet1,
                                       cex_label_gene = 0.7, cex_label_category = 0.75,
                                       cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none"))}
        }
        p <- plot_grid(p2)
        return(p)
      }
    }
  })
  
  multi_enrich12 <- reactive({
    if(!is.null(input$Gene_set8)){
      data <- multi_enrich_viewer13()
      group <- input$multi_whichGroup2_2
      if(is.null(data) || is.null(group)){
        return(NULL)
      }else{
        data2 <- dplyr::filter(data, Group == group)
        if(input$Gene_set8 != "MSigDB Hallmark" && input$Gene_set8 != "Transcription factor targets"){
          if(input$Gene_set8 == "KEGG"){
            kk1 <- enrichKEGG(data2$ENTREZID, organism =org_code6(),
                              pvalueCutoff = 0.05)
          }
          if(input$Gene_set8 == "GO"){
            kk1 <- enrichGO(data2$ENTREZID, OrgDb = org6(),
                            pvalueCutoff = 0.05)
          }
        }else{
          H_t2g <- multi_Hallmark_set3()
          kk1 <- try(enricher(data2$ENTREZID, TERM2GENE=H_t2g, pvalueCutoff = 0.05))
          if (class(kk1) == "try-error") kk1 <- NA
        }
        if(length(as.data.frame(kk1)$ID) == 0){
          cnet1 <- NULL
        } else {
          cnet1 <- setReadable(kk1, org6(), 'ENTREZID')
        }
        if (length(as.data.frame(cnet1)$ID) == 0) {
          p2 <- NULL
        } else{
          p2 <- try(as.grob(cnetplot(cnet1,
                                     cex_label_gene = 0.7, cex_label_category = 0.75,
                                     cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none")))
          if(length(class(p2)) == 1){
            if(class(p2) == "try-error") p2 <- NULL
          }else{p2 <- as.grob(cnetplot(cnet1,
                                       cex_label_gene = 0.7, cex_label_category = 0.75,
                                       cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none"))}
        }
        p <- plot_grid(p2)
        return(p)
      }
    }
  })
  
  output$multi_enrichment4 <- renderPlot({
    if(!is.null(input$Gene_set7)){
      if(is.null(multi_enrich_input3())){
        return(NULL)
      }else{
        if(input$Species6 != "not selected"){
          withProgress(message = "cnet plot",{
            p <- multi_enrich2()
            print(p)
            incProgress(1)
          })
        }else return(NULL)
      }
    }
  })
  
  output$multi_enrichment6 <- renderPlot({
    if(!is.null(input$Gene_set8)){
      if(is.null(multi_enrich_input4())){
        return(NULL)
      }else{
        if(input$Species6 != "not selected"){
          withProgress(message = "cnet plot",{
            p <- multi_enrich12()
            print(p)
            incProgress(1)
          })
        }else return(NULL)
      }
    }
  })
  
  output$download_multi_cluster_enrichment = downloadHandler(
    filename = function(){
      paste(download_multi_overview_dir(), paste0(input$Gene_set7,"_divisive_enrichment.pdf"), sep="_")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$Gene_set7 != "MSigDB Hallmark" && input$Gene_set7 != "Transcription factor targets"){
          p1 <- multi_enrich_keggGO()
        }else{
          p1 <- multi_enrich_H()
        }
        pdf(file, height = 6, width = 8)
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
        if(input$Gene_set8 != "MSigDB Hallmark" && input$Gene_set8 != "Transcription factor targets"){
          p1 <- multi_enrich_keggGO2()
        }else{
          p1 <- multi_enrich_H2()
        }
        pdf(file, height = 6, width = 8)
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
        pdf(file, height = 6, width = 6)
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
        pdf(file, height = 6, width = 6)
        print(p)
        dev.off()
        incProgress(1)
      })
    }
  )
  
  output$multi_enrichment_result <- DT::renderDataTable({
    as.data.frame(multi_enrich_viewer2())
  })
  
  output$multi_enrichment_result2 <- DT::renderDataTable({
    as.data.frame(multi_enrich_viewer12())
  })
  
  output$download_multi_enrichment_table = downloadHandler(
    filename = function() {
      paste(download_multi_overview_dir(), paste0(input$Gene_set7,"_divisive_enrichment.txt"), sep="_")
    },
    content = function(file){write.table(as.data.frame(multi_enrich_viewer2()), file, row.names = F, sep = "\t", quote = F)}
  )
  
  output$download_multi_enrichment_table2 = downloadHandler(
    filename = function() {
      paste(download_multi_overview_dir(), paste0(input$Gene_set8,"_kmeans_enrichment.txt"), sep="_")
    },
    content = function(file){write.table(as.data.frame(multi_enrich_viewer12()), file, row.names = F, sep = "\t", quote = F)}
  )
  
  #Multi PCA ------------------------------------------------------------------------------
  multiPCAdata <- reactive({
    if(is.null(multi_d_row_count_matrix())){
      return(NULL)
    }else{
      data <- multi_deg_norm_count()
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
  })
  
  multi_heat <- reactive({
    data <- as.data.frame(multi_deg_norm_count())
    meta <- multi_metadata()
    dds <- multi_dds()
    if(is.null(dds)){
      return(NULL)
    }else{
      withProgress(message = "DEG heatmap",{
        if (input$multi_data_file_type == "Row1"){
          collist <- gsub("\\_.+$", "", colnames(data))
          meta <- data.frame(condition = factor(collist))
        }else meta <- data.frame(condition=factor(meta[,1]), type=factor(meta[,2]))
        res <- results(dds)
        
        collist <- gsub("\\_.+$", "", colnames(data))
        data <- dplyr::filter(data, apply(data,1,mean) > input$basemean6)
        res <- as.data.frame(res)
        data2 <- merge(res,data, by=0)
        res <- data2[,1:7]
        rownames(res) <- res$Row.names
        
        sig_res_LRT <- res %>% 
          data.frame() %>%
          rownames_to_column(var="gene") %>%
          as_tibble() %>%
          filter(padj < input$fdr6)
        data <- data %>% 
          data.frame() %>%
          rownames_to_column(var="gene")
        data <- merge(data, sig_res_LRT, by = "gene")
        rownames(data) <- data$gene
        data <- data[,2:(1+length(collist))]
        if(length(rownames(data)) == 0){
          ht <- NULL
        }else{
        data.z <- genescale(data, axis=1, method="Z")
        ht <- as.grob(Heatmap(data.z, name = "z-score", column_order = colnames(data.z),
                              clustering_method_columns = 'ward.D2',
                              show_row_names = F, show_row_dend = T))
        return(ht)
        }
      })
    }
  })
  
  multi_pca_plot <- reactive({
    ht <- multi_heat()
    data <- multi_deg_norm_count()
    withProgress(message = "PCA",{
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
      if(is.null(ht)){
        p2 <- gridExtra::grid.arrange(gridExtra::arrangeGrob(g1, g2, g3, ncol = 1), ncol = 2)
      }else{
      p2 <- gridExtra::grid.arrange(gridExtra::arrangeGrob(g1, g2, g3, ncol = 1), ht, ncol = 2)
      }
      return(p2)
    })
  })
  
  output$download_multi_PCA = downloadHandler(
    filename = function(){
      paste0(download_multi_overview_dir(), "_PCA-MDS-dendrogram.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        pdf(file, height = 9, width = 7)
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
    multiPCAdata()
  })
  
  output$download_multi_PCA_table = downloadHandler(
    filename = function() {
      paste0(download_multi_overview_dir(), "_PCA_table.txt")
    },
    content = function(file){write.table(multiPCAdata(), file, row.names = T, sep = "\t", quote = F)}
  )
  
  # 3 conditions ------------------------------------------------------------------------------
  org2 <- reactive({
    if(input$Species2 != "not selected"){
      switch (input$Species2,
              "Mus musculus" = org <- org.Mm.eg.db,
              "Homo sapiens" = org <- org.Hs.eg.db,
              "Rattus norvegicus" = org <- org.Rn.eg.db,
              "Xenopus laevis" = org <- org.Xl.eg.db,
              "Drosophila melanogaster" = org <- org.Dm.eg.db,
              "Caenorhabditis elegans" = org <- org.Ce.eg.db)
      return(org)
    }
  })
  org_code2 <- reactive({
    if(input$Species2 != "not selected"){
      switch (input$Species2,
              "Mus musculus" = org_code <- "mmu",
              "Homo sapiens" = org_code <- "hsa",
              "Rattus norvegicus" = org_code <- "rno",
              "Xenopus laevis" = org_code <- "xla",
              "Drosophila melanogaster" = org_code <- "dme",
              "Caenorhabditis elegans" = org_code <- "cel")
      return(org_code)
    }
  })

    row_count_matrix2 <- reactive({
    if (input$data_file_type2 == "Row3"){
      tmp <- input$file4$datapath
      if(is.null(input$file4) && input$goButton2 > 0 )  tmp = "data/example4.txt"
      if(is.null(tmp)) {
        return(NULL)
      }else{
        if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
        if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
        if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
        return(df)
      }
    }else{
      tmp <- input$file5$datapath
      if(is.null(input$file5) && input$goButton2 == 0) return(NULL)
      if(is.null(input$file5) && input$goButton2 > 0 )  tmp = "data/example2.csv"
      if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
      if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
      if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
      return(df)
    }
  })
  metadata2 <- reactive({
    if (input$data_file_type2 == "Row3"){
      return(NULL)
    }else{
    tmp <- input$file6$datapath
    if(is.null(input$file6) && input$goButton2 == 0) return(NULL)
    if(is.null(input$file6) && input$goButton2 > 0 )  tmp = "data/example6.csv"
    if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
    if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
    if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
    rownames(df) <- gsub("-",".",rownames(df))
    return(df)
    }
  })
  norm_count_matrix2 <- reactive({
    data <- input$norm_file2$datapath
    if(is.null(input$norm_file2)) {
      return(NULL)
    }else{
    if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
    if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
    if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
    if(input$Species2 != "not selected"){
      if(str_detect(rownames(df)[1], "ENS")){
        rownames(df) < gsub("\\..*","", rownames(df))
      }
    }
    return(df)
    }
  })
  
  d_row_count_matrix2 <- reactive({
    row <- row_count_matrix2()
    if (input$data_file_type2 == "Row3"){
      if(is.null(row)) {
        return(NULL)
      }else{
        return(row)
      }
    }else{
      meta <- metadata2()
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
    withProgress(message = "EBSeq multiple comparison test takes a few minutes",{
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
    if(is.null(norm_count_matrix2())){
    count <- d_row_count_matrix2()
    if(is.null(count)){
      return(NULL)
    }else{
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
    }else{
      return(norm_count_matrix2())
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
    content = function(file){write.table(data_3degcount2_1(), file, row.names = T, sep = "\t", quote = F)}
  )
  output$download_3cond_DEG_table2 = downloadHandler(
    filename = function() {paste(download_cond3_dir(),"DEG_result2.txt", sep = "-")},
    content = function(file){write.table(data_3degcount2_2(), file, row.names = T, sep = "\t", quote = F)}
  )
  output$download_3cond_DEG_table3 = downloadHandler(
    filename = function() {paste(download_cond3_dir(),"DEG_result3.txt", sep = "-")},
    content = function(file){write.table(data_3degcount2_3(), file, row.names = T, sep = "\t", quote = F)}
  )


  #3conditions DEG vis------------------------

  #3conditions DEG_1------------------------
  data_3degcount1_1 <- reactive({
    data <- deg_norm_count2()
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
    result_Condm <- deg_result2_condmean()
    result_FDR <- deg_result2()
    if(str_detect(rownames(data)[1], "ENS")){
      if(length(grep("SYMBOL", colnames(data))) != 0){
        data <- data[, - which(colnames(data) == "SYMBOL")]
      }
    }
    specific = collist[1]
      FC_xlab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[2], ")"))
      FC_ylab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[3], ")"))
      result_Condm$FC_x <- log2((result_Condm$C1 + 0.01)/(result_Condm$C2 + 0.01))
      result_Condm$FC_y <- log2((result_Condm$C1 + 0.01)/(result_Condm$C3 + 0.01))
      Pattern1 <- "Pattern4"
      Pattern2 <- "Pattern5"
      result_FDR$FDR <- 1 - result_FDR$PPDE
      result <- merge(result_Condm, result_FDR, by=0)
      data$Row.names <- rownames(data)
      data2 <- merge(data, result, by="Row.names")
      result <- dplyr::filter(data2, apply(data2[,2:(Cond_1 + Cond_2 + Cond_3)],1,mean) > input$basemean2)
      sig <- rep(3, nrow(result))
      sig[which(result$FDR <= input$fdr2 & result$FC_x < log2(1/input$fc2) & result$FC_y < log2(1/input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 2
      sig[which(result$FDR <= input$fdr2 & result$FC_x > log2(input$fc2) & result$FC_y > log2(input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 1
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
  })

  data_3degcount2_1 <- reactive({
    data3 <- data_3degcount1_1()
    if(is.null(data3)){
      return(NULL)
    }else{
    if(length(unique(data3$sig)) == 1){
      data4 <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
      return(data4)
    } else {
      data4 <- dplyr::filter(data3, sig != "NS")
      if(str_detect(data4$Row.names[1], "ENS")){
        if(input$Species2 != "not selected"){
          my.symbols <- data4$Row.names
          gene_IDs<-AnnotationDbi::select(org2(),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data4 <- merge(data4, gene_IDs, by="Row.names")
        }
      }else{
        if(input$Species2 != "not selected"){
          my.symbols <- data4$Row.names
          gene_IDs<-AnnotationDbi::select(org2(),keys = my.symbols,
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
  })

  #3conditions scatter + heatmap_1

  cond3_scatter1_plot <- reactive({
    data <- deg_norm_count2()
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
    result_Condm <- deg_result2_condmean()
    result_FDR <- deg_result2()
    if(str_detect(rownames(data)[1], "ENS")){
      if(length(grep("SYMBOL", colnames(data))) != 0){
        data <- data[, - which(colnames(data) == "SYMBOL")]
      }
    }
    specific = collist[1]
    FC_xlab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[2], ")"))
    FC_ylab <- paste0(paste0(paste0("Log2(", collist[1]) ,"/"), paste0(collist[3], ")"))
    result_Condm$FC_x <- log2((result_Condm$C1 + 0.01)/(result_Condm$C2 + 0.01))
    result_Condm$FC_y <- log2((result_Condm$C1 + 0.01)/(result_Condm$C3 + 0.01))
    Pattern1 <- "Pattern4"
    Pattern2 <- "Pattern5"
    result_FDR$FDR <- 1 - result_FDR$PPDE
    result <- merge(result_Condm, result_FDR, by=0)
    data$Row.names <- rownames(data)
    data2 <- merge(data, result, by="Row.names")
    result <- dplyr::filter(data2, apply(data2[,2:(Cond_1 + Cond_2 + Cond_3)],1,mean) > input$basemean2)
    sig <- rep(3, nrow(result))
    sig[which(result$FDR <= input$fdr2 & result$FC_x < log2(1/input$fc2) & result$FC_y < log2(1/input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 2
    sig[which(result$FDR <= input$fdr2 & result$FC_x > log2(input$fc2) & result$FC_y > log2(input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 1
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
    labs_data <- subset(complete_data, padj <= input$fdr2 & Row.names !=
                          "" & (FC_x)*(FC_y) >= log2(input$fc2))
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
    p <- p  + geom_hline(yintercept = c(-log2(input$fc2), log2(input$fc2)), linetype = c(2, 2), color = c("black", "black"))+
      geom_vline(xintercept = c(-log2(input$fc2), log2(input$fc2)),linetype = c(2, 2), color = c("black", "black"))
    p <- p +
      theme_bw()+ scale_color_manual(values = col)+
      theme(legend.position = "top" , legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 10),
            axis.text.y= ggplot2::element_text(size = 10),
            text = ggplot2::element_text(size = 10),
            title = ggplot2::element_text(size = 10)) +
      xlab(FC_xlab) + ylab(FC_ylab)

    data4 <- data_3degcount2_1()
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
      paste0(download_cond3_dir(), "_scatter1.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        pdf(file, height = 6, width = 10)
        print(cond3_scatter1_plot())
        dev.off()
      })
    }
  )
  #3conditions DEG_2------------------------
  data_3degcount1_2 <- reactive({
    data <- deg_norm_count2()
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
    result_Condm <- deg_result2_condmean()
    result_FDR <- deg_result2()
    if(str_detect(rownames(data)[1], "ENS")){
      if(length(grep("SYMBOL", colnames(data))) != 0){
        data <- data[, - which(colnames(data) == "SYMBOL")]
      }
    }
    specific = collist[2]
    FC_xlab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[1], ")"))
    FC_ylab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[3], ")"))
    result_Condm$FC_x <- log2((result_Condm$C2 + 0.01)/(result_Condm$C1 + 0.01))
    result_Condm$FC_y <- log2((result_Condm$C2 + 0.01)/(result_Condm$C3 + 0.01))
    Pattern1 <- "Pattern4"
    Pattern2 <- "Pattern5"
    result_FDR$FDR <- 1 - result_FDR$PPDE
    result <- merge(result_Condm, result_FDR, by=0)
    data$Row.names <- rownames(data)
    data2 <- merge(data, result, by="Row.names")
    result <- dplyr::filter(data2, apply(data2[,2:(Cond_1 + Cond_2 + Cond_3)],1,mean) > input$basemean2)
    sig <- rep(3, nrow(result))
    sig[which(result$FDR <= input$fdr2 & result$FC_x < log2(1/input$fc2) & result$FC_y < log2(1/input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 2
    sig[which(result$FDR <= input$fdr2 & result$FC_x > log2(input$fc2) & result$FC_y > log2(input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 1
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
  })

  data_3degcount2_2 <- reactive({
    data3 <- data_3degcount1_2()
    if(is.null(data3)){
      return(NULL)
    }else{
    if(length(unique(data3$sig)) == 1){
      data4 <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
      return(data4)
    } else {
      data4 <- dplyr::filter(data3, sig != "NS")
      if(str_detect(data4$Row.names[1], "ENS")){
        if(input$Species2 != "not selected"){
          my.symbols <- data4$Row.names
          gene_IDs<-AnnotationDbi::select(org2(),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data4 <- merge(data4, gene_IDs, by="Row.names")
        }
      }else{
        if(input$Species2 != "not selected"){
          my.symbols <- data4$Row.names
          gene_IDs<-AnnotationDbi::select(org2(),keys = my.symbols,
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
  })

  #3conditions scatter + heatmap_2
  cond3_scatter2_plot <- reactive({
    data <- deg_norm_count2()
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
    result_Condm <- deg_result2_condmean()
    result_FDR <- deg_result2()
    if(str_detect(rownames(data)[1], "ENS")){
      if(length(grep("SYMBOL", colnames(data))) != 0){
        data <- data[, - which(colnames(data) == "SYMBOL")]
      }
    }
    specific = collist[2]
    FC_xlab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[1], ")"))
    FC_ylab <- paste0(paste0(paste0("Log2(", collist[2]) ,"/"), paste0(collist[3], ")"))
    result_Condm$FC_x <- log2((result_Condm$C2 + 0.01)/(result_Condm$C1 + 0.01))
    result_Condm$FC_y <- log2((result_Condm$C2 + 0.01)/(result_Condm$C3 + 0.01))
    Pattern1 <- "Pattern4"
    Pattern2 <- "Pattern5"
    result_FDR$FDR <- 1 - result_FDR$PPDE
    result <- merge(result_Condm, result_FDR, by=0)
    data$Row.names <- rownames(data)
    data2 <- merge(data, result, by="Row.names")
    result <- dplyr::filter(data2, apply(data2[,2:(Cond_1 + Cond_2 + Cond_3)],1,mean) > input$basemean2)
    sig <- rep(3, nrow(result))
    sig[which(result$FDR <= input$fdr2 & result$FC_x < log2(1/input$fc2) & result$FC_y < log2(1/input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 2
    sig[which(result$FDR <= input$fdr2 & result$FC_x > log2(input$fc2) & result$FC_y > log2(input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 1
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
    labs_data <- subset(complete_data, padj <= input$fdr2 & Row.names !=
                          "" & (FC_x)*(FC_y) >= log2(input$fc2))
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
    p <- p  + geom_hline(yintercept = c(-log2(input$fc2), log2(input$fc2)), linetype = c(2, 2), color = c("black", "black"))+
      geom_vline(xintercept = c(-log2(input$fc2), log2(input$fc2)),linetype = c(2, 2), color = c("black", "black"))
    p <- p +
      theme_bw()+ scale_color_manual(values = col)+
      theme(legend.position = "top" , legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 10),
            axis.text.y= ggplot2::element_text(size = 10),
            text = ggplot2::element_text(size = 10),
            title = ggplot2::element_text(size = 10)) +
      xlab(FC_xlab) + ylab(FC_ylab)

    data4 <- data_3degcount2_2()
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

    p<- plot_grid(p, ht, rel_widths = c(2, 1))
    return(p)
    }
  })

  output$scatter_2 <- renderPlot({
    if(is.null(deg_norm_count2())){
      return(NULL)
    }else{
    print(cond3_scatter2_plot())
    }
  })

  output$download_3cond_scatter2 = downloadHandler(
    filename = function(){
      paste0(download_cond3_dir(), "_scatter2.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        pdf(file, height = 6, width = 10)
        print(cond3_scatter2_plot())
        dev.off()
      })
    }
  )

  #3conditions DEG_3------------------------
  data_3degcount1_3 <- reactive({
    data <- deg_norm_count2()
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
    result_Condm <- deg_result2_condmean()
    result_FDR <- deg_result2()
    if(str_detect(rownames(data)[1], "ENS")){
      if(length(grep("SYMBOL", colnames(data))) != 0){
        data <- data[, - which(colnames(data) == "SYMBOL")]
      }
    }
    specific = collist[3]
    FC_xlab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[1], ")"))
    FC_ylab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[2], ")"))
    result_Condm$FC_x <- log2((result_Condm$C3 + 0.01)/(result_Condm$C1 + 0.01))
    result_Condm$FC_y <- log2((result_Condm$C3 + 0.01)/(result_Condm$C2 + 0.01))
    Pattern1 <- "Pattern4"
    Pattern2 <- "Pattern5"
    result_FDR$FDR <- 1 - result_FDR$PPDE
    result <- merge(result_Condm, result_FDR, by=0)
    data$Row.names <- rownames(data)
    data2 <- merge(data, result, by="Row.names")
    result <- dplyr::filter(data2, apply(data2[,2:(Cond_1 + Cond_2 + Cond_3)],1,mean) > input$basemean2)
    sig <- rep(3, nrow(result))
    sig[which(result$FDR <= input$fdr2 & result$FC_x < log2(1/input$fc2) & result$FC_y < log2(1/input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 2
    sig[which(result$FDR <= input$fdr2 & result$FC_x > log2(input$fc2) & result$FC_y > log2(input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 1
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
  })

  data_3degcount2_3 <- reactive({
    data3 <- data_3degcount1_3()
    if(is.null(data3)){
      return(NULL)
    }else{
    if(length(unique(data3$sig)) == 1){
      data4 <- data.frame(matrix(rep(NA, 6), nrow=1))[numeric(0), ]
      return(data4)
    } else {
      data4 <- dplyr::filter(data3, sig != "NS")
      if(str_detect(data4$Row.names[1], "ENS")){
        if(input$Species2 != "not selected"){
          my.symbols <- data4$Row.names
          gene_IDs<-AnnotationDbi::select(org2(),keys = my.symbols,
                                          keytype = "ENSEMBL",
                                          columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
          colnames(gene_IDs) <- c("Row.names","SYMBOL", "ENTREZID")
          gene_IDs <- gene_IDs %>% distinct(Row.names, .keep_all = T)
          data4 <- merge(data4, gene_IDs, by="Row.names")
        }
      }else{
        if(input$Species2 != "not selected"){
          my.symbols <- data4$Row.names
          gene_IDs<-AnnotationDbi::select(org2(),keys = my.symbols,
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
  })

  #3conditions scatter + heatmap_3
  cond3_scatter3_plot <- reactive({
    data <- deg_norm_count2()
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
    result_Condm <- deg_result2_condmean()
    result_FDR <- deg_result2()
    if(str_detect(rownames(data)[1], "ENS")){
      if(length(grep("SYMBOL", colnames(data))) != 0){
        data <- data[, - which(colnames(data) == "SYMBOL")]
      }
    }
    specific = collist[3]
    FC_xlab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[1], ")"))
    FC_ylab <- paste0(paste0(paste0("Log2(", collist[3]) ,"/"), paste0(collist[2], ")"))
    result_Condm$FC_x <- log2((result_Condm$C3 + 0.01)/(result_Condm$C1 + 0.01))
    result_Condm$FC_y <- log2((result_Condm$C3 + 0.01)/(result_Condm$C2 + 0.01))
    Pattern1 <- "Pattern4"
    Pattern2 <- "Pattern5"
    result_FDR$FDR <- 1 - result_FDR$PPDE
    result <- merge(result_Condm, result_FDR, by=0)
    data$Row.names <- rownames(data)
    data2 <- merge(data, result, by="Row.names")
    result <- dplyr::filter(data2, apply(data2[,2:(Cond_1 + Cond_2 + Cond_3)],1,mean) > input$basemean2)
    sig <- rep(3, nrow(result))
    sig[which(result$FDR <= input$fdr2 & result$FC_x < log2(1/input$fc2) & result$FC_y < log2(1/input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 2
    sig[which(result$FDR <= input$fdr2 & result$FC_x > log2(input$fc2) & result$FC_y > log2(input$fc2) & (result$MAP == Pattern1 | result$MAP == Pattern2))] = 1
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
    labs_data <- subset(complete_data, padj <= input$fdr2 & Row.names !=
                          "" & (FC_x)*(FC_y) >= log2(input$fc2))
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
    p <- p  + geom_hline(yintercept = c(-log2(input$fc2), log2(input$fc2)), linetype = c(2, 2), color = c("black", "black"))+
      geom_vline(xintercept = c(-log2(input$fc2), log2(input$fc2)),linetype = c(2, 2), color = c("black", "black"))
    p <- p +
      theme_bw()+ scale_color_manual(values = col)+
      theme(legend.position = "top" , legend.title = element_blank(),
            axis.text.x= ggplot2::element_text(size = 10),
            axis.text.y= ggplot2::element_text(size = 10),
            text = ggplot2::element_text(size = 10),
            title = ggplot2::element_text(size = 10)) +
      xlab(FC_xlab) + ylab(FC_ylab)

    data4 <- data_3degcount2_3()
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
  })

  output$scatter_3 <- renderPlot({
    if(is.null(deg_norm_count2())){
      return(NULL)
    }else{
    print(cond3_scatter3_plot())
    }
  })

  output$download_3cond_scatter3 = downloadHandler(
    filename = function(){
      paste0(download_cond3_dir(), "_scatter3.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        pdf(file, height = 6, width = 10)
        print(cond3_scatter3_plot())
        dev.off()
      })
    }
  )


  #3conditions PCA--------------
  PCA3data <- reactive({
    data <- deg_norm_count2()
    if(is.null(data)){
      return(NULL)
    }else{
    if(length(grep("SYMBOL", colnames(data))) != 0){
      data <- data[, - which(colnames(data) == "SYMBOL")]
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
  })

  cond3_pca <- reactive({
    data <- deg_norm_count2()
    if(is.null(data)){
      return(NULL)
    }else{
    if(length(grep("SYMBOL", colnames(data))) != 0){
      data <- data[, - which(colnames(data) == "SYMBOL")]
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
    p <- plot_grid(g1, g2, g3, nrow = 1)
    }
  })

  output$PCA2 <- renderPlot({
    print(cond3_pca())
  })
  output$PCA3_data <- DT::renderDataTable({
    PCA3data()
  })

  output$download_3cond_PCA = downloadHandler(
    filename = function(){
      paste0(download_cond3_dir(), "_PCA-MDS-dendrogram.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        pdf(file, height = 3.5, width = 9)
        print(cond3_pca_plot())
        dev.off()
      })
    }
  )
  #3conditions GOI------------------------------------------------------
  GOI_list2 <- reactive({
    withProgress(message = "Preparing GOI list",{
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
      incProgress(1)
    })
  })

  output$GOI2 <- renderUI({
    if(is.null(d_row_count_matrix2())){
      return(NULL)
    }else{
      withProgress(message = "Preparing GOI list",{
      selectizeInput("GOI2", "genes of interest (GOI)", c(GOI_list2()),multiple = TRUE, options = list(delimiter = " ", create = T))
        })
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
        ht <- Heatmap(data.z, name = "z-score",column_order = colnames(data.z),
                      clustering_method_columns = 'ward.D2',
                      show_row_names = T, show_row_dend = F)
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
                    title = element_text(size = 10),text = element_text(size = 20))
            + scale_fill_manual(values=c("gray", "#4dc4ff", "#ff8082")))
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
        if ((length(rowlist) > 81) && (length(rowlist) <= 200))
        {pdf_hsize <- 22.5
        pdf_wsize <- 22.5}
        if ((length(rowlist) > 64) && (length(rowlist) <= 81))
        {pdf_hsize <- 20.25
        pdf_wsize <- 20.25}
        if ((length(rowlist) > 49) && (length(rowlist) <= 64))
        {pdf_hsize <- 18
        pdf_wsize <- 18}
        if ((length(rowlist) > 36) && (length(rowlist) <= 49))
        {pdf_hsize <- 15.75
        pdf_wsize <- 15.75}
        if ((length(rowlist) > 25) && (length(rowlist) <= 36))
        {pdf_hsize <- 13.5
        pdf_wsize <- 13.5}
        if ((length(rowlist) > 16) && (length(rowlist) <= 25))
        {pdf_hsize <- 11.5
        pdf_wsize <- 11.5}
        if ((length(rowlist) > 12) && (length(rowlist) <= 16))
        {pdf_hsize <- 9
        pdf_wsize <- 9}
        if ((length(rowlist) > 9) && (length(rowlist) <= 12))
        {pdf_hsize <- 7.5
        pdf_wsize <- 9}
        if ((length(rowlist) > 6) && (length(rowlist) <= 9))
        {pdf_hsize <- 7.5
        pdf_wsize <- 6.75}
        if ((length(rowlist) > 4) && (length(rowlist) <= 6))
        {pdf_hsize <- 6
        pdf_wsize <- 9}
        if (length(rowlist) == 4)
        {pdf_hsize <- 6
        pdf_wsize <- 6}
        if (length(rowlist) == 3)
        {pdf_hsize <- 3
        pdf_wsize <- 9}
        if (length(rowlist) == 2)
        {pdf_hsize <- 3
        pdf_wsize <- 6}
        if (length(rowlist) == 1)
        {pdf_hsize <- 3
        pdf_wsize <- 3}
        if (length(rowlist) > 200)
        {pdf_hsize <- 30
        pdf_wsize <- 30}
        pdf(file, height = pdf_hsize, width = pdf_wsize)
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
        pdf(file, height = 10, width = 7)
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
  enrichment3_1_1 <- reactive({
    if(!is.null(input$Gene_set2) && input$Species2 != "not selected"){
    data4 <- data_3degcount2_1()
    data3 <- data_3degcount1_1()
    cnet_list <- list()
    if(is.null(data4)){
      return(NULL)
    }else{
      if(input$Gene_set2 != "MSigDB Hallmark" && input$Gene_set2 != "Transcription factor targets"){
        if(input$Gene_set2 == "KEGG"){
          withProgress(message = "KEGG dotplot",{
          formula_res <- try(compareCluster(ENTREZID~sig, data=data4,fun="enrichKEGG", organism =org_code2()), silent = T)
          incProgress()
          })
          }
        if(input$Gene_set2 == "GO"){
          withProgress(message = "GO dotplot",{
          formula_res <- try(compareCluster(ENTREZID~sig, data=data4,fun="enrichGO", OrgDb=org2()), silent = T)
          incProgress()
          })
          }
        if (class(formula_res) == "try-error") formula_res <- NA
        return(formula_res)
      }else{
        withProgress(message = "dotplot",{
        H_t2g <- Hallmark_cond3()
        for (name in unique(data3$sig)) {
          if (name != "NS"){
        em <- enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
        if (length(as.data.frame(em)$ID) == 0) {
          cnet1 <- NULL
        } else {
          cnet1 <- as.data.frame(setReadable(em, org2(), 'ENTREZID'))
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
  })

  keggEnrichment2_1 <- reactive({
    if(!is.null(input$Gene_set2) && input$Species2 != "not selected"){
    data4 <- data_3degcount2_1()
    data3 <- data_3degcount1_1()
    if(is.null(data4)){
      return(NULL)
    }else{
    cnet_list <- list()
    cnet_list2 <- list()
    if(input$Gene_set2 != "MSigDB Hallmark" && input$Gene_set2 != "Transcription factor targets"){
          formula_res <- enrichment3_1_1()
          if ((length(as.data.frame(formula_res)$Description) == 0) ||
              length(which(!is.na(unique(as.data.frame(formula_res)$qvalue)))) == 0) {
            d <- NULL
          } else{
            formula_res <- clusterProfiler.dplyr::filter(formula_res, !is.na(qvalue))
            d <- as.grob(dotplot(formula_res, showCategory=5, color ="qvalue" ,font.size=10))
          }
        for (name in unique(data3$sig)) {
          if (name != "NS"){
            if(input$Gene_set2 == "KEGG"){
              withProgress(message = "KEGG cnet plot",{
              kk1 <- enrichKEGG(data4$ENTREZID[data4$sig == name], organism =org_code2())
              incProgress()
              })
            }
            if(input$Gene_set2 == "GO"){
              withProgress(message = "GO cnet plot",{
              kk1 <- enrichGO(data4$ENTREZID[data4$sig == name], OrgDb=org2())
              incProgress()
              })
            }
              if (is.null(kk1)) {
                cnet1 <- NULL
              } else cnet1 <- setReadable(kk1, org2(), 'ENTREZID')
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
            
            H_t2g <- Hallmark_cond3()
        for (name in unique(data3$sig)) {
          if (name != "NS"){
            em <- enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) == 0) {
              cnet1 <- NULL
            } else {
              cnet1 <- as.data.frame(setReadable(em, org2(), 'ENTREZID'))
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
          d <- as.grob(ggplot(data, aes(x = Cluster,y=reorder(Description, GeneRatio)))+
                          geom_point(aes(color=qvalue,size=GeneRatio)) +
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
            } else cnet1 <- setReadable(em, org2(), 'ENTREZID')
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
    })
  
  output$keggenrichment2_1 <- renderPlot({
    if(!is.null(input$Gene_set2) && input$Species2 != "not selected"){
      p <- keggEnrichment2_1()
      if(!is.null(p)){
        print(p) 
      }
    }
  })

  #3conditions enrichment_2 ------------------------------------------------------------------------------
    Hallmark_cond3 <- reactive({
      if(input$Species2 != "not selected"){
      switch (input$Species2,
              "Mus musculus" = species <- "Mus musculus",
              "Homo sapiens" = species <- "Homo sapiens",
              "Rattus norvegicus" = species <- "Rattus norvegicus",
              "Xenopus laevis" = species <- "Xenopus laevis",
              "Drosophila melanogaster" = species <- "Drosophila melanogaster",
              "Caenorhabditis elegans" = species <- "Caenorhabditis elegans")
        if(input$Gene_set2 == "MSigDB Hallmark"){
      H_t2g <- msigdbr(species = species, category = "H") %>%
        dplyr::select(gs_name, entrez_gene)
      }
      if(input$Gene_set2 == "Transcription factor targets"){
        H_t2g <- msigdbr(species = species, category = "C3")
        H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "TFT:GTRD" | gs_subcat == "TFT:TFT_Legacy") %>%
          dplyr::select(gs_name, entrez_gene)
      }
      return(H_t2g)
      }else return(NULL)
    })

    enrichment3_2_1 <- reactive({
      if(!is.null(input$Gene_set2) && input$Species2 != "not selected"){
    data4 <- data_3degcount2_2()
    data3 <- data_3degcount1_2()
    cnet_list <- list()
    if(is.null(data4)){
      return(NULL)
    }else{
      if(input$Gene_set2 != "MSigDB Hallmark" && input$Gene_set2 != "Transcription factor targets"){
        if(input$Gene_set2 == "KEGG"){
          formula_res <- try(compareCluster(ENTREZID~sig, data=data4,fun="enrichKEGG", organism =org_code2()), silent = T)
        }
        if(input$Gene_set2 == "GO"){
          formula_res <- try(compareCluster(ENTREZID~sig, data=data4,fun="enrichGO", OrgDb=org2()), silent = T)
        }
        if (class(formula_res) == "try-error") formula_res <- NA
        return(formula_res)
      }else{
        H_t2g <- Hallmark_cond3()
        for (name in unique(data3$sig)) {
          if (name != "NS"){
            em <- enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) == 0) {
              cnet1 <- NULL
            } else {
              cnet1 <- as.data.frame(setReadable(em, org2(), 'ENTREZID'))
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
      }
    }
    }
  })
    keggEnrichment2_2 <- reactive({
      if(!is.null(input$Gene_set2) && input$Species2 != "not selected"){
    data4 <- data_3degcount2_2()
    data3 <- data_3degcount1_2()
    if(is.null(data4)){
      return(NULL)
    }else{
    cnet_list <- list()
    cnet_list2 <- list()
    if(input$Gene_set2 != "MSigDB Hallmark" && input$Gene_set2 != "Transcription factor targets"){
        formula_res <- enrichment3_2_1()
        if ((length(as.data.frame(formula_res)$Description) == 0) ||
            length(which(!is.na(unique(as.data.frame(formula_res)$qvalue)))) == 0) {
          d <- NULL
        } else{
          formula_res <- clusterProfiler.dplyr::filter(formula_res, !is.na(qvalue))
          d <- as.grob(dotplot(formula_res, showCategory=5, color ="qvalue" ,font.size=10))
        }
        for (name in unique(data3$sig)) {
          if (name != "NS"){
            if(input$Gene_set2 == "KEGG"){
              kk1 <- enrichKEGG(data4$ENTREZID[data4$sig == name], organism =org_code2())
            }
            if(input$Gene_set2 == "GO"){
              kk1 <- enrichGO(data4$ENTREZID[data4$sig == name], OrgDb=org2())
            }
            if (is.null(kk1)) {
              cnet1 <- NULL
            } else cnet1 <- setReadable(kk1, org2(), 'ENTREZID')
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
        H_t2g <- Hallmark_cond3()
        for (name in unique(data3$sig)) {
          if (name != "NS"){
            em <- enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) == 0) {
              cnet1 <- NULL
            } else {
              cnet1 <- as.data.frame(setReadable(em, org2(), 'ENTREZID'))
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
          d <- as.grob(ggplot(data, aes(x = Cluster,y=reorder(Description, GeneRatio)))+
                         geom_point(aes(color=qvalue,size=GeneRatio)) +
                         scale_color_continuous(low="red", high="blue",
                                                guide=guide_colorbar(reverse=TRUE)) +
                         scale_size(range=c(3, 8))+ theme_dose(font.size=8)+ylab(NULL))
        }}

        for (name in unique(data3$sig)) {
          if (name != "NS"){
            em <- enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) == 0) {
              cnet1 <- NULL
            } else cnet1 <- setReadable(em, org2(), 'ENTREZID')
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
      return(p)
    }
    }else return(NULL)
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
  enrichment3_3_1 <- reactive({
    if(!is.null(input$Gene_set2) && input$Species2 != "not selected"){
    data4 <- data_3degcount2_3()
    data3 <- data_3degcount1_3()
    if(is.null(data4)){
      return(NULL)
    }else{
    cnet_list <- list()
    if(input$Gene_set2 != "MSigDB Hallmark" && input$Gene_set2 != "Transcription factor targets"){
        if(input$Gene_set2 == "KEGG"){
          formula_res <- try(compareCluster(ENTREZID~sig, data=data4,fun="enrichKEGG", organism =org_code2()), silent = T)
        }
        if(input$Gene_set2 == "GO"){
          formula_res <- try(compareCluster(ENTREZID~sig, data=data4,fun="enrichGO", OrgDb=org2()), silent = T)
        }
        if (class(formula_res) == "try-error") formula_res <- NA
        return(formula_res)
      }else{
        H_t2g <- Hallmark_cond3()
        for (name in unique(data3$sig)) {
          if (name != "NS"){
            em <- enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) == 0) {
              cnet1 <- NULL
            } else {
              cnet1 <- as.data.frame(setReadable(em, org2(), 'ENTREZID'))
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
      }
    }
    }
  })
  keggEnrichment2_3 <- reactive({
    if(!is.null(input$Gene_set2) && input$Species2 != "not selected"){
    data4 <- data_3degcount2_3()
    data3 <- data_3degcount1_3()
    if(is.null(data4)){
      return(NULL)
    }else{
    cnet_list <- list()
    cnet_list2 <- list()
    if(input$Gene_set2 != "MSigDB Hallmark" && input$Gene_set2 != "Transcription factor targets"){
        formula_res <- enrichment3_3_1()
        if(!is.null(formula_res)){
        if ((length(as.data.frame(formula_res)$Description) == 0) ||
            length(which(!is.na(unique(as.data.frame(formula_res)$qvalue))))==0) {
          d <- NULL
        } else{
          formula_res <- clusterProfiler.dplyr::filter(formula_res, !is.na(qvalue))
          d <- as.grob(dotplot(formula_res, showCategory=5, color ="qvalue" ,font.size=10))
        }}else d <- NULL
        for (name in unique(data3$sig)) {
          if (name != "NS"){
            if(input$Gene_set2 == "KEGG"){
              kk1 <- enrichKEGG(data4$ENTREZID[data4$sig == name], organism =org_code2())
            }
            if(input$Gene_set2 == "GO"){
              kk1 <- enrichGO(data4$ENTREZID[data4$sig == name], OrgDb=org2())
            }
            if (is.null(kk1)) {
              cnet1 <- NULL
            } else cnet1 <- setReadable(kk1, org2(), 'ENTREZID')
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
        H_t2g <- Hallmark_cond3()
        for (name in unique(data3$sig)) {
          if (name != "NS"){
            em <- enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) == 0) {
              cnet1 <- NULL
            } else {
              cnet1 <- as.data.frame(setReadable(em, org2(), 'ENTREZID'))
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
          d <- as.grob(ggplot(data, aes(x = Cluster,y=reorder(Description, GeneRatio)))+
                         geom_point(aes(color=qvalue,size=GeneRatio)) +
                         scale_color_continuous(low="red", high="blue",
                                                guide=guide_colorbar(reverse=TRUE)) +
                         scale_size(range=c(3, 8))+ theme_dose(font.size=8)+ylab(NULL))
        }}

        for (name in unique(data3$sig)) {
          if (name != "NS"){
            em <- enricher(data4$ENTREZID[data4$sig == name], TERM2GENE=H_t2g, pvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) == 0) {
              cnet1 <- NULL
            } else cnet1 <- setReadable(em, org2(), 'ENTREZID')
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
      return(p)
    }
    }else return(NULL)
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
    selectInput('Gene_set2', 'Gene Set', c("KEGG", "GO", "MSigDB Hallmark", "Transcription factor targets"))
  })

  output$enrichment3_result_1 <- DT::renderDataTable({
    as.data.frame(enrichment3_1_1())
  })
  output$enrichment3_result_2 <- DT::renderDataTable({
    as.data.frame(enrichment3_2_1())
  })
  output$enrichment3_result_3 <- DT::renderDataTable({
    as.data.frame(enrichment3_3_1())
  })

  output$download_cond3_enrichment_table1 = downloadHandler(
    filename = function() {
      paste(download_cond3_dir(), paste0(input$Gene_set2,"-enrichment1.txt"), sep="_")
    },
    content = function(file){write.table(as.data.frame(enrichment3_1_1()), file, row.names = F, sep = "\t", quote = F)}
  )
  output$download_cond3_enrichment_table2 = downloadHandler(
    filename = function() {
      paste(download_cond3_dir(), paste0(input$Gene_set2,"-enrichment2.txt"), sep="_")
    },
    content = function(file){write.table(as.data.frame(enrichment3_2_1()), file, row.names = F, sep = "\t", quote = F)}
  )
  output$download_cond3_enrichment_table3 = downloadHandler(
    filename = function() {
      paste(download_cond3_dir(), paste0(input$Gene_set2,"-enrichment3.txt"), sep="_")
    },
    content = function(file){write.table(as.data.frame(enrichment3_3_1()), file, row.names = F, sep = "\t", quote = F)}
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
        pdf(file, height = 12, width = 15)
        print(plot_grid(p1, p2, p3, nrow =3))
        dev.off()
        incProgress(1)
      })
    }
  )
  
  #normalized count analysis
  #norm_count_input-------------------
  org3 <- reactive({
    if(input$Species3 != "not selected"){
      switch (input$Species3,
              "Mus musculus" = org <- org.Mm.eg.db,
              "Homo sapiens" = org <- org.Hs.eg.db,
              "Rattus norvegicus" = org <- org.Rn.eg.db,
              "Xenopus laevis" = org <- org.Xl.eg.db,
              "Drosophila melanogaster" = org <- org.Dm.eg.db,
              "Caenorhabditis elegans" = org <- org.Ce.eg.db)
      return(org)
    }
  })
  org_code3 <- reactive({
    if(input$Species3 != "not selected"){
      switch (input$Species3,
              "Mus musculus" = org_code <- "mmu",
              "Homo sapiens" = org_code <- "hsa",
              "Rattus norvegicus" = org_code <- "rno",
              "Xenopus laevis" = org_code <- "xla",
              "Drosophila melanogaster" = org_code <- "dme",
              "Caenorhabditis elegans" = org_code <- "cel")
      return(org_code)
    }
  })

  norm_count_input <- reactive({
    withProgress(message = "Importing normalized count matrix, please wait",{
      if (input$data_file_type3 == "Row5"){
        tmp <- input$file7$datapath
        if(is.null(input$file7) && input$goButton3 > 0 )  tmp = "data/example7.txt"
        if(is.null(tmp)) {
          return(NULL)
        }else{
          if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
          if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
          if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
          return(df)
        }
      }else{
        tmp <- input$file8$datapath
        if(is.null(input$file8) && input$goButton3 == 0) return(NULL)
        if(is.null(input$file8) && input$goButton3 > 0 )  tmp = "data/example7.txt"
        if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
        if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
        if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
        return(df)
      }
      incProgress(1)
    })
  })
  norm_metadata <- reactive({
    if (input$data_file_type3 == "Row5"){
      return(NULL)
    }else{
      tmp <- input$file9$datapath
      if(is.null(input$file9) && input$goButton3 == 0) return(NULL)
      if(is.null(input$file9) && input$goButton3 > 0 )  tmp = "data/example9.csv"
      if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
      if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
      if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
      rownames(df) <- gsub("-",".",rownames(df))
      return(df)
    }
  })
  gene_list <- reactive({
    data <- input$file10$datapath
    if(is.null(input$file10) && input$goButton3 == 0) return(NULL)
    if(is.null(input$file10) && input$goButton3 > 0 )  data = "data/enrich_example.txt"
    if(tools::file_ext(data) == "xlsx") df <- read.xls(data, header=TRUE, row.names = 1)
    if(tools::file_ext(data) == "csv") df <- read.csv(data, header=TRUE, sep = ",", row.names = 1)
    if(tools::file_ext(data) == "txt") df <- read.table(data, header=TRUE, sep = "\t", row.names = 1)
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
          if(is.null(gene_list)){
            return(row)
          }else{
            row <- merge(row,gene_list, by=0)
            rownames(row) <- row$Row.names
            row <- row[,-1]
            row <- row[, - which(colnames(row) == "gene")]
            return(row)
          }
        }
      }else{
        meta <- norm_metadata()
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
  normPCAdata <- reactive({
    data <- d_norm_count_matrix_cutofff()
    if(is.null(data)){
      return(NULL)
    }else{
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
  })

  norm_pca_plot <- reactive({
    data <- d_norm_count_matrix_cutofff()
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
  })

  output$download_norm_PCA = downloadHandler(
    filename = function(){
      paste0(download_norm_dir(), "PCA-MDS-dendrogram.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        pdf(file, height = 3.5, width = 9)
        print(norm_pca_plot())
        dev.off()
      })
    }
  )

  output$norm_PCA <- renderPlot({
    withProgress(message = "Clustering",{
    if(is.null(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      print(norm_pca_plot())
    }
    })
  })

  output$norm_PCA_table <- DT::renderDataTable({
    normPCAdata()
  })
  output$d_norm_count_cutoff <- DT::renderDataTable({
    d_norm_count_cutoff_uniqueID()
  })

  output$download_norm_pca_table = downloadHandler(
    filename = function() {
      paste0(download_norm_dir(), "PCA_table.txt")
    },
    content = function(file){write.table(normPCAdata(), file, row.names = T, sep = "\t", quote = F)}
  )

  norm_heat <- reactive({
    if(input$normHeat != "OFF"){
    data <- d_norm_count_matrix_cutofff()
    withProgress(message = "Heatmap",{
      if(is.null(data)){
        return(NULL)
      }else{
        data.z <-genescale(data, axis=1, method="Z")
        data.z <- na.omit(data.z)
        ht <- Heatmap(data.z, name = "z-score",column_order = colnames(data.z),
                      clustering_method_columns = 'ward.D2',
                      show_row_names = F, show_row_dend = F)
        return(ht)
      }
    })
    }else return(NULL)
  })

  output$norm_heatmap <- renderPlot({
    data <- d_norm_count_matrix_cutofff()
    withProgress(message = "Heatmap",{
    if(is.null(data)){
      return(NULL)
    }else{
      print(norm_heat())
    }
    })
  })

  output$download_norm_heatmap = downloadHandler(
    filename = function(){
      paste0(download_norm_dir(), "heatmap.pdf")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        pdf(file, height = 9, width = 3.5)
        print(norm_heat())
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
    withProgress(message = "Preparing GOI list",{
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
  })

  output$GOI3 <- renderUI({
    if(is.null(d_norm_count_matrix_cutofff())){
      return(NULL)
    }else{
      withProgress(message = "Preparing GOI list",{
        selectizeInput("GOI3", "genes of interest (GOI)", c(GOI_list3()),multiple = TRUE, options = list(delimiter = " ", create = T))
         })
    }
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
      ht <- Heatmap(data.z, name = "z-score",column_order = colnames(data.z),
                    clustering_method_columns = 'ward.D2',
                    show_row_names = T, show_row_dend = F)
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
    count <- d_norm_count_cutoff_uniqueID()
    data <- norm_GOIcount()
    if(is.null(data)){
      p <- NULL
    }else{
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
        if ((length(rowlist) > 81) && (length(rowlist) <= 200))
        {pdf_hsize <- 22.5
        pdf_wsize <- 22.5}
        if ((length(rowlist) > 64) && (length(rowlist) <= 81))
        {pdf_hsize <- 20.25
        pdf_wsize <- 20.25}
        if ((length(rowlist) > 49) && (length(rowlist) <= 64))
        {pdf_hsize <- 18
        pdf_wsize <- 18}
        if ((length(rowlist) > 36) && (length(rowlist) <= 49))
        {pdf_hsize <- 15.75
        pdf_wsize <- 15.75}
        if ((length(rowlist) > 25) && (length(rowlist) <= 36))
        {pdf_hsize <- 13.5
        pdf_wsize <- 13.5}
        if ((length(rowlist) > 16) && (length(rowlist) <= 25))
        {pdf_hsize <- 11.5
        pdf_wsize <- 11.5}
        if ((length(rowlist) > 12) && (length(rowlist) <= 16))
        {pdf_hsize <- 9
        pdf_wsize <- 9}
        if ((length(rowlist) > 9) && (length(rowlist) <= 12))
        {pdf_hsize <- 7.5
        pdf_wsize <- 9}
        if ((length(rowlist) > 6) && (length(rowlist) <= 9))
        {pdf_hsize <- 7.5
        pdf_wsize <- 6.75}
        if ((length(rowlist) > 4) && (length(rowlist) <= 6))
        {pdf_hsize <- 6
        pdf_wsize <- 9}
        if (length(rowlist) == 4)
        {pdf_hsize <- 6
        pdf_wsize <- 6}
        if (length(rowlist) == 3)
        {pdf_hsize <- 3
        pdf_wsize <- 9}
        if (length(rowlist) == 2)
        {pdf_hsize <- 3
        pdf_wsize <- 6}
        if (length(rowlist) == 1)
        {pdf_hsize <- 3
        pdf_wsize <- 3}
        if (length(rowlist) > 200)
        {pdf_hsize <- 30
        pdf_wsize <- 30}
        pdf(file, height = pdf_hsize, width = pdf_wsize)
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
        pdf(file, height = 10, width = 7)
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
                      show_row_names = F)
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

  output$norm_kmeans_count_table <- DT::renderDataTable({
    norm_kmeans_cluster()
  })

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
        pdf(file, height = 10, width = 7)
        print(norm_kmeans())
        dev.off()
        incProgress(1)
      })
    }
  )

  #venn diagram ------------------

  output$select_file1 <- renderUI({
    selectInput("selectfile1", "gene_list", choices = c("not selected", overlap_list()), selected = "not selected",multiple = F)
  })
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
        df[["day0"]] <-  c(rownames(read.table("data/example11.txt",header = T, row.names = 1,sep="\t")))
        df[["day1"]] <- c(rownames(read.table("data/example12.txt",header = T, row.names = 1,sep="\t")))
        df[["day5"]] <- c(rownames(read.table("data/example13.txt",header = T, row.names = 1,sep="\t")))
        return(df)
      }
      return(NULL)
    }else{
    for(nr in 1:length(input$files[, 1])){
      df <- read.table(file = input$files[[nr, 'datapath']], header = T ,sep="\t", row.names = 1)
      name <- c(name, gsub("\\..+$", "", input$files[nr,]$name))
      upload[[nr]] <- c(rownames(df))
    }
    names(upload) <- name
    return(upload)
    }
  })

  output$venn <- renderPlot({
    if(is.null(files_table())){
      return(NULL)
    }else{
    gene_list <- files_table()
    venn::venn(gene_list, ilabels = TRUE, zcolor = "style", ilcs = 1.0, sncs = 1.0 )
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

  count_for_venn <- reactive({
    tmp <- input$file_for_venn$datapath
    if(is.null(input$file_for_venn) && input$goButton_venn == 0){
      return(NULL)
    }else{
      if(is.null(input$file_for_venn) && input$goButton_venn > 0 )  tmp = "data/example7.txt"
      if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE, row.names = 1)
      if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
      if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
      return(df)
    }
  })

  overlap_extract <- reactive({
    cluster_file <- overlap_table2()
    rownames(cluster_file) <- cluster_file$Gene
    data <- count_for_venn()
    if(is.null(data) || is.null(cluster_file)){
      return(NULL)
    }else{
      cluster_name <- input$selectfile1
      clusterCount <- dplyr::filter(cluster_file, Group == cluster_name)
      clusterCount <- merge(clusterCount, data, by=0)
      clusterCount <- clusterCount[,-2:-3]
      rownames(clusterCount) <- clusterCount$Row.names
      clusterCount <- clusterCount[,-1]
      return(clusterCount)
    }
  })


  output$venn_result <- renderDataTable({
    overlap_table2()
  })

  output$intersection_count <- renderDataTable({
    overlap_extract()
  })

  output$download_vennplot = downloadHandler(
    filename ="venn_diagram.pdf",
    content = function(file){
      if(is.null(files_table())){
        return(NULL)
      }else{
      withProgress(message = "Preparing download",{
        gene_list <- files_table()
        pdf(file, height = 3, width = 3)
        print((venn::venn(gene_list, ilabels = TRUE, zcolor = "style", ilcs = 1.0, sncs = 1.0 )))
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
        df["day0"] <- list(read.table("data/day0.txt",header = T, row.names = 1))
        df["day1"] <- list(read.table("data/day1.txt",header = T, row.names = 1))
        df["day5"] <- list(read.table("data/day5.txt",header = T, row.names = 1))
        return(df)
      }
      return(NULL)
    }else{
      upload = list()
      name = c()
      for(nr in 1:length(input$countfiles[, 1])){
        df <- read.table(file = input$countfiles[[nr, 'datapath']], header = T ,sep="\t", row.names = 1)
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
      colnames(base) <- gsub("(.y)$", "", colnames(base))
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
      colnames(base_z) <- gsub("(.y)", "", colnames(base_z))

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
                        top_annotation = HeatmapAnnotation(condition = cond))
        }else{
          ht <- Heatmap(base_z, name = "z-score",
                        clustering_method_columns = 'ward.D2',
                        cluster_row_slices = T, show_row_names = F,
                        top_annotation = HeatmapAnnotation(condition = cond))
        }
        incProgress(1)
        return(draw(ht))
      })
    }
    }
  })

  output$intheatmap <- renderPlot({
    print(integrated_heatmap())
  })
  output$integrated_count_table <- renderDataTable({
    integrated_count()
  })
  output$integrated_count_z_table <- renderDataTable({
    integrated_count_z()
  })


  output$download_venn_result = downloadHandler(
    filename ="venn_result.txt",
    content = function(file){write.table(overlap_table2(), file, row.names = F, sep = "\t", quote = F)}
  )

  output$download_intersection_count_table = downloadHandler(
    filename = function(){paste0(paste("intersection",input$selectfile1, sep="_"), paste0(gsub("\\..+$", "", input$file_for_venn), ".txt"))},
    content = function(file){write.table(overlap_extract(), file, row.names = F, sep = "\t", quote = F)}
  )

  output$download_integrated_count_table = downloadHandler(
    filename ="integrated_count_table.txt",
    content = function(file){write.table(integrated_count(), file, row.names = T, sep = "\t", quote = F)}
  )
  output$download_integrated_z_count_table = downloadHandler(
    filename ="integrated_zscored_count_table.txt",
    content = function(file){write.table(integrated_count_z(), file, row.names = T, sep = "\t", quote = F)}
  )
  output$download_integrated_heatmap = downloadHandler(
    filename ="integrated_heatmap.pdf",
    content = function(file){
      withProgress(message = "Preparing download",{
        pdf(file, height = 8, width = 8)
        print(integrated_heatmap())
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
    if(input$Species4 != "not selected"){
      if(input$Gene_set3 != "MSigDB Hallmark" &&  input$Gene_set3 != "Transcription factor targets"){
        return(NULL)
      }else{
        switch (input$Species4,
                "Mus musculus" = species <- "Mus musculus",
                "Homo sapiens" = species <- "Homo sapiens",
                "Rattus norvegicus" = species <- "Rattus norvegicus",
                "Xenopus laevis" = species <- "Xenopus laevis",
                "Drosophila melanogaster" = species <- "Drosophila melanogaster",
                "Caenorhabditis elegans" = species <- "Caenorhabditis elegans")
        if(input$Gene_set3 == "MSigDB Hallmark"){
        H_t2g <- msigdbr(species = species, category = "H") %>%
          dplyr::select(gs_name, entrez_gene)
        }
        if(input$Gene_set3 == "Transcription factor targets"){
          H_t2g <- msigdbr(species = species, category = "C3")
          H_t2g <- H_t2g %>% dplyr::filter(gs_subcat == "TFT:GTRD" | gs_subcat == "TFT:TFT_Legacy") %>%
            dplyr::select(gs_name, entrez_gene)
        }
        return(H_t2g)
      }
    }else return(NULL)
  })
  org4 <- reactive({
    if(input$Species4 != "not selected"){
      switch (input$Species4,
              "Mus musculus" = org <- org.Mm.eg.db,
              "Homo sapiens" = org <- org.Hs.eg.db,
              "Rattus norvegicus" = org <- org.Rn.eg.db,
              "Xenopus laevis" = org <- org.Xl.eg.db,
              "Drosophila melanogaster" = org <- org.Dm.eg.db,
              "Caenorhabditis elegans" = org <- org.Ce.eg.db)
      return(org)
    }
  })
  org_code4 <- reactive({
    if(input$Species4 != "not selected"){
      switch (input$Species4,
              "Mus musculus" = org_code <- "mmu",
              "Homo sapiens" = org_code <- "hsa",
              "Rattus norvegicus" = org_code <- "rno",
              "Xenopus laevis" = org_code <- "xla",
              "Drosophila melanogaster" = org_code <- "dme",
              "Caenorhabditis elegans" = org_code <- "cel")
      return(org_code)
    }
  })

  enrich_input <- reactive({
    tmp <- input$enrich_data_file$datapath
    if(is.null(input$enrich_data_file) && input$goButton4 == 0){
      return(NULL)
    }else{
      if(is.null(input$enrich_data_file) && input$goButton4 > 0 )  tmp = "data/enrich_example.txt"
      if(tools::file_ext(tmp) == "xlsx") df <- read.xls(tmp, header=TRUE)
      if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",")
      if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t")
      return(df)
    }
  })

  output$enrichment_input <- DT::renderDataTable({
    as.data.frame(enrich_input())
  })


  enrich_viewer1 <- reactive({
    data <- enrich_input()
    if(is.null(data) || input$Species4 == "not selected"){
      return(NULL)
    }else{
      df <- data.frame(GeneID = data[,1], Group = data[,2])
      my.symbols <- df$GeneID
      if(str_detect(df$GeneID[1], "ENS")){
        gene_IDs<-AnnotationDbi::select(org4(),keys = my.symbols,
                                        keytype = "ENSEMBL",
                                        columns = c("ENSEMBL","SYMBOL", "ENTREZID"))
        colnames(gene_IDs) <- c("GeneID","SYMBOL", "ENTREZID")
      }else{
        gene_IDs <- AnnotationDbi::select(org4(), keys = my.symbols,
                                          keytype = "SYMBOL",
                                          columns = c("ENTREZID", "SYMBOL"))
        colnames(gene_IDs) <- c("GeneID","ENTREZID")
      }
      gene_IDs <- gene_IDs %>% distinct(GeneID, .keep_all = T)
      data <- merge(df, gene_IDs, by="GeneID")
      return(data)
    }
  })

  enrich_viewer2 <- reactive({
    if(!is.null(input$Gene_set3)){
    data3 <- enrich_viewer1()
    if(is.null(data3)){
      return(NULL)
    }else{
      if(input$Gene_set3 != "MSigDB Hallmark" && input$Gene_set3 != "Transcription factor targets"){
        if(input$Gene_set3 == "KEGG"){
          withProgress(message = "KEGG enrichment analysis",{
            formula_res <- try(compareCluster(ENTREZID~Group, data=data3,
                                              fun="enrichKEGG", organism=org_code4()), silent = T)
            incProgress(1)
          })
        }
        if(input$Gene_set3 == "GO"){
          withProgress(message = "GO enrichment analysis",{
            formula_res <- try(compareCluster(ENTREZID~Group, data=data3,
                                              fun="enrichGO", OrgDb=org4()), silent =T)
            incProgress(1)
          })
        }
        if (class(formula_res) == "try-error") {
          formula_res <- NULL
        }else{
          formula_res <-setReadable(formula_res, org4(), 'ENTREZID')
          return(formula_res)
        }
      }else{
        withProgress(message = "enrichment analysis",{
          H_t2g <- Hallmark_enrich()
          df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
          colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
          for (name in unique(data3$Group)) {
            em <- enricher(data3$ENTREZID[data3$Group == name], TERM2GENE=H_t2g, qvalueCutoff = 0.05)
            if (length(as.data.frame(em)$ID) != 0) {
              if(length(colnames(as.data.frame(em))) == 9){
                cnet1 <- as.data.frame(setReadable(em, org4(), 'ENTREZID'))
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
  })

  # enrichment plot ------------------------------------------------------------------------------
  enrich_keggGO <- reactive({
    if(!is.null(input$Gene_set3)){
    formula_res <- enrich_viewer2()
    if(is.null(formula_res)){
      return(NULL)
    }else{
      if(input$Gene_set3 != "MSigDB Hallmark" && input$Gene_set3 != "Transcription factor targets"){
        if ((length(as.data.frame(formula_res)$Description) == 0) ||
            length(which(!is.na(unique(as.data.frame(formula_res)$qvalue)))) == 0) {
          p1 <- NULL
        } else{
          formula_res <- clusterProfiler.dplyr::filter(formula_res, !is.na(qvalue))
          p1 <- as.grob(dotplot(formula_res, color ="qvalue", font.size = 10))
        }
        p <- plot_grid(p1)
        return(p)
      }
    }
    }
  })

  enrich_H <- reactive({
    if(!is.null(input$Gene_set3)){
    data3 <- enrich_viewer1()
    if(is.null(data3)){
      return(NULL)
    }else{
      if(input$Gene_set3 != "MSigDB Hallmark" && input$Gene_set3 != "Transcription factor targets"){
        return(NULL)
      }else{
        H_t2g <- Hallmark_enrich()
        df <- data.frame(matrix(rep(NA, 10), nrow=1))[numeric(0), ]
        colnames(df) <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", " qvalue", "geneID", "Count", "Group")
        for (name in unique(data3$Group)) {
          em <- enricher(data3$ENTREZID[data3$Group == name], TERM2GENE=H_t2g, qvalueCutoff = 0.05)
          if (length(as.data.frame(em)$ID) != 0) {
            if(length(colnames(as.data.frame(em))) == 9){
              cnet1 <- as.data.frame(setReadable(em, org4(), 'ENTREZID'))
              cnet1$Group <- name
              cnet1 <- cnet1[sort(cnet1$pvalue, decreasing = F, index=T)$ix,]
              if (length(cnet1$pvalue) > 5){
                cnet1 <- cnet1[1:5,]
              }
              df <- rbind(df, cnet1)
            }
          }
        }
        if ((length(df$Description) == 0) || length(which(!is.na(unique(df$qvalue)))) == 0) {
          p1 <- NULL
        } else{
          df["Description"] <- lapply(df["Description"], gsub, pattern="HALLMARK_", replacement = "") 
          df$GeneRatio <- parse_ratio(df$GeneRatio)
          df <- dplyr::filter(df, !is.na(qvalue))
          df <- dplyr::mutate(df, x = paste0(Group, eval(parse(text = "GeneRatio"))))
          df$x <- gsub(":","", df$x)
          df <- dplyr::arrange(df, x)
          idx <- order(df[["x"]], decreasing = FALSE)
          df$Description <- factor(df$Description,
                                   levels=rev(unique(df$Description[idx])))
          p1 <- as.grob(ggplot(df, aes(x = Group,y= Description,color=qvalue,size=GeneRatio))+
                          geom_point() +
                          scale_color_continuous(low="red", high="blue",
                                                 guide=guide_colorbar(reverse=TRUE)) +
                          scale_size(range=c(3, 8))+ theme_dose(font.size=8)+ylab(NULL))
          p <- plot_grid(p1)
          return(p)
        }
      }}
    }
  })

  output$enrichment3 <- renderPlot({
    if(!is.null(input$Gene_set3)){
      if(is.null(enrich_viewer2())){
        return(NULL)
      }else{
        withProgress(message = "Plot results",{
          if(input$Species4 != "not selected"){
            if(input$Gene_set3 != "MSigDB Hallmark" && input$Gene_set3 != "Transcription factor targets"){
              print(enrich_keggGO())
            }else{
              print(enrich_H())
            }
          }
          incProgress(1)
        })
      }
    }
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
    if(!is.null(input$Gene_set3)){
    data <- enrich_viewer1()
    group <- input$which_group
    if(is.null(data) || is.null(group)){
      return(NULL)
    }else{
      data2 <- dplyr::filter(data, Group == group)
      if(input$Gene_set3 != "MSigDB Hallmark" && input$Gene_set3 != "Transcription factor targets"){
        if(input$Gene_set3 == "KEGG"){
          kk1 <- enrichKEGG(data2$ENTREZID, organism =org_code4(),
                            qvalueCutoff = 0.05)
        }
        if(input$Gene_set3 == "GO"){
          kk1 <- enrichGO(data2$ENTREZID, OrgDb = org4(),
                          qvalueCutoff = 0.05)
        }
      }else{
        H_t2g <- Hallmark_enrich()
        kk1 <- try(enricher(data2$ENTREZID, TERM2GENE=H_t2g, qvalueCutoff = 0.05))
        if (class(kk1) == "try-error") kk1 <- NA
      }
      if(length(as.data.frame(kk1)$ID) == 0){
        cnet1 <- NULL
      } else {
        cnet1 <- setReadable(kk1, org4(), 'ENTREZID')
      }
      if (length(as.data.frame(cnet1)$ID) == 0) {
        p2 <- NULL
      } else{
        p2 <- try(as.grob(cnetplot(cnet1,
                                   cex_label_gene = 0.7, cex_label_category = 0.75,showCategory = 5,
                                   cex_category = 0.75, colorEdge = TRUE)+ guides(edge_color = "none")))
        if(length(class(p2)) == 1){
          if(class(p2) == "try-error") p2 <- NULL
        }
      }
      p <- plot_grid(p2)
      return(p)
    }
    }
  })

  output$enrichment4 <- renderPlot({
    if(!is.null(input$Gene_set3)){
      if(is.null(enrich_input())){
        return(NULL)
      }else{
        if(input$Species4 != "not selected"){
          withProgress(message = "cnet plot",{
            p <- enrich2()
            print(p)
            incProgress(1)
          })
        }else return(NULL)
      }
    }
  })

  output$download_enrichment = downloadHandler(
    filename = function(){
      paste(gsub("\\..+$", "", input$enrich_data_file), paste0(input$Gene_set3,".pdf"), sep ="-")
    },
    content = function(file) {
      withProgress(message = "Preparing download",{
        if(input$Gene_set3 != "MSigDB Hallmark" && input$Gene_set3 != "Transcription factor targets"){
          p1 <- enrich_keggGO()
        }else{
          p1 <- enrich_H()
        }
        pdf(file, height = 6, width = 8)
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
        pdf(file, height = 6, width = 6)
        print(plot_grid(p1))
        dev.off()
        incProgress(1)
      })
    }
  )

  output$Gene_set3 <- renderUI({
    selectInput('Gene_set3', 'Gene Set', c("KEGG", "GO", "MSigDB Hallmark", "Transcription factor targets"))
  })

  output$enrichment_result <- DT::renderDataTable({
    as.data.frame(enrich_viewer2())
  })

  output$download_enrichment_table = downloadHandler(
    filename = function() {
      paste(gsub("\\..+$", "", input$enrich_data_file), paste0(input$Gene_set3,"_table.txt"), sep ="-")
    },
    content = function(file){write.table(as.data.frame(enrich_viewer2()), file, row.names = F, sep = "\t", quote = F)}
  )

  #volcano navi------------------------------------------------------
  org5 <- reactive({
    if(input$Species5 != "not selected"){
      switch (input$Species5,
              "mouse" = org <- org.Mm.eg.db,
              "human" = org <- org.Hs.eg.db,
              "rat" = org <- org.Rn.eg.db,
              "fly" = org <- org.Dm.eg.db,
              "worm" = org <- org.Ce.eg.db)
      return(org)
    }
  })
  org_code5 <- reactive({
    if(input$Species5 != "not selected"){
      switch (input$Species5,
              "mouse" = org_code <- "mmu",
              "human" = org_code <- "hsa",
              "rat" = org_code <- "rno",
              "fly" = org_code <- "dme",
              "worm" = org_code <- "cel")
      return(org_code)
    }
  })
  
  degresult <- reactive({
    withProgress(message = "Importing normalized count matrix, please wait",{
      tmp <- input$deg_file1$datapath
      if(is.null(input$deg_file1) && input$goButton5 > 0 )  tmp = "data/DEGexample.txt"
      if(is.null(tmp)) {
        return(NULL)
      }else{
        if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
        if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
        return(df)
      }
      incProgress(1)
    })
  })
  norm_count_input_for_deg <- reactive({
    withProgress(message = "Importing normalized count matrix, please wait",{
      tmp <- input$deg_file2$datapath
      if(is.null(input$deg_file2) && input$goButton5 > 0 )  tmp = "data/day0.txt"
      if(is.null(tmp)) {
        return(NULL)
      }else{
        if(tools::file_ext(tmp) == "csv") df <- read.csv(tmp, header=TRUE, sep = ",", row.names = 1)
        if(tools::file_ext(tmp) == "txt") df <- read.table(tmp, header=TRUE, sep = "\t", row.names = 1)
        return(df)
      }
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
    withProgress(message = "Preparing GOI list",{
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
      withProgress(message = "Preparing GOI list",{
        selectizeInput("degGOI", "genes of interest (GOI)", c(GOI_DEG()),multiple = TRUE, options = list(delimiter = " ", create = T))
      })
    }
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
              axis.text.x= ggplot2::element_text(size = 10),
              axis.text.y= ggplot2::element_text(size = 10),
              text = ggplot2::element_text(size = 10),
              title = ggplot2::element_text(size = 10)) +
        xlab("log2 fold change") + ylab("-log10(padj)") +
        xlim(input$deg_xrange)+
        ylim(c(0, input$deg_yrange))
      if(!is.null(label_data)) {
        if(str_detect(rownames(data)[1], "ENS")){
          if(input$Species5 != "not selected"){
            v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
            v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Unique_ID),
                                              box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1)
          }else{
            v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
            v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Row.names),
                                              box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1)
          }
        }else{
          v <- v + geom_point(data=dplyr::filter(data, color == "GOI"),color="green", size=1)
          v <- v + ggrepel::geom_text_repel(data = dplyr::filter(data, color == "GOI"), mapping = aes(label = Row.names),
                                            box.padding = unit(0.35, "lines"), point.padding = unit(0.3,"lines"), force = 1)
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
        pdf(file, height = 5, width = 5)
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
      ht <- Heatmap(data.z, name = "z-score",column_order = colnames(data.z),
                    clustering_method_columns = 'ward.D2',
                    show_row_names = T, show_row_dend = F)
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
                             add = "jitter", panel.labs = NULL,
                             xlab = FALSE, ylab = "Normalized_count", ylim = c(0, NA))
      p <- (facet(p, facet.by = "Row.names",
                  panel.labs.background = list(fill = "transparent", color = "transparent"),
                  scales = "free", short.panel.labs = T)+
              theme(axis.text.x= element_text(size = 0),
                    axis.text.y= element_text(size = 10),
                    panel.background = element_rect(fill = "transparent", size = 0.5),
                    title = element_text(size = 10),text = element_text(size = 20)))
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
        if ((length(rowlist) > 81) && (length(rowlist) <= 200))
        {pdf_hsize <- 22.5
        pdf_wsize <- 22.5}
        if ((length(rowlist) > 64) && (length(rowlist) <= 81))
        {pdf_hsize <- 20.25
        pdf_wsize <- 20.25}
        if ((length(rowlist) > 49) && (length(rowlist) <= 64))
        {pdf_hsize <- 18
        pdf_wsize <- 18}
        if ((length(rowlist) > 36) && (length(rowlist) <= 49))
        {pdf_hsize <- 15.75
        pdf_wsize <- 15.75}
        if ((length(rowlist) > 25) && (length(rowlist) <= 36))
        {pdf_hsize <- 13.5
        pdf_wsize <- 13.5}
        if ((length(rowlist) > 16) && (length(rowlist) <= 25))
        {pdf_hsize <- 11.5
        pdf_wsize <- 11.5}
        if ((length(rowlist) > 12) && (length(rowlist) <= 16))
        {pdf_hsize <- 9
        pdf_wsize <- 9}
        if ((length(rowlist) > 9) && (length(rowlist) <= 12))
        {pdf_hsize <- 7.5
        pdf_wsize <- 9}
        if ((length(rowlist) > 6) && (length(rowlist) <= 9))
        {pdf_hsize <- 7.5
        pdf_wsize <- 6.75}
        if ((length(rowlist) > 4) && (length(rowlist) <= 6))
        {pdf_hsize <- 6
        pdf_wsize <- 9}
        if (length(rowlist) == 4)
        {pdf_hsize <- 6
        pdf_wsize <- 6}
        if (length(rowlist) == 3)
        {pdf_hsize <- 3
        pdf_wsize <- 9}
        if (length(rowlist) == 2)
        {pdf_hsize <- 3
        pdf_wsize <- 6}
        if (length(rowlist) == 1)
        {pdf_hsize <- 3
        pdf_wsize <- 3}
        if (length(rowlist) > 200)
        {pdf_hsize <- 30
        pdf_wsize <- 30}
        pdf(file, height = pdf_hsize, width = pdf_wsize)
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
        pdf(file, height = 10, width = 7)
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
}







# Run the app
shinyApp(ui, server)
