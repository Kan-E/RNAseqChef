popoverTempate <- 
  '<div class="popover popover-lg" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content"></div></div>'

shinyUI(
  fluidPage(
    tags$head(
        includeScript("navAppend.js")
      ),
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
      footer=p(hr(),p("Need help? Create an issue on", a("Github", href = "https://github.com/Kan-E/RNAseqChef/issues"), 
                      "or", a("contact us", href = "kaneto@kumamoto-u.ac.jp"),".",align="center",width=4)
               ),
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
                 column(12,br(),
                        column(6,br(),
                               h4(strong("Pair-wise DEG")),
                               "Detects and visualizes differentially expressed genes",
                               img(src="Pair-wise_DEG.png", width = 600,height = 300),br(),hr(),
                               h4(strong("3 conditions DEG")),
                               "Detects and visualizes differentially expressed genes by EBSeq multi-comparison analysis",
                               img(src="3cond_DEG.png", width = 600,height = 375),br(),hr(),
                               h4(strong("Multi DEG")),
                               "Detects and visualizes differentially expressed genes by DESeq2 LRT following clustering analysis",
                               img(src="Multi DEG.png", width = 600,height = 682)),
                        column(6,hr(),
                               h4(strong("Venn diagram")),
                               "Detects and visualizes the overlap between DEGs from multiple datasets",
                               img(src="Venn.png", width = 600,height = 400),br(),hr(),
                               h4(strong("Normalized count analysis")),
                               "identifies similar samples and gene expression patterns by clustering methods",
                               img(src="Normalized.png", width = 550,height = 325),hr(),
                               h4(strong("Enrichment viewer")),
                               "determines and visualizes biological functions and promoter motifs of gene sets of interest",
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
                                                span("Select a raw count matrix file"),
                                                span(icon("info-circle"), id = "icon1", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("txt", "csv", "xlsx"),
                                              multiple = FALSE,
                                              width = "80%"),
                                    bsPopover("icon1", "Count matrix format (txt, csv, or xlsx):", 
                                              content=paste(strong("The replication number"), "is represented by", strong("the underline"),".<br>", strong("Do not use it for anything else"),".<br><br>",
                                                            img(src="input_format1.png", width = 400,height = 250)), 
                                              placement = "right",options = list(container = "body")),
                   ),
                   conditionalPanel(condition="input.data_file_type=='Row2'",
                                    fileInput("file1",
                                              strong(
                                                span("Select a raw count matrix file"),
                                                span(icon("info-circle"), id = "icon2", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("txt", "csv", "xlsx"),
                                              multiple = FALSE,
                                              width = "80%"),
                                    fileInput("file2",
                                              "Select a metadata file to define samples for the following analysis",
                                              accept = c("txt", "csv","xlsx"),
                                              multiple = FALSE,
                                              width = "80%"),
                                    bsPopover("icon2", "Metadata format (txt, csv, or xlsx):", 
                                              content=paste("The first column is", strong("the sample name"), "used in the raw count data.<br>", 
                                                            "The second column is", strong("the corresponding sample name"), "that matches the sample name in the first column.<br><br>",
                                                            img(src="input_format2.png", width = 400,height = 400)),
                                              placement = "right",options = list(container = "body")),
                   ),
                   conditionalPanel(condition="input.data_file_type=='Row11'",
                                    fileInput("file11",
                                              strong(
                                                span("Select raw count matrix files"),
                                                span(icon("info-circle"), id = "icon3", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("txt", "csv", "xlsx"),
                                              multiple = TRUE,
                                              width = "80%"),
                                    bsPopover("icon3", "Count matrix format (txt, csv, or xlsx):", 
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
                     column(6, selectInput("Species", "Species", species_list, selected = "not selected"))),
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
                             accept = c("txt", "csv", "xlsx"),
                             multiple = TRUE,
                             width = "80%"),
                   bsPopover("icon4", "Option: Normalized count file (txt, csv, or xlsx):", 
                             content="You can use a normalized count data (e.g. TPM count) for basemean cutoff and boxplot.", 
                             placement = "right",options = list(container = "body")),
                   strong(span("Output plot size setting for pdf (0: default)"),
                          span(icon("info-circle"), id = "pair_pdf_icon", 
                               options = list(template = popoverTempate))),
                   fluidRow(
                     column(5, numericInput("pair_pdf_height", "pdf_height", value = 0, min = 0)),
                     column(5, numericInput("pair_pdf_width", "pdf_width", value = 0, min = 0))
                   ),
                   bsPopover("pair_pdf_icon", "Output plot size setting for pdf (default: 0): ", 
                             content=paste("You can adjust the plot size by using", strong('pdf_height'), "and", strong('pdf_width'), "parameters.<br>", 
                                           "Default size: <br>",strong("Clustering:"), "height = 3.5, width = 9<br>", strong("MA-plot:"), "height = 4, width = 7 <br>",strong("Volcano plot:"), "height = 5, width = 5 <br>",pdfSize_for_GOI,
                                           strong("Enrichment analysis:"), "height = 10, width = 12 <br>"),trigger = "click"), 
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
                   conditionalPanel(condition="input.data_file_type!='Row11'",
                                    fluidRow(column(7),
                                             column(5, downloadButton("download_pair_report", "Download summary"),
                                                    tags$head(tags$style("#download_pair_report{color: red;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                                                              tags$style("
          body {
            padding: 0 !important;
          }"
                                                              )),)
                                    )
                   ),
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
                                column(12, downloadButton("download_pair_PCA", "Download clustering analysis"),
                                       textOutput("not_cond2"),
                                       tags$head(tags$style("#not_cond2{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")))
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
                                column(4, htmlOutput("volcano_x"), htmlOutput("GOIreset_pair")),
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
                                                span("Select a raw count matrix file"),
                                                span(icon("info-circle"), id = "icon5", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("txt", "csv", "xlsx"),
                                              multiple = FALSE,
                                              width = "80%"),
                                    bsPopover("icon5", "Count matrix format (txt, csv, or xlsx):", 
                                              content=paste(strong("The replication number"), "is represented by", strong("the underline"),".<br>", strong("Do not use it for anything else"),".<br><br>",
                                                            img(src="input_format1.png", width = 400,height = 250)), 
                                              placement = "right",options = list(container = "body")),
                   ),
                   conditionalPanel(condition="input.data_file_type2=='Row4'",
                                    fileInput("file5",
                                              strong(
                                                span("Select a raw count matrix file"),
                                                span(icon("info-circle"), id = "icon6", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("txt", "csv", "xlsx"),
                                              multiple = FALSE,
                                              width = "80%"),
                                    fileInput("file6",
                                              label = "Select a metadata file to define samples for the following analysis",
                                              accept = c("txt", "csv", "xlsx"),
                                              multiple = FALSE,
                                              width = "80%"),
                                    bsPopover("icon6", "Metadata format (txt, csv, or xlsx):", 
                                              content=paste("The first column is", strong("the sample name"), "used in the raw count data.<br>", 
                                                            "The second column is", strong("the corresponding sample name"), "that matches the sample name in the first column.<br><br>",
                                                            img(src="input_format2.png", width = 400,height = 400)),
                                              placement = "right",options = list(container = "body")),
                   ),
                   fluidRow(
                     column(6, selectInput("Species2", "Species", species_list, selected = "not selected"))),
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
                             accept = c("txt", "csv", "xlsx"),
                             multiple = FALSE,
                             width = "80%"),
                   bsPopover("icon7", "Option: Normalized count file (txt, csv, or xlsx):", 
                             content="You can use a normalized count data (e.g. TPM count) for basemean cutoff and boxplot.", 
                             placement = "right",options = list(container = "body")),
                   strong(span("Output plot size setting for pdf (0: default)"),
                          span(icon("info-circle"), id = "cond3_pdf_icon", 
                               options = list(template = popoverTempate))),
                   fluidRow(
                     column(5, numericInput("cond3_pdf_height", "pdf_height", value = 0, min = 0)),
                     column(5, numericInput("cond3_pdf_width", "pdf_width", value = 0, min = 0))
                   ),
                   bsPopover("cond3_pdf_icon", "Output plot size setting for pdf (default: 0): ", 
                             content=paste("You can adjust the plot size by using", strong('pdf_height'), "and", strong('pdf_width'), "parameters.<br>", 
                                           "Default size: <br>",strong("Clustering:"), "height = 3.5, width = 9<br>", strong("scatter plot:"), "height = 6, width = 10 <br>",pdfSize_for_GOI,
                                           strong("Enrichment analysis:"), "height = 12, width = 15 <br>"),trigger = "click"), 
                   actionButton("goButton2", "example data (mouse)"),
                   tags$head(tags$style("#goButton{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                             tags$style("
          body {
            padding: 0 !important;
          }"
                             )),
                   fluidRow(column(7),
                            column(5, downloadButton("download_3cond_report", "Download summary"),
                                   tags$head(tags$style("#download_3cond_report{color: red;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                                             tags$style("
          body {
            padding: 0 !important;
          }"
                                             )))
                   )
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
                                                           column(4, downloadButton("download_cond3_d_row_count", "Download defined raw count"))
                                                         ),
                                                         dataTableOutput('D_Row_count_matrix2')
                                         )
                              )
                     ),
                     tabPanel("Result overview",
                              fluidRow(
                                column(12, downloadButton("download_3cond_PCA", "Download clustering analysis"),
                                       textOutput("not_cond3"),
                                       tags$head(tags$style("#not_cond3{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")))
                              ),
                              plotOutput("PCA2"),
                              fluidRow(
                                column(4, downloadButton("download_3cond_scatter1", "Download DEG overview"))
                              ),
                              plotOutput("scatter_1"),
                              plotOutput("scatter_2"),
                              plotOutput("scatter_3"),
                              bsCollapse(id="input_cond3_scatter1",multiple = TRUE,
                                         bsCollapsePanel(title="DEG signature 1:",
                                                         value="Scatter_plot1_panel",
                                                         fluidRow(
                                                           downloadButton("download_3cond_DEG_table1", "Download DEG signature 1")),
                                                         dataTableOutput("DEG_result2_1")
                                                         ),
                                         bsCollapsePanel(title="DEG signature 2:",
                                                         value="Scatter_plot2_panel",
                                                         fluidRow(
                                                           downloadButton("download_3cond_DEG_table2", "Download DEG signature 2")),
                                                         dataTableOutput("DEG_result2_2")
                                         ),
                                         bsCollapsePanel(title="DEG signature 3:",
                                                         value="Scatter_plot3_panel",
                                                         fluidRow(
                                                           downloadButton("download_3cond_DEG_table3", "Download DEG signature 3")),
                                                         dataTableOutput("DEG_result2_3")
                                         )
                                         ),
                              bsCollapse(id="input_collapse_3_DEG",open="cond3_result_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Result:",
                                                         value="cond3_result_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_cond3_result", "Download result"))
                                                         ),
                                                         dataTableOutput("cond3_result")
                                         ),
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
                                column(4, htmlOutput("GOI2"), htmlOutput("cond3_GOI_pair")),
                                column(4, htmlOutput("cond3_xrange"), htmlOutput("GOIreset_cond3")),
                                column(4, htmlOutput("cond3_yrange"))
                              ),
                              fluidRow(
                                column(8, plotOutput("cond3_GOIscatter")),
                                column(4, plotOutput("cond3_GOIheatmap"))
                              ),
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
                                                span("Select a raw count matrix file"),
                                                span(icon("info-circle"), id = "icon8", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("txt", "csv", "xlsx"),
                                              multiple = FALSE,
                                              width = "80%"),
                                    bsPopover("icon8", "Count matrix format (txt, csv, or xlsx):", 
                                              content=paste(strong("The replication number"), "is represented by", strong("the underline"),".<br>", strong("Do not use it for anything else"),".<br><br>",
                                                            img(src="input_format1.png", width = 400,height = 250)), 
                                              placement = "right",options = list(container = "body")),
                   ),
                   conditionalPanel(condition="input.multi_data_file_type=='Row2'",
                                    fileInput("multi_file2",
                                              strong(
                                                span("Select a raw count matrix file"),
                                                span(icon("info-circle"), id = "icon9", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("txt", "csv", "xlsx"),
                                              multiple = FALSE,
                                              width = "80%"),
                                    fileInput("multi_file3",
                                              label = "Select a metadata file to define samples for the following analysis",
                                              accept = c("txt", "csv", "xlsx"),
                                              multiple = FALSE,
                                              width = "80%"),
                                    bsPopover("icon9", "Metadata format (txt, csv, or xlsx):", 
                                              content=paste("The first column is", strong("the sample name"), "used in the raw count data.<br>", 
                                                            "The second column is", strong("the corresponding sample name"), "that matches the sample name in the first column.<br><br>",
                                                            img(src="input_format_multi.png", width = 400,height = 400)),
                                              placement = "right",options = list(container = "body")),
                   ),
                   fluidRow(column(6,  selectInput("FDR_method6", "FDR method", c("BH", "Qvalue", "IHW"), selected = "BH")),
                            column(6, selectInput("Species6", "Species", species_list, selected = "not selected"))),
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
                             accept = c("txt", "csv", "xlsx"),
                             multiple = TRUE,
                             width = "80%"),
                   bsPopover("icon10", "Option: Normalized count file (txt, csv, or xlsx):", 
                             content="You can use a normalized count data (e.g. TPM count) for basemean cutoff and boxplot.", 
                             placement = "right",options = list(container = "body")),
                   strong(span("Output plot size setting for pdf (0: default)"),
                          span(icon("info-circle"), id = "multi_pdf_icon", 
                               options = list(template = popoverTempate))),
                   fluidRow(
                     column(5, numericInput("multi_pdf_height", "pdf_height", value = 0, min = 0)),
                     column(5, numericInput("multi_pdf_width", "pdf_width", value = 0, min = 0))
                   ),
                   bsPopover("multi_pdf_icon", "Output plot size setting for pdf (default: 0): ", 
                             content=paste("You can adjust the plot size by using", strong('pdf_height'), "and", strong('pdf_width'), "parameters.<br>", 
                                           "Default size: <br>",strong("Clustering:"), "height = 9, width = 7<br>", strong("UMAP:"), "height = 3.5, width = 4.7 <br>",pdfSize_for_GOI,
                                           strong("Enrichment analysis:"), "height = 6, width = 8 <br>",strong("cnet plot:"), "height = 6, width = 6 <br>",strong("GSEA:"), "height = 5, width = 7 <br>"),trigger = "click"), 
                   actionButton("goButton6", "example data (mouse)"),
                   tags$head(tags$style("#goButton{color: black;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                             tags$style("
          body {
            padding: 0 !important;
          }"
                             )),
                   fluidRow(column(7),
                            column(5, downloadButton("download_Multi_report", "Download summary"),
                                   tags$head(tags$style("#download_Multi_report{color: red;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"),
                                             tags$style("
          body {
            padding: 0 !important;
          }"
                                             )),)
                   )
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
                              plotOutput("multi_PCA"),
                              fluidRow(
                                column(6, htmlOutput("multi_umap_n"),
                                       downloadButton("download_multi_umap", "Download umap"),
                                       textOutput("multi_umap_error"),
                                       tags$head(tags$style("#multi_umap_error{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")))
                                ),
                              div(
                                plotOutput("multi_umap", height = "100%"),
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
                              bsCollapse(id="multi_collapse_panel1",open="multi_deg_pattern_count_panel",multiple = TRUE,
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
                                              DTOutput("multi_pattern1_count")
                              ),
                              bsCollapsePanel(title= p(span("Boxplot"),span(icon("info-circle"), id = "icon_multi_boxplot", 
                                                                                 options = list(template = popoverTempate))),
                                              bsPopover("icon_multi_boxplot", "Boxplot of GOI:", 
                                                        content=paste("Please select genes in", strong("DEG_pattern_normalized_count_data"),".<br><br>",
                                                                      img(src="multi_GOIboxplot.png", width = 450,height = 640)), 
                                                        placement = "right",options = list(container = "body")),
                                              value="multi_deg_pattern_boxplot_panel",
                                              fluidRow(
                                                column(4, downloadButton("download_deg_pattern_boxplot", "Download boxplot"))
                                              ),
                                              plotOutput("multi_pattern_boxplot")
                              ),
                              bsCollapsePanel(title="Enrichment analysis:",
                                              value="multi_deg_pattern_enrichment_panel",
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
                              )
                              )
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
                              bsCollapse(id="multi_collapse_panel2",open="multi_deg_kmeans_pattern_count_panel",multiple = TRUE,
                              bsCollapsePanel(title="kmeans_result:",
                                              value="multi_deg_kmeans_pattern_panel",
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
                                              DTOutput("multi_pattern2_count")
                              ),
                              bsCollapsePanel(title= p(span("Boxplot"),span(icon("info-circle"), id = "icon_multi_boxplot2", 
                                                                            options = list(template = popoverTempate))),
                                              bsPopover("icon_multi_boxplot2", "Boxplot of GOI:", 
                                                        content=paste("Please select genes in", strong("cluster_normalized_count_data"),".<br><br>",
                                                                      img(src="multi_GOIboxplot2.png", width = 450,height = 640)), 
                                                        placement = "right",options = list(container = "body")),
                                              value="multi_deg_kmeans_boxplot_panel",
                                              fluidRow(
                                                column(4, downloadButton("download_deg_kmeans_boxplot", "Download boxplot"))
                                              ),
                                              plotOutput("multi_kmeans_GOIboxplot")
                              ),
                              bsCollapsePanel(title="Enrichment analysis:",
                                              value="multi_deg_kmeans_pattern_enrichment_panel",
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
                              )
                              )
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
                     strong(
                       span("Select gene list files"),
                       span(icon("info-circle"), id = "icon_venn", 
                            options = list(template = popoverTempate))
                     ),
                     multiple = TRUE,
                     accept = c("txt", "csv", "xlsx")
                   ),
                   bsPopover("icon_venn", "Gene list files (txt, csv, or xlsx):", 
                             content= paste("The first column is", strong("gene name"), ".<br>", 
                             "The second and subsequent columns do not affect the analysis.<br>", 
                             "The maximum number of uploads is",strong("7 files"), ".<br><br>", 
                             img(src="venn_input.png", width = 100,height = 300)),
                             placement = "right",options = list(container = "body")),
                   h4("Input for integrated heatmap"),
                   fileInput(
                     inputId = "countfiles",
                     strong(
                       span("Select normalized count files"),
                       span(icon("info-circle"), id = "icon_venn2", 
                            options = list(template = popoverTempate))
                     ),
                     multiple = TRUE,
                     accept = c("txt", "csv", "xlsx")
                   ),
                   bsPopover("icon_venn2", "Count matrix format (txt, csv, or xlsx):", 
                             content=paste(strong("The replication number"), "is represented by", strong("the underline"),".<br>", strong("Do not use it for anything else"),".<br><br>",
                                           img(src="input_format1.png", width = 400,height = 250)), 
                             placement = "right",options = list(container = "body")),
                   fluidRow(
                     column(6, selectInput("Species7", "Species", species_list, selected = "not selected")),
                     column(6, selectInput(
                       inputId = "pre_zscoring",
                       strong(
                         span("Option: Pre-zscoring"),
                         span(icon("info-circle"), id = "icon_venn3", 
                              options = list(template = popoverTempate))
                       ),
                       multiple = FALSE,choices = c("TRUE", "FALSE"), selected = "TRUE"))
                   ),
                   bsPopover("icon_venn3", "Option: Pre-zscoring:", 
                             content=paste("If True, each count data is z-scored before integrating multiple count data.<br><br>",
                                           img(src="pre-zscoring.png", width = 450,height = 200)), 
                             placement = "right",options = list(container = "body")),
                   strong(span("Output plot size setting for pdf (0: default)"),
                          span(icon("info-circle"), id = "venn_pdf_icon", 
                               options = list(template = popoverTempate))),
                   fluidRow(
                     column(5, numericInput("venn_pdf_height", "pdf_height", value = 0, min = 0)),
                     column(5, numericInput("venn_pdf_width", "pdf_width", value = 0, min = 0))
                   ),
                   bsPopover("venn_pdf_icon", "Output plot size setting for pdf (default: 0): ", 
                             content=paste("You can adjust the plot size by using", strong('pdf_height'), "and", strong('pdf_width'), "parameters.<br>", 
                                           "Default size: <br>",strong("Venn diagram:"), "height = 3, width = 3<br>", strong("Integrated heatmap:"), "height = 8, width = 8 <br>",
                                           strong("Enrichment analysis:"), "height = 6, width = 8 <br>",strong("cnet plot:"), "height = 6, width = 6 <br>"),trigger = "click"), 
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
                     tabPanel("Venn diagram analysis",
                              fluidRow(
                                column(4, downloadButton("download_vennplot", "Download venn diagram"))
                              ),
                              plotOutput("venn"),
                              downloadButton("download_venn_result", "Download venn result"),
                              dataTableOutput("venn_result")
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
                                                         DTOutput("integrated_count_table")
                                         ),
                                         bsCollapsePanel(title="integrated zscored normalized count:",
                                                         value="integrated_z_count_panel",
                                                         downloadButton("download_integrated_z_count_table", "Download integrated zscored count table"),
                                                         dataTableOutput("integrated_count_z_table")
                                         ),
                                         bsCollapsePanel(title= p(span("Boxplot"),span(icon("info-circle"), id = "icon_venn_boxplot", 
                                                options = list(template = popoverTempate))),
                                                bsPopover("icon_venn_boxplot", "Boxplot of GOI:", 
                                                          content=paste("Please select genes in", strong("integrated normalized count table"),".<br><br>",
                                                                        img(src="venn_boxplot.png", width = 450,height = 640)), 
                                                          placement = "right",options = list(container = "body")),
                                                         value="integrated_count_box_panel",
                                                         downloadButton("download_GOIbox_venn", "Download boxplot"),
                                                         plotOutput("GOIbox_venn")
                                         )
                              ),
                     ),
                     tabPanel("Enrichment analysis",
                              fluidRow(
                                column(4, htmlOutput("venn_whichGroup1")),
                                column(4, htmlOutput("Gene_set9")),
                                column(4, downloadButton("download_venn_cluster_enrichment", "Download dot plot"))
                              ),
                              fluidRow(
                                column(4, textOutput("venn_Spe"),
                                       tags$head(tags$style("#venn_Spe{color: red;
                                         font-size: 20px;
                                         font-style: bold;
                                      }")))
                              ),
                              plotOutput("venn_enrichment1"),
                              fluidRow(
                                column(4, htmlOutput("venn_whichGroup2")),
                                column(4, downloadButton("download_venn_enrichment_cnet", "Download cnet plot"))
                              ),
                              plotOutput("venn_enrichment2"),
                              fluidRow(
                                column(4, downloadButton("download_venn_enrichment_table", "Download enrichment result"))
                              ),
                              dataTableOutput('venn_enrichment_result')
                     ),
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
                                                span("Select a normalized count matrix file"),
                                                span(icon("info-circle"), id = "icon11", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("txt", "csv", "xlsx"),
                                              multiple = FALSE,
                                              width = "80%"),
                                    bsPopover("icon11", "Count matrix format (txt, csv, or xlsx):", 
                                              content=paste(strong("The replication number"), "is represented by", strong("the underline"),".<br>", strong("Do not use it for anything else"),".<br><br>",
                                                            img(src="input_format1.png", width = 400,height = 250)), 
                                              placement = "right",options = list(container = "body")),
                   ),
                   conditionalPanel(condition="input.data_file_type3=='Row6'",
                                    fileInput("file8",
                                              strong(
                                                span("Select a normalized count matrix file"),
                                                span(icon("info-circle"), id = "icon12", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("txt", "csv", "xlsx"),
                                              multiple = FALSE,
                                              width = "80%"),
                                    fileInput("file9",
                                              label = "Select a metadata file to define samples for the following analysis",
                                              accept = c("txt", "csv", "xlsx"),
                                              multiple = FALSE,
                                              width = "80%"),
                                    bsPopover("icon12", "Metadata format (txt, csv, or xlsx):", 
                                              content=paste("The first column is", strong("the sample name"), "used in the raw count data.<br>", 
                                                            "The second column is", strong("the corresponding sample name"), "that matches the sample name in the first column.<br><br>",
                                                            img(src="input_format2.png", width = 400,height = 400)), 
                                              placement = "right",options = list(container = "body")),
                   ),
                   fluidRow(
                     column(6, selectInput("Species3", "Species", species_list, selected = "not selected"))),
                   h4("Filter option 1:"),
                   fileInput("file10",
                             label = "Select a gene list file for gene extraction",
                             accept = c("txt", "csv", "xlsx"),
                             multiple = FALSE,
                             width = "80%"),
                   h4("Filter option 2:"),
                   fluidRow(
                     column(4, numericInput("basemean3", "Basemean", min   = 0, max   = NA, value = 0),
                     )
                   ),
                   strong(span("Output plot size setting for pdf (0: default)"),
                          span(icon("info-circle"), id = "norm_pdf_icon", 
                               options = list(template = popoverTempate))),
                   fluidRow(
                     column(5, numericInput("norm_pdf_height", "pdf_height", value = 0, min = 0)),
                     column(5, numericInput("norm_pdf_width", "pdf_width", value = 0, min = 0))
                   ),
                   bsPopover("norm_pdf_icon", "Output plot size setting for pdf (default: 0): ", 
                             content=paste("You can adjust the plot size by using", strong('pdf_height'), "and", strong('pdf_width'), "parameters.<br>", 
                                           "Default size: <br>",strong("Clustering:"), "height = 3.5, width = 9<br>", 
                                           strong("UMAP:"), "height = 3.5, width = 4.7 <br>",pdfSize_for_GOI),trigger = "click"), 
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
                              fluidRow(
                                column(3, downloadButton("download_norm_PCA", "Download Clustering"))
                              ),
                              plotOutput("norm_PCA"),
                              fluidRow(
                                column(6, htmlOutput("norm_umap_n"),
                                       downloadButton("download_norm_umap", "Download umap"),
                                       textOutput("norm_umap_error"),
                                       tags$head(tags$style("#norm_umap_error{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")))
                              ),
                              div(
                                plotOutput("norm_umap", height = "100%"),
                                style = "height: calc(100vh  - 100px)"
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
                                column(4, htmlOutput("GOI3"), htmlOutput("GOIreset_norm")),
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
                              bsCollapse(id="norm_kmeans_collapse_panel",open="norm_kmeans_count",multiple = TRUE,
                              bsCollapsePanel(title="k-means clustering result:",
                                              value="norm_kmeans_count",
                                              downloadButton("download_norm_kmeans_cluster", "Download k-means clustering result"),
                                              DTOutput("norm_kmeans_count_table")
                              ),
                              bsCollapsePanel(title= p(span("Boxplot"),span(icon("info-circle"), id = "icon_norm_kmeans_boxplot", 
                                                                            options = list(template = popoverTempate))),
                                              bsPopover("icon_norm_kmeans_boxplot", "Boxplot of GOI:", 
                                                        content=paste("Please select genes in", strong("k-means clustering result"),".<br><br>",
                                                                      img(src="norm_GOIboxplot.png", width = 450,height = 640)), 
                                                        placement = "right",options = list(container = "body")),
                                              value="kmeans_box_panel",
                                              downloadButton("download_norm_kmeans_box", "Download boxplot"),
                                              plotOutput("norm_kmeans_box")
                              ),
                              bsCollapsePanel(title="cluster count data:",
                                              value="norm_kmeans_extract_count",
                                              fluidRow(
                                                column(4, htmlOutput("norm_select_kmean"))
                                              ),
                                              downloadButton("download_norm_kmeans_extract_count", "Download cluster count data"),
                                              DTOutput("norm_kmeans_extract_table")
                              )
                              )
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
                   bsPopover("icon13", "Gene list format (txt, csv, or xlsx): ", 
                             content=paste(strong("Single file upload:"),
                                           "<br>First column must be", strong("gene name"), "(Gene symbol or ENSEMBL ID).<br>", "Second column must be", strong("group or cluster name"),".<br>",
                                           "You can use result files of venn diagram analysis and k-means clustering as input.<br><br>",
                                           strong("Multiple files upload:"),
                                           "<br>The first column is", strong("gene name"), ".<br>", 
                                           "The second and subsequent columns do not affect the analysis.<br>", 
                                           "File names are used as", strong("group names"),".<br><br>",
                                           img(src="input_format_enrich.png", width = 250,height = 400)), 
                             placement = "right",options = list(container = "body")),
                   fluidRow(
                     column(6, selectInput("Species4", "Species", species_list, selected = "not selected"))),
                   sliderInput("enrich_showCategory", "Most significant pathways",
                              min = 1, max = 20, value = 5,step = 1),
                   strong(span("Output plot size setting for pdf (0: default)"),
                      span(icon("info-circle"), id = "enrich_pdf_icon", 
                           options = list(template = popoverTempate))),
                   fluidRow(
                     column(5, numericInput("enrich_pdf_height", "pdf_height", value = 0, min = 0)),
                     column(5, numericInput("enrich_pdf_width", "pdf_width", value = 0, min = 0))
                   ),
                   bsPopover("enrich_pdf_icon", "Output plot size setting for pdf (default: 0): ", 
                             content=paste("You can adjust the plot size by using", strong('pdf_height'), "and", strong('pdf_width'), "parameters.<br>", 
                                           "Default size: <br>","Dotplot:", "height = 5, width = 6.5 <br>", "cnet plot:","height = 6, width = 6 <br><br>"), trigger = "click"), 
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
                     tabPanel("Input gene list",
                              dataTableOutput('enrichment_input')
                     ),
                     tabPanel("Enrichment analysis",
                              fluidRow(
                                column(4, textOutput("Spe3"),
                                       tags$head(tags$style("#Spe3{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")), downloadButton("download_enrichment", "Download dot plot")),
                                column(4, htmlOutput("Gene_set3")),
                                column(4, htmlOutput("Custom_input"))
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
                     ),
                     tabPanel("Promoter motif analysis",
                              fluidRow(
                                column(4, textOutput("motif_Spe"),
                                       tags$head(tags$style("#motif_Spe{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                                column(3, downloadButton("download_motif_plot", "Download motif plot"))
                              ),
                              fluidRow(
                                column(4, htmlOutput("promoter_upstream")),
                                column(4, htmlOutput("promoter_downstream")),
                                column(4, htmlOutput("promoter_padj"))
                                       ),
                              fluidRow(
                                column(4, actionButton("motifButton", "Start"),
                                       tags$head(tags$style("#motifButton{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"),
                                                 tags$style("
          body {
            padding: 0 !important;
          }"
                                                 ))
                                )
                              ),
                              textOutput("motif_warning"),
                              tags$head(tags$style("#motif_warning{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              plotOutput("motif_plot"),
                              bsCollapse(id="Promoter_motif_collapse_panel",open="motif_result_table",multiple = TRUE,
                                         bsCollapsePanel(title="Motif table:",
                                                         value="motif_result_table",
                                                         downloadButton("download_motif_table", "Download motif enrichment result"),
                                                         DTOutput('motif_result')
                                         ),
                                         bsCollapsePanel(title= p(span("Motif region"),span(icon("info-circle"), id = "icon_promoter_motif_region", 
                                                                                            options = list(template = popoverTempate))),
                                                         bsPopover("icon_promoter_motif_region", "Motif region:", 
                                                                   content=paste("Please select genes in", strong("k-means clustering result"),".<br><br>",
                                                                                 img(src="enrich_motif.png", width = 450,height = 600)), 
                                                                   placement = "right",options = list(container = "body")),
                                                         value="Promoter_motif_region_panel",
                                                         downloadButton("download_promoter_motif_region", "Download motif region"),
                                                         dataTableOutput("promoter_motif_region_table")
                                         )
                              )
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
                                        accept = c("txt", "csv", "xlsx"),
                                        multiple = FALSE,
                                        width = "80%"),
                              bsPopover("icon14", "File format (txt, csv, or xlsx): ", 
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
                                        accept = c("txt", "csv", "xlsx"),
                                        multiple = FALSE,
                                        width = "80%"),
                              fluidRow(
                                column(6, selectInput("Species5", "Species", c("not selected", "human", "mouse", "rat", "fly", "worm"), selected = "not selected"))
                              ),
                              fluidRow(
                                column(4, numericInput("fc4", "Fold Change", min   = 1, max   = NA, value = 2)),
                                column(4, numericInput("fdr4", "FDR", min   = 0, max   = 1, value = 0.05))
                              ),
                              strong(span("Output plot size setting for pdf (0: default)"),
                                     span(icon("info-circle"), id = "volcano_pdf_icon", 
                                          options = list(template = popoverTempate))),
                              fluidRow(
                                column(5, numericInput("volcano_pdf_height", "pdf_height", value = 0, min = 0)),
                                column(5, numericInput("volcano_pdf_width", "pdf_width", value = 0, min = 0))
                              ),
                              bsPopover("volcano_pdf_icon", "Output plot size setting for pdf (default: 0): ", 
                                        content=paste("You can adjust the plot size by using", strong('pdf_height'), "and", strong('pdf_width'), "parameters.<br>", 
                                                      "Default size: <br>","Volcano plot:", "height = 5, width = 5 <br>",pdfSize_for_GOI), trigger = "click"), 
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
                                           column(4, htmlOutput("deg_volcano_x"),htmlOutput("GOIreset_deg")),
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
                 tabPanel("DoRothEA regulon",
                          sidebarLayout(
                            # DoRothEA gene list---------------------------------
                            sidebarPanel(
                              selectInput("dorothea_Species", "Species", c("Homo sapiens", "Mus musculus"), selected = "Homo sapiens"),
                              selectInput("dorothea_TFtype", "TF type", c("activator", "repressor"), selected = "activator"),
                              radioButtons('dorothea_confidence','Confidence:',
                                           c('recommend (A, B, and C)'="recommend",
                                             'All (A, B, C, D, and E)'="All"
                                           ),selected = "recommend"),
                              selectizeInput("dorothea_target_set", label="Target name for TF search", choices = ''),
                              selectizeInput("dorothea_tf_set", label="TF name for target search", choices = '')
                              #sidebarPanel
                            ),
                            
                            # Main Panel -------------------------------------
                            mainPanel(
                              bsCollapse(id="dorothea_collapse_panel",open="dorothea_tf_panel",multiple = TRUE,
                                         bsCollapsePanel(title="TF search:",
                                                         value="dorothea_tf_panel",
                                                         downloadButton("download_dorothea_tf_list", "Download"),
                                                         dataTableOutput("dorothea_tf_list")
                                         ),
                                         bsCollapsePanel(title="Target search:",
                                                         value="dorothea_target_panel",
                                                         downloadButton("download_dorothea_target_list", "Download"),
                                                         dataTableOutput("dorothea_target_list")
                                         )
                              )
                            )
                          ) #sidebarLayout
                 ),
                 tabPanel("Reference",
                          fluidRow(
                            column(12,
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
                                   "- Konopka T (2022). _umap: Uniform Manifold Approximation and Projection_. R package version 0.2.8.0, <https://CRAN.R-project.org/package=umap>.", br(),
                                   "- T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141",br(),
                                   "- Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015 31(4):608-609",br(),
                                   "- Dolgalev I (2022). _msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format_. R
  package version 7.5.1, <https://CRAN.R-project.org/package=msigdbr>.", br(),
                                   "- Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. 'Benchmark and integration of resources for the estimation of human
  transcription factor activities.' Genome Research. 2019. DOI: 10.1101/gr.240663.118.", br(),
                                   "- Hervé Pagès, Marc Carlson, Seth Falcon and Nianhua Li (2020). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.52.0. https://bioconductor.org/packages/AnnotationDbi",br(),
                                   "- Marc Carlson (2020). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.12.0.",br(),
                                   "- Marc Carlson (2020). org.Mm.eg.db: Genome wide annotation for Mouse. R package version 3.12.0.",br(),
                                   "- Marc Carlson (2020). org.Rn.eg.db: Genome wide annotation for Rat. R package version 3.12.0.",br(),
                                   "- Marc Carlson (2020). org.Xl.eg.db: Genome wide annotation for Xenopus. R package version 3.12.0.",br(),
                                   "- Marc Carlson (2020). org.Dm.eg.db: Genome wide annotation for Fly. R package version 3.12.0.",br(),
                                   "- Marc Carlson (2020). org.Ce.eg.db: Genome wide annotation for Worm. R package version 3.12.0.",br(),
                                   "- Marc Carlson (2022). org.Bt.eg.db: Genome wide annotation for Bovine. R package version 3.15.0.",br(),
                                   "- Marc Carlson (2022). org.Cf.eg.db: Genome wide annotation for Canine. R package version 3.15.0.",br(),
                                   "- Marc Carlson (2022). org.Dr.eg.db: Genome wide annotation for Zebrafish. R package version 3.15.0.",br(),
                                   "- Marc Carlson (2022). org.Gg.eg.db: Genome wide annotation for Chicken. R package version 3.15.0.",br(),
                                   "- Marc Carlson (2022). org.Mmu.eg.db: Genome wide annotation for Rhesus. R package version 3.15.0.",br(),
                                   "- Marc Carlson (2022). org.Pt.eg.db: Genome wide annotation for Chimp. R package version 3.15.0.",br(),
                                   "- Marc Carlson (2022). org.Sc.sgd.db: Genome wide annotation for Yeast. R package version 3.15.0.",br(),
                                   "- Marc Carlson (2022). org.Ss.eg.db: Genome wide annotation for Pig. R package version 3.15.0.",br(),
                                   "- Morgan M, Shepherd L (2022). AnnotationHub: Client to access AnnotationHub resources. R package version 3.4.0.",br(),
                                   "- R. Gentleman, V. Carey, W. Huber and F. Hahne (2021). genefilter: methods for filtering genes from high-throughput experiments. R package version 1.72.1.",br(),
                                   "- Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.",br(),
                                   "- H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.", br(),
                                   "- Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.4.0. https://CRAN.R-project.org/package=ggpubr",br(),
                                   "- Adrian Dusa (2021). venn: Draw Venn Diagrams. R package version 1.10. https://CRAN.R-project.org/package=venn",br(),
                                   "- Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.7. https://CRAN.R-project.org/package=dplyr",br(),
                                   "- Hadley Wickham (2021). tidyr: Tidy Messy Data. R package version 1.1.3. https://CRAN.R-project.org/package=tidyr",br(),
                                   "- Machlab D, Burger L, Soneson C, Rijli FM, Schübeler D, Stadler MB. monaLisa: an R/Bioconductor package for identifying regulatory motifs. Bioinformatics (2022).",br(),
                                   "- Lawrence M, Huber W, Pag\`es H, Aboyoun P, Carlson M, et al. (2013) Software for Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118. doi:10.1371/journal.pcbi.1003118",br(),
                                   "- Morgan M, Wang J, Obenchain V, Lang M, Thompson R, Turaga N (2022). _BiocParallel:
  Bioconductor facilities for parallel evaluation_. R package version 1.30.3,
  <https://github.com/Bioconductor/BiocParallel>.",br(),
                                   "- Morgan M, Obenchain V, Hester J, Pagès H (2022). _SummarizedExperiment:
  SummarizedExperiment container_. R package version 1.26.1,
  <https://bioconductor.org/packages/SummarizedExperiment>.",br(),
                                   "- Baranasic D (2020). _JASPAR2020: Data package for JASPAR database (version 2020)_. R
  package version 0.99.10, <http://jaspar.genereg.net/>.",br(),
                                   "- Team BC, Maintainer BP (2019). _TxDb.Mmusculus.UCSC.mm10.knownGene: Annotation package
  for TxDb object(s)_. R package version 3.10.0.",br(),
                                   "- Team TBD (2021). _BSgenome.Mmusculus.UCSC.mm10: Full genome sequences for Mus musculus
  (UCSC version mm10, based on GRCm38.p6)_. R package version 1.4.3.",br(),
                                   "- Carlson M, Maintainer BP (2015). _TxDb.Hsapiens.UCSC.hg19.knownGene: Annotation package
  for TxDb object(s)_. R package version 3.2.2.",br(),
                                   "- Team TBD (2020). _BSgenome.Hsapiens.UCSC.hg19: Full genome sequences for Homo sapiens
  (UCSC version hg19, based on GRCh37.p13)_. R package version 1.4.3.",br(),
                                   "Tan, G., and Lenhard, B. (2016). TFBSTools: an R/bioconductor package for transcription factor
  binding site analysis. Bioinformatics 32, 1555-1556.",br(),
                                   )
                          )
                 ),
                 tabPanel("Change log",
                          fluidRow(
                            column(12,
                                   h2("Log:"),
                                   h4("v1.0.5 (2023/4/24)"),
                                   strong("・Add 'download summary' buttons in the setting panel for 'Pair-wise DEG', '3 conditions DEG', and 'Multi conditions DEG'."),br(),
                                   strong("・Add new species (Xenopus laevis and Arabidopsis thaliana) for KEGG and GO analysis."),br(),
                                   strong("・Improve the 'start button' for motif analysis in Enrichment viewer."),br(),
                                   strong("・Improve the 'condition' color of the integrated heatmap in Venn diagram."),br(),
                                   strong("・Fix the issue of column name shifting in the output table data.(2023/5/10)"),br(),
                                   strong("・Display a warning message when inappropriate data is uploaded in Pair-wise DEG and 3 conditions DEG.(2023/5/11)"),br()
                            )
                          )
                 )
      )
    )
  ))
