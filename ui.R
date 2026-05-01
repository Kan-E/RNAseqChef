popoverTempate <- 
  '<div class="popover popover-lg" role="tooltip"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content"></div></div>'

hiddenFixedRadio <- function(inputId, selected) {
  tags$div(
    style = "display:none;",
    radioButtons(inputId, NULL, setNames(selected, selected), selected = selected, inline = TRUE)
  )
}

fixedSetting <- function(label, value) {
  tags$div(
    style = "margin-bottom: 10px;",
    tags$strong(paste0(label, ": ")),
    value
  )
}

shinyUI(
  fluidPage(
    tags$head(
      tags$style(HTML("
      .shiny-output-error-validation {
        color: #ff0000;
        font-weight: bold;
      }
    "))
    ),
    
    tags$head(
        includeScript("navAppend.js")
      ),
    tags$head(tags$style(HTML("
      #pair_pdf_exact_preview_image img,
      #shared_pdf_exact_preview_image img {
        max-width: 100%;
        max-height: 65vh;
        width: auto;
        height: auto;
        display: block;
        margin: 0 auto;
      }
    "))),
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
                      "or", a("contact us", href = "mailto:omicschef@kumamoto-u.ac.jp"),".",align="center",width=4)
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
                 column(12,
                        br(),
                        h4("Current version v1.1.7 (2026/5/1)"),
                        tags$ul(
                          tags$li("Added ssGSEA analysis to Multi DEG."),
                          tags$li("Added boxplots to compare gene expression by group in k-means clustering of Normalized count analysis."),
                          tags$li("Improved app performance with lighter and more optimized processing."),
                          tags$li("Fixed bugs in the Venn diagram workflow."),
                          tags$li("Added popup PDF preview for plot downloads across all modules."),
                          tags$li("Removed promoter motif analysis from the web version. Motif analysis is available only in the Docker version.")
                        ),
                        "See the details from 'More -> Change log'",
                        h4("Publication"),
                        "Etoh K. & Nakao M. A web-based integrative transcriptome analysis, RNAseqChef, uncovers cell/tissue type-dependent action of sulforaphane. JBC, 299(6), 104810 (2023)", 
                        a("https://doi.org/10.1016/j.jbc.2023.104810",href = "https://doi.org/10.1016/j.jbc.2023.104810"),
                        br(),
                        h4("Citation"),
                        tags$div(
                          style = paste(
                            "display:flex;",
                            "align-items:center;",
                            "justify-content:space-between;",
                            "gap:16px;",
                            "flex-wrap:wrap;",
                            "margin-top:8px;"
                          ),
                          tags$div(HTML('
                            <span class="__dimensions_badge_embed__" data-doi="10.1016/j.jbc.2023.104810" data-style="small_circle"></span>
                            <script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>
                          ')),
                          tags$a(
                            href = "https://x.com/kanetoh1",
                            target = "_blank",
                            rel = "noopener noreferrer",
                            style = paste(
                              "display:inline-flex;",
                              "align-items:center;",
                              "gap:10px;",
                              "margin-left:auto;",
                              "padding:10px 18px;",
                              "background:#000;",
                              "color:#fff;",
                              "border-radius:999px;",
                              "text-decoration:none;",
                              "font-weight:700;",
                              "font-size:13px;",
                              "line-height:1;",
                              "box-shadow:0 2px 8px rgba(0,0,0,0.18);"
                            ),
                            tags$span(HTML('
                              <svg width="20" height="20" viewBox="0 0 24 24" aria-hidden="true">
                                <path fill="#ffffff" d="M18.901 1.153h3.68l-8.04 9.188L24 22.847h-7.406l-5.8-7.584-6.639 7.584H.474l8.6-9.83L0 1.153h7.594l5.243 6.932 6.064-6.932zm-1.291 19.492h2.039L6.486 3.24H4.298z"/>
                              </svg>
                            ')),
                            tags$span(
                              "Follow on X for update announcement",
                              style = "white-space:nowrap;"
                            )
                          )
                        ),
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
                               "determines and visualizes biological functions of gene sets of interest",
                               br(),
                               tags$small("Promoter motif analysis is available only in the Docker version and has been removed from the web version."),
                               br(),
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
                                c('Count_matrix'="Row1",
                                  'Option: Count_matrix + Metadata'="Row2",
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
                   hiddenFixedRadio("Level_pair", "gene_level"),
                   radioButtons('DEG_method','DEG analysis method:',
                                c('DESeq2'="DESeq2",
                                  'EBSeq'="EBSeq",
                                  'edgeR'="edgeR"
                                ),selected = "DESeq2"),
                   conditionalPanel(condition=c("input.DEG_method=='limma' || input.DEG_method=='edgeR'"),
                                    fluidRow(
                                      column(6, radioButtons("pair_prefilterON","Pre-filtering (to remove low count)",
                                                             c('ON'="ON",
                                                               'OFF'="OFF"
                                                             ), selected = "ON")),
                                      
                                      column(6, conditionalPanel(condition=c("input.pair_prefilterON=='ON'"),
                                             numericInput("pair_prefilter", "Minimum count required for at least some samples", value = 10),
                                             numericInput("pair_prefilterTotal", "Minimum total count required", value = 15)
                                             )
                                      )
                                      )
                   ),
                   conditionalPanel(condition=c("input.DEG_method=='DESeq2' || input.DEG_method=='edgeR'"),
                                    selectInput("FDR_method", "FDR method", c("BH", "Qvalue", "IHW"), selected = "BH")
                   ),
                   conditionalPanel(condition="input.DEG_method=='limma'",
                                    fluidRow(
                                      column(6, radioButtons("limma_trend","Trend",
                                                             c('TRUE'=TRUE,
                                                               'FALSE'=FALSE
                                                             ), selected = TRUE)),
                                      column(6, radioButtons("regression_mode","Regression",
                                                             c('least squares'=FALSE,
                                                               'robust'=TRUE
                                                             ), selected = FALSE))
                                    ),
                                    radioButtons("cutoff_limma", "parameter for cut-off (fdr or pval)", c('fdr'="fdr",'pval'="pval"), selected = "fdr")
                   ),
                   fluidRow(
                     column(6, selectInput("Species", "Species", species_list, selected = "not selected")),
                     conditionalPanel(condition=c("input.Species != 'not selected' && input.Species != 'Homo sapiens' &&
                   input.Species != 'Mus musculus' && input.Species != 'Rattus norvegicus' &&
                   input.Species != 'Drosophila melanogaster' && input.Species != 'Caenorhabditis elegans' &&
                   input.Species != 'Bos taurus' && input.Species != 'Canis lupus familiaris' &&
                   input.Species != 'Danio rerio' && input.Species != 'Gallus gallus' &&
                   input.Species != 'Macaca mulatta' && input.Species != 'Pan troglodytes' &&
                   input.Species != 'Saccharomyces cerevisiae' && input.Species != 'Sus scrofa' &&
                   input.Species != 'Xenopus laevis' && input.Species != 'Arabidopsis thaliana'"),
                                      column(6, selectInput("Ortholog", 
                                                            strong(
                                                              span("Ortholog"),
                                                              span(icon("info-circle"), id = "Ortholog_pair", 
                                                                   options = list(template = popoverTempate))
                                                            ),
                                                            orgDb_list, selected = "Mus musculus"),
                                             bsPopover("Ortholog_pair", "Ortholog for the pathway analysis of non-model organisms", 
                                                       content=paste(img(src="non-model organism.png", width = 500,height = 800)), 
                                                       placement = "right",options = list(container = "body"))
                                             ),
                     )
                   ),
                   conditionalPanel(condition=c("input.Species != 'not selected' && input.Species != 'Homo sapiens' &&
                   input.Species != 'Mus musculus' && input.Species != 'Rattus norvegicus' &&
                   input.Species != 'Drosophila melanogaster' && input.Species != 'Caenorhabditis elegans' &&
                   input.Species != 'Bos taurus' && input.Species != 'Canis lupus familiaris' &&
                   input.Species != 'Danio rerio' && input.Species != 'Gallus gallus' &&
                   input.Species != 'Macaca mulatta' && input.Species != 'Pan troglodytes' &&
                   input.Species != 'Saccharomyces cerevisiae' && input.Species != 'Sus scrofa' &&
                   input.Species != 'Xenopus laevis' && input.Species != 'Arabidopsis thaliana' || input.Level_pair != 'gene_level'"),
                                    column(12, selectInput("Biomart_archive", "Biomart host", ensembl_archive))
                   ),
                   h4("Cut-off conditions:"),
                   fluidRow(
                     column(4, numericInput("fc", "Fold Change", min   = 1, max   = NA, value = 2)),
                     column(4, numericInput("fdr", "FDR", min   = 0, max   = 1, value = 0.05)),
                     column(4, numericInput("basemean", "Basemean", min   = 0, max   = NA, value = 0))
                   ),
                   conditionalPanel(condition="input.DEG_method!='limma'",
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
                             content=paste0("You can use a normalized count data (e.g. TPM count) for basemean cutoff and boxplot.<br>",
                                            strong("The column names of the normalized count data must match those of the uploaded raw count data.")), 
                             placement = "right",options = list(container = "body"))
                   ),
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
                                                              )))
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
                     type = "tabs",id = "pair_tabs",
                     tabPanel("Input Data", value = "pair_input_data_tab",
                              bsCollapse(id="input_collapse_panel",open="Row_count_panel",multiple = FALSE,
                                         bsCollapsePanel(title="Count_matrix:",
                                                         value="Row_count_panel",
                                                         dataTableOutput('Row_count_matrix')
                                         ),
                                         bsCollapsePanel(title="Metadata:",
                                                         value="Metadata_panel",
                                                         dataTableOutput('Metadata')
                                         ),
                                         bsCollapsePanel(title="Defined_count_matrix:",
                                                         value="D_row_count_matrix_panel",
                                                         conditionalPanel(condition="input.data_file_type!='Row11'",
                                                                          selectizeInput("sample_order", "Select samples:", choices = "", multiple = T),
                                                                          textOutput("not_cond2_pair"),
                                                                          tags$head(tags$style("#not_cond2_pair{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))
                                                         ),
                                                         tags$div(
                                                           style = "display:none;",
                                                           radioButtons("paired_sample", "Paired-sample?", c('No' = "No"), selected = "No")
                                                         ),
                                                         htmlOutput("paired_sample_file"),
                                                         dataTableOutput('paired_table'),
                                                         fluidRow(
                                                           column(4, downloadButton("download_pair_d_row_count", "Download defined raw count"))
                                                         ),
                                                         dataTableOutput('D_Row_count_matrix')
                                         )
                              )
                     ),
                     tabPanel("Result overview", value = "pair_result_overview_tab",
                              fluidRow(
                                column(6, actionButton("preview_pair_pca", "Download clustering analysis"),
                                       textOutput("not_cond2"),
                                       tags$head(tags$style("#not_cond2{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                                column(6, selectInput("PCA_legend","Label",c("Label","Legend"),selected = "Label"))
                              ),
                              plotOutput("PCA"),
                              fluidRow(
                                column(4, actionButton("preview_pair_ma", "Download MA-plot"))
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
                     tabPanel("GOI profiling", value = "pair_goi_tab",
                              fluidRow(
                                column(4, actionButton("preview_pair_volcano", "Download volcano plot / MA plot")),
                                column(4, actionButton("preview_pair_goi_heat", "Download heatmap"))
                              ),
                              fluidRow(
                                column(4, selectInput("GOI_plot_select","Plot type",c("Volcano plot","MA plot"),
                                                      selected = "Volcano plot",multiple = F),
                                       htmlOutput("GOI_color_type")),
                                column(4, htmlOutput("volcano_x"),
                                       actionButton("GOIreset_pair", "GOI reset")),
                                column(4, htmlOutput("volcano_y"))
                              ),
                              conditionalPanel(condition="input.GOI_color_type=='pathway'",
                                               column(10,selectInput("GOI_color_pathway1","Select a gene set",choices = ""),
                                                      selectInput("GOI_color_pathway2","",choices = ""))             
                              ),
                              selectizeInput("GOI", "genes of interest (GOI)", choices = NULL, multiple = TRUE,
                                             options = list(delimiter = " ", create = TRUE, plugins = list("remove_button"), persist = FALSE)),
                              htmlOutput("uniqueID_cut"),
                              fluidRow(
                                column(8, 
                                       plotOutput("volcano1",
                                                  brush = "plot1_brush"
                                                  )
                                       ),
                                column(4, plotOutput("GOIheatmap"))
                              ),
                              htmlOutput("paired_for_GOItype"),
                              div(
                                plotOutput("GOIbox", height = "100%"),
                                style = "height: calc(100vh  - 100px)"
                              ),
                              fluidRow(
                                column(4, actionButton("preview_pair_goi_box", "Download boxplot"))
                              )),
                     tabPanel("Enrichment analysis", value = "pair_enrichment_tab",
                              fluidRow(
                                column(4, textOutput("Spe1"),
                                       tags$head(tags$style("#Spe1{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                                column(4, selectInput("Gene_set", "Gene Set", gene_set_list)),
                                column(4, htmlOutput("Custom_input_pair")),
                                column(4, actionButton("preview_pair_enrichment", "Download"))
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
                     ),
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
                                  'Option: Raw_count_matrix + Metadata'="Row4",
                                  'Option: Recode.Rdata'="RowRecode_cond3"
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
                   conditionalPanel(condition="input.data_file_type2=='RowRecode_cond3'",
                                    fileInput("file_recode_cond3",
                                              strong(
                                                span("Select a Recode.Rdata"),
                                                span(icon("info-circle"), id = "icon_recode", 
                                                     options = list(template = popoverTempate))
                                              ),
                                              accept = c("Rdata"),
                                              multiple = FALSE,
                                              width = "80%"),
                                    bsPopover("icon_recode", "Recode.Rdata:", 
                                              content=paste("Using this mode, you can save time to EBSeq analysis.<br>",
                                                            "The recode.Rdata file from the previous analysis using '3 conditions DEG' is needed.<br>", 
                                                            "You can obtain a recode.Rdata file by clicking 'Download summary' buttom after the analysis with 3 conditions DEG"
                                                            ), 
                                              placement = "right",options = list(container = "body")),
                   ),
                   radioButtons("EBSeq_mode","EBSeq",
                                c('v2 (recommend)'=TRUE,
                                  'v1 (<= RNAseqChef v.1.1.0)'=FALSE),
                                selected = TRUE),
                   hiddenFixedRadio("Level_cond3", "gene_level"),
                   fluidRow(
                     column(6, selectInput("Species2", "Species", species_list, selected = "not selected")),
                     conditionalPanel(condition=c("input.Species2 != 'not selected' && input.Species2 != 'Homo sapiens' &&
                   input.Species2 != 'Mus musculus' && input.Species2 != 'Rattus norvegicus' &&
                   input.Species2 != 'Drosophila melanogaster' && input.Species2 != 'Caenorhabditis elegans' &&
                   input.Species2 != 'Bos taurus' && input.Species2 != 'Canis lupus familiaris' &&
                   input.Species2 != 'Danio rerio' && input.Species2 != 'Gallus gallus' &&
                   input.Species2 != 'Macaca mulatta' && input.Species2 != 'Pan troglodytes' &&
                   input.Species2 != 'Saccharomyces cerevisiae' && input.Species2 != 'Sus scrofa' &&
                   input.Species2 != 'Xenopus laevis' && input.Species2 != 'Arabidopsis thaliana'"),
                                      column(6, selectInput("Ortholog2", strong(
                                        span("Ortholog"),
                                        span(icon("info-circle"), id = "Ortholog_cond3", 
                                             options = list(template = popoverTempate))
                                      ), orgDb_list, selected = "Mus musculus"),
                                      bsPopover("Ortholog_cond3", "Ortholog for the pathway analysis of non-model organisms", 
                                                content=paste(img(src="non-model organism.png", width = 500,height = 800)), 
                                                placement = "right",options = list(container = "body"))
                                      ),
                                      column(12, selectInput("Biomart_archive2", "Biomart host", ensembl_archive)))
                     ),
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
                             content=paste0("You can use a normalized count data (e.g. TPM count) for basemean cutoff and boxplot.<br>",
                                            strong("The column names of the normalized count data must match those of the uploaded raw count data.")), 
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
                                                         conditionalPanel(condition="input.data_file_type2 != 'RowRecode_cond3'",
                                                                          selectizeInput("sample_order_cond3", "Select samples:", choices = "", multiple = T),
                                                                          textOutput("not_cond3_select"),
                                                                          tags$head(tags$style("#not_cond3_select{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))
                                                         ),
                                                         fluidRow(
                                                           column(4, downloadButton("download_cond3_d_row_count", "Download defined raw count"))
                                                         ),
                                                         dataTableOutput('D_Row_count_matrix2')
                                         )
                              )
                     ),
                     tabPanel("Result overview",
                              fluidRow(
                                column(6, actionButton("preview_cond3_pca", "Download clustering analysis"),
                                       textOutput("not_cond3"),
                                       tags$head(tags$style("#not_cond3{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                                column(6, selectInput("PCA_legend_cond3","Label",c("Label","Legend"),selected = "Label"))
                              ),
                              plotOutput("PCA2"),
                              fluidRow(
                                column(4, actionButton("preview_cond3_overview", "Download DEG overview"))
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
                                column(4, actionButton("preview_cond3_goi_scatter", "Download scatter plot")),
                                column(4, actionButton("preview_cond3_goi_heat", "Download heatmap"))
                              ),
                              fluidRow(
                                column(4,
                                       htmlOutput("cond3_GOI_color_type")),
                                column(4,
                                       htmlOutput("cond3_xrange"),
                                       actionButton("GOIreset_cond3", "GOI reset")),
                                column(4, htmlOutput("cond3_yrange"))
                              ),
                              conditionalPanel(condition="input.cond3_GOI_color_type=='pathway'",
                                               column(6,selectInput("cond3_GOI_color_pathway1","Select a gene set",choices = ""),
                                                      selectInput("cond3_GOI_color_pathway2","",choices = ""))             
                              ),
                              selectizeInput("GOI2", "genes of interest (GOI)", choices = NULL, multiple = TRUE,
                                             options = list(delimiter = " ", create = TRUE, plugins = list("remove_button"), persist = FALSE)),
                              htmlOutput("cond3_GOI_pair"),
                              htmlOutput("cond3_uniqueID_cut"),
                              fluidRow(
                                column(8, plotOutput("cond3_GOIscatter",
                                                     brush = "plot1_brush_cond3")),
                                column(4, plotOutput("cond3_GOIheatmap"))
                              ),
                              div(
                                plotOutput("cond3_GOIboxplot", height = "100%"),
                                style = "height: calc(100vh  - 100px)"
                              ),
                              fluidRow(
                                column(4, actionButton("preview_cond3_goi_box", "Download boxplot"))
                              )
                     ),
                     tabPanel("Enrichment analysis",
                              fluidRow(
                                column(4, textOutput("Spe2"),
                                       tags$head(tags$style("#Spe2{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                                column(4, selectInput("Gene_set2", "Gene Set", gene_set_list)),
                                column(4, actionButton("preview_cond3_enrich", "Download"))
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
                                c('Count_matrix (One-factor multi-condition)'="Row1",
                                  'Count_matrix + metadata (Two-factor multi-condition)'="Row2"
                                ),selected = "Row1"),
                   conditionalPanel(condition="input.multi_data_file_type=='Row1'",
                                    fileInput("multi_file1",
                                              strong(
                                                span("Select a count matrix file"),
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
                                                span("Select a count matrix file"),
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
                   hiddenFixedRadio("Level_multi", "gene_level"),
                   fixedSetting("DEG analysis method", "DESeq2 for raw count"),
                   hiddenFixedRadio("DEG_method_multi", "DESeq2"),
                   conditionalPanel(condition="input.DEG_method_multi=='limma'",
                                    fluidRow(
                                      column(6, radioButtons("multi_prefilterON","Pre-filtering (to remove low count)",
                                                             c('ON'="ON",
                                                               'OFF'="OFF"
                                                             ), selected = "ON")),
                                      column(6, conditionalPanel(condition=c("input.multi_prefilterON=='ON'"),
                                             numericInput("multi_prefilter", "Minimum count required for at least some samples", value = 0),
                                             numericInput("multi_prefilterTotal", "Minimum total count required", value = 0))
                                      )
                                    ),
                                    fluidRow(
                                      column(6, radioButtons("limma_trend_multi","Trend",
                                                             c('TRUE'=TRUE,
                                                               'FALSE'=FALSE
                                                             ), selected = TRUE)),
                                      column(6, radioButtons("regression_mode_multi","Regression",
                                                             c('least squares'=FALSE,
                                                               'robust'=TRUE
                                                             ), selected = FALSE))
                                    )
                   ),
                   conditionalPanel(condition="input.DEG_method_multi=='DESeq2'",
                                    fluidRow(
                                      column(6,  selectInput("FDR_method6", "FDR method", c("BH", "Qvalue", "IHW"), selected = "BH"))
                                    )
                   ),
                   fluidRow(
                            column(6, selectInput("Species6", "Species", species_list, selected = "not selected")),
                            conditionalPanel(condition=c("input.Species6 != 'not selected' && input.Species6 != 'Homo sapiens' &&
                   input.Species6 != 'Mus musculus' && input.Species6 != 'Rattus norvegicus' &&
                   input.Species6 != 'Drosophila melanogaster' && input.Species6 != 'Caenorhabditis elegans' &&
                   input.Species6 != 'Bos taurus' && input.Species6 != 'Canis lupus familiaris' &&
                   input.Species6 != 'Danio rerio' && input.Species6 != 'Gallus gallus' &&
                   input.Species6 != 'Macaca mulatta' && input.Species6 != 'Pan troglodytes' &&
                   input.Species6 != 'Saccharomyces cerevisiae' && input.Species6 != 'Sus scrofa' &&
                   input.Species6 != 'Xenopus laevis' && input.Species6 != 'Arabidopsis thaliana'"),
                                             column(6, selectInput("Ortholog6", strong(
                                               span("Ortholog"),
                                               span(icon("info-circle"), id = "Ortholog_multi", 
                                                    options = list(template = popoverTempate))
                                             ), orgDb_list, selected = "Mus musculus"),
                                             bsPopover("Ortholog_multi", "Ortholog for the pathway analysis of non-model organisms", 
                                                       content=paste(img(src="non-model organism.png", width = 500,height = 800)), 
                                                       placement = "right",options = list(container = "body"))),
                                             column(12, selectInput("Biomart_archive6", "Biomart host", ensembl_archive)))
                            ),
                   h4("Cut-off conditions:"),
                   fluidRow(
                     column(4, numericInput("fc6", "Fold Change", min   = 0, max   = NA, value = 1.5)),
                     column(4, numericInput("fdr6", "FDR", min   = 0, max   = 1, value = 0.05)),
                     column(4, numericInput("basemean6", "Basemean", min   = 0, max   = NA, value = 0))
                   ),
                   conditionalPanel(condition="input.DEG_method_multi!='limma'",
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
                             placement = "right",options = list(container = "body"))
                   ),
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
                    id = "multi_main_tab",
                    type = "tabs",
                    tabPanel("Input Data", value = "multi_input_data_tab",
                              bsCollapse(id="multi_input_collapse_panel",open="multi_Row_count_panel",multiple = TRUE,
                                         bsCollapsePanel(title="Count_matrix:",
                                                         value="multi_Row_count_panel",
                                                         dataTableOutput('multi_Row_count_matrix')
                                         ),
                                         bsCollapsePanel(title="Metadata:",
                                                         value="multi_Metadata_panel",
                                                         dataTableOutput('multi_Metadata')
                                         ),
                                         bsCollapsePanel(title="Defined_count_matrix:",
                                                         value="multi_d_Row_count_panel",
                                                         dataTableOutput('multi_d_Row_count_matrix')
                                         )
                              )
                     ),
                    tabPanel("Result overview", value = "multi_result_overview_tab",
                              fluidRow(
                                column(4, actionButton("preview_multi_pca", "Download clustering analysis")),
                                column(6, selectInput("PCA_legend_multi","Label",c("Label","Legend"),selected = "Label"))
                              ),
                              plotOutput("multi_PCA"),
                              fluidRow(
                                column(6, htmlOutput("multi_umap_n"),
                                       actionButton("preview_multi_umap", "Download umap"),
                                       textOutput("multi_umap_error"),
                                       tags$head(tags$style("#multi_umap_error{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                                column(6, selectInput("multi_umap_label","Label",c("Label","None"),selected = "Label"))
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
                    tabPanel("Divisive clustering", value = "multi_divisive_clustering_tab",
                              fluidRow(
                                column(6, htmlOutput("selectFC")),
                                column(6, htmlOutput("selectFC_2")),
                              ),
                              textOutput("multi_DEG_total1"),
                              tags$head(tags$style("#multi_DEG_total1{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }")),
                              htmlOutput("topP"),
                              fluidRow(
                                column(3, actionButton("start_multi", "Start"),
                                       tags$head(tags$style("#start_multi{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"),
                                                 tags$style("
          body {
            padding: 0 !important;
          }"
                                                 ))),
                              ),
                              fluidRow(
                                column(4, actionButton("preview_multi_div_box", "Download boxplots"))
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
                                                column(4, selectInput("multi_selectfile1", "cluster_list",
                                                                      choices = character(0),
                                                                      selected = NULL,
                                                                      multiple = FALSE)),
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
                                                column(4, actionButton("preview_multi_div_goi_box", "Download boxplot"))
                                              ),
                                              plotOutput("multi_pattern_boxplot")
                              ),
                              bsCollapsePanel(title="Enrichment analysis:",
                                              value="multi_deg_pattern_enrichment_panel",
                                              fluidRow(
                                                column(4, htmlOutput("multi_whichGroup1_1")),
                                                column(4, selectInput("Gene_set7", "Gene Set", gene_set_list)),
                                                column(4, actionButton("preview_multi_div_enrich", "Download dot plot"))
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
                                                column(4, actionButton("preview_multi_div_cnet", "Download cnet plot"))
                                              ),
                                              plotOutput("multi_enrichment4"),
                                              fluidRow(
                                                column(4, downloadButton("download_multi_enrichment_table", "Download enrichment result"))
                                              ),
                                              dataTableOutput('multi_enrichment_result')
                              )
                              )
                     ),
                    tabPanel("k-means clustering", value = "multi_kmeans_clustering_tab",
                              fluidRow(
                                column(6, htmlOutput("selectFC2")),
                                column(6, htmlOutput("selectFC2_2"))
                              ),
                              textOutput("multi_DEG_total2"),
                              tags$head(tags$style("#multi_DEG_total2{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }")),
                              fluidRow(
                                column(6, htmlOutput("topP2"))
                              ),
                              fluidRow(
                                column(4, htmlOutput("multi_kmeans_num"),
                                       actionButton("kmeans_start_multi", "Start"),
                                       textOutput("multi_kmeans_progress_text"),
                                       tags$head(tags$style("#kmeans_start_multi{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }
                                 #multi_kmeans_progress_text{
                                 color: #d9534f;
                                 font-size: 14px;
                                 font-weight: bold;
                                 white-space: pre-line;
                                 margin-top: 8px;
                                 }"),
                                                 tags$style("
          body {
            padding: 0 !important;
          }"
                                                 )),
                                       actionButton("preview_multi_kmeans_heat", "Download heatmap"),
                                       actionButton("preview_multi_kmeans_box", "Download boxplots")),
                                column(8, htmlOutput("kmeans_order_multi"),
                                       plotOutput("multi_kmeans_heatmap"))
                              ),
                              div(
                                plotOutput("multi_kmeans_boxplot", height = "100%"),
                                style = "height: calc(100vh  - 100px)"
                              ),
                              bsCollapse(id="multi_collapse_panel2",open="multi_deg_kmeans_pattern_panel",multiple = TRUE,
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
                                                column(4, actionButton("preview_multi_kmeans_goi_box", "Download boxplot"))
                                              ),
                                              plotOutput("multi_kmeans_GOIboxplot")
                              ),
                              bsCollapsePanel(title="Enrichment analysis:",
                                              value="multi_deg_kmeans_pattern_enrichment_panel",
                                              fluidRow(
                                                column(4, htmlOutput("multi_whichGroup2_1")),
                                                column(4, selectInput("Gene_set8", "Gene Set", gene_set_list)),
                                                column(4, actionButton("preview_multi_kmeans_enrich", "Download dot plot"))
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
                                                column(4, actionButton("preview_multi_kmeans_cnet", "Download cnet plot"))
                                              ),
                                              plotOutput("multi_enrichment6"),
                                              fluidRow(
                                                column(4, downloadButton("download_multi_enrichment_table2", "Download enrichment result"))
                                              ),
                                              dataTableOutput('multi_enrichment_result2')
                              )
                              )
                     ),
                    tabPanel("GSEA", value = "multi_gsea_tab",
                              fluidRow(
                                column(3, selectizeInput("selectEnrich_pair", "Select a pair for GSEA", choices = NULL,
                                                         multiple = TRUE, options = list(maxItems = 2))),
                                column(4, selectInput("Gene_set6", "Gene Set", choices = "")),
                                column(4, actionButton("preview_multi_gsea", "Download"))
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
                     ),
                    tabPanel("ssGSEA", value = "multi_ssgsea_tab",
                              fluidRow(
                                column(4, selectInput("Gene_set_ssGSEA", "Gene Set", choices = "", selected = "")),
                              column(4, htmlOutput("Custom_input_ssGSEA"))
                              ),
                              fluidRow(
                                column(4, textOutput("multi_Spe1_ssGSEA"),
                                       tags$head(tags$style("#multi_Spe1_ssGSEA{color: red;
                                         font-size: 20px;
                                         font-style: bold;
                                      }")))
                              ),
                              bsCollapse(id="input_collapse_multi_ssGSEA",open="ssGSEA_score_panel",multiple = TRUE,
                                         bsCollapsePanel(title="ssGSEA score:",
                                                         value="ssGSEA_score_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_multi_ssGSEA_score", "Download"))
                                                         ),
                                                         dataTableOutput('multi_ssGSEA_score')
                                         ),
                                         bsCollapsePanel(title="Differential analysis:",
                                                         value="ssGSEA_limma_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_multi_ssGSEA_limma", "Download"))
                                                         ),
                                                         dataTableOutput('multi_ssGSEA_limma')
                                         ),
                                         bsCollapsePanel(title="Differential pathways:",
                                                         value="ssGSEA_limma_panel",
                                                         fluidRow(
                                                           column(4, downloadButton("download_multi_ssGSEA_limma_dp", "Download"))
                                                         ),
                                                         dataTableOutput('multi_ssGSEA_limma_dp')
                                         )
                              ),
                              fluidRow(
                                column(4, actionButton("preview_multi_ssgsea_heat", "Download heatmap"))
                              ),
                              fluidRow(
                                column(4, 
                                       numericInput("ssGSEA_fdr","FDR cut-off value",value = 0.01,max = 1,min = 0,step = 0.0001),
                                       radioButtons('GOI_type_multi_ssGSEA','Pathways:',
                                                    c('Select all differential pathways'="ALL",
                                                      'Custom'="custom"
                                                    ),selected = "custom"),
                                       htmlOutput("GOI_type_multi_ssGSEA_all"),
                                       htmlOutput("GOI_type_multi_ssGSEA_custom"),
                                       conditionalPanel(
                                         condition = "input.GOI_type_multi_ssGSEA=='ALL' || (input.GOI_type_multi_ssGSEA=='custom' && input.GOI_type_multi_ssGSEA_custom=='Manual')",
                                         selectizeInput("GOI_multi_ssGSEA", "Pathways of interest", choices = NULL, multiple = TRUE,
                                                        options = list(delimiter = " ", create = TRUE, plugins = list("remove_button"), persist = FALSE)),
                                         actionButton("GOIreset_multi_ssGSEA", "reset")
                                       )),
                                column(8, plotOutput("multi_ssGSEA_GOIheatmap"))
                              ),
                              fluidRow(
                                column(4, htmlOutput("statistics_multi_ssGSEA")),
                                column(4, selectInput('PlotType_multi_ssGSEA', 'PlotType', c("Boxplot", "Barplot", "Violin plot","Errorplot"))),
                                column(4, actionButton("preview_multi_ssgsea_box", "Download boxplot"))
                              ),
                              div(
                                plotOutput("multi_ssGSEA_GOIboxplot", height = "100%"),
                                style = "height: calc(100vh  - 100px)"
                              ),
                              fluidRow(
                                column(4, downloadButton("download_multi_ssGSEA_statisics", "Download table"))
                              ),
                              dataTableOutput("statistical_table_multi_ssGSEA"),
                              bsCollapse(id="input_collapse_multi_ssGSEA_dorothea",multiple = TRUE,
                                         bsCollapsePanel(title = "Correlation between ssGSEA score and expression of pathway genes",
                                                         value="ssGSEA_contribute_panel",
                                                         selectizeInput("selectssGSEA_contribute_pathway", "Select pathway of interest",
                                                                        choices = NULL, multiple = FALSE,
                                                                        options = list(create = FALSE, persist = FALSE)),
                                                         fluidRow(
                                                           column(4, actionButton("preview_multi_ssgsea_contrib", "Download"))
                                                         ),
                                                         htmlOutput("multi_ssGSEA_contribute_note"),
                                                         div(
                                                           plotOutput("multi_ssGSEA_contribute", height = "100%"),
                                                           style = "height: calc(100vh  - 100px)"
                                                         ),
                                                         downloadButton("download_multi_ssGSEA_contribute_table", "Download"),
                                                         dataTableOutput('multi_ssGSEA_contribute_table')
                                         ),
                                         bsCollapsePanel(title = "Correlation between ssGSEA score and TF expression",
                                                         value="ssGSEA_dorothea_panel",
                                                         fluidRow(
                                                           column(4, actionButton("preview_multi_ssgsea_dorothea", "Download"))
                                                         ),
                                                         plotOutput('multi_ssGSEA_dorothea'),
                                                         downloadButton("download_multi_ssGSEA_dorothea_table", "Download"),
                                                         dataTableOutput('multi_ssGSEA_dorothea_table')
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
                             img(src="venn_input.png", width = 400,height = 300)),
                             placement = "right",options = list(container = "body")),
                   hiddenFixedRadio("Level_venn", "gene_level"),
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
                     conditionalPanel(condition=c("input.Species7 != 'not selected' && input.Species6 != 'Homo sapiens' &&
                   input.Species7 != 'Mus musculus' && input.Species7 != 'Rattus norvegicus' &&
                   input.Species7 != 'Drosophila melanogaster' && input.Species7 != 'Caenorhabditis elegans' &&
                   input.Species7 != 'Bos taurus' && input.Species7 != 'Canis lupus familiaris' &&
                   input.Species7 != 'Danio rerio' && input.Species7 != 'Gallus gallus' &&
                   input.Species7 != 'Macaca mulatta' && input.Species7 != 'Pan troglodytes' &&
                   input.Species7 != 'Saccharomyces cerevisiae' && input.Species7 != 'Sus scrofa' &&
                   input.Species7 != 'Xenopus laevis' && input.Species7 != 'Arabidopsis thaliana'"),
                     column(6, selectInput("Ortholog7", strong(
                       span("Ortholog"),
                       span(icon("info-circle"), id = "Ortholog_venn", 
                            options = list(template = popoverTempate))
                     ), orgDb_list, selected = "Mus musculus"),
                     bsPopover("Ortholog_venn", "Ortholog for the pathway analysis of non-model organisms", 
                               content=paste(img(src="non-model organism.png", width = 500,height = 800)), 
                               placement = "right",options = list(container = "body"))),
                     column(12, selectInput("Biomart_archive7", "Biomart host", ensembl_archive)))
                   ),
                   fluidRow(
                   column(6, selectInput(
                     inputId = "pre_zscoring",
                     strong(
                       span("Option: Pre-zscoring"),
                       span(icon("info-circle"), id = "icon_venn3", 
                            options = list(template = popoverTempate))
                     ),
                     multiple = FALSE,choices = c("TRUE", "FALSE"), selected = "TRUE"))),
                   bsPopover("icon_venn3", "Option: Pre-zscoring:", 
                             content=paste("If True, each count data is z-scored before integrating multiple count data.<br><br>",
                                           img(src="pre-zscoring.png", width = 450,height = 200)), 
                             placement = "right",options = list(container = "body")),
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
                                column(4, actionButton("preview_venn_diagram", "Download venn diagram"))
                              ),
                              fluidRow(
                                column(4,radioButtons("venn_type","venn_type",c("default"="default","eulerr package"="eulerr"),"eulerr")),
                                column(4,htmlOutput("eulerr_label"))
                              ),
                              plotOutput("venn"),
                              downloadButton("download_venn_result", "Download venn result"),
                              dataTableOutput("venn_result")
                              ),
                     tabPanel("integrated_heatmap",
                              fluidRow(
                                column(4, htmlOutput("select_file2")),
                                column(4, actionButton("preview_venn_heat", "download integrated heatmap"))
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
                                                         actionButton("preview_venn_box", "Download boxplot"),
                                                         plotOutput("GOIbox_venn")
                                         )
                              ),
                     ),
                     tabPanel("Enrichment analysis",
                              fluidRow(
                                column(4, htmlOutput("venn_whichGroup1")),
                                column(4, selectInput("Gene_set9", "Gene Set", gene_set_list)),
                                column(4, actionButton("preview_venn_enrich", "Download dot plot"))
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
                                column(4, actionButton("preview_venn_cnet", "Download cnet plot"))
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
                   hiddenFixedRadio("Level_norm", "gene_level"),
                   fluidRow(
                     column(6, selectInput("Species3", "Species", species_list, selected = "not selected")),
                     conditionalPanel(condition=c("input.Species3 != 'not selected' && input.Species3 != 'Homo sapiens' &&
                   input.Species3 != 'Mus musculus' && input.Species3 != 'Rattus norvegicus' &&
                   input.Species3 != 'Drosophila melanogaster' && input.Species3 != 'Caenorhabditis elegans' &&
                   input.Species3 != 'Bos taurus' && input.Species3 != 'Canis lupus familiaris' &&
                   input.Species3 != 'Danio rerio' && input.Species3 != 'Gallus gallus' &&
                   input.Species3 != 'Macaca mulatta' && input.Species3 != 'Pan troglodytes' &&
                   input.Species3 != 'Saccharomyces cerevisiae' && input.Species3 != 'Sus scrofa' &&
                   input.Species3 != 'Xenopus laevis' && input.Species3 != 'Arabidopsis thaliana'"),
                     column(6, selectInput("Ortholog3", strong(
                       span("Ortholog"),
                       span(icon("info-circle"), id = "Ortholog_norm", 
                            options = list(template = popoverTempate))
                     ), orgDb_list, selected = "Mus musculus"),
                     bsPopover("Ortholog_norm", "Ortholog for the pathway analysis of non-model organisms", 
                               content=paste(img(src="non-model organism.png", width = 500,height = 800)), 
                               placement = "right",options = list(container = "body"))),
                     column(12, selectInput("Biomart_archive3", "Biomart host", ensembl_archive)))
                     ),
                   fluidRow(
                     column(6, h4("Filter option 1:"),
                            radioButtons("norm_filter","Gene set filter",
                                c("not selected"="not selected","Gene set"="Gene set","cusom"="custom"))),
                   conditionalPanel(condition="input.norm_filter=='Gene set'",
                                    column(6,selectInput("gene_set_forFilter","Select a gene set for gene extraction",choices = ""),
                                           selectInput("gene_set_forFilter2","",choices = ""))
                                    
                                    
                                    
                   ),
                   conditionalPanel(condition="input.norm_filter=='custom'",
                                    column(6,fileInput("file10",
                                              label = "Select a gene list file for gene extraction",
                                              accept = c("txt", "csv", "xlsx"),
                                              multiple = FALSE,
                                              width = "80%"))
                   )
                   ),
                   h4("Filter option 2:"),
                   fluidRow(
                     column(4, numericInput("fc3", "Fold Change", min   = 1, max   = NA, value = 1)),
                     column(4, numericInput("basemean3", "Basemean", min   = 0, max   = NA, value = 0),
                     )
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
                                                         selectizeInput("sample_order_norm", "Select samples:", choices = "", multiple = T),
                                                         fluidRow(
                                                           column(4, downloadButton("download_d_norm_count", "Download defined normalized count"))
                                                         ),
                                                         dataTableOutput('d_norm_count')
                                         )
                              )
                     ),
                     tabPanel("Clustering",
                              fluidRow(
                                column(3, actionButton("preview_norm_pca", "Download Clustering")),
                                column(6, selectInput("PCA_legend_norm","Label",c("Label","Legend"),selected = "Label"))
                              ),
                              plotOutput("norm_PCA"),
                              fluidRow(
                                column(6, htmlOutput("norm_umap_n"),
                                       actionButton("preview_norm_umap", "Download umap"),
                                       textOutput("norm_umap_error"),
                                       tags$head(tags$style("#norm_umap_error{color: red;
                                 font-size: 20px;
            font-style: bold;
            }"))),
                                column(6, selectInput("norm_umap_label","Label",c("Label","None"),selected = "Label"))
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
                                column(4, actionButton("preview_norm_goi_heat", "Download heatmap"))
                              ),
                              fluidRow(
                                column(4, 
                                       radioButtons("normGOI_filter_on","Filter (basemean, fold change)",c("ON" = "ON","OFF"="OFF"),selected = "OFF"),
                                       conditionalPanel(condition="input.normGOI_filter_on=='ON'",
                                       htmlOutput("selectFC_normGOI"),
                                       textOutput("filtered_regionGOI"),
                                       tags$head(tags$style("#filtered_regionGOI{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"))
                                       ),
                                       radioButtons('GOI3_type','Genes:',
                                                    c('Select all genes'="ALL",
                                                      'Custom'="custom"
                                                    ),selected = "custom"),
                                       selectizeInput("GOI3", "genes of interest (GOI)", choices = NULL, multiple = TRUE,
                                                      options = list(delimiter = " ", create = TRUE, plugins = list("remove_button"), persist = FALSE)),
                                       actionButton("GOIreset_norm", "GOI reset"),
                                       htmlOutput("norm_uniqueID_cut")),
                                column(8, plotOutput("norm_GOIheatmap"))
                              ),
                              fluidRow(
                                column(3, radioButtons('PlotType', 'PlotType', c("Boxplot"="Boxplot", 
                                                                                 "Barplot"="Barplot",
                                                                                 "Violin plot"="Violin plot",
                                                                                 "Errorplot"="Errorplot"),selected = "Boxplot")),
                                column(3, htmlOutput("Color_design_norm")),
                                column(3, htmlOutput("Color_norm")),
                                column(3, htmlOutput("Color_rev_norm"))
                              ),
                              fluidRow(
                                column(3, numericInput("norm_ymin","min value of y-axis",value = 0)),
                                column(3, textInput("norm_ylab","title of y-axis","Normalized count"))
                              ),
                              htmlOutput("statistics"),
                              actionButton("preview_norm_goi_box", "Download boxplot"),
                              div(
                                plotOutput("norm_GOIboxplot", height = "100%"),
                                style = "height: calc(100vh  - 100px)"
                              ),
                              downloadButton("download_statisics", "Download table"),
                              dataTableOutput("statistical_table")
                     ),
                     tabPanel("Correlation analysis",
                              fluidRow(
                                column(4, actionButton("preview_norm_corr", "Download correlation plot"))
                              ),
                              fluidRow(
                                column(4,
                                       radioButtons('corr_mode','Mode:',
                                                    c('Screening'="corr_mode1",
                                                      'Selected pair'="corr_mode2"
                                                    ),selected = "corr_mode2"),
                                       selectInput("corr_statistics","Statistics",c("spearman","pearson"),
                                                   selected = "spearman",multiple = F),
                                       selectizeInput("GOI_x", "Select GOI", choices = NULL, multiple = FALSE,
                                                      options = list(create = FALSE, persist = FALSE)),
                                       conditionalPanel(
                                         condition="input.corr_mode=='corr_mode2'",
                                         selectizeInput("GOI_y", "Select GOI (y_axis)", choices = NULL, multiple = FALSE,
                                                        options = list(create = FALSE, persist = FALSE)),
                                         selectizeInput("corr_color","Color", choices = NULL, multiple = FALSE,
                                                        options = list(create = FALSE, persist = FALSE))
                                       ),
                                       conditionalPanel(condition="input.corr_mode=='corr_mode1'",
                                       actionButton("corr_start", "Start"),
                                       tags$head(tags$style("#corr_start{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"),
                                                 tags$style("
          body {
            padding: 0 !important;
          }"
                                                 ))
                                       ),
                                       htmlOutput("norm_corr_cutoff"),
                                       htmlOutput("corr_fdr"),
                                       htmlOutput("corr_rho")
                              ),
                              column(8, plotOutput("norm_corrplot",brush = "plot1_brush_corr"))
                              ),
                              column(4, downloadButton("download_statisics_corrplot", "Download table")),
                              dataTableOutput("statistical_table_corrplot"),
                              fluidRow(
                              column(4, htmlOutput("norm_corr_selected_list"),
                                     conditionalPanel(
                                       condition="input.corr_mode=='corr_mode1'",
                                       selectizeInput("corr_color_selected","Color", choices = NULL, multiple = FALSE,
                                                      options = list(create = FALSE, persist = FALSE))
                                     ),
                                     actionButton("preview_norm_corr_selected", "Download correlation plot (all genes from 'select GOI')"))
                              ),
                              plotOutput("norm_corrplot_selected")
                     ),
                     tabPanel("k-means clustering",
                              fluidRow(
                                column(4, htmlOutput("selectFC_norm")),
                                column(4, htmlOutput("selectFC_norm2"))
                                ),
                              textOutput("filtered_region"),
                              tags$head(tags$style("#filtered_region{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }")),
                              fluidRow(
                                column(4, htmlOutput("norm_kmeans_num"),
                                       htmlOutput("kmeans_cv"),
                                       actionButton("kmeans_start", "Start"),
                                       textOutput("norm_kmeans_progress_text"),
                                       tags$head(tags$style("#kmeans_start{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }
                                 #norm_kmeans_progress_text{
                                 color: #d9534f;
                                 font-size: 14px;
                                 font-weight: bold;
                                 white-space: pre-line;
                                 margin-top: 8px;
                                 }"),
                                                 tags$style("
          body {
            padding: 0 !important;
          }"
                                                 ))),
                                column(8, actionButton("preview_norm_kmeans_heat", "Download heatmap"),
                                       actionButton("preview_norm_kmeans_box", "Download boxplots"),
                                       htmlOutput("kmeans_order"),
                                       plotOutput("norm_kmeans_heatmap")),
                              ),
                              div(
                                plotOutput("norm_kmeans_boxplot", height = "100%"),
                                style = "height: calc(100vh  - 100px)"
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
                                              actionButton("preview_norm_kmeans_goi_box", "Download boxplot"),
                                              plotOutput("norm_kmeans_box")
                              ),
                              bsCollapsePanel(title="cluster count data:",
                                              value="norm_kmeans_extract_count",
                                              fluidRow(
                                                column(4, htmlOutput("norm_select_file2"))
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
                               span("Select gene list files (txt, csv, xlsx)"),
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
                                           img(src="input_format_enrich.png", width = 400,height = 250)), 
                             placement = "right",options = list(container = "body")),
                   hiddenFixedRadio("Level_enrich", "gene_level"),
                   fluidRow(
                     column(6, selectInput("Species4", "Species", species_list, selected = "not selected")),
                     conditionalPanel(condition=c("input.Species4 != 'not selected' && input.Species4 != 'Homo sapiens' &&
                   input.Species4 != 'Mus musculus' && input.Species4 != 'Rattus norvegicus' &&
                   input.Species4 != 'Drosophila melanogaster' && input.Species4 != 'Caenorhabditis elegans' &&
                   input.Species4 != 'Bos taurus' && input.Species4 != 'Canis lupus familiaris' &&
                   input.Species4 != 'Danio rerio' && input.Species4 != 'Gallus gallus' &&
                   input.Species4 != 'Macaca mulatta' && input.Species4 != 'Pan troglodytes' &&
                   input.Species4 != 'Saccharomyces cerevisiae' && input.Species4 != 'Sus scrofa' &&
                   input.Species4 != 'Xenopus laevis' && input.Species4 != 'Arabidopsis thaliana'"),
                     column(6, selectInput("Ortholog4", strong(
                       span("Ortholog"),
                       span(icon("info-circle"), id = "Ortholog_enrich", 
                            options = list(template = popoverTempate))
                     ), orgDb_list, selected = "Mus musculus"),
                     bsPopover("Ortholog_enrich", "Ortholog for the pathway analysis of non-model organisms", 
                               content=paste(img(src="non-model organism.png", width = 500,height = 800)), 
                               placement = "right",options = list(container = "body"))),
                     column(12, selectInput("Biomart_archive4", "Biomart host", ensembl_archive)))
                     ),
                   sliderInput("enrich_showCategory", "Most significant pathways",
                              min = 1, max = 20, value = 5,step = 1),
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
                              htmlOutput("pre_enrich_input_choice"),
                              htmlOutput("enrich_input_choice"),
                              dataTableOutput('enrichment_input')
                     ),
                     tabPanel("Enrichment analysis",
                              fluidRow(
                                column(4, textOutput("Spe3"),
                                       tags$head(tags$style("#Spe3{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")), actionButton("preview_enrich_dot", "Download dot plot")),
                                column(4, selectInput("Gene_set3", "Gene Set", gene_set_list)),
                                column(4, htmlOutput("Custom_input"))
                              ),
                              plotOutput("enrichment3"),
                              fluidRow(
                                column(4, htmlOutput("whichGroup")),
                                column(4, actionButton("preview_enrich_cnet", "Download cnet plot"))
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
                                        accept = c("txt", "csv", "xlsx"),
                                        multiple = FALSE,
                                        width = "80%"),
                              bsPopover("icon14", "File format (txt, csv, or xlsx): ", 
                                        content=paste("First column must be gene name (Gene symbol or ENSEMBL ID).<br>", 
                                                      "The file must contain", strong("log2FoldChange"), "and", strong("padj"), "columns.<br>",
                                                      "You can use a pair-wise DEG result file as input.<br><br>", 
                                                      img(src="input_format_volcano.png", width = 480,height = 230)), 
                                        placement = "right",options = list(container = "body")),
                              hiddenFixedRadio("Level_volcano", "gene_level"),
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
                                column(6, selectInput("Species5", "Species", species_list, selected = "not selected")),
                                conditionalPanel(condition=c("input.Species5 != 'not selected' && input.Species5 != 'Homo sapiens' &&
                   input.Species5 != 'Mus musculus' && input.Species5 != 'Rattus norvegicus' &&
                   input.Species5 != 'Drosophila melanogaster' && input.Species5 != 'Caenorhabditis elegans' &&
                   input.Species5 != 'Bos taurus' && input.Species5 != 'Canis lupus familiaris' &&
                   input.Species5 != 'Danio rerio' && input.Species5 != 'Gallus gallus' &&
                   input.Species5 != 'Macaca mulatta' && input.Species5 != 'Pan troglodytes' &&
                   input.Species5 != 'Saccharomyces cerevisiae' && 
                   input.Species5 != 'Xenopus laevis' && input.Species5 != 'Arabidopsis thaliana'"),
                                                 column(6, selectInput("Ortholog5", strong(
                                                   span("Ortholog"),
                                                   span(icon("info-circle"), id = "Ortholog_volcano", 
                                                        options = list(template = popoverTempate))
                                                 ), orgDb_list, selected = "Mus musculus"),
                                                 bsPopover("Ortholog_volcano", "Ortholog for the pathway analysis of non-model organisms", 
                                                           content=paste(img(src="non-model organism.png", width = 500,height = 800)), 
                                                           placement = "right",options = list(container = "body"))),
                                column(12, selectInput("Biomart_archive5", "Biomart host", ensembl_archive)))
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
                                           column(4, actionButton("preview_volcano_plot", "Download volcano plot")),
                                           column(4, actionButton("preview_volcano_heat", "Download heatmap"))
                                         ),
                                         fluidRow(
                                           column(4, htmlOutput("deg_GOI_color_type")),
                                           column(4, htmlOutput("deg_volcano_x"),actionButton("GOIreset_deg", "GOI reset")),
                                           column(4, htmlOutput("deg_volcano_y"))
                                         ),
                                         conditionalPanel(condition="input.deg_GOI_color_type=='pathway'",
                                                          column(6,selectInput("deg_GOI_color_pathway1","Select a gene set",choices = ""),
                                                                 selectInput("deg_GOI_color_pathway2","",choices = ""))             
                                         ),
                                         selectizeInput("degGOI", "genes of interest (GOI)", choices = NULL, multiple = TRUE,
                                                        options = list(delimiter = " ", create = TRUE, plugins = list("remove_button"), persist = FALSE)),
                                         htmlOutput("deg_uniqueID_cut"),
                                         fluidRow(
                                           column(8, plotOutput("deg_volcano1",brush = "plot1_brush_volcano")),
                                           column(4, plotOutput("deg_GOIheatmap"))
                                         ),
                                         div(
                                           plotOutput("deg_GOIbox", height = "100%"),
                                           style = "height: calc(100vh  - 100px)"
                                         ),
                                         fluidRow(
                                           column(4, actionButton("preview_volcano_box", "Download boxplot"))
                                         ),
                                         bsCollapse(id="DEGlist_volcano_collapse_panel",open="up_panel",multiple = FALSE,
                                                    bsCollapsePanel(title="up (red) list:",
                                                                    value="up_panel",
                                                                    downloadButton("download_uplist_volcano", "Download"),
                                                                    dataTableOutput('uplist_volcano')
                                                    ),
                                                    bsCollapsePanel(title="down (blue) list:",
                                                                    value="down_panel",
                                                                    downloadButton("download_downlist_volcano", "Download"),
                                                                    dataTableOutput('downlist_volcano')
                                                    )
                                         )
                                         )
                              )
                            )
                          ) #sidebarLayout
                 ),
                 tabPanel("ENSEMBL ID to SYMBOL",
                          sidebarLayout(
                            # ENSEMBL ID to SYMBOL---------------------------------
                            sidebarPanel(
                              fileInput("data_file_ens",
                                        strong(
                                          span("Select gene list files (txt, csv, xlsx)"),
                                          span(icon("info-circle"), id = "icon_venn", 
                                               options = list(template = popoverTempate))
                                        ),
                                        accept = c("txt", "csv", "xlsx"),
                                        multiple = FALSE,
                                        width = "80%"),
                              bsPopover("icon_venn", "Gene list files (txt, csv, or xlsx):", 
                                        content= paste("The first column is", strong("gene name"), ".<br>", 
                                                       "The second and subsequent columns do not affect the analysis.<br>", 
                                                       img(src="venn_input.png", width = 400,height = 300)),
                                        placement = "right",options = list(container = "body")),
                              fluidRow(
                                column(6, selectInput("Species_ens", "Species", species_list, selected = "not selected")),
                                conditionalPanel(condition=c("input.Species_ens != 'not selected' && input.Species_ens != 'Homo sapiens' &&
                   input.Species_ens != 'Mus musculus' && input.Species_ens != 'Rattus norvegicus' &&
                   input.Species_ens != 'Drosophila melanogaster' && input.Species_ens != 'Caenorhabditis elegans' &&
                   input.Species_ens != 'Bos taurus' && input.Species_ens != 'Canis lupus familiaris' &&
                   input.Species_ens != 'Danio rerio' && input.Species_ens != 'Gallus gallus' &&
                   input.Species_ens != 'Macaca mulatta' && input.Species_ens != 'Pan troglodytes' &&
                   input.Species_ens != 'Saccharomyces cerevisiae' && input.Species_ens != 'Sus scrofa' &&
                   input.Species_ens != 'Xenopus laevis' && input.Species_ens != 'Arabidopsis thaliana'"),
                                                 column(6, selectInput("Ortholog_ens", strong(
                                                   span("Ortholog"),
                                                   span(icon("info-circle"), id = "Ortholog_enrich", 
                                                        options = list(template = popoverTempate))
                                                 ), orgDb_list, selected = "Mus musculus"),
                                                 bsPopover("Ortholog_enrich", "Ortholog for the pathway analysis of non-model organisms", 
                                                           content=paste(img(src="non-model organism.png", width = 500,height = 800)), 
                                                           placement = "right",options = list(container = "body"))),
                                                 column(12, selectInput("Biomart_archive_ens", "Biomart host", ensembl_archive)))
                              ),
                              actionButton("goButton_ens", "example data (mouse)"),
                              tags$head(tags$style("#goButton_ens{color: black;
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
                              h4("Input file"),
                              dataTableOutput('input_ens'),
                              h4("Result file"),
                              downloadButton("download_ens2symbol", "Download"),
                              textOutput("Spe_ens"),
                              tags$head(tags$style("#Spe_ens{color: red;
                                 font-size: 20px;
            font-style: bold;
            }")),
                              dataTableOutput('input_ens2symbol'),
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
                 tabPanel("Reference", value = "reference",
                          uiOutput("reference_tab")
                 ),
                 tabPanel("Change log",value = "log",
                          uiOutput("change_log_tab")
                 ),
                 tabPanel("SessionInfo",
                          fluidRow(
                            column(12,h2("SessionInfo:"),verbatimTextOutput("sessionInfo"))
                          )
                 )
      )
    )
  ))
