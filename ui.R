library(shiny)
library(shinythemes)

ui <- navbarPage(title = "Workflow Hub for Automated Metagenomic Exploration",
                 theme = shinytheme("superhero"),
                       tabPanel("Home",
                                  fluidRow(column(7, 
                                                  wellPanel(h3("Methods Workflow"),
                                                            img(src='Figure_1.png', height = 647.7, width = 500))),
                                           column(5,
                                                  mainPanel(h3("Resources"),
                                                            textOutput("resource_text1"),
                                                            uiOutput("samp_url"),
                                                            htmlOutput("resource_text2"),
                                                            uiOutput("hmp_url"),
                                                            htmlOutput("humann2wham_text"),
                                                            uiOutput("humann2wham_url"),
                                                            h3("External Information"),
                                                            textOutput("ext_text"),
                                                            uiOutput("ruggles_url"),
                                                            uiOutput("github"),
                                                            uiOutput("github_url"),
                                                            uiOutput("citation"))))
                                ),
                       tabPanel("Upload",
                                fluidPage(
                                  sidebarLayout(
                                    sidebarPanel(uiOutput("file_selector"),
                                      tags$hr(),
                                      #checkboxInput('header', 'Header', TRUE),
                                      radioButtons('input_type', "Input Source",
                                                   choiceNames = c("Biobakery", "EBI"),
                                                   choiceValues = c("Biobakery", "EBI")),
                                      checkboxInput("testme", "Try a Sample Dataset!", 
                                                    value = FALSE)
                                    ),
                                    mainPanel(
                                      fluidRow(
                                        column(7, sliderInput("filter_level", "Variance Filter",
                                                   min = 0, max = 1, value = 0.75))),
                                      fluidRow(
                                        column(7, uiOutput("filter_message1"))),
                                      tags$head(tags$style(
                                        "#filter_message1{font-size: 18px}")),
                                      fluidRow(column(7, uiOutput("filter_message2"))),
                                      tags$head(tags$style(
                                        "#filter_message2{color: #df691a; font-size: 18px}")),
                                      uiOutput("preview_shower")
                                    )
                                  ))
                                
                       ),
                        tabPanel("Groups",
                                 fluidPage(
                                   titlePanel("Assign each sample to a Group"),
                                   sidebarLayout(
                                     sidebarPanel(
                                       numericInput("numInputs", "Select Number of Groups", 1, min = 1, max = 10),
                                       uiOutput("group_pre")
                                     ),
                                    mainPanel(
                                   fluidPage(
                                     h4("Group Selection"),
                                     textOutput("tutorialGroup"),
                                     tags$head(tags$style(
                                       "#tutorialGroup{color: red; font-size: 18px}")),
                                     textOutput("group_warning"),
                                     tags$head(tags$style(
                                       "#group_warning{color: red; font-size: 18px}")),
                                     fluidRow(uiOutput("key0")),
                                     downloadButton("legend_download", "Download Legend"),
                                     # place to hold dynamic inputs
                                     uiOutput("inputGroup"))
                                   )
                                 )
                        )),
                        tabPanel("Explore Your Data",
                                 tabsetPanel(
                                 tabPanel("Explore Taxa",
                                          mainPanel(
                                               fluidRow(column(5, sliderInput("taxaDims", 
                                                                     "Taxa Level", value = 6, 
                                                                     min=6, max = 7, step = 1)),
                                                        column(5, sliderInput("taxa_limit",
                                                                              "Taxa Proportions",
                                                                              1, min = 0, max = 1, step = 0.05,
                                                                              value = c(0.75,1)))),
                                               fluidRow(column(5, uiOutput("ex_delimiter")),
                                                        column(5, uiOutput("prop_exp"))),
                                               tags$head(tags$style(
                                                 "#ex_delimiter {color:#df691a; font-size:18px}")),
                                               tags$head(tags$style(
                                                 "#prop_exp {color:#df691a; font-size:18px}")),
                                               #fluidRow(uiOutput("TaxaDimExp")),
                                               fluidPage(
                                                 fluidRow(column(12, uiOutput("key1"))),
                                                 fluidRow(column(12, plotlyOutput("species_explore")))),
                                               fluidRow(downloadButton("species_download", "Download Plot"), 
                                                        downloadButton("species_legend_download", "Download Legend"),
                                                        downloadButton("species_raw_data", "Download Table")),
                                               #fluidPage(
                                               fluidRow(uiOutput("da_taxa_ui")),
                                               width = 11)),
                                               #fluidRow(tableOutput("da_taxa_labs")),
                                               #width = 11)),
                                 tabPanel("Explore Features",
                                          mainPanel(
                                            textOutput("instructions"),
                                            fluidRow(uiOutput("feat_selectors")),
                                            fluidRow(textOutput("curr_select_exp"),
                                                     tags$head(tags$style("#curr_select_exp{font-size: 20px}"))),
                                            fluidRow(
                                              column(6, fluidPage(
                                                uiOutput("key4"),
                                                plotlyOutput("gene_da"))),
                                              column(6, fluidPage(
                                                uiOutput("key5"),
                                                plotlyOutput("Plot3")))),
                                            fluidRow(column(6, fluidPage(downloadButton("gene_explore_download", 
                                                                              "Download Heatmap"))),
                                                     column(6, fluidPage(downloadButton("gene_explore_taxa_download", 
                                                                              "Download Plot"),
                                                            downloadButton("gene_explore_taxa_legend_download", 
                                                                           "Download Legend")))),
                                            fluidRow(column(6, uiOutput("da_feat_stat_ui"))),
                                            width = 12))
                        )),
                 tabPanel("Query Your Data",
                          tabsetPanel(
                            tabPanel("Feature Search",
                                     fluidRow(column(4,
                                                     selectizeInput('acc_list', choices=NULL,
                                                                    label = h3("Select a Feature"),
                                                                    multiple = TRUE))#,
                                              #column(10, checkboxInput('xy_switch', label = 'group by Gene Family?')
                                              ),
                                     hr(),
                                     fluidPage(
                                       fluidRow(textOutput("curr_select_search"),
                                                tags$head(tags$style("#curr_select_search{font-size: 20px}"))),
                                       fluidRow(column(6, fluidPage(uiOutput("key7"),
                                                                    plotlyOutput("plot1", height = '500px'))),
                                                column(6, fluidPage(uiOutput("key6"),
                                                                    plotlyOutput("spec_select", height = '500px')))
                                                ),
                                       fluidRow(column(6, fluidPage(downloadButton("sel_explore_download", 
                                                                                   "Download Heatmap"))),
                                                column(6, fluidPage(downloadButton("sel_explore_taxa_download", 
                                                                                   "Download Plot"),
                                                                    downloadButton("sel_explore_taxa_legend_download", 
                                                                                   "Download Legend")))),
                                       fluidRow(column(6, uiOutput("select_feat_stat_ui")))
                                       )
                            ),
    
                            tabPanel("Correlation", selectizeInput('sig_select', choices=NULL,
                                                                   label = h3("Begin by selecting two Features"),
                                                                   multiple = TRUE),
                                     fluidRow(column(4, verbatimTextOutput("sig_message"))),
                                     fluidPage(column(10, plotOutput("corr_plot"))),
                                     fluidRow(downloadButton("corr_download", "Download Plot")),
                                     fluidPage(column(10, uiOutput("group_corrs"))),
                                     fluidRow(uiOutput("group_download")),
                                     fluidPage(column(10, tableOutput("corr_labels"))),
                                     fluidRow(downloadButton("corr_table_download", "Download Labels")))
                            )
                          ),
                        tabPanel(title = loadingLogo("https://www.youtube.com/watch?v=pIgZ7gMze7A", "wham_logo_trans.png",
                                                            'wham_grey_inf.gif', height = 135, width = 280)),
                 tags$head(tags$style(
                   '.navbar { font-size: 18px}',
                   '.navbar-brand {font-size:32px}')),
                 tags$style(HTML(".control-label {font-size:20px; font-weight:normal}")),
                 tags$style(HTML(".radio label {font-size:20px; font-weight:normal}")),
                 tags$style(HTML(".checkbox label {font-size:20px; font-weight:normal}")),
                 tags$style(HTML("label {font-size:20px; font-weight:normal}")),
                 tags$style(HTML(".btn {font-size:20px; font-weight:normal}")),
                 tags$style(HTML(".shiny-output-error-validation {font-size:18px;color:red}")),
                 tags$head(tags$style("*{ font-family: Helvetica; }"))
                                                            
)