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
                                      radioButtons('sep', 'Separator',
                                                   c(Tab='\t'),
                                                   '\t'),
                                      checkboxInput("testme", "Try a Sample Dataset!", 
                                                    value = FALSE)
                                    ),
                                    mainPanel(
                                      numericInput("filter_level", "Enter Filter Level",
                                                   min = 0, max = 1, value = 0.9),
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
                                   fluidRow(
                                     h4("Group Selection"),
                                     textOutput("tutorialGroup"),
                                     tags$head(tags$style(
                                       "#tutorialGroup{color: red; font-size: 18px}")),
                                     textOutput("group_warning"),
                                     tags$head(tags$style(
                                       "#group_warning{color: red; font-size: 18px}")),
                                     # place to hold dynamic inputs
                                     uiOutput("inputGroup"))
                                   )
                                 )
                        )),
                        tabPanel("Explore Your Data",
                                 tabsetPanel(
                                 tabPanel("Explore Taxa",
                                          mainPanel(
                                               fluidRow(column(4, sliderInput("taxaDims", 
                                                                     "Taxa Level",
                                                                     5, min = 1, max = 8)),
                                                        column(4, sliderInput("taxa_limit",
                                                                              "Taxa Proportions",
                                                                              1, min = 0, max = 1, step = 0.05,
                                                                              value = c(0.75,1)))),
                                               fluidRow(column(6, uiOutput("ex_delimiter"))),
                                               tags$head(tags$style(
                                                 "#ex_delimiter {color:white; font-size:16px}")),
                                               fluidRow(uiOutput("TaxaDimExp")),
                                               fluidPage(
                                                 fluidRow(column(12, plotOutput("key1", height = '50px', width = '100%'))),
                                                 fluidRow(column(12, plotlyOutput("species_explore")))),
                                               fluidRow(downloadButton("species_download", "Download Plot"), 
                                                        downloadButton("species_legend_download", "Download Legend"),
                                                        downloadButton("species_raw_data", "Download Table")),
                                               fluidRow(uiOutput("taxa_selectors")),
                                               #fluidPage(
                                               fluidRow(textOutput("curr_taxa_exp")),
                                               fluidRow(column(12, plotOutput("key2", height = '50px', width = "100%"))),
                                               fluidRow(column(12, uiOutput("da_taxa_heat_UI"))),
                                               fluidRow(column(12, uiOutput("da_taxa_stat_ui"))),
                                               width = 11)),
                                               #fluidRow(tableOutput("da_taxa_labs")),
                                               #width = 11)),
                                 tabPanel("Explore Genes",
                                          mainPanel(
                                            textOutput("instructions"),
                                            fluidRow(uiOutput("feat_selectors")),
                                            fluidRow(textOutput("curr_select_exp")),
                                            fluidRow(
                                              column(6, fluidPage(
                                                plotOutput("key4", height = '50px'),
                                                plotlyOutput("gene_da"))),
                                              column(6, fluidPage(
                                                plotOutput("key5", height = '50px'),
                                                plotlyOutput("Plot3")))),
                                            fluidRow(column(6, fluidPage(uiOutput("da_feat_stat_ui")))),
                                            #fluidPage(
                                             # fluidRow(column(8, plotOutput("key3", height = '40px'))),
                                            #  fluidRow(column(8, plotlyOutput("gene_explore")))),
                                            fluidRow(downloadButton("gene_explore_download", "Download Plot")),
                                            width = 12))
                        )),
                 tabPanel("Query Your Data",
                          tabsetPanel(
                            tabPanel("Gene Search",
                                     fluidRow(column(4,
                                                     selectizeInput('acc_list', choices=NULL,
                                                                    label = h3("Select a Feature"),
                                                                    multiple = TRUE))#,
                                              #column(10, checkboxInput('xy_switch', label = 'group by Gene Family?')
                                              ),
                                     hr(),
                                     fluidPage(
                                       fluidRow(textOutput("curr_select_search")),
                                       fluidRow(column(7, fluidPage(plotOutput("key7", height = '50px'),
                                                                    plotlyOutput("plot1", height = '500px'))),
                                                #column(7,plotOutput("plot1", height = '500px', click = 'plot_click'))
                                                column(5, fluidPage(plotOutput("key6", height = '50px'),
                                                                    plotlyOutput("spec_select", height = '450px')))
                                                ),
                                       fluidRow(column(7, fluidPage(uiOutput("select_feat_stat_ui"))))
                                       ),
                                     fluidRow(downloadButton("expression_download", "Download Plot")),
                                     #h5("Pairwise T-test results are performed here if at least 2 groups are selected"),
                                     fluidRow(column(8,uiOutput("exp_heat"))),
                                     fluidRow(uiOutput("exp_heat_download")),
                                     fluidRow(column(6, align= "center", tableOutput("exp_table"))),
                                     fluidRow(downloadButton("expression_table_download", "Download Table"))
                            ),
    
                            tabPanel("Correlation", selectizeInput('sig_select', choices=NULL,
                                                                   label = h3("Begin by selecting two gene families of interest"),
                                                                   multiple = TRUE),
                                     fluidPage(column(4, verbatimTextOutput("sig_message"))),
                                     fluidPage(column(10, plotOutput("corr_plot"))),
                                     fluidRow(downloadButton("corr_download", "Download Plot")),
                                     fluidPage(column(10, uiOutput("group_corrs"))),
                                     fluidPage(uiOutput("group_download")),
                                     fluidPage(column(10, tableOutput("corr_labels"))),
                                     fluidPage(downloadButton("corr_table_download", "Download Labels")))
                            )
                          ),
                        tabPanel(title = loadingLogo("https://www.youtube.com/watch?v=pIgZ7gMze7A", "wham_logo_trans.png",
                                                            'wham_grey_inf.gif', height = 135, width = 280)),
                 tags$head(tags$style(
                   '.navbar { font-size: 18px}',
                   '.navbar-brand {font-size:32px}')),
                 tags$style(HTML(".control-label {font-size:20px; font-weight:normal}")),
                 tags$style(HTML(".radio label {font-size:18px; font-weight:normal}")),
                 tags$style(HTML(".checkbox label {font-size:20px; font-weight:normal}")),
                 tags$style(HTML("label {font-size:20px; font-weight:normal}")),
                 tags$style(HTML(".btn {font-size:18px; font-weight:normal}")),
                 tags$head(tags$style("*{ font-family: Helvetica; }"))
                                                            
)