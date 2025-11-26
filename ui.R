#' *------------------------------------------------------*
#' Author: Julien Chevreau
#' Mail: julien.chevreau(at)univ-rouen.fr
#' Date: 08/10/25
#' Affiliation: University of Rouen Normandie
#' Code: UI part of my differential expression analysis project
#' 
#' ProTip: On RStudio, quickly navigate through sections thanks to 
#' document outline (top right and bottom left of editor)
#' 
#' NB: Listen to Chris Christodolou, Christopher Larkin and 
#' Gentle love
#' *------------------------------------------------------*
#' TODO : Patch la sidebar à gauche avec le texte qui fait n'imp
#' Faire un about et un home
#' Patch le titre en noir -> blanc

srcfile("global.R") # source libraries

  # Personal color theme ####
  shyning <- create_theme(
    bs4dash_vars(
      # Light theme colors
      navbar_light_color = "#F44336",       
      navbar_light_active_color = "#FF0000",
      navbar_light_hover_color = "#FF5252",
      # Dark theme colors
      sidebar_dark_bg = "#121212",          
      sidebar_dark_color = "#E0E0E0",       
      sidebar_dark_hover_color = "#FF1744", 
      sidebar_dark_active_color = "#FF5252" 
    ),
    bs4dash_layout(
      main_bg = "#000000"          # Black background
    ),
    
    # Redefine status colors 
    bs4dash_status(
      primary = "#D32F2F",
      danger = "#B71C1C", 
      success = "#388E3C",   
      warning = "#F57C00",
      info = "#1976D2"    
    ),
    bs4dash_color(
      gray_900 = "#0D0D0D", 
      white = "#FFFFFF",
      red = "#C62828"
    ),
    bs4dash_yiq(
      contrasted_threshold = 150,  # high contrast for readability
      text_dark = "#E0E0E0",       # light grey on dark background
      text_light = "#FFFFFF"       # white for high contrasts
    )
  )
  
  #function to set titles and text in white
  white_h1 = function(text) {h1(text, style = "color: white;")}
  white_h2 = function(text) {h2(text, style = "color:white")}
  white_p = function(text) {p(text, style = "color:white")}
  
  
  # Dashboard main visible content ####
  dashboardPage(
    freshTheme = shyning, # Theme to use (custom; see above)
    dark = NULL, # Remove dark mode switch button
    help = NULL, # Remove the help button
    
    ## Header ####
    dashboardHeader(title = "DontBShyning"),
  
  ## Sidebar ####
  dashboardSidebar( 
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")), # Homepage 
      fileInput( # Button to upload input data
        inputId = "user_file", 
        label="Upload your file",
        buttonLabel = "Browse...",
        placeholder = "No file selected",
        accept = c("text/csv","text/plain"),
        multiple = FALSE),
      selectInput( # List to select the data to display
        inputId = "selectChoice",
        label = "The choice is yours !",
        choices =  list(`Metroidvania` = list("Hollow Knight", "Ori", "Guacamelee"),
                        `Roguelike` = list("Hades", "Dead cells", "Binding of Isaac"),
                        `Deckbuilding` = list("Monster Train", "Slay the Spire", "Vault of the Void")),
        multiple = FALSE
      ),
      # Other tabs
      menuItem("Whole Data Inspection", tabName = "wholeData", icon = icon("eye")),
      menuItem("Go Term Enrichment", tabName = "GO", icon = icon("money-bill-1")),
      menuItem("Pathway Enrichment", tabName = "pathway", icon = icon("flask-vial")),
      menuItem("About", tabName = "about", icon = icon("address-card"))
    )
  ),
  
  ## Body ####
  dashboardBody(
    tabItems(
      ### Premier tab content ####
      tabItem(tabName = "home",
              white_h1("Home sweet home !"),
              fluidRow( # Zone de plot et sélection
                box(title = "Histogram",
                    width = 7,
                    status = "maroon", 
                    background = NULL,
                    plotlyOutput("volcanoPlot", height = 250),
                    p("You can download the plot in PNG in the plot utils")
                  ),
                box(
                  title="Seuils",
                  status = "maroon",
                  width = 5,
                  uiOutput("slider_log2FC"),
                  uiOutput("slider_padj")
                )
              ),
              fluidRow(
                box(
                  title = "Data Tables",
                  width = 12,
                  status = "maroon",
                  tabsetPanel(
                    id = "table_panels",
                    tabPanel( # Panel of user data (formatted)
                      title ="Your data",
                      width = 12,
                      status = "maroon", 
                      DTOutput('user_table')
                    ),
                    tabPanel( # Panel of selected genes with associated table
                      title ="Selected genes",
                      width = 12, 
                      status = "maroon", 
                      HTML("You can select genes in 'Your Data' table or use plotly selection tools. </br>
                        Selection reset might not show on the plot itself.</br>"),
                      downloadButton('downloadSelected', 'Download CSV of selected genes', icon=icon("download")),
                      actionButton("reset_selection", "Reset selection", icon = icon("arrow-rotate-right")),
                      DTOutput('user_selected_table')
                   ),
                   tabPanel( # Panel of over-expressed genes with associated table
                     title = "Over-expressed Genes",
                     width = 12, 
                     status = "maroon", 
                     downloadButton('downloadOver', 'Download CSV of over-expressed genes', icon=icon("download")),
                     DTOutput("over_expressed_table")),
                   tabPanel( # Panel of under-expressed genes with associated table
                     "Under-expressed Genes",
                     width = 12, 
                     status = "maroon", 
                     downloadButton('downloadUnder', 'Download CSV of under-expressed genes', icon=icon("download")),
                     DTOutput("under_expressed_table"))
                 )
                )
              )
      ),
      ### Second tab content ####
      tabItem(tabName = "wholeData",
              white_h2("Look at all this data ! Too dark to be seen, unfortunately.")
      ),
      ### Troisième tab content ####
      tabItem(tabName = "GO",
              white_h2("Go, GO Power rangers !")
      ),
      ### Quatrième tab content ####
      tabItem(tabName = "pathway",
              white_h2("I know the path ! You know the way !")
      ),
      ### Cinquième tab content ####
      tabItem(tabName = "about",
              white_h2("Don't Be Shyning"),
              white_p("This name is a play on word as I confuse tend to write shyni instead of shiny from the sound of it."),
              fluidRow(
                box(
                  title = "The project",
                  status = "maroon",
                  p("The goal of this application is to process RNA-seq data to analyse gene expression
                    and show nice plots, eventually to find biologically relevant clues"),
                  p("Truth is, this project objective up to now was to try on Rshiny.")),
                box(
                  title = "How to use",
                  status = "maroon",
                  p("For now, only the Home page has content you can upload a csv file from the sidebar menu.
                  Mandatory columns are 'GeneName', 'ID', 'baseMean', 'log2FC', 'pval', 'padj'
                  Once this is done, a volcano plot should show on the left panel of Home page.
                  The sliders in the right panels allow you to change the thresholds shown on the graph.
                  Finally, in the bottom panel you will find different tables to download.")
                  ),
                box(
                  title = "App status",
                  status = "maroon",
                  p("The last part of the project involves group development and there
                    is not so much chance for my app to be the template so, honestly, I won't put much more effort in this app."),
                  p("The experience was not really to my liking, mostly frustrating, so I'll have fun coding elsewhere. 
                    I find the frontend/backend mix weird and unintuitive in Rshiny."),
                  p("Still, if you wish to use this code, here's the Github link :"),
                  a("DontBShyning",
                      href = "https://github.com/Knight-JChev/DontBeShyning",
                      target = "_blank"
                  ),
                  
                  )
                )# End of fluidRow
              ) # End of tabItem
      ) # End of tabItems
    ) # End of dashboardBody
  ) # End of dashboardPage
