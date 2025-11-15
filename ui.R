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
  ),
  # Special html code to specify h1 headers colors between light and dark modes
  adminlte_css = "
    h1 {
      color: #FFFFFF !important;  /* Light mode; Light text color because background is dark */
    }
    .dark-mode h1 {
      color: #FFFFFF !important;  /* Dark mode; Light text color */
    }
  "
)


# Dashboard main visible content ####
dashboardPage(
  freshTheme = shyning, # Theme to use (custom; see above)
  
  ## Header ####
  dashboardHeader(title = "DontBShyning"),
  
  
  ## Sidebar ####
  dashboardSidebar( 
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")), # Homepage
      fileInput( # Button to upload input data
        inputId = "user_File", 
        label="Télécharger un fichier",
        buttonLabel = "Parcourir...",
        placeholder = "Aucun fichier sélectionné",
        accept = c("text/csv","text/plain"),
        multiple = FALSE),
      selectInput( # List to select the data to display
        inputId = "selectChoice",
        label = "A vous le choix !",
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
              h1("Bienvenue au manoir !"),
              fluidRow( # Zone de plot et sélection
                box(title = "Histogram",
                    width = 7,
                    status = "maroon", 
                    solidHeader = T,
                    background = NULL,
                    plotlyOutput("volcanoPlot", height = 250)
                  ),
                box(
                  title="Seuils",
                  width = 5,
                  uiOutput("slider_log2FC"),
                  uiOutput("slider_padj"),
                  downloadButton('tableData', 'Download Data')
                )
              ),
              fluidRow(
                box(DTOutput('user_table'), width = 12, title ="Tableau",
                    status = "maroon", 
                    solidHeader = T))
      ),
      ### Second tab content ####
      tabItem(tabName = "wholeData",
              h2("Overview des données")
      ),
      ### Troisième tab content ####
      tabItem(tabName = "GO",
              h2("Go, GO Power rangers !")
      ),
      ### Quatrième tab content ####
      tabItem(tabName = "pathway",
              h2("I know the path ! You know the way !")
      ),
      ### Cinquième tab content ####
      tabItem(tabName = "about",
              h2("Shall be me")
      )
    )
  )
)
