#*------------------------------------------------------*
#' Author: Julien Chevreau
#' Mail: julien.chevreau(at)univ-rouen.fr
#' Date: 08/10/25
#' Affiliation: University of Rouen Normandie
#' Code: Server part of my differential expression analysis project
#' NB: Listen to All Good Things, STARSET and Alpine Universe
#' *------------------------------------------------------*
#' TODO : edit the dowload handler; 

srcfile("global.R") # source libraries

server <- function(input, output) {

  # Check and load user data ####
  user_data <- reactive({
    req(input$user_File)  # Halt execution in absence of input file
    
    user_file <- input$user_File
    ext <- tools::file_ext(user_file$datapath)
    ## Checking extension
    if (ext != "csv"){
      showNotification("Provided file extension is not csv. Make sure it is well formatted.", 
                       duration=5, closeButton = T, type = "warning")
    }
    ## Read the csv file if checks passed
    read.csv(user_file$datapath, header = T, sep=",")
    })
  
  ## Shortcuts to maxima and minima
  top_log_value = reactive({ # Max Log2 Fold change value
    ceiling(max(abs(user_data()$log2FC)))
  })
  
  max_padj = reactive({ # Max Negative log Padjusted value
    ceiling(max(-10*log(user_data()$padj)))
  })
  
  min_padj = reactive({ # Min Negative log Padjusted value
    floor(min(-10*log(user_data()$padj)))
  })
  
  # DataTable ####
  output$user_table = renderDT(datatable(user_data(),filter = "top", selection = "multiple"))
  
  ## Download table data 
  output$tableData = downloadHandler(filename = "iris.csv", content = iris)
  
  ## Personnalized sliders ####
  ## Use dataframe values to scale sliders
  output$slider_log2FC = renderUI({ # Log2FC slider
    sliderInput(
      inputId = "slider_log2FC",
      label = "log2FC threshold",
      min = 0,
      max = top_log_value(), # Dynamic max value of slider
      value = top_log_value()/2,
      step = 0.1
    )
  })
  
  output$slider_padj = renderUI({ # Adjusted Pvalue slider
    sliderInput(
      inputId = "slider_padj",
      label = "P-value threshold",
      min = 0,
      max = max_padj(), # Dynamic max value of slider
      value = max_padj()/2,
      step = 1
    )
  })
  
  ## Volcano plot ####
  output$volcanoPlot <- renderPlotly({
    req(user_data()) # Require data to proceed
    
    ### Prepare dataframe ####
    local_dataframe = user_data() # Load as df
    local_dataframe$padj <- as.numeric(as.character(local_dataframe$padj)) # Force numeric interpretation
    local_dataframe$log2FC <- as.numeric(as.character(local_dataframe$log2FC)) # Force numeric interpretation
    local_dataframe$negLogPadj = -10*log(local_dataframe$padj) # Precompute negative log of padj
    
    ### Define expression groups ####
    local_dataframe$groups = ifelse( # Over expressed genes
        local_dataframe$negLogPadj>input$slider_padj & local_dataframe$log2FC>input$slider_log2FC, "Sur-exprimé",
      ifelse ( # Under expressed genes
        local_dataframe$negLogPadj>input$slider_padj & local_dataframe$log2FC< -input$slider_log2FC, "Sous-exprimé",
        "Non significatif") # Genes under thresholds
      )
    # Assign colors to groups
    group_colors <- c(
      "Sur-exprimé" = "red",
      "Sous-exprimé" = "black",
      "Non significatif" = "gray"
    )
    
    ### Split plot traces by groups ####
    plot_ly() %>% add_trace( # Over expressed genes
      data = local_dataframe[local_dataframe$group == "Sur-exprimé",], # Rows chosen
      x = ~log2FC, y = ~negLogPadj, # Columns chosen
      type = "scatter", mode = "markers", # Type of plot
      marker = list(color = group_colors["Sur-exprimé"]), # Color of points
      name = "Sur-exprimé", # Name of variable in legend
      text = ~paste0( # Text to show on hover
        "<b>Gene:</b> ", GeneName, "<br>",
        "log2FC: ", round(log2FC, 2), "<br>",
        "-10log(Padj): ", round(-10*log(padj), 2), "<br>"
      ),
      hoverinfo = "text"
    ) %>% add_trace( #Under expressed genes
        data = local_dataframe[local_dataframe$group == "Sous-exprimé",],
        x = ~log2FC, y = ~negLogPadj,
        type = "scatter", mode = "markers",
        marker = list(color = group_colors["Sous-exprimé"]),
        name = "Sous-exprimé",
        text = ~paste0(
          "<b>Gene:</b> ", GeneName, "<br>",
          "log2FC: ", round(log2FC, 2), "<br>",
          "-10log(Padj): ", round(-10*log(padj), 2), "<br>"
        ),
        hoverinfo = "text"
    ) %>% add_trace( # Genes under thresholds
        data = local_dataframe[local_dataframe$group == "Non significatif",],
        x = ~log2FC, y = ~negLogPadj,
        type = "scatter", mode = "markers",
        marker = list(color = group_colors["Non significatif"]),
        name = "Non significatif",
        text = ~paste0(
          "<b>Gene:</b> ", GeneName, "<br>",
          "log2FC: ", round(log2FC, 2), "<br>",
          "-10log(Padj): ", round(-10*log(padj), 2), "<br>",
          "<i>Ne dépasse pas les seuils</i>"),
        hoverinfo = "text"
      ) %>%
      ### Change layout and add threshold lines ####
      layout(
        title = "Dynamic volcano plot",
        xaxis = list(title = "Log2 Fold Change"),
        yaxis = list(title = "-10.Log(Padj)"),
        shapes = list(
          # Log2FC Vertical lines
          list( # Rightmost line
            type = "line",
            x0 = input$slider_log2FC, x1 = input$slider_log2FC,
            y0 = min_padj(), y1=max_padj(),
            line = list(color = "firebrick", dash = "dash")
          ),
          list( #Leftmost line
            type = "line",
            x0 = -input$slider_log2FC, x1 = -input$slider_log2FC,
            y0 = min_padj(), y1=max_padj(),
            line = list(color = "firebrick", dash = "dash")
          ),
          # Padj horizontal line
          list(
            type = "line",
              x0 = floor(min(local_dataframe$log2FC)), x1 = top_log_value(),
            y0 = input$slider_padj, y1 = input$slider_padj,
            line = list(color = "indianred", dash = "dash")
          )
        )
      )
  })
  
  # Plot-Table interactions ####
  observe({
    # Retrieve selected points on plot
    selected_points <- event_data("plotly_selected")
    # Reset if nothing is selected
    if (is.null(selected_points)) {
      table_proxy <- dataTableProxy("table")
      selectRows(table_proxy, NULL)
      return()
    }
    # Find corresponging indices in the table
    selected_rows <- selected_points$pointNumber + 1  # +1 because R is 1-indexed
    table_proxy <- dataTableProxy("table")
    selectRows(table_proxy, selected_rows) # Select rows
  })
}
