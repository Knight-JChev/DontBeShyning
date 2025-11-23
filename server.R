#*------------------------------------------------------*
#' Author: Julien Chevreau
#' Mail: julien.chevreau(at)univ-rouen.fr
#' Date: 08/10/25
#' Affiliation: University of Rouen Normandie
#' Code: Server part of my differential expression analysis project
#' NB: Listen to All Good Things, STARSET and Alpine Universe
#' *------------------------------------------------------*
#' TODO : edit the dowload handler;
#' 

srcfile("global.R") # source libraries

server <- function(input, output) {
  
  # Check and load user data ####
  user_data <- reactiveVal(NULL) # Initialise reactive value variable to sotre uploaded data
  required_cols <- c("GeneName", "ID", "baseMean", "log2FC", "pval", "padj") # Mandatory cols
  # What happens when data is loaded
  observeEvent(input$user_file, {
    req(input$user_file) # Require upload
    
    # Safely try to read file
    df <- tryCatch(
      read.csv(input$user_file$datapath, header = TRUE),
      error = function(e) NULL
    )
    
    # Send popup if not readable
    if (is.null(df)) { 
      showModal(modalDialog( # Popup
        title = "Erreur de fichier",
        "Impossible de lire le fichier. Vérifiez qu'il s'agit bien d'un CSV.",
        easyClose = TRUE # Click anywhere to close
      ))
      return(NULL)
    }
    
    # Retrieve missing cols
    missing_cols <- setdiff(required_cols, colnames(df))
    
    # Stop loading if any cols is missing
    if (length(missing_cols) > 0) {
      showModal(modalDialog( # Popup
        title = "Colonnes manquantes",
        paste(
          "Le fichier ne contient pas les colonnes suivantes :",
          paste(missing_cols, collapse = ", ")
        ),
        easyClose = TRUE
      ))
      
      # Block file usage
      return(NULL)
    }
    # ️Load data if all checks were passed
    user_data(df)
    
  }) # End of observeEvent
  
  ## Shortcuts to maxima and minima
  top_log_value = reactive({ # Max Log2 Fold change value
    req(user_data())
    ceiling(max(abs(user_data()$log2FC)))
  })
  max_padj = reactive({ # Max Negative log Padjusted value
    req(user_data())
    ceiling(max(-10*log(user_data()$padj)))
  })
  min_padj = reactive({ # Min Negative log Padjusted value
    req(user_data())
    floor(min(-10*log(user_data()$padj)))
  })
  
  # DataTable ####
  output$user_table = renderDT(datatable(user_data(),filter = "top", selection = "multiple"))
  
  ## Download table data 
  output$tableData = downloadHandler(filename = "iris.csv", content = iris)
  
  ## Personnalized sliders ####
  ## Use dataframe values to scale sliders
  output$slider_log2FC = renderUI({ # Log2FC slider
    req(user_data())
    
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
    req(user_data())
      
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
    plot_ly(source ="volcSource") %>% add_trace( # Over expressed genes
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
  ## Table to plot interaction ####
  observeEvent(input$user_table_rows_selected, {
    volcanoProxy = plotlyProxy("volcanoPlot")
    # Unselect everything
    plotlyProxyInvoke(volcanoProxy, "restyle", list(opacity = 1), list(0))
    
    if (!is.null(input$user_table_rows_selected) &&
        length(input$user_table_rows_selected) > 0) { 
      # When row is selected in the table, update the plot
      plotlyProxyInvoke(
        volcanoProxy,
        "restyle",
        list(marker = list(color = ifelse(
          seq_len(nrow(user_data())) %in% input$user_table_rows_selected,
          "red",
          "blue"
        ))),
        list(0)
      )
    }
  })
  
  
  ## Plot to table interaction ####
  observe({
    # Retrieve selected points on plot
    selected_points <- event_data("plotly_selected", source = "volcSource")
    # Reset if nothing is selected
    if (is.null(selected_points)) {
      table_proxy <- dataTableProxy("user_table")
      selectRows(table_proxy, NULL)
      return()
    }
    # Find corresponging indices in the table
    selected_rows <- selected_points$pointNumber + 1  # +1 because R is 1-indexed
    table_proxy <- dataTableProxy("user_table")
    selectRows(table_proxy, selected_rows) # Select rows
  })
  
  ## Reactive table of selected values ####
  output$user_selected_table <- renderDT({
    req(input$user_table_rows_selected)  # Need selected rows
    isolate({
      if (!is.null(input$user_table_rows_selected) && length(input$user_table_rows_selected) > 0) {
        selected_rows <- user_data()[input$user_table_rows_selected, ] # Save selected rows
      } else { # Empty data if nothing selected
        selected_rows <- data.frame(GeneName = character(0), log2FC = numeric(0), padj = numeric(0))
      }
      datatable(selected_rows, options = list(pageLength = 5)) # Return the datatable
    })
  })
  
  ## Reactive table of significant values ####
  ## Reactive table of selected values ####
  output$user_significant_table <- renderDT({
    req(user_data())  # Need data
    
    ### Prepare dataframe ####
    local_dataframe = user_data() # Load as df
    local_dataframe$padj <- as.numeric(as.character(local_dataframe$padj)) # Force numeric interpretation
    local_dataframe$log2FC <- as.numeric(as.character(local_dataframe$log2FC)) # Force numeric interpretation
    local_dataframe$negLogPadj = -10*log(local_dataframe$padj) # Precompute negative log of padj
    
  })
}
user_data[user_data$padj<1e-30 & user_data$log2FC< -2,]
