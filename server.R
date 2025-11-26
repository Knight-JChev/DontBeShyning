#*------------------------------------------------------*
#' Author: Julien Chevreau
#' Mail: julien.chevreau(at)univ-rouen.fr
#' Date: 08/10/25
#' Affiliation: University of Rouen Normandie
#' Code: Server part of my differential expression analysis project
#' NB: Listen to All Good Things, STARSET and Alpine Universe
#' *------------------------------------------------------*
#' TODO : Patch la vérification de fichier
#' Patch le réinitialiser

srcfile("global.R") # source libraries

server = function(input, output, session) {
  
  # Check and load user data ####
  user_data = reactiveVal(NULL) # Initialise reactive value variable to sotre uploaded data
  required_cols = c("GeneName", "ID", "baseMean", "log2FC", "pval", "padj") # Mandatory cols
  # What happens when data is loaded
  observeEvent(input$user_file, {
    req(input$user_file) # Require upload
    
    # Safely try to read file
    df = tryCatch(
      read.csv(input$user_file$datapath, header = TRUE),
      error = function(e) NULL
    )
    
    # Send popup if not readable
    if (is.null(df)) { 
      showModal(modalDialog( # Popup
        title = "File type error",
        "Can't read file. Check it is CSV-formatted.",
        easyClose = TRUE # Click anywhere to close
      ))
      return(NULL)
    }
    
    # Retrieve missing cols
    missing_cols = setdiff(required_cols, colnames(df))
    
    # Stop loading if any cols is missing
    if (length(missing_cols) > 0) {
      showModal(modalDialog( # Popup
        title = "Missing cols",
        paste(
          "The file does not have the following cols : \n",
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
  
  ## Shortcuts to maxima and minima ####
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
      label = "Adjusted P-value threshold",
      min = 0,
      max = max_padj(), # Dynamic max value of slider
      value = max_padj()/2,
      step = 1
    )
  })
  
  user_formatted = reactive({
    req(user_data()) # Require data to proceed
    
    ### Prepare dataframe ####
    local_dataframe = user_data() # Load as df
    local_dataframe$padj = as.numeric(as.character(local_dataframe$padj)) # Force numeric interpretation
    local_dataframe$log2FC = as.numeric(as.character(local_dataframe$log2FC)) # Force numeric interpretation
    local_dataframe$negLogPadj = -10*log(local_dataframe$padj) # Precompute negative log of padj
    
    ### Define expression groups ####
    local_dataframe$groups = ifelse( # Over expressed genes
      local_dataframe$negLogPadj>input$slider_padj & local_dataframe$log2FC>input$slider_log2FC, "Over-expressed",
      ifelse ( # Under expressed genes
        local_dataframe$negLogPadj>input$slider_padj & local_dataframe$log2FC< -input$slider_log2FC, "Under-expressed",
        "Under thresholds") # Genes under thresholds
    )
    local_dataframe
  })
  
  ## Volcano plot ####
  volcanoPlot = reactive({
    req(user_formatted())
    local_dataframe = user_formatted()
    
    # Assign colors to groups
    group_colors = c(
      "Over-expressed" = "red",
      "Under-expressed" = "black",
      "Under thresholds" = "gray"
    )
    
    ### Split plot traces by groups ####
    plot_ly(source ="volcSource") %>% add_trace( # Over expressed genes
      data = local_dataframe[local_dataframe$group == "Over-expressed",], # Rows chosen
      x = ~log2FC, y = ~negLogPadj, # Columns chosen
      type = "scatter", mode = "markers", # Type of plot
      marker = list(color = group_colors["Over-expressed"]), # Color of points
      name = "Over-expressed", # Name of variable in legend
      text = ~paste0( # Text to show on hover
        "<b>Gene:</b> ", GeneName, "<br>",
        "log2FC: ", round(log2FC, 2), "<br>",
        "-10log(Padj): ", round(-10*log(padj), 2), "<br>"
      ),
      hoverinfo = "text"
    ) %>% add_trace( #Under expressed genes
      data = local_dataframe[local_dataframe$group == "Under-expressed",],
      x = ~log2FC, y = ~negLogPadj,
      type = "scatter", mode = "markers",
      marker = list(color = group_colors["Under-expressed"]),
      name = "Under-expressed",
      text = ~paste0(
        "<b>Gene:</b> ", GeneName, "<br>",
        "log2FC: ", round(log2FC, 2), "<br>",
        "-10log(Padj): ", round(-10*log(padj), 2), "<br>"
      ),
      hoverinfo = "text"
    ) %>% add_trace( # Genes under thresholds
      data = local_dataframe[local_dataframe$group == "Under thresholds",],
      x = ~log2FC, y = ~negLogPadj,
      type = "scatter", mode = "markers",
      marker = list(color = group_colors["Under thresholds"]),
      name = "Under thresholds",
      text = ~paste0(
        "<b>Gene:</b> ", GeneName, "<br>",
        "log2FC: ", round(log2FC, 2), "<br>",
        "-10log(Padj): ", round(-10*log(padj), 2), "<br>",
        "<i>Under thresholds</i>"),
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
  # Send volcano plot to UI
  output$volcanoPlot = renderPlotly({volcanoPlot()})
  
  # Plot-Table interactions ####
  default_colors = reactive({
    req(user_formatted())
    
    group_colors = c("Over-expressed" = "red", "Under-expressed" = "black", "Under thresholds" = "gray")
    colors = group_colors[user_formatted()$groups]
    colors
  })
  
  ## Plot to table interaction ####
  observe({
    # Retrieve selected points on plot
    selected_points = event_data("plotly_selected", source = "volcSource")
    # Reset if nothing is selected
    if (is.null(selected_points)) {
      table_proxy = dataTableProxy("user_table")
      selectRows(table_proxy, NULL)
      return()
    }
    # Find corresponging indices in the table
    selected_rows = selected_points$pointNumber + 1  # +1 because R is 1-indexed
    table_proxy = dataTableProxy("user_table")
    selectRows(table_proxy, selected_rows) # Select rows
  })
  
  # UI DataTables ####
  # User Table formatted ####
  output$user_table = renderDT(datatable(user_formatted(),filter = "top", selection = "multiple"))
  
  ## Reactive table of selected values ####
  output$user_selected_table = renderDT({
    req(input$user_table_rows_selected)  # Need selected rows
    isolate({
      if (!is.null(input$user_table_rows_selected) && length(input$user_table_rows_selected) > 0) {
        selected_rows = user_formatted()[input$user_table_rows_selected, ] # Save selected rows
      } else { # Empty data if nothing selected
        selected_rows = data.frame(GeneName = character(0), log2FC = numeric(0), padj = numeric(0))
      }
      datatable(selected_rows, selection = "none", options = list(pageLength = 10), filter = "top") # Return the datatable
    })
  })
  
  ### Reactive table of over/under expressed genes ####
  output$over_expressed_table = renderDT({
    req(user_formatted()) # require formatted data
    local_dataframe = user_formatted()
    
    over_expressed = user_data()[local_dataframe$negLogPadj > input$slider_padj & local_dataframe$log2FC > input$slider_log2FC, ]
    datatable(over_expressed, options = list(pageLength = 10), selection = "none", filter = "none")
  })
  
  output$under_expressed_table = renderDT({ # Same on under-expressed genes
    req(user_formatted()) # require formatted data
    local_dataframe = user_formatted()
    
    under_expressed = user_data()[local_dataframe$negLogPadj > input$slider_padj & local_dataframe$log2FC < -input$slider_log2FC, ]
    datatable(under_expressed, options = list(pageLength = 10), selection = "none", filter = "none")
  })
  
  ## Download buttons events for each table ####
  # Selected genes
  output$downloadSelected = downloadHandler(filename = "selected_genes.csv", 
                                            content =  function(file) { # select the proper table to write
                                              req(input$user_table_rows_selected)
                                              selected_rows = user_formatted()[input$user_table_rows_selected, ]
                                              write.csv(selected_rows, file, row.names = FALSE)
                                            },
                                            contentType = "text/csv")
  # Over-Expressed genes
  output$downloadOver = downloadHandler(filename = "over_expressed_genes.csv",
                                        content = function(file) {
                                          req(user_formatted())
                                          selected_rows = user_formatted()[input$user_table_rows_selected, ]
                                          write.csv(user_formatted()[user_formatted()$groups=="Over-expressed",],
                                                    file, row.names = FALSE)
                                        },
                                        contentType = "text/csv")
  # Under-Expressed genes
  output$downloadUnder = downloadHandler(filename = "under_expressed_genes.csv",
                                         content = function(file) {
                                           req(user_formatted())
                                           selected_rows = user_formatted()[input$user_table_rows_selected, ]
                                           write.csv(user_formatted()[user_formatted()$groups=="Under-expressed",],
                                                     file, row.names = FALSE)
                                         },
                                         contentType = "text/csv")

  
  # UI plot ####
  ## Reset selection ####
  observeEvent(input$reset_selection, {
    # Reset selction in table
    table_proxy = dataTableProxy("user_table")
    selectRows(table_proxy, NULL)
    
    # Reset colors in plots
    volcanoProxy = plotlyProxy("volcanoPlot", session)
    plotlyProxyInvoke(volcanoProxy, "restyle", list(list(selectedpoints = NULL)))
    plotlyProxyInvoke(volcanoProxy, "restyle", list(list("marker.opacity" = 1)))
  })
}
