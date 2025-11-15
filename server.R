#' *------------------------------------------------------*
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
  
  top_log_value = reactive({
    ceiling(max(abs(user_data()$log2FC)))
  })
  
  max_pvalue = reactive({
    ceiling(max(-10*log(user_data()$p_value)))
  })
  
  min_pvalue = reactive({
    floor(min(-10*log(user_data()$p_value)))
  })
  
  # DataTable ####
  output$user_table = renderDT(datatable(user_data(),filter = "top", selection = "multiple"))
  
  ## Download table data 
  output$tableData = downloadHandler(filename = "iris.csv", content = iris)
  
  ## Personnalized sliders
  ## Use dataframe values to scale sliders
  output$slider_LFC = renderUI({ 
    sliderInput(
      inputId = "slider_Log2FC",
      label = "Log2FC threshold",
      min = 0,
      max = top_log_value(),
      value = top_log_value()
    )
  })
  
  output$slider_Pval = renderUI({
    sliderInput(
      inputId = "slider_Pvalue",
      label = "P-value threshold",
      min = 0,
      max = ceiling(max(user_data()$p_value)),
      value = ifelse(max_pvalue()>0.05, 0.05, max_pvalue()/2)
    )
  })
  
  ## Volcano plot ####
  output$volcanoPlot <- renderPlotly({
    req(user_data())
    plot_ly(
      data = user_data(),
      x = ~log2FC,  # remplacer par vos colonnes
      y = -10*log(user_data()$p_value),
      type = "scatter",
      mode = "markers",
      text = ~paste("<b>Gene:</b>", Gene, 
                    '<br>Log2FC:', round(user_data()$log2FC, digits = 2),
                    '<br>P-value:', round(user_data()$pvalue, digits = 2))
    ) %>%
      layout(
        title = "Dynamic volcano plot",
        shapes = list(
          # P values vertical lines
          list(
            type = "line",
            x0 = input$slider_Log2FC, x1 = input$slider_Log2FC,
            y0 = min_pvalue(), y1=max_pvalue(),
            line = list(color = "black", dash = "dash")
          ),
          # Ligne horizontale du seuil Y
          list(
            type = "line",
            x0 = floor(min(user_data()$log2FC)), x1 = top_log_value(),
            y0 = input$slider_Pvalue, y1 = input$slider_Pvalue,
            line = list(color = "black", dash = "dash")
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
