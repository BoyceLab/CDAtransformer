library(shiny)
library(xml2)
library(jsonlite)
library(dplyr)
library(tidyr)
library(DT)
library(devtools)

# Increase the maximum upload size to 30MB
options(shiny.maxRequestSize = 30 * 1024^2)

# Function to parse CDA document and extract relevant elements
parse_cda <- function(file) {
  file_content <- tryCatch({
    readLines(file, warn = FALSE)
  }, error = function(e) {
    stop("Error reading file: ", e)
  })

  xml_content <- paste(file_content, collapse = "\n")

  doc <- tryCatch({
    read_xml(xml_content)
  }, error = function(e) {
    stop("Error parsing XML: ", e)
  })

  elements <- tryCatch({
    xml_find_all(doc, "//*[local-name()='ClinicalDocument']//*")
  }, error = function(e) {
    stop("Error finding elements: ", e)
  })

  if (length(elements) == 0) {
    stop("No elements found in the ClinicalDocument section.")
  }

  element_list <- lapply(elements, function(node) {
    name <- xml_name(node)
    value <- xml_text(node)
    list(name = name, value = value)
  })

  element_list <- do.call(rbind, lapply(element_list, function(x) as.data.frame(x, stringsAsFactors = FALSE)))

  return(element_list)
}

# Function to recursively flatten a list
flatten_list <- function(x, name = NULL) {
  if (is.atomic(x)) {
    return(data.frame(name = name, value = as.character(x), stringsAsFactors = FALSE))
  } else if (is.list(x)) {
    out <- do.call(rbind, lapply(names(x), function(n) {
      flatten_list(x[[n]], paste0(name, if (!is.null(name)) ".", n))
    }))
    return(out)
  }
}

# Function to parse multiple FHIR resources and flatten them
parse_fhir <- function(file) {
  file_content <- tryCatch({
    readLines(file, warn = FALSE)
  }, error = function(e) {
    stop("Error reading file: ", e)
  })

  all_flattened <- lapply(file_content, function(json_content) {
    fhir_data <- tryCatch({
      fromJSON(json_content)
    }, error = function(e) {
      stop("Error parsing JSON: ", e)
    })

    flatten_list(fhir_data)
  })

  combined_flattened <- do.call(rbind, all_flattened)

  return(combined_flattened)
}

# Function to transpose data frame
transpose_data <- function(data) {
  t_data <- t(data)
  colnames(t_data) <- t_data[1,]
  t_data <- t_data[-1, , drop = FALSE]
  return(as.data.frame(t_data, stringsAsFactors = FALSE))
}

# Define UI
ui <- fluidPage(
  titlePanel("Boyce Lab C-CDA and FHIR Data Converter"),
  tags$head(
    tags$style(HTML("
      .data-table th, .data-table td { padding: 10px; }
      .data-table tr:nth-child(even) { background-color: #f2f2f2; }
      .data-table th { background-color: #4CAF50; color: white; }
      .section-spacing { margin-top: 20px; margin-bottom: 20px; }
      .dataTables_wrapper .dataTables_scrollBody { overflow: visible; }
      .highlight { background-color: yellow; }
      .dataTables_wrapper .dataTables_scroll { padding-left: 0px; }
      .dataTables_wrapper .dataTables_scrollHead { margin-bottom: -1px; }
      .shiny-output-error { color: red; }
      .dt-left { text-align: left; }
      .dt-right { text-align: right; }
      .dt-center { text-align: center; }
      .dt-top { vertical-align: top; }
      .dt-bottom { vertical-align: bottom; }
      .dataTables_wrapper .dataTables_scrollBody { overflow: auto; }
      .main-panel .dataTables_wrapper { margin-left: -20px; margin-right: -20px; } /* Reduce side margins to remove extra space */
    "))
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("file", "Choose CDA or FHIR Text File", accept = c(".txt")),
      radioButtons("fileType", "File Type", choices = c("CDA", "FHIR")),
      actionButton("parse", "Parse Document"),
      downloadButton("downloadData", "Download CSV")
    ),
    mainPanel(
      width = 9,
      div(class = "main-panel", DTOutput("parsedData"))
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  parsed_data <- reactiveVal()

  observeEvent(input$parse, {
    req(input$file)
    data <- tryCatch({
      if (input$fileType == "CDA") {
        parse_cda(input$file$datapath)
      } else if (input$fileType == "FHIR") {
        parse_fhir(input$file$datapath)
      }
    }, error = function(e) {
      data.frame(name = "Error", value = as.character(e))
    })

    parsed_data(data)
  })

  output$parsedData <- renderDT({
    req(parsed_data())
    datatable(
      parsed_data()[, c("value", "name")], # Switch order of columns
      options = list(
        pageLength = 100, # Default to 100 rows
        autoWidth = TRUE,
        searchHighlight = TRUE,
        scrollX = TRUE,
        ordering = FALSE, # Disable column ordering
        rowCallback = JS('function(row, data, index) {
                            $(row).find("td:eq(0)").css("text-align", "left");
                            $(row).find("td:eq(1)").css("text-align", "right");
                         }'),
        columnDefs = list(
          list(targets = 0, className = 'dt-left'), # Value column
          list(targets = 1, className = 'dt-right') # Name column
        )
      ),
      class = "display compact",
      rownames = FALSE # Hide row numbers
    )
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("parsed_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      data <- parsed_data()
      t_data <- transpose_data(data)
      write.csv(t_data, file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
