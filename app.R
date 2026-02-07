# ===========================
# Boyce Lab C-CDA and FHIR Data Converter 
#
# Supports:
# - C-CDA (CDA XML):
#     * .txt (single document)
#     * .csv table export with columns: patient_id, doc   (doc = full XML per row)
# - FHIR:
#     * .txt (each non-empty line is a JSON object: Bundle or Resource)
#     * .csv table export with columns: patient_id, record_body (record_body = JSON per row)
#
# Export formats:
#   (1) Long (recommended)  -> tidy rows; preserves patient_id (from CSVs)
#   (2) Wide (1 row)        -> collapses duplicates then transposes; keeps one patient_id column
#   (3) FHIR split ZIP      -> one CSV per resourceType; preserves patient_id
# ===========================

options(shiny.maxRequestSize = 250 * 1024^2)  # 250 MB (adjust as needed)

library(shiny)
library(xml2)
library(jsonlite)
library(dplyr)
library(tidyr)
library(DT)
library(tibble)
library(zip)      # install.packages("zip") if needed
library(readr)
library(tools)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------
# Helpers
# ---------------------------

clean_long <- function(df) {
  df %>%
    mutate(
      path  = trimws(path),
      value = trimws(value)
    ) %>%
    filter(!is.na(value), value != "")
}

# Flatten nested lists into (path, value) rows, preserving indices (good for repeated fields)
flatten_list <- function(x, path = NULL) {
  if (is.null(x)) return(NULL)
  
  if (is.atomic(x)) {
    if (length(x) == 0) return(NULL)
    
    if (length(x) > 1) {
      out <- lapply(seq_along(x), function(i) {
        data.frame(
          path  = paste0(path, "[", i, "]"),
          value = as.character(x[[i]]),
          stringsAsFactors = FALSE
        )
      })
      return(do.call(rbind, out))
    }
    
    return(data.frame(path = path, value = as.character(x), stringsAsFactors = FALSE))
  }
  
  if (is.list(x)) {
    nms <- names(x)
    
    if (is.null(nms)) {
      out <- lapply(seq_along(x), function(i) {
        flatten_list(x[[i]], paste0(path, "[", i, "]"))
      })
      out <- Filter(Negate(is.null), out)
      if (!length(out)) return(NULL)
      return(do.call(rbind, out))
    }
    
    out <- lapply(nms, function(n) {
      flatten_list(x[[n]], if (is.null(path)) n else paste0(path, ".", n))
    })
    out <- Filter(Negate(is.null), out)
    if (!length(out)) return(NULL)
    return(do.call(rbind, out))
  }
  
  NULL
}

# Transpose (name/value) into one-row wide
transpose_data <- function(data) {
  t_data <- t(data)
  colnames(t_data) <- t_data[1, ]
  t_data <- t_data[-1, , drop = FALSE]
  as.data.frame(t_data, stringsAsFactors = FALSE)
}

# ---------------------------
# C-CDA parsing
# ---------------------------

# Parse a single C-CDA XML string into long rows (returns NULL if parsing fails)
parse_ccda_xml_string <- function(xml_string, doc_id = "doc1", patient_id = NA_character_, row_index = NA_integer_) {
  if (is.na(xml_string) || !nzchar(trimws(xml_string))) return(NULL)
  
  doc <- tryCatch(read_xml(xml_string), error = function(e) NULL)
  if (is.null(doc)) return(NULL)
  
  # Reduce redundancy: keep only nodes with non-whitespace text
  nodes <- xml_find_all(doc, "//*[local-name()='ClinicalDocument']//*[normalize-space(text())]")
  if (length(nodes) == 0) return(NULL)
  
  tibble(
    doc_id = doc_id,
    source = "C-CDA",
    patient_id = as.character(patient_id),
    resource_type = NA_character_,
    row_index = row_index,
    bundle_entry_index = NA_integer_,
    resource_id = NA_character_,
    path = xml_name(nodes),
    value = xml_text(nodes)
  )
}

# Parse C-CDA from:
# - .txt (single XML document)
# - .csv table export with columns: patient_id, doc
parse_ccda <- function(file, doc_id = "doc1") {
  ext <- tolower(file_ext(file))
  
  if (ext == "csv") {
    df <- readr::read_csv(file, show_col_types = FALSE, progress = FALSE)
    
    if (!"doc" %in% names(df)) stop("C-CDA CSV must contain a column named 'doc' (XML).")
    if (!"patient_id" %in% names(df)) stop("C-CDA CSV must contain a column named 'patient_id'.")
    
    out <- lapply(seq_len(nrow(df)), function(i) {
      parse_ccda_xml_string(
        xml_string = df$doc[[i]],
        doc_id = doc_id,
        patient_id = df$patient_id[[i]],
        row_index = i
      )
    })
    
    out <- do.call(rbind, Filter(Negate(is.null), out))
    if (is.null(out) || nrow(out) == 0) {
      return(tibble(
        doc_id = character(),
        source = character(),
        patient_id = character(),
        resource_type = character(),
        row_index = integer(),
        bundle_entry_index = integer(),
        resource_id = character(),
        path = character(),
        value = character()
      ))
    }
    return(as_tibble(out))
  }
  
  # .txt single document
  xml_content <- paste(readLines(file, warn = FALSE), collapse = "\n")
  out <- parse_ccda_xml_string(xml_content, doc_id = doc_id, patient_id = NA_character_, row_index = NA_integer_)
  if (is.null(out)) {
    return(tibble(
      doc_id = character(),
      source = character(),
      patient_id = character(),
      resource_type = character(),
      row_index = integer(),
      bundle_entry_index = integer(),
      resource_id = character(),
      path = character(),
      value = character()
    ))
  }
  out
}

# ---------------------------
# FHIR parsing
# ---------------------------

# Parse FHIR from a vector of JSON strings (+ optional patient_id vector)
parse_fhir_from_json_vector <- function(json_vec, patient_id_vec = NULL, doc_id = "doc1") {
  all_rows <- list()
  
  process_resource <- function(res, row_index, bundle_entry_index, patient_id) {
    rt  <- res$resourceType %||% NA_character_
    rid <- res$id %||% NA_character_
    
    flat <- flatten_list(res, path = NULL)
    if (is.null(flat) || nrow(flat) == 0) return(NULL)
    
    flat$doc_id <- doc_id
    flat$source <- "FHIR"
    flat$patient_id <- patient_id %||% NA_character_
    flat$resource_type <- rt
    flat$row_index <- row_index
    flat$bundle_entry_index <- bundle_entry_index
    flat$resource_id <- rid
    flat
  }
  
  for (i in seq_along(json_vec)) {
    js <- json_vec[[i]]
    if (is.na(js) || !nzchar(trimws(js))) next
    
    pid <- NA_character_
    if (!is.null(patient_id_vec) && length(patient_id_vec) >= i) {
      pid <- as.character(patient_id_vec[[i]])
    }
    
    fhir_data <- tryCatch(fromJSON(js, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(fhir_data) || !is.list(fhir_data)) next
    
    # Bundle-like if it has entry[]
    is_bundle_like <- !is.null(fhir_data$entry) && is.list(fhir_data$entry)
    
    if (is_bundle_like) {
      entries <- fhir_data$entry
      for (j in seq_along(entries)) {
        res <- entries[[j]]$resource
        # fallback: sometimes entry itself is a resource
        if (is.null(res) && !is.null(entries[[j]]$resourceType)) res <- entries[[j]]
        
        if (!is.null(res) && !is.null(res$resourceType)) {
          all_rows[[length(all_rows) + 1]] <- process_resource(
            res,
            row_index = i,
            bundle_entry_index = j,
            patient_id = pid
          )
        }
      }
    } else {
      # Single resource
      if (!is.null(fhir_data$resourceType)) {
        all_rows[[length(all_rows) + 1]] <- process_resource(
          fhir_data,
          row_index = i,
          bundle_entry_index = NA_integer_,
          patient_id = pid
        )
      }
    }
  }
  
  out <- do.call(rbind, Filter(Negate(is.null), all_rows))
  if (is.null(out) || nrow(out) == 0) {
    return(tibble(
      doc_id = character(),
      source = character(),
      patient_id = character(),
      resource_type = character(),
      row_index = integer(),
      bundle_entry_index = integer(),
      resource_id = character(),
      path = character(),
      value = character()
    ))
  }
  
  out <- out[, c("doc_id","source","patient_id","resource_type","row_index",
                 "bundle_entry_index","resource_id","path","value")]
  as_tibble(out)
}

# Parse FHIR from:
# - .csv with columns: record_body + patient_id
# - .txt where each non-empty line is JSON (no patient_id)
parse_fhir <- function(file, doc_id = "doc1") {
  ext <- tolower(file_ext(file))
  
  if (ext == "csv") {
    df <- readr::read_csv(file, show_col_types = FALSE, progress = FALSE)
    
    if (!"record_body" %in% names(df)) stop("FHIR CSV must contain a column named 'record_body'.")
    if (!"patient_id" %in% names(df))  stop("FHIR CSV must contain a column named 'patient_id'.")
    
    json_vec <- as.character(df$record_body)
    pid_vec  <- df$patient_id
    return(parse_fhir_from_json_vector(json_vec, patient_id_vec = pid_vec, doc_id = doc_id))
  }
  
  # .txt fallback
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  parse_fhir_from_json_vector(lines, patient_id_vec = NULL, doc_id = doc_id)
}

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Boyce Lab C-CDA and FHIR Data Converter"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput(
        "file",
        "Choose C-CDA (.txt or .csv) or FHIR (.txt or .csv)",
        accept = c(".txt", ".csv")
      ),
      
      # Label as C-CDA, keep internal value "CDA" to avoid breaking logic if you prefer.
      radioButtons("fileType", "File Type", choices = c("C-CDA" = "CCDA", "FHIR" = "FHIR")),
      
      radioButtons(
        "exportFormat",
        "Export Format",
        choices = c(
          "Long (recommended)" = "long",
          "Wide (1 row)" = "wide",
          "FHIR split (ZIP: one CSV per resource type)" = "fhir_zip"
        ),
        selected = "long"
      ),
      
      actionButton("parse", "Parse Document"),
      downloadButton("downloadData", "Download")
    ),
    mainPanel(
      width = 9,
      DTOutput("parsedData")
    )
  )
)

# ---------------------------
# Server
# ---------------------------
server <- function(input, output, session) {
  parsed_data <- reactiveVal()
  
  observeEvent(input$parse, {
    req(input$file)
    
    doc_id <- paste0("doc-", format(Sys.time(), "%Y%m%d-%H%M%S"))
    
    data <- tryCatch({
      if (input$fileType == "CCDA") {
        parse_ccda(input$file$datapath, doc_id = doc_id)
      } else {
        parse_fhir(input$file$datapath, doc_id = doc_id)
      }
    }, error = function(e) {
      tibble(
        doc_id = doc_id,
        source = input$fileType,
        patient_id = NA_character_,
        resource_type = NA_character_,
        row_index = NA_integer_,
        bundle_entry_index = NA_integer_,
        resource_id = NA_character_,
        path = "Error",
        value = as.character(e)
      )
    })
    
    parsed_data(clean_long(data))
  })
  
  output$parsedData <- renderDT({
    req(parsed_data())
    df <- parsed_data()
    
    # Show key provenance columns; resource_type/resource_id relevant mostly for FHIR
    datatable(
      df[, c("patient_id","value","path","resource_type","row_index","bundle_entry_index","resource_id")],
      options = list(
        pageLength = 100,
        autoWidth = TRUE,
        searchHighlight = TRUE,
        scrollX = TRUE,
        ordering = FALSE
      ),
      class = "display compact",
      rownames = FALSE
    )
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      base <- paste0("parsed_data-", Sys.Date(), "-", input$exportFormat)
      if (input$exportFormat == "fhir_zip") paste0(base, ".zip") else paste0(base, ".csv")
    },
    content = function(file) {
      df <- parsed_data()
      req(df)
      
      # ---- Long (default) ----
      if (input$exportFormat == "long") {
        write.csv(df, file, row.names = FALSE)
        return()
      }
      
      # ---- Wide (1 row) ----
      if (input$exportFormat == "wide") {
        # Keep one patient_id value (if present) in the wide output
        pid <- df$patient_id[!is.na(df$patient_id)][1] %||% NA_character_
        
        wide_pairs <- df %>%
          select(path, value) %>%
          group_by(path) %>%
          summarise(value = paste(unique(value), collapse = " | "), .groups = "drop") %>%
          rename(name = path)
        
        t_data <- transpose_data(wide_pairs)
        t_data <- cbind(patient_id = pid, t_data)
        write.csv(t_data, file, row.names = FALSE)
        return()
      }
      
      # ---- FHIR split ZIP (one CSV per resource_type) ----
      if (input$exportFormat == "fhir_zip") {
        fhir <- df %>% filter(source == "FHIR")
        if (nrow(fhir) == 0) {
          writeLines("FHIR split export selected, but current parsed data is not FHIR.", con = file)
          return()
        }
        
        fhir <- fhir %>%
          mutate(resource_type = if_else(
            is.na(resource_type) | trimws(resource_type) == "",
            "Unknown",
            resource_type
          ))
        
        out_dir <- tempfile("fhir_split_")
        dir.create(out_dir)
        
        by_rt <- split(fhir, fhir$resource_type)
        csv_paths <- character(0)
        
        for (rt in names(by_rt)) {
          safe_rt <- gsub("[^A-Za-z0-9_\\-]+", "_", rt)
          csv_path <- file.path(out_dir, paste0(safe_rt, ".csv"))
          write.csv(by_rt[[rt]], csv_path, row.names = FALSE)
          csv_paths <- c(csv_paths, csv_path)
        }
        
        zip::zipr(zipfile = file, files = csv_paths)
        return()
      }
    }
  )
}

shinyApp(ui = ui, server = server)
