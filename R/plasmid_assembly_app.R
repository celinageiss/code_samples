#' ==============
#' MoClo Assembly
#' ==============
#  Celina Geiss <celina.geiss@stud.uni-heidelberg.de>
#' 
#' This app allows you to select different pieces for each module and assemble
#' a MoClo plasmid. Multiple selections are possible for each piece, so that the 
#' output will include all possible permutations from the selection. 
#' The result of each assembly is a fasta file (single assemblies) or multi-fasta 
#' file, where the respective headers correspond to the selected pieces.
#' 
#' The app also allows you to search the database (see tab 'Browse database') in 
#' an interactive way.


# Package imports
library(shiny)
library(dplyr)

# Load GG sites and extarct as character vector
GG_sites <- readRDS(here::here("results", "GG_sites.rds")) %>% pull(GG_site)
# Load MoClo database from file (NOTE: update file name if database was changed)
moclo_db <- readRDS(here::here("results", "MoClo_database_simple_pieces.rds")) %>% 
  dplyr::arrange(match(right_GG, GG_sites))  # reorder rows according to GG position

pieces <- moclo_db %>% dplyr::pull(piece) %>% unique()

# Function to extract name of selected pieces
extract_pieces <- function(selected_piece) {
  moclo_db %>% 
    dplyr::filter(piece == selected_piece) %>% 
    dplyr::pull(name)
}


#### Define UI ####

ui <- pageWithSidebar(
  
  headerPanel("MoClo assembly"),
  
  sidebarPanel(
    # Select pieces for assembly
    # Note: use selectizeInput instead of selectInput to enable type search and multiple selections
    selectizeInput(inputId = "backbone", 
                   label = "Backbone (BB)", 
                   choices = c(Choose = "",  # Placeholder for empty selection
                               extract_pieces("BB")),  # extract names of all "BB" pieces from moclo_db
                   multiple = TRUE  # multiple selction is possible
    ),
    selectizeInput(inputId = "insulator", 
                   label = "Insulator (B-2)", 
                   choices = c(Choose = "", extract_pieces("B-2")),
                   multiple = TRUE
    ),
    selectizeInput(inputId = "promoter", 
                   label = "Promoter (1-3)", 
                   choices = c(Choose = "", extract_pieces("1-3")),
                   multiple = TRUE
    ),
    selectizeInput(inputId = "UTR5", 
                   label = "5' UTR (2-4)", 
                   choices = c(Choose = "", extract_pieces("2-4")),
                   multiple = TRUE
    ),
    selectizeInput(inputId = "reporter", 
                   label = "Reporter gene (3-5)", 
                   choices = c(Choose = "", extract_pieces("3-5")),
                   multiple = TRUE
    ),
    selectizeInput(inputId = "linker", 
                   label = "Linker (4-6)", 
                   choices = c(Choose = "", extract_pieces("4-6")),
                   multiple = TRUE
    ),
    selectizeInput(inputId = "selection_marker", 
                   label = "Selection marker (5-7)", 
                   choices = c(Choose = "", extract_pieces("5-7")),
                   multiple = TRUE
    ),
    selectizeInput(inputId = "SMAR", 
                   label = "S/MAR (6-8)", 
                   choices = c(Choose = "", extract_pieces("6-8")),
                   multiple = TRUE
    ),
    selectizeInput(inputId = "UTR3", 
                   label = "3' UTR (7-9)", 
                   choices = c(Choose = "", extract_pieces("7-9")),
                   multiple = TRUE
    ),
    selectizeInput(inputId = "polyA", 
                   label = "polyA tail (8-B)", 
                   choices = c(Choose = "", extract_pieces("8-B")),
                   multiple = TRUE
    ),
    
    # Download button for fasta file
    downloadButton(outputId = "download", label = "Download fasta")
  ),
  
  mainPanel(
    tabsetPanel(
      # Display resulting assembly
      tabPanel("Result",
               br(),
               tableOutput("result")
      ),
      # Show full MoClo Database
      tabPanel("Browse database",
               titlePanel("MoClo database"),
               dataTableOutput("moclo_db")
      )
    )
  )
)


#### Define SERVER ####

server <- function(input, output, session) {
  
  # Filter selection interactively from moclo_db
  selected_pieces <- reactive({  
    selection <- reactiveValuesToList(input) %>%  # inputs as list
      unlist()  # necessary to unpack multi-selections
    
    moclo_db %>% 
      dplyr::filter(name %in% selection)  # extract full entries from moclo_db
  })
  
  # List of assembled plasmids
  all_plasmids <- reactive({
    selected_pieces() %>% 
      split(.$piece) %>%  # split into named list by piece
      purrr::map(., ~pull(., name)) %>%  # extract names
      dplyr::bind_cols() %>%  # make table from pieces
      expand.grid() %>%  # permute all combinations (result is dataframe with factors)
      dplyr::distinct() %>%  # remove duplicates
      dplyr::as_tibble() %>%  # convert from dataframe to tibble
      dplyr::mutate_if(is.factor, as.character) %>%  # convert columns back to character
      dplyr::rowwise() %>%  # group by rows
      dplyr::group_split() %>%  # split assemblies into list (each row is one entry)
      purrr::map(., ~filter(moclo_db, name %in% .))  # filter resp. combination from moclo_db
  })
  
  # Create fasta output
  fasta <- reactive({
    purrr::map(all_plasmids(), function(plasmid) {
      res <- plasmid %>% 
        dplyr::mutate(left_GG_sequence = paste0(left_GG, sequence))  # join left GG with sequence
      
      seq <- res %>%
        dplyr::pull(left_GG_sequence) %>%  # extract column as character vector
        stringr::str_c(collapse = "")  # collapse all sequences to one character vector
      
      header <- paste0("> pMoClo ",  # fasta header
                       stringr::str_c(res$name, collapse = " - "),  # pieces
                       " (", nchar(seq), " bp)")  # total length
      
      fasta <- paste(header, seq, sep = "\n")  # join header an sequence with newline in between
      return(fasta)
    })
  })
  
  # Collapse single fasta files to mfasta
  mfasta <- reactive({
    stringr::str_c(fasta(), collapse = "\n\n")
  })
  
  
  #### OUTPUTS ####
  
  # Display table with selection
  output$result <- renderTable(
    selected_pieces() %>% 
      dplyr::mutate(sequence = paste0(stringr::str_sub(sequence, 1, 20), "...")) %>%  # trim displayed sequence to 20 bases
      dplyr::select(piece, name, length, left_GG, sequence, right_GG)  # reorder columns
  )
  
  # Download fasta
  output$download <-  downloadHandler(
    filename = function(){
      if_else(length(fasta()) == 1, 
              # One plasmid: selected pieces as file name
              paste0("pMoClo ",
                     stringr::str_c(selected_pieces()$name, collapse = " - "),
                     ".fa"),
              # Multi-assemblies: unique filename from date-time
              paste0("pMoClo multi ",
                     format(Sys.time(), "%Y%m%d%-%H%M%S"),
                     ".fa")
      )},
    content = function(file) {
      readr::write_file(mfasta(), file)
    }
  )
  
  # Display moclo_db as interactive data table
  output$moclo_db <- renderDataTable(
    moclo_db %>% 
      dplyr::mutate(sequence = paste0(stringr::str_sub(sequence, 1, 30), "..."))  # trim displayed sequence to 30 bases
  )
  
}


#### Call APP ####
shinyApp(ui, server)