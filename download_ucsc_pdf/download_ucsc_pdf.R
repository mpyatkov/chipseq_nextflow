library(httr)
library(rvest)
library(qpdf)
library(readxl)
library(tools)
library(ggpubr)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(stringr)
library(purrr)
library(dplyr)
library(argparser)

ParseArguments <- function() {
  p <- arg_parser('Downloading pdf files from public UCSC browser session')
  p <- add_argument(p,'--session', default="ucsc_session_name", help="Session name")
  p <- add_argument(p,'--user', default="Djw", help="Username")
  p <- add_argument(p,'--db', default="mm9", help="UCSC database, ex. mm9, mm10, hg19,...")
  p <- add_argument(p,"--zoom", default = "3" , help = "zoom factor for each region")
  p <- add_argument(p, "--input_file", default = NA, help = "input bed (without header) / xlsx (with header) file")
  return(parse_args(p))
}

argv <- ParseArguments()

DEBUG <- F
if (DEBUG) {
  argv$session <- "G235_DAR_ATAC_Combined_GroupScaled"
  argv$user <- "Djw"
  argv$db <- "mm9"
  argv$zoom <- "1.5"
  argv$input_file <- "/projectnb/wax-dk/max/G235_DAR_ATAC/RESULTS_G235_DAR_ATAC/summary/manorm2_diffreps_overlap_top25/Female_AAV_Cre_ATAC_vs_Female_Null_ATAC_top25_DOWNREGULATED.xlsx"
}


# download_pdf function
# example: download_pdf("chr1:176696608-176697356", "test.pdf")
download_pdf <- function(session, url, outname) {
  
  # receive pdf link
  pdflink <- session %>% 
    session_jump_to(url) %>% 
    read_html() %>%
    html_nodes("#pdfLink") %>%
    html_attr("href") %>% 
    str_replace(., "..", "https://genome.ucsc.edu/")
  
  final_pdf_link <- session %>% 
    session_jump_to(pdflink) %>% 
    read_html() %>% 
    html_node(xpath = "/html/body/div[5]/ul[1]/li[1]/a") %>% 
    html_attr("href") %>%
    str_replace(., "..", "https://genome.ucsc.edu")
  
  # pdf_link
  tryCatch({
    download.file(final_pdf_link, destfile = outname)
  }, error = function(err) {
    print(paste0("ERROR: Can't process the following pdf file: ", final_pdf_link, ", output file: ", outname))
  })
}

## returns new (zoomed) coordinates of fragment
calculate_zoom_factor <- function(start, end, zoom){
  
  frag_len <-  floor((end-start)/2)
  mid_point <-  start+frag_len
  new_start <- mid_point - zoom*frag_len
  new_end <- mid_point + zoom*frag_len
  
  list(start = floor(new_start),
       end = floor(new_end))
}

## read and save bed/xls as table
export_table_as_pdf <- function(file_path, outdir, add_annotations = TRUE){
  ## for bed file without colnames
  ## for xlsx should be header provided
  
  ## landscape mode
  OUTPUT_HEIGHT <- 8.5
  OUTPUT_WIDTH <- 11
  CHUNK_SIZE <- 30
  
  ext <- tools::file_ext(file_path)
  file_data <- NULL
  
  if (ext == "xlsx" || ext == "xls") {
    file_data <- read_excel(file_path, col_names = T) 
  } else {
    file_data <- read_tsv(file_path, col_names = F) 
  }
  
  ## if annotations are not required make portrait mode
  ## change default number of chunks
  ## show only 3 first column
  if(!add_annotations){
    file_data <- file_data %>% 
      select(chr = 1, start = 2, end = 3)
    
    ## portrait mode
    OUTPUT_HEIGHT <- 11
    OUTPUT_WIDTH <- 8.5
    CHUNK_SIZE <- 40
  }
  
  ## split data.frame if it has more then 40 rows
  
  n <- nrow(file_data)
  r <- rep(1:ceiling(n/CHUNK_SIZE),each=CHUNK_SIZE)[1:n]
  dlist <- split(file_data,r)
  
  map(dlist, function(chunk_table){
    cowplot::plot_grid(ggtexttable(chunk_table, rows = NULL, theme = ttheme("minimal", base_size = 8)))
  }) %>% 
    marrangeGrob(nrow =1, ncol=1) %>% 
    ggsave(str_c(outdir, "/00000_coordinates.pdf"), plot = ., width = OUTPUT_WIDTH, height = OUTPUT_HEIGHT)
  
}

## read bed/xls/xlsx files
read_data <- function(file_path, pdfdir, zoom){
  ext <- tools::file_ext(file_path)
  file_data <- NULL
  
  if (ext == "xlsx" || ext == "xls") {
    file_data <- read_excel(file_path, col_names = T) %>% 
      select(chr = 1, start = 2, end = 3)
  } else {
    file_data <- read_tsv(file_path, col_names = F) %>% 
      select(chr = 1, start = 2, end = 3)  
  }
  
  file_data <- pmap_dfr(file_data, function(chr, start, end){
    new_coords <- calculate_zoom_factor(start, end, zoom)
    tibble(chr = chr, start = new_coords$start, end = new_coords$end)
  }) %>% 
    mutate(id = row_number()) %>% 
    rowwise() %>% 
    mutate(coords = paste0(chr,":",start,"-",end),
           outname = paste0(pdfdir,
                            paste(str_pad(id,width = 3, pad = "0"), 
                                  chr, 
                                  start, 
                                  end, 
                                  sep = "_", collapse = ""), 
                            ".pdf")) %>% 
    select(coords, outname)
  file_data
}

init_local <- function(login, session_name, db) {
  sessionUrl <- paste0("https://genome.ucsc.edu/s/", login,"/", session_name)
  session <- session(URLencode(sessionUrl))
  
  main_url <- "https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A26664573%2D26693533&hgsid=966765435_KMLSgA0F83EB9IJZmt7zQ5uPhnOr"  
  main_url <- str_replace(main_url,"mm9",db)
  main_url <- str_extract(main_url, "(?<=^)(.+)(?<=position=)")
  
  list(session=session, main_url=main_url)
}



## MAIN

# tmp dir for pdf
pdfdir <- tempfile("tmp", tmpdir = "./")
if (!dir.exists(pdfdir)) {
  dir.create(pdfdir)  
}
pdfdir <- normalizePath(pdfdir)

# create name for combined pdf ex. tmp111111_combined.pdf
combined_name <- str_c(argv$session,"__",file_path_sans_ext(basename(argv$input_file)),"_ucsc.pdf")

## create pdf with table
export_table_as_pdf(argv$input_file, pdfdir, add_annotations = TRUE)

# load bed file
bed <- read_data(argv$input_file, pdfdir = paste0(pdfdir,"/"), as.numeric(argv$zoom))

## init_params <- init(input$login, input$password, str_trim(input$session), input$db)  
init_params <- init_local(argv$user, argv$session, argv$db)

n <- nrow(bed)

for (i in 1:n) {
  # Increment the progress bar, and update the detail text.
  correct_url <- paste0(init_params$main_url, as.character(bed[i, 1]))
  download_pdf(init_params$session, correct_url, as.character(bed[i, 2]))
  Sys.sleep(2)
}


# combine all files
pdffiles <- sort(list.files(pdfdir, pattern = "pdf", full.names = T))
qpdf::pdf_combine(input = pdffiles, output = combined_name)

# remove tmp* directory with pdf files
unlink(pdfdir, recursive = TRUE)


