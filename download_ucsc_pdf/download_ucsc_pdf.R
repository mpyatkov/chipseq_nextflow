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
library(readr)
library(argparser)

ParseArguments <- function() {
  p <- arg_parser('Download UCSC browser PDFs by coordinates (requires login)')
  p <- add_argument(p,'--session', help="UCSC session name")
  p <- add_argument(p,'--user', help="UCSC username")
  p <- add_argument(p,'--password', help="UCSC password")
  p <- add_argument(p,'--db', default="mm9", help="UCSC database, e.g., mm9, mm10, hg19")
  p <- add_argument(p,"--zoom", default = "3" , help = "zoom factor for each region")
  p <- add_argument(p, "--input_file", help = "input BED (no header) or XLS/XLSX (with header)")
  argv <- parse_args(p)

  # Fail fast if missing required params
  required <- c("session","user","password","input_file")
  missing <- required[map_lgl(required, ~ is.null(argv[[.x]]) || is.na(argv[[.x]]) || argv[[.x]]=="")]
  if (length(missing) > 0) {
    stop(sprintf("Missing required argument(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  argv
}

argv <- ParseArguments()

download_pdf <- function(sess, url, outname) {
  pdflink <- sess %>%
    session_jump_to(url) %>%
    read_html() %>%
    html_nodes("#pdfLink") %>%
    html_attr("href") %>%
    str_replace(., "..", "https://genome.ucsc.edu/")

  final_pdf_link <- sess %>%
    session_jump_to(pdflink) %>%
    read_html() %>%
    html_node(xpath = "/html/body/div[5]/ul[1]/li[1]/a") %>%
    html_attr("href") %>%
    str_replace(., "..", "https://genome.ucsc.edu")

  tryCatch({
    download.file(final_pdf_link, destfile = outname, mode = "wb", quiet = TRUE)
  }, error = function(err) {
    message(sprintf("ERROR: Can't process pdf: %s -> %s", final_pdf_link, outname))
  })
}

calculate_zoom_factor <- function(start, end, zoom){
  frag_len <- floor((end - start)/2)
  mid_point <- start + frag_len
  new_start <- mid_point - zoom*frag_len
  new_end <- mid_point + zoom*frag_len
  list(start = floor(new_start), end = floor(new_end))
}

export_table_as_pdf <- function(file_path, outdir, title_text = "", add_annotations = TRUE){
  OUTPUT_HEIGHT <- 8.5; OUTPUT_WIDTH <- 16; CHUNK_SIZE <- 30
  ext <- tools::file_ext(file_path)
  file_data <- if (ext %in% c("xlsx","xls")) read_excel(file_path, col_names = TRUE) else read_tsv(file_path, col_names = FALSE)

  if (!add_annotations){
    file_data <- file_data %>% select(chr = 1, start = 2, end = 3)
    OUTPUT_HEIGHT <- 16; OUTPUT_WIDTH <- 8.5; CHUNK_SIZE <- 40
  }

  n <- nrow(file_data)
  r <- rep(1:ceiling(n/CHUNK_SIZE), each = CHUNK_SIZE)[1:n]
  dlist <- split(file_data, r)

  map(dlist, function(chunk_table){
    tab <- ggtexttable(chunk_table, rows = NULL,
                       theme = ttheme(base_size = 8, padding = unit(c(15, 3), "mm"))) %>%
      tab_add_title(text = title_text, face = "bold", size = 8, padding = unit(1, "line"))
    cowplot::plot_grid(tab)
  }) %>%
    marrangeGrob(nrow = 1, ncol = 1) %>%
    ggsave(str_c(outdir, "/00000_coordinates.pdf"), plot = ., width = OUTPUT_WIDTH, height = OUTPUT_HEIGHT)
}

read_data <- function(file_path, pdfdir, zoom){
  ext <- tools::file_ext(file_path)
  file_data <- if (ext %in% c("xlsx","xls")) {
    read_excel(file_path, col_names = TRUE) %>% select(chr = 1, start = 2, end = 3)
  } else {
    read_tsv(file_path, col_names = FALSE) %>% select(chr = 1, start = 2, end = 3)
  }

  file_data <- pmap_dfr(file_data, function(chr, start, end){
    z <- calculate_zoom_factor(start, end, zoom)
    tibble(chr = chr, start = z$start, end = z$end)
  }) %>%
    mutate(id = row_number()) %>%
    rowwise() %>%
    mutate(coords = paste0(chr,":",start,"-",end),
           outname = paste0(pdfdir,
                            paste(str_pad(id, width = 3, pad = "0"), chr, start, end,
                                  sep = "_", collapse = ""), ".pdf")) %>%
    ungroup() %>%
    select(coords, outname)
  file_data
}

init <- function(login, password, session_name, db) {
  sessionUrl <- paste0("https://genome.ucsc.edu/s/", login, "/", session_name)
  main_url <- "https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A26664573%2D26693533&hgsid=966765435_KMLSgA0F83EB9IJZmt7zQ5uPhnOr"
  main_url <- str_replace(main_url, "mm9", db)
  main_url <- str_extract(main_url, "(?<=^)(.+)(?<=position=)")

  loginUrl <- "https://genome.ucsc.edu/cgi-bin/hgLogin?hgLogin.do.displayLoginPage=1"
  sess <- session(loginUrl)

  params_list <- list(hgLogin_userName = login, hgLogin_password = password)
  form <- html_form(sess)[[1]] %>% html_form_set(!!!params_list)
  sess <- session_submit(sess, form = form, submit = "hgLogin.do.displayLogin")

  # Fail fast: attempt to reach the session page; if login page detected, stop.
  sess <- sess %>% session_jump_to(URLencode(sessionUrl))
  page <- read_html(sess)
  if (length(html_nodes(page, "form[action*='hgLogin']")) > 0) {
    stop("Login failed: UCSC credentials rejected or session inaccessible.", call. = FALSE)
  }
  list(session = sess, main_url = main_url)
}

# MAIN
pdfdir <- tempfile("tmp", tmpdir = "./"); if (!dir.exists(pdfdir)) dir.create(pdfdir); pdfdir <- normalizePath(pdfdir)
combined_name <- str_c(argv$session, "__", file_path_sans_ext(basename(argv$input_file)), "_ucsc.pdf")

title_text <- str_glue("Session: {argv$session}\nDB: {argv$db}\nZoom: {argv$zoom}x")
export_table_as_pdf(argv$input_file, pdfdir, title_text = title_text, add_annotations = TRUE)

bed <- read_data(argv$input_file, pdfdir = paste0(pdfdir,"/"), as.numeric(argv$zoom))
init_params <- init(argv$user, argv$password, str_trim(argv$session), argv$db)

n <- nrow(bed)
for (i in seq_len(n)) {
  correct_url <- paste0(init_params$main_url, as.character(bed[i, 1]))
  download_pdf(init_params$session, correct_url, as.character(bed[i, 2]))
  Sys.sleep(3)
}

pdffiles <- sort(list.files(pdfdir, pattern = "pdf", full.names = TRUE))
qpdf::pdf_combine(input = pdffiles, output = combined_name)

unlink(pdfdir, recursive = TRUE)
