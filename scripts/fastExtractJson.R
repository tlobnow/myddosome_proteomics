pacman::p_load(dplyr, tidyr, stringr, fs, jsonlite, purrr, utils, data.table)

source(file = ifelse(exists("https://raw.githubusercontent.com/tlobnow/myddosome_proteomics/main/scripts/functions.R"), 
                     yes =  "https://raw.githubusercontent.com/tlobnow/myddosome_proteomics/main/scripts/functions.R",
                     no  =  "~/Documents/Github/myddosome_proteomics/scripts/functions.R"))

# run_extraction <- function(LOC, MAIN = NULL, SUMMARY_FOLDER = NULL, ADD_2_EXISTING_DF = F, EXISTING_DF = NULL) {
#   
#   if (is.null(MAIN)) {
#     MAIN <- paste0("/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/", LOC, "/")
#   }
#   
#   if (is.null(SUMMARY_FOLDER)) {
#     SUMMARY_FOLDER <- "/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/SUMMARIES/"
#   }
#   
#   if (isTRUE(ADD_2_EXISTING_DF) & is.null(EXISTING_DF)) {
#     print("If you wish to add the extracted data to an existing DF, please provide the path for EXISTING_DF.")
#   }
#   
#   summaryWithRecycles = paste0(SUMMARY_FOLDER, "RECYCLES/", LOC, "_summaryWithRecycles.csv")
#   summary             = paste0(SUMMARY_FOLDER, LOC, ".csv")
#   if (file.exists(summaryWithRecycles))  {file.remove(summaryWithRecycles)}
#   if (file.exists(summary))              {file.remove(summary)}
#   
#   COUNTER <- 0
#   FILES   <- list.files(MAIN, pattern = "_x")
#   LEN     <- as.numeric(length(FILES))
#   
#   # if you wish to extract a single file, unhash below and provide file name
#   #FILES = "ARL8B_MOUSE_x1_UN93B_MOUSE_x1_TLR7_MOUSE_x1"
#   #FILE = "ARL8B_MOUSE_x1_UN93B_MOUSE_x1_TLR7_MOUSE_x1"
#  
#   # Loop through each FILE in FILES
#   for (FILE in FILES) {
#     
#     tryCatch({
#       
#       FOLDER = paste0(MAIN, FILE, "/")
#       dir_create(FOLDER, "CSV")
#       csv_file <- paste0(FOLDER, "CSV/", FILE, "_withRecycles.csv")
#       csv_file2 <- paste0(FOLDER, "CSV/", FILE, ".csv")
#       
#       if (file.exists(csv_file))  {file.remove(csv_file)}
#       if (file.exists(csv_file2)) {file.remove(csv_file2)}
#       
#       json_folder <- ifelse(dir.exists(file.path(FOLDER, "JSON")), yes = file.path(FOLDER, "JSON"), file.path(FOLDER))
#       maxJSON <- list.files(json_folder, pattern = ".json")
#       
#       if (is_empty(maxJSON)) {
#         next
#       }
#       
#       json_files <- dir_ls(json_folder, regexp = "\\.json$", recurse = TRUE)
#       
#       for (file in json_files) {
#         json <- readLines(file)
#         json <- str_replace_all(json, "Infinity", "9999")
#         writeLines(json, file)
#       }
#       
#       if (length(maxJSON) > 0) {
#         for (i in 1:length(maxJSON)) {
#           JSON <- file.path(json_folder, maxJSON[i])
#           OUT <- paste(FOLDER, "CSV", FILE, sep = "/")
#           jsonExtract(JSON = JSON, OUT = OUT, FILE = FILE)
#         }
#       }
#       
#       CSV_FILES <- c(csv_file, csv_file2)
#       
#       for (CSV_FILE in CSV_FILES) {
#         JSON_EXTRACT <- data.table::fread(CSV_FILE, header = FALSE) %>%
#           dplyr::mutate(ORIGIN = FILE) %>%
#           dplyr::rename(FILE = V1, MODEL = V2, RECYCLE = V3, TOL = V4, pLDDT = V5, pTM = V6, piTM = V7, iScore = V8, iRes = V9, iCnt = V10, FILE_MODEL = V11, NUM_CLUSTERS = V12, N_MONOMERS = V13)
#         
#         JE <- JSON_EXTRACT %>%
#           dplyr::mutate(FILE_RECYCLE = paste0(FILE_MODEL, "_RECYCLE_", RECYCLE), RANK = NA) %>%
#           dplyr::distinct(FILE_RECYCLE, .keep_all = TRUE) %>%
#           dplyr::group_by(FILE) %>%
#           dplyr::mutate(RANK = frank(desc(iScore), ties.method = "min"))
#         
#         data.table::fwrite(JE, CSV_FILE, row.names = FALSE)
#         data.table::fwrite(JE, summaryWithRecycles, row.names = FALSE, append = T)
#       }
#       
#       # Increment counter after each successful round
#       COUNTER <- COUNTER + 1
#       
#       # Print current count out of total
#       cat("\r", paste(COUNTER, "/", LEN, " done"), fill = F)
#       
#     }, error = function(e) {
#       # Handle error
#       cat("An error occurred for FILE:", FILE, "\n")
#       cat("Error message:", e$message, "\n")
#     })
#     
#   }
#   
#   bigboy      <- fread(summaryWithRecycles, fill = T) %>% distinct(RECYCLE, .keep_all = T)
#   bigboy$DATE <- ifelse(test = str_detect(bigboy$RECYCLE, pattern = "multimer"), 
#                         yes = sapply(strsplit(bigboy$RECYCLE, "_"), function(x) x[6]),
#                         no = sapply(strsplit(bigboy$RECYCLE, "_"), function(x) x[5]))
#   bigboy$DATE <- lubridate::as_date(bigboy$DATE)
#   data.table::fwrite(bigboy, summaryWithRecycles, row.names = FALSE, append = F)
#   
#   if (isTRUE(ADD_2_EXISTING_DF) &  !is.null(EXISTING_DF)) {
#     add2Summary(NEW_DF = summaryWithRecycles, EXISTING_DF = EXISTING_DF)
#   }
#   
#   smolboy <- bigboy %>% filter(!str_detect(RECYCLE, "_recycled_")) %>% distinct(RECYCLE, .keep_all = T)
#   data.table::fwrite(smolboy, summary, row.names = FALSE, append = F)
#   
#   if (isTRUE(ADD_2_EXISTING_DF) &  !is.null(EXISTING_DF)) {
#     add2Summary(NEW_DF = summary, EXISTING_DF = EXISTING_DF)
#   }
# }

# LOC = "BDLD_KAUR_ALN"
# LOC = "DHF_ALL"
# LOC = "BDLD"
# LOC = "CHIMY_T6BM"
# LOC = "CONTROLS"
# LOC = "MYD88"
# LOC = "TEMP"
# LOC = "BDLD"
LOC = "ADLD"
LOC = "EDLD"

# MAIN    = file.path("/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS", LOC, "/")
# MAIN    = file.path("/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS_2/", LOC, "/")
# MAIN    = file.path("/Users/u_lobnow/Documents/Github/", LOC, "/")
# MAIN    = file.path("/Users/u_lobnow/Documents/Github/transferGit/", LOC, "/")
MAIN    = file.path("/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/", LOC, "/")

run_extraction(LOC = LOC, MAIN = MAIN)


