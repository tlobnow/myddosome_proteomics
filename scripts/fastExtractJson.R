pacman::p_load(dplyr, tidyr, stringr, fs, jsonlite, purrr, utils, data.table)

source(file = ifelse(exists("https://raw.githubusercontent.com/tlobnow/myddosome_proteomics/main/scripts/functions.R"), 
                     yes =  "https://raw.githubusercontent.com/tlobnow/myddosome_proteomics/main/scripts/functions.R",
                     no  =  "~/Documents/Github/myddosome_proteomics/scripts/functions.R"))

# LOC = "BDLD_KAUR_ALN"
# LOC = "DHF_ALL"
# LOC = "BDLD"
# LOC = "CHIMY_T6BM"
# LOC = "CONTROLS"
# LOC = "MYD88"
LOC = "TEMP"
# MAIN    = file.path("/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS", LOC, "/")
# MAIN    = file.path("/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS_2/", LOC, "/")
# MAIN    = file.path("/Users/u_lobnow/Documents/Github/", LOC, "/")
MAIN    = file.path("/Users/u_lobnow/Documents/Github/transferGit/", LOC, "/")

SUMMARY_FOLDER = "/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/SUMMARIES/"
summaryWithRecycles = paste0(SUMMARY_FOLDER, "RECYCLES/", LOC, "_summaryWithRecycles.csv")
summary             = paste0(SUMMARY_FOLDER, LOC, ".csv")
if (file.exists(summaryWithRecycles))  {file.remove(summaryWithRecycles)}
if (file.exists(summary))              {file.remove(summary)}

COUNTER <- 0
FILES   <- list.files(MAIN, pattern = "_x")
LEN     <- as.numeric(length(FILES))

# if you wish to extract a single file, unhash below and provide file name
#FILES = "ARL8B_MOUSE_x1_UN93B_MOUSE_x1_TLR7_MOUSE_x1"
#FILE = "ARL8B_MOUSE_x1_UN93B_MOUSE_x1_TLR7_MOUSE_x1"

jsonExtract <- function(JSON, OUT, FILE) {
  json       <- fromJSON(JSON)
  MODEL      <- strsplit(json$order[[1]], "_", fixed = TRUE)[[1]][2]
  TOL        <- json$tol_values %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "TOL") 
  pLDDT      <- json$plddts %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "pLDDT") 
  pTM        <- json$ptms %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "pTM") 
  piTM       <- json$pitms %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "piTM") 
  iScore     <- json$`interface score` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iScore") 
  iRes       <- json$`interfacial residue number` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iRes") 
  iCnt       <- json$`interficial contact number` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iCnt") 
  FILE_MODEL <- paste(FILE, "MODEL", MODEL, sep = "_")
  NUM_CLUSTERS <- json$clusters[[iScore$RECYCLE[1]]]$num_clusters
  N_MONOMERS <- length(json$chains)
  
  EXTRACT <- cbind(FILE, MODEL, TOL, pLDDT, pTM, piTM, iScore, iRes, iCnt, FILE_MODEL, NUM_CLUSTERS, N_MONOMERS)
  EXTRACT <- EXTRACT[, !duplicated(colnames(EXTRACT))]
  
  write.table(EXTRACT, file = paste0(OUT,"_withRecycles.csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
  
  EXTRACT_noRecycle <- EXTRACT %>% filter(!str_detect(RECYCLE, "_recycled_"))
  write.table(EXTRACT_noRecycle, file = paste0(OUT,".csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
}

# Loop through each FILE in FILES
for (FILE in FILES) {
  
  tryCatch({
    
    FOLDER = paste0(MAIN, FILE, "/")
    dir_create(FOLDER, "CSV")
    csv_file <- paste0(FOLDER, "CSV/", FILE, "_withRecycles.csv")
    csv_file2 <- paste0(FOLDER, "CSV/", FILE, ".csv")
    
    if (file.exists(csv_file))  {file.remove(csv_file)}
    if (file.exists(csv_file2)) {file.remove(csv_file2)}
    
    json_folder <- ifelse(dir.exists(file.path(FOLDER, "JSON")), yes = file.path(FOLDER, "JSON"), file.path(FOLDER))
    maxJSON <- list.files(json_folder, pattern = ".json")
    
    if (is_empty(maxJSON)) {
      next
    }
    
    json_files <- dir_ls(json_folder, regexp = "\\.json$", recurse = TRUE)
    
    for (file in json_files) {
      json <- readLines(file)
      json <- str_replace_all(json, "Infinity", "9999")
      writeLines(json, file)
    }
    
    if (length(maxJSON) > 0) {
      for (i in 1:length(maxJSON)) {
        JSON <- file.path(json_folder, maxJSON[i])
        OUT <- paste(FOLDER, "CSV", FILE, sep = "/")
        jsonExtract(JSON = JSON, OUT = OUT, FILE = FILE)
      }
    }
    
    CSV_FILES <- c(csv_file, csv_file2)
    
    for (CSV_FILE in CSV_FILES) {
      JSON_EXTRACT <- data.table::fread(CSV_FILE, header = FALSE) %>%
        dplyr::mutate(ORIGIN = FILE) %>%
        dplyr::rename(FILE = V1, MODEL = V2, RECYCLE = V3, TOL = V4, pLDDT = V5, pTM = V6, piTM = V7, iScore = V8, iRes = V9, iCnt = V10, FILE_MODEL = V11, NUM_CLUSTERS = V12, N_MONOMERS = V13)
      
      JE <- JSON_EXTRACT %>%
        dplyr::mutate(FILE_RECYCLE = paste0(FILE_MODEL, "_RECYCLE_", RECYCLE), RANK = NA) %>%
        dplyr::distinct(FILE_RECYCLE, .keep_all = TRUE) %>%
        dplyr::group_by(FILE) %>%
        dplyr::mutate(RANK = frank(desc(iScore), ties.method = "min"))
      
      data.table::fwrite(JE, CSV_FILE, row.names = FALSE)
      data.table::fwrite(JE, summaryWithRecycles, row.names = FALSE, append = T)
    }
    
    # Increment counter after each successful round
    COUNTER <- COUNTER + 1
    
    # Print current count out of total
    cat("\r", paste(COUNTER, "/", LEN, " done"), fill = F)
    
  }, error = function(e) {
    # Handle error
    cat("An error occurred for FILE:", FILE, "\n")
    cat("Error message:", e$message, "\n")
  })
  
}

bigboy      <- fread(summaryWithRecycles, fill = T) %>% distinct(RECYCLE, .keep_all = T)
bigboy$DATE <- ifelse(test = str_detect(bigboy$RECYCLE, pattern = "multimer"), 
                      yes = sapply(strsplit(bigboy$RECYCLE, "_"), function(x) x[6]),
                      no = sapply(strsplit(bigboy$RECYCLE, "_"), function(x) x[5]))
bigboy$DATE <- lubridate::as_date(bigboy$DATE)
data.table::fwrite(bigboy, summaryWithRecycles, row.names = FALSE, append = F)

add2Summary(NEW_DF = summaryWithRecycles,
            EXISTING_DF = "/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/SUMMARIES/RECYCLES/MYD88_summaryWithRecycles.csv")

smolboy <- bigboy %>% filter(!str_detect(RECYCLE, "_recycled_")) %>% distinct(RECYCLE, .keep_all = T)
data.table::fwrite(smolboy, summary, row.names = FALSE, append = F)

add2Summary(NEW_DF = summary,
            EXISTING_DF = "/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/SUMMARIES/MYD88.csv")

rm(list = ls())
