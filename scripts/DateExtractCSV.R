# Required packages: tidyr, data.table, lubridate, stringr

# Source functions
source(file = ifelse(exists("https://raw.githubusercontent.com/tlobnow/myddosome_proteomics/main/scripts/functions.R"), 
                     yes =  "https://raw.githubusercontent.com/tlobnow/myddosome_proteomics/main/scripts/functions.R",
                     no  =  "~/Documents/Github/myddosome_proteomics/scripts/functions.R"))

# Load the required packages
pacman::p_load(tidyr,data.table,lubridate,stringr,dplyr)
# Set the parent directory where the subfolders are located
# parent_directory <- "/Volumes/TAYLOR-LAB/Finn/RESULTS/BDLD"
# parent_directory <- "/Volumes/TAYLOR-LAB/Finn/RESULTS/BDLD_rep2"
# parent_directory <- "/Volumes/TAYLOR-LAB/Finn/RESULTS/BDLD_rep3"
# parent_directory <- "/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS_2/MYD88/"
# parent_directory <- "~/Documents/Github/temp_storage/BDLD_x10_rep2/"
# parent_directory <- "~/Documents/Github/temp_storage/BDLD_x10_rep3/"
# parent_directory <- "~/Documents/Github/temp_storage/CHIMY_T6BM_rep2/"
# parent_directory <- "/Volumes/TAYLOR-LAB/Finn/RESULTS/CHIMY_T6BM_rep3/"
# parent_directory = "/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/BDLD_KAUR_ALN/"
# parent_directory = "/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS_2/MYD88/"
parent_directory = "/Users/u_lobnow/Documents/Github/transferGit/TEMP/"

# Get a list of all subfolders in the parent directory
subfolders <- list.dirs(parent_directory, recursive = FALSE)

# Initialize an empty list to store the data
data_list <- list()
# Loop through each subfolder
for (subfolder in subfolders) {
  # Get a list of CSV files in the subfolder
  csv_files <- list.files(path = paste0(subfolder, "/CSV/"), pattern = "\\.csv", full.names = TRUE)
  
  # Read and bind the CSV files into a single data.table
  if (length(csv_files) > 0) {
    data <- rbindlist(lapply(csv_files, fread), fill = T)
    data_list[[subfolder]] <- data
  }
}

# Combine all the data.tables into a single data.table
combined_data <- rbindlist(data_list, fill = T)

# Extract the date from the 'RECYCLE' column by splitting the string and selecting the sixth part:
# considering data structures: "model_1_multimer_v3_p1_230517_742373" or "model_4_multimer_v3_p1_230628_004801_recycled_00"
combined_data$DATE <- sapply(strsplit(combined_data$RECYCLE, "_"), function(x) x[6])

# Convert the 'DATE' column to Date format
combined_data$DATE <- lubridate::as_date(combined_data$DATE)

# Retain unique rows based on 'RECYCLE' column
combined_data <- combined_data %>% distinct(RECYCLE, .keep_all = TRUE) #%>% select(-c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13))

# Write the combined data to a summary file with recycles
write.csv(x = combined_data, file = paste0(parent_directory, "/", basename(parent_directory), "_summary_withRecycles.csv"))

# Filter out rows with '_recycled_' in the 'RECYCLE' column
combined_data <- combined_data %>% filter(!str_detect(RECYCLE, pattern = "_recycled_"))

# Write the filtered data to a summary file without recycles
write.csv(x = combined_data, file = paste0(parent_directory, "/", basename(parent_directory), "_summary.csv"))

# If you need to add to a larger df
# add2Summary(NEW_DF = paste0(parent_directory, "/", basename(parent_directory), "_summary_withRecycles.csv"), 
#             EXISTING_DF = "/path/to/existing/DF_summary_withRecycles.csv")
# 
# add2Summary(NEW_DF = paste0(parent_directory, "/", basename(parent_directory), "_summary.csv"), 
#             EXISTING_DF = "/path/to/existing/DF_summary.csv")
