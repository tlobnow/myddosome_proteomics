#          .                                                   ####  ##        #
#       ":"                               ####              ###########        #
#     ___:____     |"\/"|               ########              #######          #
#   ,'        `.    \  /                  #####                                #
#   |  O        \___/  |                                                       #
# ~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~#

### LOAD LIBRARIES
pacman::p_load(tidyverse,data.table,jsonlite,janitor,ggrepel,UniprotR,knitr,svglite,readxl,plotly,fs,stringr)

### LOAD FUNCTIONS
source(file = ifelse(exists("https://raw.githubusercontent.com/tlobnow/master_thesis/main/scripts/functions.R"), 
                     yes =  "https://raw.githubusercontent.com/tlobnow/master_thesis/main/scripts/functions.R",
                     no  =  "~/Documents/Github/myddosome_proteomics/scripts/functions.R"))

### SET MODES
JSON_XTRCT    = T
JSON_PROCESS  = T
ANNOTATE      = T
# avoid SLURM extraction (unfortunately quite error-prone..)
SLURM_XTRCT   = F
PROCESS_SLURM = F

### DEFINE PATHS
MAIN_FOLDER        = "IP_MS_2"
FOLDER             = "MYD88"
# FOLDER             = "IRAK4"
# FOLDER             = "IRAK1"
# FOLDER             = "MYD88_signif"

FILES_LOC          = "~/Documents/Github/master_thesis/"
MAIN               = ifelse(dir.exists(paste0("/Volumes/TAYLOR-LAB/Finn/RESULTS/", MAIN_FOLDER, "/")), 
                            yes =  paste0("/Volumes/TAYLOR-LAB/Finn/RESULTS/", MAIN_FOLDER, "/"),
                            no  =  "~/Documents/Github/transferGit/")

OUT                = paste0(FILES_LOC, FOLDER)
RAW_SUMMARIES_PATH = paste0(FILES_LOC, "raw_summaries/")        # preprocessed files (raw extracted json, slurms)
SUMMARIES_PATH     = paste0(FILES_LOC, "summaries/")            # processed files
TAXA_PATH          = paste0(FILES_LOC, "taxa_lists/")
ANNOTATED_PATH     = paste0(FILES_LOC, "summaries_annotated/")  # added taxa info

### EXTRACT ALL JSON FILES CONTAINED IN PROVIDED FOLDER
if (JSON_XTRCT == T) {source(paste0(FILES_LOC, "scripts/JSON_XTRCT.R"))}

### PROCESS JSON SUMMARY FILE
if (JSON_PROCESS == T) {source(paste0(FILES_LOC, "scripts/JSON_PROCESS.R"))}

### ANNOTATE JSON SUMMARY FILE
if (ANNOTATE == T) {
  source(paste0(FILES_LOC, "scripts/ANNOTATE_1.R"))
} else source(paste0(FILES_LOC, "scripts/ANNOTATE_2.R"))





