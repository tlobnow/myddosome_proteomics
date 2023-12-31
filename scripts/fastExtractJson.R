source(file = ifelse(exists("https://raw.githubusercontent.com/tlobnow/myddosome_proteomics/main/scripts/functions.R"), 
                     yes =  "https://raw.githubusercontent.com/tlobnow/myddosome_proteomics/main/scripts/functions.R",
                     no  =  "~/Documents/Github/myddosome_proteomics/scripts/functions.R"))

# LOC = "BDLD_KAUR_ALN"
# LOC = "DHF_ALL"
# LOC = "BDLD"
# LOC = "CHIMY_T6BM"
LOC = "CONTROLS"
# LOC = "MYD88"
# LOC = "TEMP"
# LOC = "BDLD"
# LOC = "ADLD"
# LOC = "EDLD"


# MAIN    = file.path("/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS", LOC, "/")
# MAIN    = file.path("/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS_2/", LOC, "/")
# MAIN    = file.path("/Users/u_lobnow/Documents/Github/", LOC, "/")
# MAIN    = file.path("/Users/u_lobnow/Documents/Github/transferGit/", LOC, "/")
# MAIN    = file.path("/Volumes/TAYLOR-LAB/Finn/CURATED_RESULTS/", LOC, "/")

run_extraction(LOC = LOC)

plot_alphafold_results(LOC = LOC)# filter_by = "EDLD_18_Q9N9X2_GEOCY_Geodia_cydonium_x6")

plot_alphafold_results(LOC = LOC, 
                       pattern = "DD",
                       best_only = T)
