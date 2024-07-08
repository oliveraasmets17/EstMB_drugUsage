

# Command line arguments ----
#------------------------------------------------#
#                                                #
#          SET UP COMMAND LINE ARGUMENTS         # 
#                                                #
#------------------------------------------------#

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)              
# Expected input:
# Model_nr  

# NB! For testing - remove
# args = c("JACCARD")

cat(args, sep = "\n")


# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 


# Declare distance matrix name 
distance_algo = args[1]






# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load packages
library("dplyr")
library("vegan")
library("stringr")

# Read distance matrix
if (distance_algo == "CLR"){
  distance_matrix_raw <- readRDS("RData/Interim/CLR_dist_mOTUs_filtered.rds")
} else if (distance_algo == "BRAY"){
  distance_matrix_raw <- readRDS("RData/Interim/Bray_dist_mOTUs_filtered.rds")
} else if (distance_algo == "JACCARD"){
  distance_matrix_raw <- readRDS("RData/Interim/Jaccard_dist_mOTUs_filtered.rds")
}


# Clean the phenotype data
disease_factors <- readRDS("RData/Interim/Disease_factors_analyzed.rds")

phenotype_data <- readRDS("RData/Interim/Data_master.rds") %>% 
  dplyr::select(skood, BMI, Age_at_MBsample, gender, all_of(disease_factors))


# Drug usage data
currentlyUsed_ATC3 <- readRDS("RData/Interim/Data_medications_currentlyUsed_ATC3.rds") 
currentlyUsed_ATC4 <- readRDS("RData/Interim/Data_medications_currentlyUsed_ATC4.rds") 
currentlyUsed_ATC5 <- readRDS("RData/Interim/Data_medications_currentlyUsed_ATC5.rds") 

drugs_ATC3 <- readRDS("RData/Interim/Medications_currentlyUsed_n20_ATC3.rds") 
drugs_ATC4 <- readRDS("RData/Interim/Medications_currentlyUsed_n20_ATC4.rds") 
drugs_ATC5 <- readRDS("RData/Interim/Medications_currentlyUsed_n20_ATC5.rds") 

cumulative_drugUsage_ATC3 <- readRDS("RData/Interim/Cumulative_usage_ATC3.rds") %>% dplyr::select(skood, ATC3, n_prescription) %>% tidyr::spread(ATC3, n_prescription, fill = 0)
cumulative_drugUsage_ATC4 <- readRDS("RData/Interim/Cumulative_usage_ATC4.rds") %>% dplyr::select(skood, ATC4, n_prescription) %>% tidyr::spread(ATC4, n_prescription, fill = 0)
cumulative_drugUsage_ATC5 <- readRDS("RData/Interim/Cumulative_usage_ATC5.rds") %>% dplyr::select(skood, ATC5, n_prescription) %>% tidyr::spread(ATC5, n_prescription, fill = 0)



# Merge data for long-term usage
Cumulative_drugUsage <- phenotype_data %>% 
  dplyr::select(skood) %>% 
  dplyr::full_join(cumulative_drugUsage_ATC3[ , c("skood", drugs_ATC3)], by = "skood") %>% 
  dplyr::full_join(cumulative_drugUsage_ATC4[ , c("skood", drugs_ATC4)], by = "skood") %>% 
  dplyr::full_join(cumulative_drugUsage_ATC5[ , c("skood", drugs_ATC5)], by = "skood") %>% 
  tidyr::gather(key, value, -skood) %>% 
  dplyr::filter(str_detect(key, "<NA>") == F) %>% 
  dplyr::mutate(value = ifelse(is.na(value), 0, value)) %>% 
  tidyr::spread(key, value)

colnames(Cumulative_drugUsage) <- c("skood", paste("cumulative_", colnames(Cumulative_drugUsage[ ,2:ncol(Cumulative_drugUsage)]), sep = ""))



# Create binary dataset for long-term usage
Cumulative_drugUsage_binary <- phenotype_data %>% 
  dplyr::select(skood) %>% 
  dplyr::full_join(cumulative_drugUsage_ATC3[ , c("skood", drugs_ATC3)], by = "skood") %>% 
  dplyr::full_join(cumulative_drugUsage_ATC4[ , c("skood", drugs_ATC4)], by = "skood") %>% 
  dplyr::full_join(cumulative_drugUsage_ATC5[ , c("skood", drugs_ATC5)], by = "skood") %>% 
  tidyr::gather(key, value, -skood) %>% 
  dplyr::filter(str_detect(key, "<NA>") == F) %>% 
  dplyr::mutate(value = ifelse(is.na(value), 0, value), 
                value = ifelse(value > 0, 1, 0)) %>% 
  tidyr::spread(key, value)

colnames(Cumulative_drugUsage_binary) <- c("skood", paste("binary_", colnames(Cumulative_drugUsage_binary[ ,2:ncol(Cumulative_drugUsage_binary)]), sep = ""))


# Merge the data
PERMANOVA_data <- phenotype_data %>% 
  dplyr::left_join(currentlyUsed_ATC3[ , c("skood", drugs_ATC3)], by = "skood") %>% 
  dplyr::left_join(currentlyUsed_ATC4[ , c("skood", drugs_ATC4)], by = "skood") %>% 
  dplyr::left_join(currentlyUsed_ATC5[ , c("skood", drugs_ATC5)], by = "skood") %>% 
  dplyr::left_join(Cumulative_drugUsage, by = "skood") %>% 
  dplyr::left_join(Cumulative_drugUsage_binary, by = "skood") %>% 
  dplyr::filter(skood %in% rownames(as.matrix(distance_matrix_raw)))






# Function that run PERMANOVA ----
#------------------------------------------------#
#                                                #
#       PERMANOVA FOR COMMUNITY DIFFERENCES      # 
#                                                #
#------------------------------------------------#

# Dist to matrix
distance_matrix <- distance_matrix_raw %>% 
  as.matrix()


# Define the function
run_PERMANOVA = function(factor_chosen){
  
  # Subset complete observations for the given factor
  metadata = PERMANOVA_data %>%
    dplyr::mutate(analysis_factor = PERMANOVA_data %>% dplyr::pull(factor_chosen)) %>%
    dplyr::select(skood, analysis_factor) %>%
    dplyr::filter(complete.cases(.)) %>% 
    as.data.frame()
 
  
  # Subset distance matrix - keep rows and cols without missing values
  distance_matrix_subset <- distance_matrix[metadata$skood, metadata$skood]
  
  
  # Create formulas for the analysis
  formula = as.formula("distance_matrix_subset ~ analysis_factor")

  
  # Run PERMANOVA
  set.seed(1)
  permanova <- vegan::adonis(formula, 
                             data = metadata,
                             permutations = 10000, 
                             parallel = 10)
  
  # Clean output
  permanova_output <- permanova$aov.tab %>% 
    as.data.frame() %>%
    dplyr::filter(str_detect(rownames(.), coll("analysis_factor")) == TRUE) %>%
    dplyr::mutate(analysis_factor = factor_chosen) %>%
    dplyr::select(analysis_factor, everything())
  
  return(permanova_output) 
  
}
 




# Apply the function for numerous factors ----
#------------------------------------------------#
#                                                #
#            ANALYZE MULTIPLE FACTORS            # 
#                                                #
#------------------------------------------------#

# Run PERMANOVA on all the codes chosen
PERMANOVA_summary_output <- data.frame()

run_all_PERMANOVA <- TRUE
if (run_all_PERMANOVA == TRUE){
  for (i in setdiff(colnames(PERMANOVA_data), "skood")){
    output_run = run_PERMANOVA(factor_chosen = i)
    PERMANOVA_summary_output = dplyr::bind_rows(PERMANOVA_summary_output, output_run)
  }
  
  # Save the results
  output_file_name <- paste("RData/Results/PERMANOVA_results_", distance_algo, ".rds", sep = "")
  saveRDS(PERMANOVA_summary_output, file = output_file_name)
}



