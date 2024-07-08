



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
# args = c("5", "ATC4", "PA")

cat(args, sep = "\n")


# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 


# Declare model number 
model_nr = as.numeric(args[1])


# Declare ATC level
ATC_level = args[2]


# Declare data transformation
MB_preprocessing = args[3]







# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load packages
library("readr")
library("dplyr")
library("effsize")
library("ppcor")


# Drug usage data
#----------------------------#
current_usage_data_ATC3 <- readRDS("RData/Interim/Data_medications_currentlyUsed_ATC3.rds") 
current_usage_data_ATC4 <- readRDS("RData/Interim/Data_medications_currentlyUsed_ATC4.rds") 
current_usage_data_ATC5 <- readRDS("RData/Interim/Data_medications_currentlyUsed_ATC5.rds") 

drugs_considered_ATC3 <- readRDS("RData/Interim/Medications_currentlyUsed_n20_ATC3.rds")
drugs_considered_ATC4 <- readRDS("RData/Interim/Medications_currentlyUsed_n20_ATC4.rds")
drugs_considered_ATC5 <- readRDS("RData/Interim/Medications_currentlyUsed_n20_ATC5.rds")

last_used_data_ATC3 <- readRDS("RData/Interim/Drug_last_used_ATC3.rds")
last_used_data_ATC4 <- readRDS("RData/Interim/Drug_last_used_ATC4.rds")
last_used_data_ATC5 <- readRDS("RData/Interim/Drug_last_used_ATC5.rds")

cumulative_usage_ATC3 <- readRDS("RData/Interim/Cumulative_usage_ATC3.rds") %>% dplyr::rename("ATC_used" = "ATC3")
cumulative_usage_ATC4 <- readRDS("RData/Interim/Cumulative_usage_ATC4.rds") %>% dplyr::rename("ATC_used" = "ATC4")
cumulative_usage_ATC5 <- readRDS("RData/Interim/Cumulative_usage_ATC5.rds") %>% dplyr::rename("ATC_used" = "ATC5")


# Phenotype data
#----------------------------#
drugs_old <-  readRDS("RData/Interim/Medication_factors_analyzed.rds")

# Remove old drug variables from the phenotype data and add new
phenotype_data <- readRDS("RData/Interim/Data_master.rds") %>% 
  dplyr::select(-drugs_old) %>% 
  dplyr::left_join(current_usage_data_ATC3[ ,c("skood", drugs_considered_ATC3)], by = "skood") %>% 
  dplyr::left_join(current_usage_data_ATC4[ ,c("skood", drugs_considered_ATC4)], by = "skood") %>%
  dplyr::left_join(current_usage_data_ATC5[ ,c("skood", drugs_considered_ATC5)], by = "skood") %>% 
  dplyr::filter(is.na(BMI) == F)



# Microbiome data 
#----------------------------#

# Abundance data
if (MB_preprocessing == "PA"){
  countData_raw <- readRDS("RData/Interim/mOTUs_PA_prevalence10p.rds")
} else if (MB_preprocessing == "CLR"){
  countData_raw <- readRDS("RData/Interim/mOTUs_CLR_prevalence10p.rds")
}else if (MB_preprocessing == "ABS"){
  countData_raw = readRDS("RData/Interim/mOTUs2_cellcounts_filtered_prevalence10p.rds") %>% 
    tibble::column_to_rownames(var = "skood")
  
  pseudocount = min(countData_raw[countData_raw != 0])/2
  countData_raw[countData_raw == 0] <- pseudocount
  
  countData_raw = log10(countData_raw) %>% 
    tibble::rownames_to_column(var = "skood")
}

PB_ratio <- readRDS("RData/Interim/PB_ratio_mOTUs.rds") %>% 
  dplyr::select(skood, PB_ratio)

diversity_data <- readRDS("RData/Interim/mOTUs_alphaDiversity_unfiltered.rds") 



# Define antibiotics users
#----------------------------#
q_AB_users <- current_usage_data_ATC3 %>% 
  dplyr::select(skood, starts_with("J01")) %>% 
  tidyr::gather(key, value, -skood) %>% 
  dplyr::group_by(skood) %>% 
  dplyr::summarise(J01 = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(J01 > 0) %>% 
  dplyr::pull(skood)


# Define factor groups
#----------------------------#
factors_diseases <- c(readRDS("RData/Interim/Disease_factors_analyzed.rds"), "gumDiseaseDiagnosed", "seasonableAllergy", "chickenpox")
factors_procedures <- c(readRDS("RData/Interim/Procedure_factors_analyzed.rds"), "cecumRemoved", "tonsilsRemoved")
factors_dietary <- c(readRDS("RData/Interim/Dietary_factors_analyzed.rds"), "eatingHabit_name")
factors_lifestyle <- c("has_smoked", "alcohol_has_used", "doesPhysicalExercise")
factors_stool <- c("usualStoolType_category", "frequencyGutEmpting_name") 
factors_drugs <- c(drugs_considered_ATC3, drugs_considered_ATC4, drugs_considered_ATC5)




# Define which drugs to use in analysis
#----------------------------#
if (ATC_level == "ATC3"){
  drugs_considered = drugs_considered_ATC3
  last_used_data = last_used_data_ATC3
} else if (ATC_level == "ATC4"){
  drugs_considered = drugs_considered_ATC4
  last_used_data = last_used_data_ATC4
  
} else if (ATC_level == "ATC5"){
  drugs_considered = drugs_considered_ATC5
  last_used_data = last_used_data_ATC5
  
}








# Deconfounding analysis for V1 ----
#------------------------------------------------#
#                                                #
#            DECONFOUNDING ANALYSIS              # 
#                                                #
#------------------------------------------------#

if (model_nr == 1){
  
  # Read the result files
  #----------------------------------#
  V1_res_ATC3 <- readRDS(paste("RData/Results/drugUsage/V1_drugUsage_analysis_ATC3_", MB_preprocessing, ".rds", sep = ""))
  V1_res_ATC4 <- readRDS(paste("RData/Results/drugUsage/V1_drugUsage_analysis_ATC4_", MB_preprocessing, ".rds", sep = ""))
  V1_res_ATC5 <- readRDS(paste("RData/Results/drugUsage/V1_drugUsage_analysis_ATC5_", MB_preprocessing, ".rds", sep = ""))
  
  V1_res_ATC3_sens <- readRDS(paste("RData/Results/drugUsage/V1_drugUsage_analysis_ATC3_", MB_preprocessing, "_sensitivityAB.rds", sep = ""))
  V1_res_ATC4_sens <- readRDS(paste("RData/Results/drugUsage/V1_drugUsage_analysis_ATC4_", MB_preprocessing, "_sensitivityAB.rds", sep = ""))
  V1_res_ATC5_sens <- readRDS(paste("RData/Results/drugUsage/V1_drugUsage_analysis_ATC5_", MB_preprocessing, "_sensitivityAB.rds", sep = ""))
  
  
  # Deconfounding analysis results 
  #------------------------------------------#
  V1_deconf_ATC3 <- readRDS(paste("RData/Results/deconfounding/V1_deconfounding_analysis_ATC3_", MB_preprocessing, ".rds", sep = ""))
  V1_deconf_ATC4 <- readRDS(paste("RData/Results/deconfounding/V1_deconfounding_analysis_ATC4_", MB_preprocessing, ".rds", sep = ""))
  V1_deconf_ATC5 <- readRDS(paste("RData/Results/deconfounding/V1_deconfounding_analysis_ATC5_", MB_preprocessing, ".rds", sep = ""))
  
  
  # Collect naiive results ----
  #----------------------------------#
  
  # Combine results, filter significant hits
  combined_results <- dplyr::bind_rows(V1_res_ATC3, V1_res_ATC4, V1_res_ATC5) %>% 
    dplyr::mutate(group = "raw") %>% 
    dplyr::bind_rows(V1_res_ATC3_sens, V1_res_ATC4_sens, V1_res_ATC5_sens) %>% 
    dplyr::mutate(group = ifelse(is.na(group) == TRUE, "sensitivity", group)) %>% 
    dplyr::filter((substr(drug, 1, 3) == "J01" & group == "raw") | 
                    (substr(drug, 1, 3) != "J01" & group == "sensitivity")) %>% 
    dplyr::select(-one_of("group")) %>% 
    dplyr::filter(timeGroup != "after_90d") %>% 
    dplyr::mutate(feature_group = ifelse(MB_feature %in% c("shannon", "observed", "PB_ratio"), "diversity", "univariate"))
  
  # Stage 1 - hits with current usage 
  stage1_FDR <- combined_results %>% 
    dplyr::filter(timeGroup == "Current_user") %>% 
    dplyr::group_by(feature_group) %>% 
    dplyr::mutate(FDR = p.adjust(p.value, method = "BH")) %>% 
    dplyr::ungroup()
  
  stage1_hits <- stage1_FDR %>% 
    dplyr::filter(p.value <= 0.05) %>% 
    dplyr::mutate(stage1_hit = "Y")
  
  
  # Sage 2 - hits with historical usage
  stage2_FDR <- combined_results %>% 
    dplyr::filter(timeGroup != "Current_user") %>% 
    # Focus on only these hits that had a signal in stage1 - immediate usage
    dplyr::left_join(stage1_hits[ ,c("drug", "MB_feature", "stage1_hit")], by = c("drug", "MB_feature")) %>% 
    dplyr::filter(stage1_hit == "Y") %>% 
    # Adjust p-value
    dplyr::group_by(feature_group) %>% 
    dplyr::mutate(FDR = p.adjust(p.value, method = "BH")) %>% 
    dplyr::ungroup()  
  
  # Combine hits
  naive_hits <- stage1_FDR %>% 
    dplyr::bind_rows(stage2_FDR) %>% 
    dplyr::select(-one_of("stage1_hit")) %>% 
    dplyr::filter(p.value <= 0.05) %>% 
    dplyr::mutate(ATC_group = case_when(nchar(as.character(drug)) == 4 ~ "ATC3", 
                                        nchar(as.character(drug)) == 5 ~ "ATC4", 
                                        nchar(as.character(drug)) == 7 ~ "ATC5", 
                                        TRUE ~ "Something else")) %>% 
    dplyr::filter(ATC_group == ATC_level)
  

  
  # Prepare raw datasets for analysis
  # -------------------------------#
  recentUsage_data <- readRDS(paste("RData/Interim/Drug_last_used_", ATC_level, ".rds", sep = ""))
  colnames(recentUsage_data) <- c("skood", "ATC_used", "timeGroup")
  
  recentUsage_analysis_data <- phenotype_data %>% 
    dplyr::select(skood) %>% 
    dplyr::full_join(recentUsage_data, by = "skood") %>% 
    dplyr::select(skood, ATC_used) %>% 
    tidyr::expand(skood, ATC_used) %>% 
    dplyr::filter(complete.cases(.)) %>% 
    dplyr::left_join(recentUsage_data, by = c("skood", "ATC_used")) %>% 
    dplyr::mutate(timeGroup = ifelse(is.na(timeGroup), "non_user", as.character(timeGroup))) %>% 
    dplyr::select(skood, ATC_used, timeGroup) %>% 
    dplyr::mutate(help = 1) %>% 
    tidyr::spread(timeGroup, help, fill = 0) %>%
    dplyr::mutate(Current_user = ifelse(Current_user == 1, 1, 0),
                  after_90d = ifelse(Current_user == 0 & (`<1y` + `1y-2y` + `2y-3y` + `3y-4y` + `4y-5y` > 0), 1, 0),
                  after_1y = ifelse(Current_user == 0 & `<1y` == 0 & (`1y-2y` + `2y-3y` + `3y-4y` + `4y-5y` > 0), 1, 0),
                  after_2y = ifelse(Current_user == 0 & `<1y` == 0 & `1y-2y` == 0 & (`2y-3y` + `3y-4y` + `4y-5y` > 0), 1, 0),
                  after_3y = ifelse(Current_user == 0 & `<1y` == 0 & `1y-2y` == 0 & `2y-3y` == 0 & (`3y-4y` + `4y-5y` > 0), 1, 0),
                  after_4y = ifelse(Current_user == 0 & `<1y` == 0 & `1y-2y` == 0 & `2y-3y` == 0 & `3y-4y` == 0 & (`4y-5y` > 0), 1, 0)) %>% 
    dplyr::select(skood, ATC_used, non_user, Current_user, after_90d, after_1y, after_2y, after_3y, after_4y) %>% 
    tidyr::gather(key, value, -c("skood", "ATC_used", "non_user")) %>% 
    dplyr::filter(non_user == 1 | value == 1) %>% 
    dplyr::mutate(group = ifelse(value == 1 & non_user == 0, 1, 0)) %>% 
    dplyr::select(skood, ATC_used, key, group)
  
  
  
  # Run the deconfounding models
  # -------------------------------#
  
  confounding_results_df = data.frame()
  
  # Run over all drug hits
  for (i in 1:nrow(naive_hits)){
    
    run_drug = naive_hits[i, ] %>% dplyr::pull(drug) %>% as.character()
    run_bug = naive_hits[i, ] %>% dplyr::pull(MB_feature) %>% as.character()
    run_timeGroup = naive_hits[i, ] %>% dplyr::pull(timeGroup) %>% as.character()
    
    
    # Filter original data
    if (substr(run_drug, 1, 3) != "J01"){
      run_data = recentUsage_analysis_data %>% 
        dplyr::filter(!(skood %in% q_AB_users)) %>% 
        dplyr::filter(ATC_used == run_drug & key == run_timeGroup) %>% 
        dplyr::left_join(phenotype_data, by = "skood") %>% 
        dplyr::left_join(countData_raw, by = "skood") %>% 
        dplyr::left_join(PB_ratio, by = "skood") %>% 
        dplyr::left_join(diversity_data, by = "skood") 
      
      run_data = run_data %>%
        dplyr::mutate(bug = run_data %>% dplyr::pull(run_bug))
      
    } else {
      run_data = recentUsage_analysis_data  %>% 
        dplyr::filter(ATC_used == run_drug & key == run_timeGroup) %>% 
        dplyr::left_join(phenotype_data, by = "skood") %>% 
        dplyr::left_join(countData_raw, by = "skood") %>%
        dplyr::left_join(PB_ratio, by = "skood") %>% 
        dplyr::left_join(diversity_data, by = "skood") 
      
      run_data = run_data %>%
        dplyr::mutate(bug = run_data %>% dplyr::pull(run_bug)) 
    }
    
    # Define potential confounder list
    confounder_list = data.frame(drug = run_drug, 
                                 confounder = c(factors_drugs, factors_diseases, factors_procedures, factors_dietary, factors_lifestyle, factors_stool)) %>% 
      dplyr::mutate(match = case_when(substr(drug, 1, 4) == substr(confounder, 1, 4) ~ 1, 
                                      substr(drug, 1, 5) == substr(confounder, 1, 5) ~ 1, 
                                      substr(drug, 1, 7) == substr(confounder, 1, 7) ~ 1, 
                                      TRUE ~ 0)) %>% 
      dplyr::filter(match == 0) %>% 
      dplyr::pull(confounder)
    
    
    # Run over all confounders
    for (j in confounder_list){
      
      run_data_clean <- run_data %>% 
        dplyr::mutate(confounder = run_data %>% dplyr::pull(j)) %>% 
        dplyr::filter(is.na(confounder) == F)
      
      # Define the models
      f0_1 = paste("bug ~  BMI + gender + Age_at_MBsample + group", sep = "")
      f0_2 = paste("bug ~  BMI + gender + Age_at_MBsample + `", j, "`", sep = "")
      f_full = paste("bug ~  BMI + gender + Age_at_MBsample + group + `", j, "`", sep = "")
      
      # Run models
      if (MB_preprocessing == "PA" & !(run_bug %in% c("observed", "shannon", "PB_ratio"))){
        m0_1 = glm(f0_1, data = run_data_clean, family = "binomial")
        m0_2 = glm(f0_2, data = run_data_clean, family = "binomial")
        m_full = glm(f_full, data = run_data_clean, family = "binomial")
      } else{
        m0_1 = lm(f0_1, data = run_data_clean)
        m0_2 = lm(f0_2, data = run_data_clean)
        m_full = lm(f_full, data = run_data_clean)
      }
      
      
      # Compare the models
      if (MB_preprocessing == "PA" & !(run_bug %in% c("observed", "shannon", "PB_ratio"))){
        anova_ai = anova(m_full, m0_2, test = "Chisq")
        anova_bi = anova(m_full, m0_1, test = "Chisq")
        p_anova_ai = anova_ai$`Pr(>Chi)`[2]
        p_anova_bi = anova_bi$`Pr(>Chi)`[2]
      } else{
        anova_ai = anova(m_full, m0_2)
        anova_bi = anova(m_full, m0_1)
        p_anova_ai = anova_ai$`Pr(>F)`[2]
        p_anova_bi = anova_bi$`Pr(>F)`[2]
      }
      
      
      # Output results
      run_output = data.frame(drug = run_drug, 
                              confounding_factor = j,
                              MB_feature = run_bug,
                              timeGroup = run_timeGroup, 
                              group = c("A", "B"),
                              p.value = c(p_anova_ai, p_anova_bi), 
                              transformation = MB_preprocessing)
      
      confounding_results_df = dplyr::bind_rows(confounding_results_df, run_output)
      
    } # End of j cycle
  } # End of i cycle
  
  # Save the results
  filename = paste("RData/Results/deconfounding/V1_deconfounding_analysis_", ATC_level, "_", MB_preprocessing, ".rds", sep = "")
  saveRDS(confounding_results_df, filename)
} # End of RUN MODEL  














# Deconfounding analysis for V2 ----
#------------------------------------------------#
#                                                #
#            DECONFOUNDING ANALYSIS              # 
#                                                #
#------------------------------------------------#

if (model_nr == 2){
  
  # Read the result files
  #----------------------------------#
  V2_res_ATC3 <- readRDS(paste("RData/Results/drugUsage/V2_drugUsage_analysis_ATC3_", MB_preprocessing, ".rds", sep = ""))
  V2_res_ATC4 <- readRDS(paste("RData/Results/drugUsage/V2_drugUsage_analysis_ATC4_", MB_preprocessing, ".rds", sep = ""))
  V2_res_ATC5 <- readRDS(paste("RData/Results/drugUsage/V2_drugUsage_analysis_ATC5_", MB_preprocessing, ".rds", sep = ""))
  
  V2_res_ATC3_sens <- readRDS(paste("RData/Results/drugUsage/V2_drugUsage_analysis_ATC3_", MB_preprocessing, "_sensitivityAB.rds", sep = ""))
  V2_res_ATC4_sens <- readRDS(paste("RData/Results/drugUsage/V2_drugUsage_analysis_ATC4_", MB_preprocessing, "_sensitivityAB.rds", sep = ""))
  V2_res_ATC5_sens <- readRDS(paste("RData/Results/drugUsage/V2_drugUsage_analysis_ATC5_", MB_preprocessing, "_sensitivityAB.rds", sep = ""))
  
  
  # Deconfounding analysis results 
  #------------------------------------------#
  V2_deconf_ATC3 <- readRDS(paste("RData/Results/deconfounding/V2_deconfounding_analysis_ATC3_", MB_preprocessing, ".rds", sep = ""))
  V2_deconf_ATC4 <- readRDS(paste("RData/Results/deconfounding/V2_deconfounding_analysis_ATC4_", MB_preprocessing, ".rds", sep = ""))
  V2_deconf_ATC5 <- readRDS(paste("RData/Results/deconfounding/V2_deconfounding_analysis_ATC5_", MB_preprocessing, ".rds", sep = ""))
  
  
  # Collect naiive results ----
  #----------------------------------#
  
  # Combine results, filter significant hits
  combined_results <- dplyr::bind_rows(V2_res_ATC3, V2_res_ATC4, V2_res_ATC5) %>% 
    dplyr::mutate(group = "raw") %>% 
    dplyr::bind_rows(V2_res_ATC3_sens, V2_res_ATC4_sens, V2_res_ATC5_sens) %>% 
    dplyr::mutate(group = ifelse(is.na(group) == TRUE, "sensitivity", group)) %>% 
    dplyr::filter((substr(drug, 1, 3) == "J01" & group == "raw") | 
                    (substr(drug, 1, 3) != "J01" & group == "sensitivity")) %>% 
    dplyr::select(-one_of("group")) %>% 
    dplyr::filter(timeGroup != "after_90d") %>% 
    dplyr::mutate(feature_group = ifelse(MB_feature %in% c("shannon", "observed", "PB_ratio"), "diversity", "univariate"))
  
  # Stage 1 - hits with current usage 
  stage1_FDR <- combined_results %>% 
    dplyr::filter(timeGroup == "Current_user") %>% 
    dplyr::group_by(feature_group) %>% 
    dplyr::mutate(FDR = p.adjust(p.value, method = "BH")) %>% 
    dplyr::ungroup()
  
  stage1_hits <- stage1_FDR %>% 
    dplyr::filter(p.value <= 0.05) %>% 
    dplyr::mutate(stage1_hit = "Y")
  
  
  # Sage 2 - hits with historical usage
  stage2_FDR <- combined_results %>% 
    dplyr::filter(timeGroup != "Current_user") %>% 
    # Focus on only these hits that had a signal in stage1 - immediate usage
    dplyr::left_join(stage1_hits[ ,c("drug", "MB_feature", "stage1_hit")], by = c("drug", "MB_feature")) %>% 
    dplyr::filter(stage1_hit == "Y") %>% 
    # Adjust p-value
    dplyr::group_by(feature_group) %>% 
    dplyr::mutate(FDR = p.adjust(p.value, method = "BH")) %>% 
    dplyr::ungroup()  
  
  
  # Combine hits
  naive_hits <- stage1_FDR %>% 
    dplyr::bind_rows(stage2_FDR) %>% 
    dplyr::select(-one_of("stage1_hit")) %>% 
    dplyr::filter(p.value <= 0.05) %>% 
    dplyr::mutate(ATC_group = case_when(nchar(as.character(drug)) == 4 ~ "ATC3", 
                                        nchar(as.character(drug)) == 5 ~ "ATC4", 
                                        nchar(as.character(drug)) == 7 ~ "ATC5", 
                                        TRUE ~ "Something else")) %>% 
    dplyr::filter(ATC_group == ATC_level)
  

  # Prepare raw datasets for analysis
  # -------------------------------#
  recentUsage_data <- readRDS(paste("RData/Interim/Drug_last_used_", ATC_level, ".rds", sep = ""))
  colnames(recentUsage_data) <- c("skood", "ATC_used", "timeGroup")
  
  recentUsage_analysis_data <- phenotype_data %>% 
    dplyr::select(skood) %>% 
    dplyr::full_join(recentUsage_data, by = "skood") %>% 
    dplyr::select(skood, ATC_used) %>% 
    tidyr::expand(skood, ATC_used) %>% 
    dplyr::filter(complete.cases(.)) %>% 
    dplyr::left_join(recentUsage_data, by = c("skood", "ATC_used")) %>% 
    dplyr::mutate(timeGroup = ifelse(is.na(timeGroup), "non_user", as.character(timeGroup))) %>% 
    dplyr::select(skood, ATC_used, timeGroup) %>% 
    dplyr::mutate(help = 1) %>% 
    tidyr::spread(timeGroup, help, fill = 0)  %>% 
    tidyr::gather(key, value, -c("skood", "ATC_used", "non_user")) %>% 
    dplyr::filter(non_user == 1 | value == 1) %>% 
    dplyr::mutate(group = ifelse(value == 1 & non_user == 0, 1, 0)) %>% 
    dplyr::select(skood, ATC_used, key, group)

  
  
  # Run the deconfounding models
  # -------------------------------#
  
  confounding_results_df = data.frame()
  
  # Run over all drug hits
  for (i in 1:nrow(naive_hits)){
    
    run_drug = naive_hits[i, ] %>% dplyr::pull(drug) %>% as.character()
    run_bug = naive_hits[i, ] %>% dplyr::pull(MB_feature) %>% as.character()
    run_timeGroup = naive_hits[i, ] %>% dplyr::pull(timeGroup) %>% as.character()
    
    
    # Filter original data
    if (substr(run_drug, 1, 3) != "J01"){
      run_data = recentUsage_analysis_data %>% 
        dplyr::filter(!(skood %in% q_AB_users)) %>% 
        dplyr::filter(ATC_used == run_drug & key == run_timeGroup) %>% 
        dplyr::left_join(phenotype_data, by = "skood") %>% 
        dplyr::left_join(countData_raw, by = "skood") %>% 
        dplyr::left_join(PB_ratio, by = "skood") %>% 
        dplyr::left_join(diversity_data, by = "skood") 
      
      run_data = run_data %>%
        dplyr::mutate(bug = run_data %>% dplyr::pull(run_bug))

    } else {
      run_data = recentUsage_analysis_data  %>% 
        dplyr::filter(ATC_used == run_drug & key == run_timeGroup) %>% 
        dplyr::left_join(phenotype_data, by = "skood") %>% 
        dplyr::left_join(countData_raw, by = "skood") %>%
        dplyr::left_join(PB_ratio, by = "skood") %>% 
        dplyr::left_join(diversity_data, by = "skood") 
      
      run_data = run_data %>%
        dplyr::mutate(bug = run_data %>% dplyr::pull(run_bug)) 
    }
    
    # Define potential confounder list
    confounder_list = data.frame(drug = run_drug, 
                                 confounder = c(factors_drugs, factors_diseases, factors_procedures, factors_dietary, factors_lifestyle, factors_stool)) %>% 
      dplyr::mutate(match = case_when(substr(drug, 1, 4) == substr(confounder, 1, 4) ~ 1, 
                                      substr(drug, 1, 5) == substr(confounder, 1, 5) ~ 1, 
                                      substr(drug, 1, 7) == substr(confounder, 1, 7) ~ 1, 
                                      TRUE ~ 0)) %>% 
      dplyr::filter(match == 0) %>% 
      dplyr::pull(confounder)
    
    
    # Run over all confounders
    for (j in confounder_list){
      
      run_data_clean <- run_data %>% 
        dplyr::mutate(confounder = run_data %>% dplyr::pull(j)) %>% 
        dplyr::filter(is.na(confounder) == F)
      
      # Define the models
      f0_1 = paste("bug ~  BMI + gender + Age_at_MBsample + group", sep = "")
      f0_2 = paste("bug ~  BMI + gender + Age_at_MBsample + `", j, "`", sep = "")
      f_full = paste("bug ~  BMI + gender + Age_at_MBsample + group + `", j, "`", sep = "")
      
      # Run models
      if (MB_preprocessing == "PA" & !(run_bug %in% c("observed", "shannon", "PB_ratio"))){
        m0_1 = glm(f0_1, data = run_data_clean, family = "binomial")
        m0_2 = glm(f0_2, data = run_data_clean, family = "binomial")
        m_full = glm(f_full, data = run_data_clean, family = "binomial")
      } else{
        m0_1 = lm(f0_1, data = run_data_clean)
        m0_2 = lm(f0_2, data = run_data_clean)
        m_full = lm(f_full, data = run_data_clean)
      }
      

      
      # Compare the models
      if (MB_preprocessing == "PA" & !(run_bug %in% c("observed", "shannon", "PB_ratio"))){
        anova_ai = anova(m_full, m0_2, test = "Chisq")
        anova_bi = anova(m_full, m0_1, test = "Chisq")
        p_anova_ai = anova_ai$`Pr(>Chi)`[2]
        p_anova_bi = anova_bi$`Pr(>Chi)`[2]
      } else{
        anova_ai = anova(m_full, m0_2)
        anova_bi = anova(m_full, m0_1)
        p_anova_ai = anova_ai$`Pr(>F)`[2]
        p_anova_bi = anova_bi$`Pr(>F)`[2]
      }

      
      # Output results
      run_output = data.frame(drug = run_drug, 
                              confounding_factor = j,
                              MB_feature = run_bug,
                              timeGroup = run_timeGroup, 
                              group = c("A", "B"),
                              p.value = c(p_anova_ai, p_anova_bi), 
                              transformation = MB_preprocessing)
      
      confounding_results_df = dplyr::bind_rows(confounding_results_df, run_output)
      
    } # End of j cycle
  } # End of i cycle
  
  # Save the results
  filename = paste("RData/Results/deconfounding/V2_deconfounding_analysis_", ATC_level, "_", MB_preprocessing, ".rds", sep = "")
  saveRDS(confounding_results_df, filename)
} # End of RUN MODEL  



  
  
  
  


# Deconfounding analysis for diseases ----
#------------------------------------------------#
#                                                #
#            DECONFOUNDING ANALYSIS              # 
#                                                #
#------------------------------------------------#

if (model_nr == 0){
  
  # Naive hits
  disease_naive_hits <- readRDS(paste("RData/Interim/Naive_disease_associations_", MB_preprocessing, ".rds", sep = "")) %>% 
    dplyr::filter(p.value <= 0.05)
  
  
  # Prepare raw datasets for analysis
  # -------------------------------#
  cumulative_usage <- dplyr::bind_rows(cumulative_usage_ATC3, cumulative_usage_ATC4, cumulative_usage_ATC5) %>% 
    dplyr::select(skood, ATC_used, n_prescription) %>% 
    dplyr::filter(ATC_used %in% c(drugs_considered_ATC3, drugs_considered_ATC4, drugs_considered_ATC5)) %>% 
    dplyr::mutate(ATC_used = paste("Cumulative_", ATC_used, sep = ""))
  
  Cumulative_Usage_data <- phenotype_data %>% 
    dplyr::select(skood) %>% 
    dplyr::full_join(cumulative_usage, by = c("skood")) %>% 
    tidyr::expand(skood, ATC_used) %>%
    dplyr::filter(complete.cases(.)) %>% 
    dplyr::left_join(cumulative_usage, by = c("skood", "ATC_used")) %>% 
    dplyr::mutate(n_prescription = ifelse(is.na(n_prescription) == TRUE, 0, n_prescription)) %>% 
    tidyr::spread(ATC_used, n_prescription)

  disease_df = phenotype_data %>% 
    dplyr::left_join(Cumulative_Usage_data, by = "skood") %>% 
    dplyr::left_join(countData_raw, by = "skood") %>% 
    dplyr::left_join(PB_ratio, by = "skood") %>% 
    dplyr::left_join(diversity_data, by = "skood")

  
  # Run the deconfounding models
  # -------------------------------#
  
  confounding_results_df = data.frame()
  
  # Run over all drug hit
  for (i in 1:nrow(disease_naive_hits)){
    
    run_disease = disease_naive_hits[i, ] %>% dplyr::pull(disease) %>% as.character()
    run_bug = disease_naive_hits[i, ] %>% dplyr::pull(MB_feature) %>% as.character()

    
    # Prepare data
    run_data = disease_df %>% 
      dplyr::mutate(bug = disease_df %>% dplyr::pull(run_bug), 
                    disease = disease_df %>% dplyr::pull(run_disease))
    
    
    # Define potential confounder list
    confounder_list = c(factors_drugs, 
                        setdiff(factors_diseases, run_disease),
                        paste("Cumulative_", factors_drugs, sep = ""), 
                        factors_procedures, factors_dietary, 
                        factors_lifestyle, factors_stool)

    # Run over all confounders
    for (k in confounder_list){
      
      run_data_clean <- run_data %>% 
        dplyr::mutate(confounder = run_data %>% dplyr::pull(k)) %>% 
        dplyr::filter(is.na(confounder) == F)

      # Run models
      if (MB_preprocessing == "PA" & !(run_bug %in% c("observed", "shannon", "PB_ratio"))){
        m1_1 = glm(paste("bug ~ BMI + gender + Age_at_MBsample + disease", sep = ""), data = run_data_clean, family = "binomial")
        m1_2 = glm(paste("bug ~ BMI + gender + Age_at_MBsample + `", k, "`", sep = ""), data = run_data_clean, family = "binomial")
        m1_full = glm(paste("bug ~ BMI + gender + Age_at_MBsample + disease + `", k, "`", sep = ""), data = run_data_clean, family = "binomial")
        
        m1_anova_ai = anova(m1_full, m1_2, test = "Chisq")
        m1_anova_bi = anova(m1_full, m1_1, test = "Chisq")
        
        ai_p = m1_anova_ai$`Pr(>Chi)`[2]
        bi_p = m1_anova_bi$`Pr(>Chi)`[2]
      } else{
        m1_1 = lm(paste("bug ~ BMI + gender + Age_at_MBsample + disease", sep = ""), data = run_data_clean)
        m1_2 = lm(paste("bug ~ BMI + gender + Age_at_MBsample + `", k, "`", sep = ""), data = run_data_clean)
        m1_full = lm(paste("bug ~ BMI + gender + Age_at_MBsample + disease + `", k, "`", sep = ""), data = run_data_clean)
        
        m1_anova_ai = anova(m1_full, m1_2)
        m1_anova_bi = anova(m1_full, m1_1)
        
        ai_p = m1_anova_ai$`Pr(>F)`[2]
        bi_p = m1_anova_bi$`Pr(>F)`[2]
      }

      # Output results
      run_output = data.frame(disease = run_disease, 
                              confounding_factor = k,
                              MB_feature = run_bug,
                              group = c("A", "B"),
                              p.value = c(ai_p, bi_p))
      
      confounding_results_df = dplyr::bind_rows(confounding_results_df, run_output)
    } # End of k cycle
  } # End of i cycle
  
  # Save the results
  filename = paste("RData/Results/deconfounding/Disease_deconfounding_analysis_",  MB_preprocessing, ".rds", sep = "")
  saveRDS(confounding_results_df, filename)
} # End of RUN MODEL 

