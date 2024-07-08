

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load the packages
library("tidyverse")







# Run decounfounding analysis ----
#------------------------------------------------#
#                                                #
#             ANALYZE DECONFOUNDING              # 
#                                                #
#------------------------------------------------#

V1_run_deconfounding <- function(MB_preprocessing){
  
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
  V1_deconf_ATC3 <- readRDS(paste("RData/Results/deconfounding/V1_deconfounding_analysis_ATC3_", MB_preprocessing, ".rds", sep = "")) %>% 
    dplyr::mutate(p.value = ifelse(is.na(p.value) == T, 0, p.value))
  V1_deconf_ATC4 <- readRDS(paste("RData/Results/deconfounding/V1_deconfounding_analysis_ATC4_", MB_preprocessing, ".rds", sep = "")) %>% 
    dplyr::mutate(p.value = ifelse(is.na(p.value) == T, 0, p.value))
  V1_deconf_ATC5 <- readRDS(paste("RData/Results/deconfounding/V1_deconfounding_analysis_ATC5_", MB_preprocessing, ".rds", sep = "")) %>% 
    dplyr::mutate(p.value = ifelse(is.na(p.value) == T, 0, p.value))
  
  
  
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
    dplyr::filter(MB_feature != "PB_ratio") %>% 
    dplyr::mutate(feature_group = ifelse(MB_feature %in% c("shannon", "observed", "PB_ratio"), "diversity", "univariate"))
  
  # Stage 1 - hits with current usage 
  stage1_FDR <- combined_results %>% 
    dplyr::filter(timeGroup == "Current_user") %>% 
    dplyr::group_by(feature_group) %>% 
    dplyr::mutate(FDR = p.adjust(p.value, method = "BH")) %>% 
    dplyr::ungroup()
  
  stage1_hits <- stage1_FDR %>% 
    dplyr::filter(FDR <= 0.1) %>% 
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
    dplyr::select(-one_of("stage1_hit"))
  
  saveRDS(naive_hits, paste("RData/Interim/Naive_univariateResults_V1_", MB_preprocessing, ".rds", sep = ""))
  
  
  
  # Analyze the number of hits after confounding analysis ----
  #----------------------------------#
  V1_postHoc_helpDf <- dplyr::bind_rows(V1_deconf_ATC3, V1_deconf_ATC4, V1_deconf_ATC5) %>% 
    dplyr::mutate(posthoc_sign = ifelse(p.value <= 0.05, 1, 0)) %>% 
    dplyr::mutate(posthoc_sign = ifelse(is.na(posthoc_sign) & substr(confounding_factor, 1, 3) == "J01", 1, posthoc_sign)) %>% 
    dplyr::select(-one_of("p.value")) %>% 
    tidyr::spread(group, posthoc_sign) %>% 
    dplyr::mutate(help1 = ifelse((A == 1 & B == 0) | (A == 1 & B == 1), 1, 0), 
                  help2 = ifelse(A == 0 & B == 0, 1, 0), 
                  help3 = ifelse(A == 0 & B == 1, 1, 0),
                  help4 = ifelse(A == 1 & B == 1, 1, 0))
  
  V1_postHoc_assessment <- V1_postHoc_helpDf %>% 
    dplyr::group_by(drug, MB_feature, timeGroup) %>% 
    dplyr::summarise(n_confDecon = sum(help1 == 1, na.rm = T), 
                     n_ambigDecon = sum(help2 == 1, na.rm = T),
                     n_conf = sum(help3 == 1, na.rm = T), 
                     confounders = paste(unique(confounding_factor[help3 == 1]), collapse = ", "), 
                     n_both = sum(help4 == 1, na.rm = T),
                     n_total = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(confounding_status = case_when(n_confDecon == n_total ~ "Confidently deconfounded",
                                                 n_conf > 0 ~ "Confounded",
                                                 n_ambigDecon > 0 ~ "Ambiguously deconfounded",
                                                 TRUE ~ "Something else"))
  
  saveRDS(V1_postHoc_assessment, paste("RData/Results/deconfounding/V1_postHoc_assessment_", MB_preprocessing, ".rds", sep = ""))
}

V1_run_deconfounding(MB_preprocessing = "CLR")
V1_run_deconfounding(MB_preprocessing = "PA")







# V2 decounfounding  ----
#------------------------------------------------#
#                                                #
#           V2 ANALYSIS DECONFOUNDING            # 
#                                                #
#------------------------------------------------#

V2_run_deconfounding <- function(MB_preprocessing){
  
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
    dplyr::filter(MB_feature != "PB_ratio") %>% 
    dplyr::mutate(feature_group = ifelse(MB_feature %in% c("shannon", "observed", "PB_ratio"), "diversity", "univariate"))
  
  # Stage 1 - hits with current usage 
  stage1_FDR <- combined_results %>% 
    dplyr::filter(timeGroup == "Current_user") %>% 
    dplyr::group_by(feature_group) %>% 
    dplyr::mutate(FDR = p.adjust(p.value, method = "BH")) %>% 
    dplyr::ungroup()
  
  stage1_hits <- stage1_FDR %>% 
    dplyr::filter(FDR <= 0.1) %>% 
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
    dplyr::select(-one_of("stage1_hit"))
  
  saveRDS(naive_hits, paste("RData/Interim/Naive_univariateResults_V2_", MB_preprocessing, ".rds", sep = ""))
  
  
  
  # Analyze the number of hits after confounding analysis ----
  #----------------------------------#
  V2_postHoc_helpDf <- dplyr::bind_rows(V2_deconf_ATC3, V2_deconf_ATC4, V2_deconf_ATC5) %>% 
    dplyr::mutate(posthoc_sign = ifelse(p.value <= 0.05, 1, 0)) %>% 
    dplyr::mutate(posthoc_sign = ifelse(is.na(posthoc_sign) & substr(confounding_factor, 1, 3) == "J01", 1, posthoc_sign)) %>% 
    dplyr::select(-one_of("p.value")) %>% 
    tidyr::spread(group, posthoc_sign) %>% 
    dplyr::mutate(help1 = ifelse((A == 1 & B == 0) | (A == 1 & B == 1), 1, 0), 
                  help2 = ifelse(A == 0 & B == 0, 1, 0), 
                  help3 = ifelse(A == 0 & B == 1, 1, 0),
                  help4 = ifelse(A == 1 & B == 1, 1, 0)) 
  
  V2_postHoc_assessment <- V2_postHoc_helpDf %>% 
    dplyr::filter(complete.cases(.)) %>% 
    dplyr::group_by(drug, MB_feature, timeGroup) %>% 
    dplyr::summarise(n_confDecon = sum(help1 == 1), 
                     n_ambigDecon = sum(help2 == 1),
                     n_conf = sum(help3 == 1), 
                     confounders = paste(unique(confounding_factor[help3 == 1]), collapse = ", "), 
                     n_both = sum(help4 == 1),
                     n_total = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(confounding_status = case_when(n_confDecon == n_total ~ "Confidently deconfounded",
                                                 n_conf > 0 ~ "Confounded",
                                                 n_ambigDecon > 0 ~ "Ambiguously deconfounded",
                                                 TRUE ~ "Something else"))
  
  saveRDS(V2_postHoc_assessment, paste("RData/Results/deconfounding/V2_postHoc_assessment_", MB_preprocessing, ".rds", sep = ""))
}

V2_run_deconfounding(MB_preprocessing = "CLR")
V2_run_deconfounding(MB_preprocessing = "PA")











# Disease deconfounding ----
#------------------------------------------------#
#                                                #
#             SUMMARIZE THE RESULTS              # 
#                                                #
#------------------------------------------------#

# Read the result files
#----------------------------------#

# Deconfounding analysis results 
Disease_deconf <- readRDS("RData/Results/deconfounding/Disease_deconfounding_analysis_CLR.rds")

# Used drugs
currentlyUsed_ATC3_considered <- readRDS("RData/Interim/Medications_currentlyUsed_n20_ATC3.rds")
currentlyUsed_ATC4_considered <- readRDS("RData/Interim/Medications_currentlyUsed_n20_ATC4.rds")
currentlyUsed_ATC5_considered <- readRDS("RData/Interim/Medications_currentlyUsed_n20_ATC5.rds")

drugs <- c(currentlyUsed_ATC5_considered, currentlyUsed_ATC4_considered, currentlyUsed_ATC3_considered)




# Analyze the number of hits after confounding analysis ----
#----------------------------------#

# Classify signals to confounder groups based on post-hoc analysis 
Disease_postHoc_helpDf <- Disease_deconf %>% 
  dplyr::mutate(posthoc_sign = ifelse(p.value <= 0.05, 1, 0)) %>% 
  dplyr::select(-one_of("p.value")) %>% 
  tidyr::spread(group, posthoc_sign) %>% 
  dplyr::mutate(help1 = ifelse((A == 1 & B == 0) | (A == 1 & B == 1), 1, 0), 
                help2 = ifelse(A == 0 & B == 0, 1, 0), 
                help3 = ifelse(A == 0 & B == 1, 1, 0),
                help4 = ifelse(A == 1 & B == 1, 1, 0)) 

# Base
Disease_postHoc_assessment_base <- Disease_postHoc_helpDf %>% 
  dplyr::filter(str_detect(confounding_factor, "Cumulative_") == F) %>% 
  dplyr::filter(!(confounding_factor %in% drugs)) %>%
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::group_by(disease, MB_feature) %>% 
  dplyr::summarise(n_confDecon = sum(help1 == 1), 
                   n_ambigDecon = sum(help2 == 1),
                   n_conf = sum(help3 == 1), 
                   confounders = paste(unique(confounding_factor[help3 == 1]), collapse = ", "), 
                   n_both = sum(help4 == 1),
                   n_total = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(confounding_status = case_when(n_confDecon == n_total ~ "Confidently deconfounded",
                                               n_conf > 0 ~ "Confounded",
                                               n_ambigDecon > 0 ~ "Ambiguously deconfounded",
                                               TRUE ~ "Something else"))

# Base + drugs
Disease_postHoc_assessment_drugs <- Disease_postHoc_helpDf %>% 
  dplyr::filter(str_detect(confounding_factor, "Cumulative_") == F) %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::group_by(disease, MB_feature) %>% 
  dplyr::summarise(n_confDecon = sum(help1 == 1), 
                   n_ambigDecon = sum(help2 == 1),
                   n_conf = sum(help3 == 1), 
                   confounders = paste(unique(confounding_factor[help3 == 1]), collapse = ", "), 
                   n_both = sum(help4 == 1),
                   n_total = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(confounding_status = case_when(n_confDecon == n_total ~ "Confidently deconfounded",
                                               n_conf > 0 ~ "Confounded",
                                               n_ambigDecon > 0 ~ "Ambiguously deconfounded",
                                               TRUE ~ "Something else"))


# Base + drugs + cumulative
Disease_postHoc_assessment_cumulative <- Disease_postHoc_helpDf %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::group_by(disease, MB_feature) %>% 
  dplyr::summarise(n_confDecon = sum(help1 == 1), 
                   n_ambigDecon = sum(help2 == 1),
                   n_conf = sum(help3 == 1), 
                   confounders = paste(unique(confounding_factor[help3 == 1]), collapse = ", "), 
                   n_both = sum(help4 == 1),
                   n_total = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(confounding_status = case_when(n_confDecon == n_total ~ "Confidently deconfounded",
                                               n_conf > 0 ~ "Confounded",
                                               n_ambigDecon > 0 ~ "Ambiguously deconfounded",
                                               TRUE ~ "Something else"))


saveRDS(Disease_postHoc_assessment_base, paste("RData/Results/deconfounding/Disease_postHoc_assessment_", Disease_preprocessing, "_base.rds", sep = ""))
saveRDS(Disease_postHoc_assessment_drugs, paste("RData/Results/deconfounding/Disease_postHoc_assessment_", Disease_preprocessing, "_drugs.rds", sep = ""))
saveRDS(Disease_postHoc_assessment_cumulative, paste("RData/Results/deconfounding/Disease_postHoc_assessment_", Disease_preprocessing, "_cumulative.rds", sep = ""))




