



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
# args = c("1", "ATC4", "PA")

cat(args, sep = "\n")


# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 


# Declare model number 
model_nr = as.numeric(args[1])


# Declare ATC level
ATC_level = args[2]


# Declare MB preprocessing
MB_preprocessing = args[3]







# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load packages
library("effsize")
library("ppcor")
library("tidyverse")

# Read the data
current_usage_data <- readRDS(paste("RData/Interim/Data_medications_currentlyUsed_", ATC_level, ".rds", sep = "")) 
current_drugs_considered <- readRDS(paste("RData/Interim/Medications_currentlyUsed_n20_", ATC_level, ".rds", sep = ""))

phenotype_data <- readRDS("RData/Interim/Data_master.rds") 



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

# Diversity indexes based on raw data
diversity_data <- readRDS("RData/Interim/mOTUs_alphaDiversity_unfiltered.rds")


# Define antibiotics users
q_AB_users <- current_usage_data %>% 
  dplyr::select(skood, starts_with("J01")) %>% 
  tidyr::gather(key, value, -skood) %>% 
  dplyr::group_by(skood) %>% 
  dplyr::summarise(J01 = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(J01 > 0) %>% 
  dplyr::pull(skood)


# PB ratio
PB_ratio <- readRDS("RData/Interim/PB_ratio_mOTUs.rds") %>% 
  dplyr::select(skood, PB_ratio)








# Version 1  ----
#------------------------------------------------#
#                                                #
#         ANALYZING CUMULATIVE EFFECTS           # 
#                                                #
#------------------------------------------------#


# Run the analysis
#--------------------------------------------------#
if (model_nr == 1){
  
  # Categorize subject by the recent drug usage
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
  
  
  # Run over all drugs analyzed
  V1_cumulative_analysis = data.frame()
  
  for (i in current_drugs_considered){
    for (j in c("Current_user", "after_90d", "after_1y", "after_2y", "after_3y", "after_4y")){
      
      run_df <- recentUsage_analysis_data %>% 
        dplyr::filter(ATC_used == i) %>% 
        dplyr::filter(key == j) %>% 
        dplyr::left_join(countData_raw, by = "skood") %>% 
        dplyr::left_join(diversity_data, by = "skood") %>% 
        dplyr::left_join(PB_ratio, by = "skood") %>% 
        dplyr::left_join(phenotype_data[ ,c("skood", "BMI", "gender", "Age_at_MBsample")], by = "skood")
      
      if (length(table(run_df$group)) == 2){
        
        # Run over all microbial features
        for (k in c(setdiff(colnames(countData_raw), "skood"),
                    "PB_ratio", 
                    setdiff(colnames(diversity_data), "skood"))){
          
          f0 = paste("`", k, "` ~ BMI + gender + Age_at_MBsample", sep = "")
          f1 = paste("`", k, "` ~ BMI + gender + Age_at_MBsample + group", sep = "")
          
          if (MB_preprocessing == "PA" & !(k %in% c("PB_ratio", "observed", "shannon"))){
            m0 = glm(f0, data = run_df, family = "binomial")
            m1 = glm(f1, data = run_df, family = "binomial")
            lr_test = anova(m0, m1, test = "Chisq")
            p_anova = lr_test$`Pr(>Chi)`[2]
            
            effect_OR = broom::tidy(m1) %>% 
              dplyr::filter(term == "group") %>% 
              dplyr::mutate(OR = exp(estimate)) %>% 
              dplyr::pull(OR)
            
          } else{
            m0 = lm(f0, data = run_df)
            m1 = lm(f1, data = run_df)
            lr_test = anova(m0, m1)
            p_anova = lr_test$`Pr(>F)`[2]
            
            effect_OR = as.numeric(NA)
          }
          
          run_df_clean <- run_df %>% 
            dplyr::select(all_of(k), group, gender, Age_at_MBsample, BMI) %>% 
            dplyr::filter(complete.cases(.))
          
          eff_pearson = spcor.test(x = run_df_clean %>% pull(k),
                                   y = run_df_clean %>% pull(group) %>% as.numeric(), 
                                   z = run_df_clean %>% dplyr::select(gender, Age_at_MBsample, BMI), 
                                   method = "pearson")
          
          eff_cliff = cliff.delta(d = run_df %>% dplyr::pull(k), 
                                  f = factor(run_df %>% dplyr::pull("group"), levels = c(1, 0)))
          
          run_output = data.frame(drug = i, 
                                  timeGroup = j, 
                                  MB_feature = k, 
                                  p.value = p_anova, 
                                  effect_OR, 
                                  effect_pearson = eff_pearson$estimate,
                                  effect_cliff = eff_cliff$estimate, 
                                  transformation = MB_preprocessing)
          
          V1_cumulative_analysis = dplyr::bind_rows(V1_cumulative_analysis, run_output)
        }
      }
    }
  }
  filename = paste("RData/Results/drugUsage/V1_drugUsage_analysis_", ATC_level, "_", MB_preprocessing, ".rds", sep = "")
  saveRDS(V1_cumulative_analysis, filename)
}





# Run the analysis - sensitivity analysis to AB usage
#--------------------------------------------------#
if (model_nr == 12){
  
  # Categorize subject by the recent drug usage
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
    dplyr::select(skood, ATC_used, key, group) %>% 
    dplyr::filter(!(skood %in% q_AB_users))
  
  
  # Run over all drugs analyzed
  V1_cumulative_analysis2 = data.frame()
  
  for (i in current_drugs_considered){
    for (j in c("Current_user", "after_90d", "after_1y", "after_2y", "after_3y", "after_4y")){
      
      run_df <- recentUsage_analysis_data %>% 
        dplyr::filter(ATC_used == i) %>% 
        dplyr::filter(key == j) %>% 
        dplyr::left_join(countData_raw, by = "skood") %>% 
        dplyr::left_join(diversity_data, by = "skood") %>% 
        dplyr::left_join(PB_ratio, by = "skood") %>% 
        dplyr::left_join(phenotype_data[ ,c("skood", "BMI", "gender", "Age_at_MBsample")], by = "skood")
      
      if (length(table(run_df$group)) == 2){
        # Run over all microbial features
        for (k in c(setdiff(colnames(countData_raw), "skood"),
                    "PB_ratio",
                    setdiff(colnames(diversity_data), "skood"))){
          
          f0 = paste("`", k, "` ~ BMI + gender + Age_at_MBsample", sep = "")
          f1 = paste("`", k, "` ~ BMI + gender + Age_at_MBsample + group", sep = "")
          
          if (MB_preprocessing == "PA" & !(k %in% c("PB_ratio", "observed", "shannon"))){
            m0 = glm(f0, data = run_df, family = "binomial")
            m1 = glm(f1, data = run_df, family = "binomial")
            lr_test = anova(m0, m1, test = "Chisq")
            p_anova = lr_test$`Pr(>Chi)`[2]
            
            effect_OR = broom::tidy(m1) %>% 
              dplyr::filter(term == "group") %>% 
              dplyr::mutate(OR = exp(estimate)) %>% 
              dplyr::pull(OR)
            
          } else{
            m0 = lm(f0, data = run_df)
            m1 = lm(f1, data = run_df)
            lr_test = anova(m0, m1)
            p_anova = lr_test$`Pr(>F)`[2]
            
            effect_OR = as.numeric(NA)
          }
          
          run_df_clean <- run_df %>% 
            dplyr::select(all_of(k), group, gender, Age_at_MBsample, BMI) %>% 
            dplyr::filter(complete.cases(.))
          
          eff_pearson = spcor.test(x = run_df_clean %>% pull(k),
                                   y = run_df_clean %>% pull(group) %>% as.numeric(), 
                                   z = run_df_clean %>% dplyr::select(gender, Age_at_MBsample, BMI), 
                                   method = "pearson")
          
          eff_cliff = cliff.delta(d = run_df %>% dplyr::pull(k), 
                                  f = factor(run_df %>% dplyr::pull("group"), levels = c(1, 0)))
          
          run_output = data.frame(drug = i, 
                                  timeGroup = j, 
                                  MB_feature = k, 
                                  p.value = p_anova, 
                                  effect_OR,
                                  effect_pearson = eff_pearson$estimate,
                                  effect_cliff = eff_cliff$estimate, 
                                  transformation = MB_preprocessing)
          
          V1_cumulative_analysis2 = dplyr::bind_rows(V1_cumulative_analysis2, run_output)
        }
      }
    }
  }
  filename = paste("RData/Results/drugUsage/V1_drugUsage_analysis_", ATC_level, "_", MB_preprocessing, "_sensitivityAB.rds", sep = "")
  saveRDS(V1_cumulative_analysis2, filename)
}








# Version 2 (Analyse the carry-over effects by time groups)  ----
#------------------------------------------------#
#                                                #
#         ANALYZING CUMULATIVE EFFECTS           # 
#                                                #
#------------------------------------------------#

# Run the analysis
#--------------------------------------------------#
if (model_nr == 2){
  
  # Categorize subject by the recent drug usage
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
  
  
  # Run over all drugs analyzed
  V2_cumulative_analysis = data.frame()
  
  for (i in current_drugs_considered){
    for (j in c("Current_user", "<1y", "1y-2y", "2y-3y", "3y-4y", "4y-5y")){
      
      run_df <- recentUsage_analysis_data %>% 
        dplyr::filter(ATC_used == i) %>% 
        dplyr::filter(key == j) %>% 
        dplyr::left_join(countData_raw, by = "skood") %>% 
        dplyr::left_join(diversity_data, by = "skood") %>% 
        dplyr::left_join(PB_ratio, by = "skood") %>% 
        dplyr::left_join(phenotype_data[ ,c("skood", "BMI", "gender", "Age_at_MBsample")], by = "skood")
      
      if (length(table(run_df$group)) == 2){
        
        # Run over all microbial features
        for (k in c(setdiff(colnames(countData_raw), "skood"),
                    "PB_ratio",
                    setdiff(colnames(diversity_data), "skood"))){
          
          f0 = paste("`", k, "` ~ BMI + gender + Age_at_MBsample", sep = "")
          f1 = paste("`", k, "` ~ BMI + gender + Age_at_MBsample + group", sep = "")
          
          
          if (MB_preprocessing == "PA" & !(k %in% c("PB_ratio", "observed", "shannon"))){
            m0 = glm(f0, data = run_df, family = "binomial")
            m1 = glm(f1, data = run_df, family = "binomial")
            lr_test = anova(m0, m1, test = "Chisq")
            p_anova = lr_test$`Pr(>Chi)`[2]
            
            effect_OR = broom::tidy(m1) %>% 
              dplyr::filter(term == "group") %>% 
              dplyr::mutate(OR = exp(estimate)) %>% 
              dplyr::pull(OR)
            
          } else{
            m0 = lm(f0, data = run_df)
            m1 = lm(f1, data = run_df)
            lr_test = anova(m0, m1)
            p_anova = lr_test$`Pr(>F)`[2]
            
            effect_OR = as.numeric(NA)
          }
          
          run_df_clean <- run_df %>% 
            dplyr::select(all_of(k), group, gender, Age_at_MBsample, BMI) %>% 
            dplyr::filter(complete.cases(.))
          
          eff_pearson = spcor.test(x = run_df_clean %>% pull(k),
                                   y = run_df_clean %>% pull(group) %>% as.numeric(), 
                                   z = run_df_clean %>% dplyr::select(gender, Age_at_MBsample, BMI), 
                                   method = "pearson")
          
          eff_cliff = cliff.delta(d = run_df %>% dplyr::pull(k), 
                                  f = factor(run_df %>% dplyr::pull("group"), levels = c(1, 0)))
          
          run_output = data.frame(drug = i, 
                                  timeGroup = j, 
                                  MB_feature = k, 
                                  p.value = p_anova,
                                  effect_OR,
                                  effect_pearson = eff_pearson$estimate,
                                  effect_cliff = eff_cliff$estimate, 
                                  transformation = MB_preprocessing)
          
          V2_cumulative_analysis = dplyr::bind_rows(V2_cumulative_analysis, run_output)
        }
      }
    }
  }
  filename = paste("RData/Results/drugUsage/V2_drugUsage_analysis_", ATC_level, "_", MB_preprocessing, ".rds", sep = "")
  saveRDS(V2_cumulative_analysis, filename)
}







# Run the analysis (AB sensitivity)
#--------------------------------------------------#
if (model_nr == 22){
  
  # Categorize subject by the recent drug usage
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
    dplyr::select(skood, ATC_used, key, group) %>% 
    dplyr::filter(!(skood %in% q_AB_users))
  
  
  # Run over all drugs analyzed
  V22_cumulative_analysis = data.frame()
  
  for (i in current_drugs_considered){
    for (j in c("Current_user", "90d-1y", "1y-2y", "2y-3y", "3y-4y", "4y-5y")){
      
      run_df <- recentUsage_analysis_data %>% 
        dplyr::filter(ATC_used == i) %>% 
        dplyr::filter(key == j) %>% 
        dplyr::left_join(countData_raw, by = "skood") %>% 
        dplyr::left_join(diversity_data, by = "skood") %>% 
        dplyr::left_join(PB_ratio, by = "skood") %>% 
        dplyr::left_join(phenotype_data[ ,c("skood", "BMI", "gender", "Age_at_MBsample")], by = "skood")
      
      if (length(table(run_df$group)) == 2){
        
        # Run over all microbial features
        for (k in c(setdiff(colnames(countData_raw), "skood"),
                    "PB_ratio",
                    setdiff(colnames(diversity_data), "skood"))){
          
          f0 = paste("`", k, "` ~ BMI + gender + Age_at_MBsample", sep = "")
          f1 = paste("`", k, "` ~ BMI + gender + Age_at_MBsample + group", sep = "")
          
          if (MB_preprocessing == "PA" & !(k %in% c("PB_ratio", "observed", "shannon"))){
            m0 = glm(f0, data = run_df, family = "binomial")
            m1 = glm(f1, data = run_df, family = "binomial")
            lr_test = anova(m0, m1, test = "Chisq")
            p_anova = lr_test$`Pr(>Chi)`[2]
            
            effect_OR = broom::tidy(m1) %>% 
              dplyr::filter(term == "group") %>% 
              dplyr::mutate(OR = exp(estimate)) %>% 
              dplyr::pull(OR)
            
          } else{
            m0 = lm(f0, data = run_df)
            m1 = lm(f1, data = run_df)
            lr_test = anova(m0, m1)
            p_anova = lr_test$`Pr(>F)`[2]
            
            effect_OR = as.numeric(NA)
          }
          
          run_df_clean <- run_df %>% 
            dplyr::select(all_of(k), group, gender, Age_at_MBsample, BMI) %>% 
            dplyr::filter(complete.cases(.))
          
          eff_pearson = spcor.test(x = run_df_clean %>% pull(k),
                                   y = run_df_clean %>% pull(group) %>% as.numeric(), 
                                   z = run_df_clean %>% dplyr::select(gender, Age_at_MBsample, BMI), 
                                   method = "pearson")
          
          eff_cliff = cliff.delta(d = run_df %>% dplyr::pull(k), 
                                  f = factor(run_df %>% dplyr::pull("group"), levels = c(1, 0)))
          
          run_output = data.frame(drug = i, 
                                  timeGroup = j, 
                                  MB_feature = k, 
                                  p.value = p_anova, 
                                  effect_OR,
                                  effect_pearson = eff_pearson$estimate,
                                  effect_cliff = eff_cliff$estimate, 
                                  transformation = MB_preprocessing)
          
          V22_cumulative_analysis = dplyr::bind_rows(V22_cumulative_analysis, run_output)
        }
      }
    }
  }
  filename = paste("RData/Results/drugUsage/V2_drugUsage_analysis_", ATC_level, "_", MB_preprocessing, "_sensitivityAB.rds", sep = "")
  saveRDS(V22_cumulative_analysis, filename)
}






# Version 3 (AIC comparisons)  ----
#------------------------------------------------#
#                                                #
#         ANALYZING CUMULATIVE EFFECTS           # 
#                                                #
#------------------------------------------------#

# Run the analysis
#--------------------------------------------------#
if (model_nr == 3){

  # Run over all drugs analyzed
  V3_cumulative_analysis = data.frame()
  
  # Read necessary data
  cumulative_usage_data <- readRDS(paste("RData/Interim/Cumulative_usage_", ATC_level, ".rds", sep = ""))
  colnames(cumulative_usage_data)[2] <- "ATC_used"
  
  currentUsage <- current_usage_data %>% 
    tidyr::gather(ATC_used, user_current, -skood)
  
  users_5y <- cumulative_usage_data %>% 
    dplyr::mutate(user_5y = ifelse(n_prescription == 0, 0, 1)) %>% 
    dplyr::select(skood, ATC_used, user_5y, n_prescription)
  
  cumulativeMedicationUsage_analysisData = currentUsage %>% 
    dplyr::left_join(users_5y, by = c("skood", "ATC_used")) %>% 
    dplyr::mutate(user_5y = ifelse(is.na(user_5y) == TRUE, 0, user_5y),
                  n_prescription = ifelse(is.na(n_prescription) == TRUE, 0, n_prescription))

  # Run the analysis
  for (i in current_drugs_considered){
    
    # Divide subjects into 3 groups
    analysis_data = cumulativeMedicationUsage_analysisData %>% 
      dplyr::filter(ATC_used == i) %>% 
      dplyr::ungroup()  %>% 
      dplyr::left_join(countData_raw, by = "skood") %>% 
      dplyr::left_join(diversity_data, by = "skood") %>% 
      dplyr::left_join(PB_ratio, by = "skood") %>% 
      dplyr::left_join(phenotype_data[ ,c("skood", "BMI", "Age_at_MBsample", "gender")], by = "skood") 
    
    # Run over all microbial features
    for (k in c(setdiff(colnames(countData_raw), "skood"),
                "PB_ratio",
                setdiff(colnames(diversity_data), "skood"))){
      
      
      
      if (MB_preprocessing == "PA" & !(k %in% c("observed", "shannon", "PB_ratio"))){
        m1 = glm(paste("`", k, "` ~ BMI + gender + Age_at_MBsample", sep = ""), data = analysis_data, family = "binomial")
        m2 = glm(paste("`", k, "` ~ BMI + gender + Age_at_MBsample + user_current", sep = ""), data = analysis_data, family = "binomial")
        m3 = glm(paste("`", k, "` ~ BMI + gender + Age_at_MBsample + user_current + user_5y", sep = ""), data = analysis_data, family = "binomial")
        m4 = glm(paste("`", k, "` ~ BMI + gender + Age_at_MBsample + user_current + n_prescription ", sep = ""), data = analysis_data, family = "binomial")
        m5 = glm(paste("`", k, "` ~ BMI + gender + Age_at_MBsample + user_current + n_prescription + user_current*n_prescription", sep = ""), data = analysis_data, family = "binomial")
      } else{
        m1 = lm(paste("`", k, "` ~ BMI + gender + Age_at_MBsample", sep = ""), data = analysis_data)
        m2 = lm(paste("`", k, "` ~ BMI + gender + Age_at_MBsample + user_current", sep = ""), data = analysis_data)
        m3 = lm(paste("`", k, "` ~ BMI + gender + Age_at_MBsample + user_current + user_5y", sep = ""), data = analysis_data)
        m4 = lm(paste("`", k, "` ~ BMI + gender + Age_at_MBsample + user_current + n_prescription ", sep = ""), data = analysis_data)
        m5 = lm(paste("`", k, "` ~ BMI + gender + Age_at_MBsample + user_current + n_prescription + user_current*n_prescription", sep = ""), data = analysis_data)
      }

      
      if (MB_preprocessing == "PA" & !(k %in% c("observed", "shannon", "PB_ratio"))){
        a2 = anova(m2, m1, test = "Chisq")
        a3 = anova(m3, m1, test = "Chisq")
        a4 = anova(m4, m1, test = "Chisq")
        a5 = anova(m5, m1, test = "Chisq")
        
        a32 = anova(m3, m2, test = "Chisq")
        a42 = anova(m4, m2, test = "Chisq")
        a52 = anova(m5, m2, test = "Chisq")
        
        p_a2 = a2$`Pr(>Chi)`[2]
        p_a3 = a3$`Pr(>Chi)`[2]
        p_a4 = a4$`Pr(>Chi)`[2]
        p_a5 = a4$`Pr(>Chi)`[2]
        p_a32 = a32$`Pr(>Chi)`[2]
        p_a42 = a42$`Pr(>Chi)`[2]
        p_a52 = a52$`Pr(>Chi)`[2]
        
      } else{
        a2 = anova(m2, m1)
        a3 = anova(m3, m1)
        a4 = anova(m4, m1)
        a5 = anova(m5, m1)
        
        a32 = anova(m3, m2)
        a42 = anova(m4, m2)
        a52 = anova(m5, m2)
        
        p_a2 = a2$`Pr(>F)`[2]
        p_a3 = a3$`Pr(>F)`[2]
        p_a4 = a4$`Pr(>F)`[2]
        p_a5 = a4$`Pr(>F)`[2]
        p_a32 = a32$`Pr(>F)`[2]
        p_a42 = a42$`Pr(>F)`[2]
        p_a52 = a52$`Pr(>F)`[2]
      }

      run_output = data.frame(drug = i, 
                              MB_feature = k, 
                              model = c("Base", "Current_usage", "Current + 5y", "Current + cumulative", "Current + interaction"),
                              p_vs_base = c(as.numeric(NA), p_a2, p_a3, p_a4, p_a5),
                              p_vs_current = c(as.numeric(NA), as.numeric(NA), p_a32, p_a42, p_a52),
                              AIC = c(AIC(m1), AIC(m2), AIC(m3), AIC(m4), AIC(m5)),
                              BIC = c(BIC(m1), BIC(m2), BIC(m3), BIC(m4), BIC(m5)),
                              transformation = MB_preprocessing)
      
      V3_cumulative_analysis = dplyr::bind_rows(V3_cumulative_analysis, run_output)
      
    }
  }
  filename = paste("RData/Results/drugUsage/V3_drugUsage_analysis_", ATC_level, "_", MB_preprocessing, ".rds", sep = "")
  saveRDS(V3_cumulative_analysis, filename)
} 

