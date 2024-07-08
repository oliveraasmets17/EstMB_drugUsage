
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
library("stringr")
library("tidyverse")
library("ggpubr")


# Read the data
drugs_ATC4 <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Medications_currentlyUsed_n20_ATC4.rds")


# Previously analyzed data
phenotype_data <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Data_master.rds")


# ATC names 
ATC_names <- readRDS("RData/Interim/DrugUsage_ATC_names.rds")


# Drug starting period
MBAD_drugInitiation_time_ATC3 <- readRDS("RData/Interim/MBAD_drugUsage_ATC3_q2.rds") %>% dplyr::rename("ATC" = "ATC3")
MBAD_drugInitiation_time_ATC4 <- readRDS("RData/Interim/MBAD_drugUsage_ATC4_q2.rds") %>% dplyr::rename("ATC" = "ATC4")
MBAD_drugInitiation_time_ATC5 <- readRDS("RData/Interim/MBAD_drugUsage_ATC5_q2.rds") %>% dplyr::rename("ATC" = "ATC5")

MBAD_drugInitiation_time <- dplyr::bind_rows(MBAD_drugInitiation_time_ATC3, MBAD_drugInitiation_time_ATC4, MBAD_drugInitiation_time_ATC5)

MBAD_drugInitiation_summary_ATC3 <- readRDS("RData/Interim/MBAD_drugUsage_ATC3_q3.rds")
MBAD_drugInitiation_summary_ATC4 <- readRDS("RData/Interim/MBAD_drugUsage_ATC4_q3.rds")
MBAD_drugInitiation_summary_ATC5 <- readRDS("RData/Interim/MBAD_drugUsage_ATC5_q3.rds")


# Count data
EstMB_countdata_raw <- readRDS("RData/Interim/mOTUs_CLR_prevalence10p.rds")
MBAD_countdata_raw <- readRDS("RData/Interim/mOTUs_TP2_CLR_prevalence10p.rds")


# Diversity
diversity_TP1 <- readRDS("RData/Interim/mOTUs_alphaDiversity_unfiltered.rds")
diversity_TP2 <- readRDS("RData/Interim/mOTUs_TP2_alphaDiversity_unfiltered.rds")






# Preprocess count data ----
#------------------------------------------------#
#                                                #
#             PREPROCESS COUNT DATA              # 
#                                                #
#------------------------------------------------#

# Data to long format 
EstMB_countdata <- EstMB_countdata_raw %>%
  tidyr::gather(taxa, EstMB_value ,-c("skood"))

MBAD_countdata <- MBAD_countdata_raw %>%
  tidyr::gather(taxa, MBAD_value ,-c("skood"))


# Merge EstMB and MBAD data, calculate abundance delta 
Analysis_countdata <- MBAD_countdata %>% 
  dplyr::left_join(EstMB_countdata, by = c("skood", "taxa")) %>% 
  dplyr::mutate(delta_abundance = MBAD_value-EstMB_value) %>% 
  dplyr::filter(is.na(delta_abundance) == F) %>% 
  dplyr::select(skood, taxa, delta_abundance) %>% 
  tidyr::spread(taxa, delta_abundance)


# Calculate diversity difference
Diversity_df <- diversity_TP1 %>%
  dplyr::left_join(diversity_TP2, by = c("skood")) %>% 
  dplyr::mutate(delta_observed = observed_TP2 - observed,
                delta_shannon = shannon_TP2 - shannon) %>% 
  dplyr::select(skood, delta_observed, delta_shannon) %>% 
  dplyr::filter(complete.cases(.))







# Analyze drug initiation effects ----
#------------------------------------------------#
#                                                #
#                  RUN ANALYSIS                  # 
#                                                #
#------------------------------------------------#

# Define phenotype as non-users at all time points vs users at MBAD timepoint
scodes_all <- unique(MBAD_countdata_raw$skood)
drugs1 <- c(MBAD_drugInitiation_summary_ATC3$ATC3[MBAD_drugInitiation_summary_ATC3$`Initiator, active usage` >= 10], 
            MBAD_drugInitiation_summary_ATC4$ATC4[MBAD_drugInitiation_summary_ATC4$`Initiator, active usage` >= 10], 
            MBAD_drugInitiation_summary_ATC5$ATC5[MBAD_drugInitiation_summary_ATC5$`Initiator, active usage` >= 10])

shell1 <- data.frame(skood = rep(scodes_all, each = length(drugs1)), 
                    ATC = rep(drugs1, length(scodes_all)))


Analysis_df1 <- shell1 %>% 
  dplyr::left_join(MBAD_drugInitiation_time, by = c("skood", "ATC")) %>% 
  dplyr::mutate(starter_group = ifelse(is.na(After_MBAD) == TRUE, "Non-user", starter_group)) %>%
  dplyr::filter(starter_group %in% c("Non-user", "Initiator, active usage")) %>% 
  dplyr::left_join(phenotype_data[ ,c("skood", "Age_at_MBsample", "BMI", "gender")], by = "skood")

run_analysis <- FALSE
if (run_analysis == TRUE){
  
  analysis1_output <- data.frame()
  
  # Run over all drugs
  for (i in drugs1){
    
    Analysis_df_run <- Analysis_df1 %>% 
      dplyr::filter(ATC == i) %>% 
      dplyr::filter(is.na(starter_group) == FALSE) %>% 
      dplyr::select(skood, Age_at_MBsample, BMI, gender, starter_group) %>% 
      dplyr::left_join(Analysis_countdata, by = "skood") %>% 
      dplyr::left_join(Diversity_df, by = "skood") %>% 
      dplyr::filter(complete.cases(.))
    
    # Run over all bugs
    if (nrow(Analysis_df_run) > 1 & length(unique(Analysis_df_run$starter_group)) == 2 & min(table(Analysis_df_run$starter_group)) >= 2){
      for (j in c(setdiff(colnames(MBAD_countdata_raw), "skood"), "delta_observed", "delta_shannon")){
        
        eff_cliff = cliff.delta(d = Analysis_df_run %>% dplyr::pull(j), 
                                f = factor(Analysis_df_run %>% dplyr::pull("starter_group"), levels = c("Initiator, active usage", "Non-user")))
        
        eff_pearson = spcor.test(x = Analysis_df_run %>% dplyr::pull(j),
                                 y = Analysis_df_run %>% dplyr::pull(starter_group) %>% factor(levels = c("Non-user", "Initiator, active usage")) %>% as.numeric(), 
                                 z = Analysis_df_run %>% dplyr::select(gender, Age_at_MBsample, BMI), 
                                 method = "pearson")
        
        run_output = data.frame(drug = i, 
                                MB_feature = j, 
                                Delta_p.value = eff_pearson$p.value, 
                                Delta_effect_pearson = eff_pearson$estimate,
                                Delta_effect_cliff = eff_cliff$estimate, 
                                n_nonusers = sum(Analysis_df_run$starter_group == "Non-user"), 
                                n_initiators = sum(Analysis_df_run$starter_group == "Initiator, active usage"))
        
        analysis1_output = dplyr::bind_rows(analysis1_output, run_output)
      }
    }
  }
  saveRDS(analysis1_output, "RData/Results/drugUsage/MBAD_analysis1_output.rds")
} 








# Analyze drug discontinuation effects ----
#------------------------------------------------#
#                                                #
#                  RUN ANALYSIS                  # 
#                                                #
#------------------------------------------------#

# Define phenotype as non-users at all time points vs users at MBAD timepoint
scodes_all <- unique(MBAD_countdata_raw$skood)
drugs2 <- c(MBAD_drugInitiation_summary_ATC3$ATC3[MBAD_drugInitiation_summary_ATC3$`Stopper, active usage` >= 10], 
            MBAD_drugInitiation_summary_ATC4$ATC4[MBAD_drugInitiation_summary_ATC4$`Stopper, active usage` >= 10], 
            MBAD_drugInitiation_summary_ATC5$ATC5[MBAD_drugInitiation_summary_ATC5$`Stopper, active usage` >= 10])

shell2 <- data.frame(skood = rep(scodes_all, each = length(drugs2)), 
                     ATC = rep(drugs2, length(scodes_all)))

Analysis_df2 <- shell2 %>% 
  dplyr::left_join(MBAD_drugInitiation_time, by = c("skood", "ATC")) %>% 
  dplyr::mutate(starter_group = ifelse(is.na(After_MBAD) == TRUE, "Non-user", starter_group)) %>%
  dplyr::filter(starter_group %in% c("Non-user", "Stopper, active usage")) %>% 
  dplyr::left_join(phenotype_data[ ,c("skood", "Age_at_MBsample", "BMI", "gender")], by = "skood")

run_analysis <- FALSE
if (run_analysis == TRUE){
  
  analysis2_output <- data.frame()
  
  # Run over all drugs
  for (i in drugs2){
    
    Analysis_df_run <- Analysis_df2 %>% 
      dplyr::filter(ATC == i) %>% 
      dplyr::filter(is.na(starter_group) == FALSE) %>% 
      dplyr::select(skood, Age_at_MBsample, BMI, gender, starter_group) %>% 
      dplyr::left_join(Analysis_countdata, by = "skood") %>% 
      dplyr::left_join(Diversity_df, by = "skood") %>% 
      dplyr::filter(complete.cases(.))
    
    # Run over all bugs
    if (nrow(Analysis_df_run) > 1 & length(unique(Analysis_df_run$starter_group)) == 2 & min(table(Analysis_df_run$starter_group)) >= 2){
      for (j in c(setdiff(colnames(MBAD_countdata_raw), "skood"), "delta_observed", "delta_shannon")){
        
        eff_cliff = cliff.delta(d = Analysis_df_run %>% dplyr::pull(j), 
                                f = factor(Analysis_df_run %>% dplyr::pull("starter_group"), levels = c("Stopper, active usage", "Non-user")))
        
        eff_pearson = spcor.test(x = Analysis_df_run %>% dplyr::pull(j),
                                 y = Analysis_df_run %>% dplyr::pull(starter_group) %>% factor(levels = c("Non-user", "Stopper, active usage")) %>% as.numeric(), 
                                 z = Analysis_df_run %>% dplyr::select(gender, Age_at_MBsample, BMI), 
                                 method = "pearson")
        
        run_output = data.frame(drug = i, 
                                MB_feature = j, 
                                Delta_p.value = eff_pearson$p.value, 
                                Delta_effect_pearson = eff_pearson$estimate,
                                Delta_effect_cliff = eff_cliff$estimate, 
                                n_nonusers = sum(Analysis_df_run$starter_group == "Non-user"), 
                                n_initiators = sum(Analysis_df_run$starter_group == "Stopper, active usage"))
        
        analysis2_output = dplyr::bind_rows(analysis2_output, run_output)
      }
    }
  }
  saveRDS(analysis2_output, "RData/Results/drugUsage/MBAD_analysis2_output.rds")
} 







# Analyze drug initiation effects (carryover) ----
#------------------------------------------------#
#                                                #
#                  RUN ANALYSIS                  # 
#                                                #
#------------------------------------------------#

# Define phenotype as non-users at all time points vs users at MBAD timepoint
scodes_all <- unique(MBAD_countdata_raw$skood)
drugs3 <- c(MBAD_drugInitiation_summary_ATC3$ATC3[MBAD_drugInitiation_summary_ATC3$`Initiator, past usage` >= 10], 
            MBAD_drugInitiation_summary_ATC4$ATC4[MBAD_drugInitiation_summary_ATC4$`Initiator, past usage` >= 10], 
            MBAD_drugInitiation_summary_ATC5$ATC5[MBAD_drugInitiation_summary_ATC5$`Initiator, past usage` >= 10])

shell3 <- data.frame(skood = rep(scodes_all, each = length(drugs3)), 
                     ATC = rep(drugs3, length(scodes_all)))


Analysis_df3 <- shell3 %>% 
  dplyr::left_join(MBAD_drugInitiation_time, by = c("skood", "ATC")) %>% 
  dplyr::mutate(starter_group = ifelse(is.na(After_MBAD) == TRUE, "Non-user", starter_group)) %>%
  dplyr::filter(starter_group %in% c("Non-user", "Initiator, past usage")) %>% 
  dplyr::left_join(phenotype_data[ ,c("skood", "Age_at_MBsample", "BMI", "gender")], by = "skood")

run_analysis <- FALSE
if (run_analysis == TRUE){
  
  analysis3_output <- data.frame()
  
  # Run over all drugs
  for (i in drugs3){
    
    Analysis_df_run <- Analysis_df3 %>% 
      dplyr::filter(ATC == i) %>% 
      dplyr::filter(is.na(starter_group) == FALSE) %>% 
      dplyr::select(skood, Age_at_MBsample, BMI, gender, starter_group) %>% 
      dplyr::left_join(Analysis_countdata, by = "skood") %>% 
      dplyr::left_join(Diversity_df, by = "skood") %>% 
      dplyr::filter(complete.cases(.))
    
    # Run over all bugs
    if (nrow(Analysis_df_run) > 1 & length(unique(Analysis_df_run$starter_group)) == 2 & min(table(Analysis_df_run$starter_group)) >= 2){
      for (j in c(setdiff(colnames(MBAD_countdata_raw), "skood"), "delta_observed", "delta_shannon")){
        
        eff_cliff = cliff.delta(d = Analysis_df_run %>% dplyr::pull(j), 
                                f = factor(Analysis_df_run %>% dplyr::pull("starter_group"), levels = c("Initiator, past usage", "Non-user")))
        
        eff_pearson = spcor.test(x = Analysis_df_run %>% dplyr::pull(j),
                                 y = Analysis_df_run %>% dplyr::pull(starter_group) %>% factor(levels = c("Non-user", "Initiator, past usage")) %>% as.numeric(), 
                                 z = Analysis_df_run %>% dplyr::select(gender, Age_at_MBsample, BMI), 
                                 method = "pearson")
        
        run_output = data.frame(drug = i, 
                                MB_feature = j, 
                                Delta_p.value = eff_pearson$p.value, 
                                Delta_effect_pearson = eff_pearson$estimate,
                                Delta_effect_cliff = eff_cliff$estimate, 
                                n_nonusers = sum(Analysis_df_run$starter_group == "Non-user"), 
                                n_initiators = sum(Analysis_df_run$starter_group == "Initiator, past usage"))
        
        analysis3_output = dplyr::bind_rows(analysis3_output, run_output)
      }
    }
  }
  saveRDS(analysis3_output, "RData/Results/drugUsage/MBAD_analysis3_output.rds")
} 




