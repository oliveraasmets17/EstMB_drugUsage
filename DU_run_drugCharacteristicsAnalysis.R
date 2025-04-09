

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
phenotype_data <- readRDS("RData/Interim/Data_master.rds") 


# Abundance data
countData_raw <- readRDS("RData/Interim/mOTUs_CLR_prevalence10p.rds")


# Diversity indexes based on raw data
diversity_data <- readRDS("RData/Interim/mOTUs_alphaDiversity_unfiltered.rds") 

  
  




# Run analysis for drug dosage (last max dosage) ----
#------------------------------------------------#
#                                                #
#                   RUN ANALYSIS                 # 
#                                                #
#------------------------------------------------#

# Drug characteritics
drug_usage_characteristics_ATC5 <- readRDS("RData/Interim/Drug_dosage_data_ATC5.rds")

ATC5_analyzed <- readRDS("RData/Interim/Medications_currentlyUsed_n20_ATC5.rds")



# Analyze drug dosage effects
#---------------------------------
run_dosageAnalysis <- F
if (run_dosageAnalysis == T){
  
  max_dosage_res <- data.frame()
  for (i in ATC5_analyzed){
    
    run_df <- drug_usage_characteristics_ATC5 %>% 
      dplyr::filter(ATC5 == i) %>% 
      dplyr::left_join(phenotype_data[ ,c("skood", "Age_at_MBsample", "BMI", "gender")], by = "skood") %>% 
      dplyr::left_join(countData_raw, by = "skood") %>%
      dplyr::left_join(diversity_data, by = "skood") %>% 
      dplyr::filter(is.na(Last_max_conc) == F)
    
    n_dosage <- run_df %>% nrow()
    
    n_unique_dose <- length(table(run_df$Last_max_conc))
    min_unique_dosage <- min(table(run_df$Last_max_conc))
    
    if (n_dosage >= 30 & n_unique_dose != 1 & min_unique_dosage >= 5){
      for (j in setdiff(c(colnames(countData_raw),
                          colnames(diversity_data)),
                        "skood")){
        
        m0 = lm(paste("`", j, "` ~ gender + Age_at_MBsample + BMI", sep = ""), data = run_df)
        m1 = lm(paste("`", j, "` ~ gender + Age_at_MBsample + BMI + factor(Last_max_conc)", sep = ""), data = run_df)
        a1 <- anova(m0, m1)
        
        run_res <- data.frame(drug = i, 
                              bug = j, 
                              p.value_adjusted = a1$`Pr(>F)`[2], 
                              n_dosage, 
                              n_unique_dose, 
                              min_unique_dosage)
        
        max_dosage_res = dplyr::bind_rows(max_dosage_res, run_res)
      }
    }
  }
  # Final modifications
  max_dosage_res <- max_dosage_res %>% 
    dplyr::mutate(FDR_adjusted = p.adjust(p.value_adjusted, method = "BH"))
  
  saveRDS(max_dosage_res, "RData/Interim/Dosage_analysis_max.rds")
} else{
  max_dosage_res <- readRDS("RData/Interim/Dosage_analysis_max.rds")
}


