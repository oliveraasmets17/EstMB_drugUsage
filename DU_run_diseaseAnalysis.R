

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load the packages
library("tidyverse")


# Define preprocessing
MB_preprocessing <- "CLR"


# Read phenotype data
phenotype_data <- readRDS("RData/Interim/Data_master.rds") 
factors_diseases <- readRDS("RData/Interim/Disease_factors_analyzed.rds")


# Cunt data
countData_raw <- readRDS(paste("RData/Interim/mOTUs_", MB_preprocessing, "_prevalence10p.rds", sep = ""))


# PB ratio
PB_ratio <- readRDS("RData/Interim/PB_ratio_mOTUs.rds") %>% 
  dplyr::select(skood, PB_ratio)


# Diversity indexes based on raw data
diversity_data <- readRDS("RData/Interim/mOTUs_alphaDiversity_unfiltered.rds") 







# Analyze naive disease associations ----
#------------------------------------------------#
#                                                #
#              UNIVARIATE ANALYSIS               # 
#                                                #
#------------------------------------------------#

# Prepare data
analysis_data <- phenotype_data %>% 
  dplyr::left_join(countData_raw, by = "skood") %>% 
  dplyr::left_join(diversity_data, by = "skood") %>% 
  dplyr::left_join(PB_ratio, by = "skood") %>% 
  dplyr::filter(is.na(BMI) == F & is.na(Age_at_MBsample) == F)


# Run the analysis
run_disease_analysis = TRUE
if (run_disease_analysis == TRUE){
  disease_res_df = data.frame()
  for (i in factors_diseases){
    
    # Run over all microbial features
    for (j in c(setdiff(colnames(countData_raw), "skood"),
                "PB_ratio",
                setdiff(colnames(diversity_data), "skood"))){
      
      f0 = paste("`", j, "` ~ Age_at_MBsample + gender + BMI", sep = "")
      f1 = paste("`", j, "` ~ Age_at_MBsample + gender + BMI + ", i, sep = "")
      
      if (MB_preprocessing == "CLR" | j %in% c("PB_ratio", "observed", "shannon")){
        m0 = lm(f0, data = analysis_data)
        m1 = lm(f1, data = analysis_data)
        
        a1 = anova(m1, m0)
        anova_p = a1$`Pr(>F)`[2]
      } else{
        m0 = glm(f0, data = analysis_data, family = "binomial")
        m1 = glm(f1, data = analysis_data, family = "binomial")
        
        a1 = anova(m1, m0, test = "Chisq")
        anova_p = a1$`Pr(>Chi)`[2]
      }

      run_df = data.frame(disease = i, 
                          MB_feature = j, 
                          p.value = anova_p)
      
      disease_res_df = dplyr::bind_rows(disease_res_df, run_df)
      
    }
  }
  
  # Merge deconfounding analysis results
  disease_associations_df <- disease_res_df %>% 
    dplyr::group_by(disease) %>% 
    dplyr::mutate(FDR = p.adjust(p.value), method = "BH") %>% 
    dplyr::ungroup() 
  
  saveRDS(disease_associations_df, paste("RData/Interim/Naive_disease_associations_", MB_preprocessing, ".rds", sep = ""))
} else{
  disease_associations_df <- readRDS(paste("RData/Interim/Naive_disease_associations_", MB_preprocessing, ".rds", sep = ""))
}
