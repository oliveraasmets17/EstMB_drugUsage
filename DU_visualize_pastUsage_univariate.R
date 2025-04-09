

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load the packages
library("tidyverse")
library("ggpubr")
library("ggthemes")


# Taxonomy information 
taxonomy_raw <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/mOTUs_taxonomy_used.rds") %>% 
  tidyr::separate(taxonomy, into = c("k", "p", "c", "o", "f", "g", "s"), sep = "\\|", remove = FALSE) %>% 
  dplyr::arrange(k, p, c, o, f, g, s)

p_order <- taxonomy_raw %>% 
  dplyr::group_by(p) %>%
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::pull(p)


# ATC names - curate some names manually for consistence
ATC_names <- readRDS("RData/Interim/DrugUsage_ATC_names.rds")


# Phenotype data
phenotype_data <- readRDS("RData/Interim/Data_master.rds")


# Currentkly used drugs
currentlyUsed_ATC3_considered <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Medications_currentlyUsed_n20_ATC3.rds")
currentlyUsed_ATC4_considered <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Medications_currentlyUsed_n20_ATC4.rds")
currentlyUsed_ATC5_considered <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Medications_currentlyUsed_n20_ATC5.rds")


# Microbiome data
n_prescriptions_ATC4 <- readRDS("RData/Interim/Cumulative_usage_ATC4.rds")

currentUsage_ATC4 <- readRDS("RData/Interim/Data_medications_currentlyUsed_ATC4.rds") %>% 
  tidyr::gather(drug, currentUsage, -skood)









# Visualize the effects of treatment last used ----
#------------------------------------------------#
#                                                #
#            ANALYZE LONG TERM EFFECTS           # 
#                                                #
#------------------------------------------------#

f_plot_carryover <- function(carryOver_version, carryOver_preprocessing, carryOver_ATC){
  # Read the data
  #----------------------#
  
  # Combined results - output from analyzeDeconfounding.R script
  Univariate_results <- readRDS(paste("RData/Interim/Naive_univariateResults_", carryOver_version, "_", carryOver_preprocessing, ".rds", sep = ""))
  Univariate_confounding <- readRDS(paste("RData/Results/deconfounding/", carryOver_version, "_postHoc_assessment_", carryOver_preprocessing, ".rds", sep = ""))
  
  Univariate_results_merged <- Univariate_results %>% 
    dplyr::left_join(Univariate_confounding[ ,c("drug", "timeGroup", "MB_feature", "confounding_status", "confounders")], by = c("drug", "timeGroup", "MB_feature")) %>% 
    dplyr::left_join(ATC_names, by = c("drug" = "drug_code"))
  
  
  # Focus on only hits that had CS signal
  Univariate_results_direct <- Univariate_results_merged %>% 
    dplyr::filter(timeGroup == "Current_user") %>% 
    dplyr::filter(FDR <= 0.1 & confounding_status == "Confidently deconfounded") %>% 
    dplyr::select(drug, MB_feature) %>% 
    dplyr::mutate(CS_hit = "Yes")
  
  Univariate_results_clean <- Univariate_results_merged %>% 
    dplyr::left_join(Univariate_results_direct, by = c("drug", "MB_feature")) %>% 
    dplyr::filter(CS_hit == "Yes")
  
  
  # Bugs with most carryover effects
  Univariate_results_clean %>%
    dplyr::filter(FDR <= 0.1) %>% 
    dplyr::filter(timeGroup != "Current_user") %>% 
    dplyr::group_by(MB_feature) %>% 
    dplyr::summarise(n = n(), 
                     drugs = paste(drug, collapse = ",")) %>% 
    dplyr::arrange(desc(n))
  
  
  # Summarize the number of hits - current drug usage
  #----------------------#
  n_ATC = case_when(carryOver_ATC == "ATC3" ~ 4,
                    carryOver_ATC == "ATC4" ~ 5,
                    carryOver_ATC == "ATC5" ~ 7)
  
  CarryOver_analysis_df <- Univariate_results_clean %>% 
    # ATC4 level
    dplyr::filter(nchar(as.character(drug)) == n_ATC) %>% 
    dplyr::group_by(drug_name, timeGroup) %>% 
    dplyr::summarise(n = sum(FDR <= 0.1 & confounding_status == "Confidently deconfounded")) %>% 
    dplyr::ungroup()
  
  if (carryOver_version == "V1"){
    CarryOver_analysis_df = CarryOver_analysis_df %>% 
      dplyr::filter(timeGroup != "after_90d") %>% 
      dplyr::mutate(timeGroup = factor(timeGroup, levels = rev(c("Current_user", "after_1y", "after_2y", "after_3y", "after_4y")),
                                       labels = rev(c("Active usage", "Last used more than 1y ago", "Last used more than 2y ago", 
                                                      "Last used more than 3y ago", "Last used more than 4y ago"))))
  } else{
    CarryOver_analysis_df = CarryOver_analysis_df %>% 
      dplyr::filter(timeGroup != "<1y") %>% 
      dplyr::mutate(timeGroup = factor(timeGroup, levels = rev(c("Current_user", "1y-2y", "2y-3y", "3y-4y", "4y-5y")),
                                       labels = rev(c("Active usage", "1-2 years ago", "2-3 years ago", "3-4 years ago", "4-5 years ago"))))
  }
  
  # Drug order - based on CS hits
  carryover_drug_order <- CarryOver_analysis_df %>% 
    dplyr::filter(timeGroup == "Active usage") %>% 
    dplyr::arrange(desc(n)) %>%
    dplyr::filter(n >= 5) %>% 
    dplyr::pull(drug_name)
  
  saveRDS(carryover_drug_order, file = paste("RData/Interim/", "PastUsage_drugOrder_", carryOver_preprocessing, "_", carryOver_version, "_", carryOver_ATC, ".rds", sep = ""))
  
  # Plot
  p_CarryOver <- ggplot(CarryOver_analysis_df %>% dplyr::filter(drug_name %in% carryover_drug_order), 
                        aes(y = factor(drug_name, levels = rev(carryover_drug_order)), 
                            x = n, 
                            fill = timeGroup)) + 
    geom_bar(stat = "identity", position = "dodge", color = "gray10") + 
    xlab("Number of associations (FDR <= 0.1)") + 
    ylab("") + 
    scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150)) + 
    scale_fill_manual(name = "", values = c("#c994c7", "#df65b0", "#e7298a", "#ce1256", "#0570b0")) + 
    theme_classic() + 
    theme(legend.position = c(0.7, 0.1),
          legend.background = element_rect(fill = "white"),
          legend.text = element_text(size = 18),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 18),
          panel.grid.major.x = element_line())
  
  if (carryOver_preprocessing == "PA"){
    p_CarryOver = p_CarryOver + 
      scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 150)) 
  }
  
  ggsave(p_CarryOver, filename = paste("Figures/Past_usage/p_CarryOver_", carryOver_preprocessing, "_", carryOver_version, "_", carryOver_ATC, ".png", sep = ""), height = 10, width = 10)
  ggsave(p_CarryOver, filename = paste("Figures/Past_usage/p_CarryOver_", carryOver_preprocessing, "_", carryOver_version, "_", carryOver_ATC, ".pdf", sep = ""), height = 10, width = 10)
  
}


# Apply function
f_plot_carryover(carryOver_version = "V1", carryOver_preprocessing = "CLR", carryOver_ATC = "ATC4")
f_plot_carryover(carryOver_version = "V1", carryOver_preprocessing = "CLR", carryOver_ATC = "ATC3")

f_plot_carryover(carryOver_version = "V1", carryOver_preprocessing = "PA", carryOver_ATC = "ATC4")
f_plot_carryover(carryOver_version = "V1", carryOver_preprocessing = "PA", carryOver_ATC = "ATC3")






# Visualize the AIC based modelling results ----
#------------------------------------------------#
#                                                #
#                AKAIKE EVALUATION               # 
#                                                #
#------------------------------------------------#
visualize_Akaike <- function(Akaike_preprocessing, Akaike_version, Akaike_ATC){
  
  # Univariate hits
  #----------------------#
  Univariate_results <- readRDS(paste("RData/Interim/Naive_univariateResults_", Akaike_version, "_", Akaike_preprocessing, ".rds", sep = ""))
  Univariate_confounding <- readRDS(paste("RData/Results/deconfounding/", Akaike_version, "_postHoc_assessment_", Akaike_preprocessing, ".rds", sep = ""))
  
  Univariate_results_merged <- Univariate_results %>% 
    dplyr::left_join(Univariate_confounding[ ,c("drug", "timeGroup", "MB_feature", "confounding_status", "confounders")], by = c("drug", "timeGroup", "MB_feature")) %>% 
    dplyr::left_join(ATC_names, by = c("drug" = "drug_code"))
  
  # Focus on only hits that had CS signal
  Univariate_results_direct <- Univariate_results_merged %>% 
    dplyr::filter(timeGroup == "Current_user") %>% 
    dplyr::mutate(CS_hit = ifelse(FDR <= 0.1 & confounding_status == "Confidently deconfounded", "Yes", ""), 
                  CS_hit_nominal = ifelse(p.value <= 0.05, "Yes", "")) %>% 
    dplyr::select(drug, MB_feature, CS_hit, CS_hit_nominal) 

  
  # Summarize the data - pick the simpliest model with AIC - AIC_min < 2
  #----------------------#
  AIC_df <- readRDS(paste("RData/Results/drugUsage/V5_drugUsage_analysis_", Akaike_ATC, "_", Akaike_preprocessing, ".rds", sep = "")) %>% 
    dplyr::filter(!(model %in% c("Current + interaction", "Base"))) %>% 
    dplyr::group_by(drug, MB_feature) %>% 
    dplyr::mutate(min_AIC = min(AIC)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(delta_AIC = AIC - min_AIC,
                  FDR_vs_base = p.adjust(p_vs_base, method = "BH"),
                  FDR_vs_current = p.adjust(p_vs_current, method = "BH")) %>% 
    dplyr::mutate(AIC_real = ifelse(delta_AIC < 2, min_AIC, AIC),
                  model = factor(model, levels = c("Base", "Current_usage", "Current + 5y", "Current + cumulative", "Current + interaction")), 
                  help = 1) %>% 
    dplyr::filter(AIC_real == min_AIC) %>% 
    dplyr::arrange(drug, MB_feature, model) %>% 
    dplyr::group_by(drug, MB_feature) %>% 
    dplyr::mutate(counter = cumsum(help)) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(counter == 1) %>%
    dplyr::left_join(Univariate_results_direct, by = c("drug", "MB_feature"))
  
  # Summarize the number of hits by the best fitting model
  AIC_df_aggregated <- AIC_df %>% 
    dplyr::filter(CS_hit== "Yes") %>% 
    dplyr::mutate(model = ifelse(model == "Current + interaction", "Current + cumulative", as.character(model))) %>% 
    dplyr::group_by(drug, model) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::left_join(ATC_names, by = c("drug" = "drug_code")) %>% 
    dplyr::select(drug, drug_name, model, n)
  
  # Visualize
  AIC_druglist <- readRDS(file = paste("RData/Interim/", "PastUsage_drugOrder_", Akaike_preprocessing, "_", Akaike_version, "_", Akaike_ATC, ".rds", sep = ""))
  
  Akaike_plot <- ggplot(AIC_df_aggregated %>% dplyr::filter(model != "Base") %>% dplyr::filter(drug_name %in% AIC_druglist), 
                        aes(x = factor(drug_name, levels = rev(AIC_druglist)), 
                            y = n, 
                            fill = factor(model, levels = rev(c("Current_usage", "Current + 5y", "Current + cumulative")),
                                          labels = rev(c("Active usage (M1)", "Active + past usage (M2)", "Active usage + amount of usage (M3)"))))) + 
    geom_bar(stat = "identity") + 
    #scale_fill_manual(name = "", values = c("#c7d7e2", "#509aca", "#ce1256")) +  
    scale_fill_manual(name = "Best-fitting model", values = rev(c("#509aca", "#c7d7e2", "#FF6F00FF"))) +  
    theme_classic() + 
    coord_flip() + 
    ylab("Number of associations (FDR <= 0.1)") + 
    xlab("") + 
    scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150)) + 
    theme(axis.title = element_text(size = 18), 
          axis.text = element_text(size = 18),
          legend.position = c(0.66, 0.08), 
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          panel.grid.major.x = element_line())
  
  if (Akaike_preprocessing == "PA"){
    Akaike_plot = Akaike_plot + 
      scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 150)) 
  }
  
  ggsave(Akaike_plot, filename = paste("Figures/Past_usage/Akaike_plot_", Akaike_preprocessing, "_", Akaike_ATC, ".png", sep = ""), height = 10, width = 10)
  ggsave(Akaike_plot, filename = paste("Figures/Past_usage/Akaike_plot_", Akaike_preprocessing, "_", Akaike_ATC, ".pdf", sep = ""), height = 10, width = 10)
}


visualize_Akaike(Akaike_preprocessing = "CLR", Akaike_version = "V1", Akaike_ATC = "ATC4")
visualize_Akaike(Akaike_preprocessing = "CLR", Akaike_version = "V1", Akaike_ATC = "ATC3")
visualize_Akaike(Akaike_preprocessing = "PA", Akaike_version = "V1", Akaike_ATC = "ATC3")
visualize_Akaike(Akaike_preprocessing = "PA", Akaike_version = "V1", Akaike_ATC = "ATC4")








# Save Akaike results ----
#------------------------------------------------#
#                                                #
#                AKAIKE EVALUATION               # 
#                                                #
#------------------------------------------------#

# --------------------------------#
# Combined results - output from analyzeDeconfounding.R script
Univariate_results_CLR <- readRDS("RData/Interim/Naive_univariateResults_V1_CLR.rds")
Univariate_confounding_CLR <- readRDS("RData/Results/deconfounding/V1_postHoc_assessment_CLR.rds")

# Focus on only hits that had CS signal
Univariate_results_merged_CLR <- Univariate_results_CLR %>% 
  dplyr::left_join(Univariate_confounding_CLR[ ,c("drug", "timeGroup", "MB_feature", "confounding_status", "confounders")], by = c("drug", "timeGroup", "MB_feature")) %>% 
  dplyr::left_join(ATC_names, by = c("drug" = "drug_code")) %>% 
  dplyr::mutate(CS_hit = ifelse(FDR <= 0.1 & confounding_status == "Confidently deconfounded", "Yes", "")) %>% 
  dplyr::filter(timeGroup == "Current_user") %>% 
  dplyr::select(drug, MB_feature, CS_hit) 

# Summarize the data - pick the simpliest model with AIC - AIC_min < 2
AIC_CLR_df <- readRDS("RData/Results/drugUsage/V5_drugUsage_analysis_ATC3_CLR.rds") %>%
  dplyr::bind_rows(readRDS("RData/Results/drugUsage/V5_drugUsage_analysis_ATC4_CLR.rds")) %>% 
  dplyr::bind_rows(readRDS("RData/Results/drugUsage/V5_drugUsage_analysis_ATC5_CLR.rds")) %>% 
  dplyr::mutate(ATC_group = case_when(nchar(as.character(drug)) == 4 ~ "ATC3", 
                                      nchar(as.character(drug)) == 5 ~ "ATC4", 
                                      nchar(as.character(drug)) == 7 ~ "ATC5")) %>% 
  dplyr::filter(!(model %in% c("Current + interaction", "Base"))) %>% 
  dplyr::group_by(ATC_group, drug, MB_feature) %>% 
  dplyr::mutate(min_AIC = min(AIC)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(delta_AIC = AIC - min_AIC) %>% 
  dplyr::mutate(AIC_real = ifelse(delta_AIC < 2, min_AIC, AIC),
                model = factor(model, levels = c("Current_usage", "Current + 5y", "Current + cumulative"), 
                               labels = c("Active usage (M1)", "Active + past usage (M2)", "Active usage + amount of usage (M3)")), 
                help = 1) %>% 
  dplyr::arrange(ATC_group, drug, MB_feature, model) %>% 
  dplyr::group_by(ATC_group, drug, MB_feature) %>% 
  dplyr::mutate(counter1 = cumsum(AIC_real == min_AIC),
                counter2 = cumsum(counter1)) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(Univariate_results_merged_CLR, by = c("drug", "MB_feature")) %>% 
  dplyr::filter(CS_hit == "Yes") %>% 
  dplyr::left_join(ATC_names, by = c("drug" = "drug_code")) %>% 
  dplyr::mutate(Best_fit = ifelse(counter2 == 1, "Yes", "")) %>% 
  dplyr::mutate(MB_feature_group = ifelse(MB_feature %in% c("observed", "shannon"), "Diversity", "Centered log-ratio")) %>% 
  dplyr::select(ATC_group, drug_name, MB_feature_group, MB_feature, model, AIC, Best_fit) 
  




# --------------------------------#
# Combined results - output from analyzeDeconfounding.R script
Univariate_results_PA <- readRDS("RData/Interim/Naive_univariateResults_V1_PA.rds")
Univariate_confounding_PA <- readRDS("RData/Results/deconfounding/V1_postHoc_assessment_PA.rds")

# Focus on only hits that had CS signal
Univariate_results_merged_PA <- Univariate_results_PA %>% 
  dplyr::left_join(Univariate_confounding_PA[ ,c("drug", "timeGroup", "MB_feature", "confounding_status", "confounders")], by = c("drug", "timeGroup", "MB_feature")) %>% 
  dplyr::left_join(ATC_names, by = c("drug" = "drug_code")) %>% 
  dplyr::mutate(CS_hit = ifelse(FDR <= 0.1 & confounding_status == "Confidently deconfounded", "Yes", "")) %>% 
  dplyr::filter(timeGroup == "Current_user") %>% 
  dplyr::select(drug, MB_feature, CS_hit) 

# Summarize the data - pick the simpliest model with AIC - AIC_min < 2
AIC_PA_df <- readRDS("RData/Results/drugUsage/V5_drugUsage_analysis_ATC3_PA.rds") %>%
  dplyr::bind_rows(readRDS("RData/Results/drugUsage/V5_drugUsage_analysis_ATC4_PA.rds")) %>% 
  dplyr::bind_rows(readRDS("RData/Results/drugUsage/V5_drugUsage_analysis_ATC5_PA.rds")) %>% 
  dplyr::mutate(ATC_group = case_when(nchar(as.character(drug)) == 4 ~ "ATC3", 
                                      nchar(as.character(drug)) == 5 ~ "ATC4", 
                                      nchar(as.character(drug)) == 7 ~ "ATC5")) %>% 
  dplyr::filter(!(model %in% c("Current + interaction", "Base"))) %>% 
  dplyr::filter(!(MB_feature %in% c("observed", "shannon", "PB_ratio"))) %>% 
  dplyr::group_by(ATC_group, drug, MB_feature) %>% 
  dplyr::mutate(min_AIC = min(AIC)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(delta_AIC = AIC - min_AIC) %>% 
  dplyr::mutate(AIC_real = ifelse(delta_AIC < 2, min_AIC, AIC),
                model = factor(model, levels = c("Current_usage", "Current + 5y", "Current + cumulative"), 
                               labels = c("Active usage (M1)", "Active + past usage (M2)", "Active usage + amount of usage (M3)")), 
                help = 1) %>% 
  dplyr::arrange(ATC_group, drug, MB_feature, model) %>% 
  dplyr::group_by(ATC_group, drug, MB_feature) %>% 
  dplyr::mutate(counter1 = cumsum(AIC_real == min_AIC),
                counter2 = cumsum(counter1)) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(Univariate_results_merged_PA, by = c("drug", "MB_feature")) %>% 
  dplyr::filter(CS_hit == "Yes") %>% 
  dplyr::left_join(ATC_names, by = c("drug" = "drug_code")) %>% 
  dplyr::mutate(Best_fit = ifelse(counter2 == 1, "Yes", "")) %>% 
  dplyr::mutate(MB_feature_group = "Presence-absence") %>% 
  dplyr::select(ATC_group, drug_name, MB_feature_group, MB_feature, model, AIC, Best_fit) 

# Save results
AIC_df <- AIC_PA_df %>% 
  dplyr::bind_rows(AIC_CLR_df) %>%
  dplyr::arrange(ATC_group, MB_feature_group, drug_name, MB_feature, model) 
# 
xlsx::write.xlsx(x = AIC_df %>% as.data.frame(),
                 file = "Results/Extended Data Table 8. Univariate analysis, best fit (AIC).xlsx",
                 sheetName = "Table8",
                 col.names = TRUE,
                 row.names = FALSE,
                 append = FALSE)


