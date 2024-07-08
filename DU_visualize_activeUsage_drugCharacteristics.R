

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load the packages
library("tidyverse")
library("ggpubr")


# Abundance data
countData_raw <- readRDS("RData/Interim/mOTUs_CLR_prevalence10p.rds")


# Dosage analysis results
daily_dosage_res <- readRDS("RData/Interim/Dosage_analysis_daily.rds")
max_dosage_res <- readRDS("RData/Interim/Dosage_analysis_max.rds")


# Dosage information
drug_usage_characteristics_ATC5 <- readRDS("RData/Interim/Drug_dosage_data_ATC5.rds")


# Phenotype data
phenotype_data <- readRDS("RData/Interim/Data_master.rds")


# ATC names - curate some names manually for consistence
ATC_names <- readRDS("RData/Interim/DrugUsage_ATC_names.rds")


# Current users
ATC_last_used = readRDS("RData/Interim/Drug_last_used_ATC5.rds")



# Active usage results
Univariate_results <- readRDS("RData/Interim/Naive_univariateResults_V1_CLR.rds")
Univariate_confounding <- readRDS("RData/Results/deconfounding/V1_postHoc_assessment_CLR.rds")

Univariate_results_merged <- Univariate_results %>% 
  dplyr::left_join(Univariate_confounding[ ,c("drug", "timeGroup", "MB_feature", "confounding_status", "confounders")], by = c("drug", "timeGroup", "MB_feature")) %>% 
  dplyr::left_join(ATC_names, by = c("drug" = "drug_code")) %>% 
  dplyr::mutate(CS_hit = ifelse(FDR <= 0.1 & confounding_status == "Confidently deconfounded", "Yes", "")) %>% 
  dplyr::filter(timeGroup == "Current_user")









# Visualize maximum drug dosage effects ----
#------------------------------------------------#
#                                                #
#         VISUALIZE DRUG DOSAGE EFFECTS          # 
#                                                #
#------------------------------------------------#

# Define drug group
drug_group <- phenotype_data %>% 
  dplyr::select(skood) %>% 
  dplyr::full_join(ATC_last_used, by = "skood") %>% 
  tidyr::expand(skood, ATC5) %>% 
  dplyr::left_join(ATC_last_used, by = c("skood", "ATC5")) %>% 
  dplyr::mutate(timeGroup = ifelse(is.na(timeGroup), "non_user", timeGroup))


# Run over all hits
Max_dosage_hits <- max_dosage_res %>% 
  dplyr::filter(FDR_adjusted <= 0.1) %>% 
  dplyr::left_join(ATC_names, by = c("drug" = "drug_code")) %>% 
  dplyr::arrange(FDR_adjusted) %>% 
  dplyr::left_join(Univariate_results_merged[ ,c("drug", "MB_feature", "CS_hit")], by = c("drug", "bug" = "MB_feature"))

for (i in 1:nrow(Max_dosage_hits)){
  run_drug = Max_dosage_hits[i, ] %>% dplyr::pull(drug)
  run_drug_name = Max_dosage_hits[i, ] %>% dplyr::pull(drug_name)
  run_bug = Max_dosage_hits[i, ] %>% dplyr::pull(bug)
  run_CS_hit = Max_dosage_hits[i, ] %>% dplyr::pull(CS_hit)
  
  plot_bug_clean = paste(ifelse(str_detect(run_bug, "ref"), substr(run_bug, 1, nchar(run_bug)-21), substr(run_bug, 1, nchar(run_bug)-22)),
                         " [", substr(run_bug, nchar(run_bug) - 5, nchar(run_bug)-1), "]", sep = "")
  
  # Prepare data
  plotdata_help <- drug_usage_characteristics_ATC5 %>% 
    dplyr::filter(ATC5 == run_drug)
  
  dose_plot_df <- drug_group  %>% 
    dplyr::filter(ATC5 == run_drug) %>% 
    dplyr::left_join(plotdata_help, by = c("skood", "ATC5")) %>% 
    dplyr::left_join(phenotype_data[ ,c("skood", "Age_at_MBsample", "BMI", "gender")], by = "skood") %>% 
    dplyr::left_join(countData_raw, by = "skood") %>%
    dplyr::filter(timeGroup %in% c("non_user", "Current_user")) %>% 
    dplyr::mutate(Last_max_conc = case_when(timeGroup == "non_user" ~ "Non-user", 
                                            timeGroup == "Current_user" & is.na(Last_max_conc) == F ~ as.character(Last_max_conc),
                                            TRUE ~ as.character(NA))) %>% 
    dplyr::filter(is.na(Last_max_conc) == F)
  
  dose_plot_df$bug = dose_plot_df %>% dplyr::pull(run_bug)
  
  
  # Define dose levels 
  dose_levels <- c("Non-user", as.character(sort(as.numeric(setdiff(unique(dose_plot_df$Last_max_conc), "Non-user")))))

  my_comparisons <- list(c(dose_levels[1], "Low dose"), c(dose_levels[1], "High dose"))
  
  # Plot
  plot <- ggplot(dose_plot_df, aes(x = factor(Last_max_conc, levels = dose_levels, labels = c("Non-user", "Low dose", "High dose")), 
                                   y = bug, 
                                   fill = factor(Last_max_conc, levels = dose_levels, labels = c("Non-user", "Low dose", "High dose")))) + 
    geom_boxplot() + 
    stat_compare_means(comparisons = my_comparisons, size = 5) + 
    scale_fill_manual(values = c("#1f78b4", "#fec44f", "#fe9929", "#ec7014", "#cc4c02"), guide = F) + 
    xlab(paste(run_drug_name, " dosage", sep = "")) + 
    ylab(plot_bug_clean) + 
    theme_bw() + 
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18))
  
  # Save output
  run_bug = str_replace(run_bug, ":", "")
  plot_name = paste(run_bug, "_", run_drug, ".png", sep = "")
  plot_name_pdf = paste(run_bug, "_", run_drug, ".pdf", sep = "")
  
  ggsave(plot, file = paste("Figures/Max_daily_dosage/", plot_name, sep = ""), width = 8, height = 6)
  ggsave(plot, file = paste("Figures/Max_daily_dosage/", plot_name_pdf, sep = ""), width = 8, height = 6)
  
}



# Save the dosage results
max_dosage_res_output <- max_dosage_res %>% 
  dplyr::left_join(Univariate_results_merged[ ,c("drug", "MB_feature", "CS_hit")], by = c("drug", "bug" = "MB_feature")) %>% 
  dplyr::left_join(ATC_names, by = c("drug" = "drug_code")) %>% 
  dplyr::filter(p.value_adjusted <= 0.05) %>% 
  dplyr::rename("Q1_hit" = "CS_hit") %>% 
  dplyr::mutate(ATC_group = "ATC5", 
                timeGroup = "Active usage") %>% 
  dplyr::rename("MB_feature" = "bug", 
                "p.value" = "p.value_adjusted", 
                "FDR" = "FDR_adjusted") %>%
  dplyr::select(ATC_group, drug_name, MB_feature, timeGroup, p.value, FDR, Q1_hit)

xlsx::write.xlsx(x = max_dosage_res_output,
                 file = "Results/Extended Data Table 6. Drug dosage analysis.xlsx",
                 sheetName = "Table6",
                 col.names = TRUE,
                 row.names = FALSE, 
                 append = FALSE)



