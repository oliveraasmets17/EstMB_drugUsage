  
  
  
  
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
taxonomy <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/mOTUs_taxonomy_used.rds") %>% 
  tidyr::separate(taxonomy, into = c("k", "p", "c", "o", "f", "g", "s"), sep = "\\|", remove = FALSE) %>% 
  dplyr::arrange(k, p, c, o, f, g, s)

p_order <- taxonomy %>% 
  dplyr::group_by(p) %>%
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::pull(p)
  
  
# ATC names - curate some names manually for consistence
ATC_names <- readRDS("RData/Interim/DrugUsage_ATC_names.rds")

  
# Abundance data
countData_raw <- readRDS("RData/Interim/mOTUs_CLR_prevalence10p.rds")
  
  
# Phenotype data
phenotype_data <- readRDS("RData/Interim/Data_master.rds")
  
  
# Currentkly used drugs
currentlyUsed_ATC3_considered <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Medications_currentlyUsed_n20_ATC3.rds")
currentlyUsed_ATC4_considered <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Medications_currentlyUsed_n20_ATC4.rds")
currentlyUsed_ATC5_considered <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Medications_currentlyUsed_n20_ATC5.rds")

# Factor names
factor_names <- readRDS("RData/Interim/Factor_groups_names.rds")






# Summarize and visualize PERMANOVA results ----
#------------------------------------------------#
#                                                #
#                    PERMANOVA                   # 
#                                                #
#------------------------------------------------#

f_plot_permanova = function(PERMANOVA_preprocessing, PERMANOVA_ATC){
  
  # PERMANOVA results
  PERMANOVA_results_raw <- readRDS(paste("RData/Results/PERMANOVA_results_", PERMANOVA_preprocessing, ".rds", sep = "")) %>% dplyr::mutate(MB_preprocessing = PERMANOVA_preprocessing)
  
  
  # Summarize the number of hits
  #----------------------#
  PERMANOVA_results <- PERMANOVA_results_raw %>% 
    dplyr::filter(analysis_factor %in% c(currentlyUsed_ATC3_considered, currentlyUsed_ATC4_considered, currentlyUsed_ATC5_considered)) %>% 
    dplyr::left_join(ATC_names, by = c("analysis_factor" = "drug_code")) %>% 
    dplyr::mutate(ATC_group = case_when(nchar(analysis_factor) == 4 ~ "ATC3",
                                        nchar(analysis_factor) == 5 ~ "ATC4",
                                        nchar(analysis_factor) == 7 ~ "ATC5",
                                        TRUE ~ "Something else")) %>% 
    dplyr::rename("p.value" = "Pr(>F)") %>% 
    dplyr::group_by(ATC_group, MB_preprocessing) %>% 
    dplyr::mutate(FDR = p.adjust(p.value, method = "BH")) %>% 
    dplyr::ungroup() %>%
    dplyr::select(ATC_group, MB_preprocessing, analysis_factor, drug_name, everything()) %>% 
    dplyr::arrange(ATC_group, desc(R2)) %>% 
    dplyr::mutate(significance = case_when(FDR <= 0.1 ~ "FDR <= 0.1", 
                                           p.value <= 0.05 ~ "Nominal", 
                                           TRUE ~ "NS"))
  
  
  # Total number of hits by factor group
  PERMANOVA_results %>% 
    dplyr::group_by(ATC_group, MB_preprocessing) %>% 
    dplyr::summarise(n = n(), 
                     n_nominal = sum(p.value <= 0.05), 
                     n_FDR = sum(FDR <= 0.1)) %>%
    dplyr::ungroup()
  
  
  
  
  # Visualize the results
  #----------------------#
  PERMANOVA_factorLevels = PERMANOVA_results %>%
    dplyr::filter(FDR <= 0.1) %>% 
    dplyr::filter(ATC_group %in% PERMANOVA_ATC) %>% 
    dplyr::filter(MB_preprocessing == PERMANOVA_preprocessing) %>% 
    dplyr::arrange(R2) %>% 
    dplyr::pull(drug_name)
  
  plot_PERMANOVA <- ggplot(PERMANOVA_results %>% dplyr::filter(drug_name %in% PERMANOVA_factorLevels) %>% dplyr::filter(MB_preprocessing == PERMANOVA_preprocessing) , 
                           aes(x = factor(drug_name, levels = rev(PERMANOVA_factorLevels)), 
                               y = R2*100, 
                               fill = ifelse(FDR <= 0.1, "FDR <= 0.1", "NS"))) + 
    geom_bar(stat = "identity", color = "gray30") + 
    scale_fill_manual(values = c("slategray4", "gray80"), name = "") + 
    theme_classic() + 
    xlab("") + 
    ylab("R2 (%)") + 
    #coord_flip() + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 10), 
          axis.text.y = element_text(size = 14), 
          axis.title = element_text(size = 18),
          legend.background = element_rect(color = "black"), 
          legend.text = element_text(size = 12), 
          legend.position = "none")
  
  plot_PERMANOVA
  
  ggsave(plot_PERMANOVA, filename = paste("Figures/Active_usage/PERMANOVA_currentUsage_", PERMANOVA_preprocessing, "_", PERMANOVA_ATC, ".png", sep = ""), height = 6, width = 20)
  ggsave(plot_PERMANOVA, filename = paste("Figures/Active_usage/PERMANOVA_currentUsage_", PERMANOVA_preprocessing, "_", PERMANOVA_ATC, ".pdf", sep = ""), height = 6, width = 20)
  
}

f_plot_permanova(PERMANOVA_preprocessing = "CLR", PERMANOVA_ATC = "ATC4")





# Save PERMANOVA results
# -------------------------------
PERMANOVA_results_CLR <- readRDS("RData/Results/PERMANOVA_results_CLR.rds") %>% 
  dplyr::filter(analysis_factor %in% c(currentlyUsed_ATC3_considered, currentlyUsed_ATC4_considered, currentlyUsed_ATC5_considered)) %>% 
  dplyr::left_join(ATC_names, by = c("analysis_factor" = "drug_code")) %>% 
  dplyr::mutate(ATC_group = case_when(nchar(analysis_factor) == 4 ~ "ATC3",
                                      nchar(analysis_factor) == 5 ~ "ATC4",
                                      nchar(analysis_factor) == 7 ~ "ATC5",
                                      TRUE ~ "Something else")) %>% 
  dplyr::rename("p.value" = "Pr(>F)") %>% 
  dplyr::group_by(ATC_group) %>% 
  dplyr::mutate(FDR = p.adjust(p.value, method = "BH"), 
                R2 = R2*100) %>% 
  dplyr::ungroup() %>%
  dplyr::select(ATC_group, drug_name, F.Model, R2, p.value, FDR) %>%
  dplyr::arrange(ATC_group, desc(R2)) 

xlsx::write.xlsx(x = PERMANOVA_results_CLR %>% as.data.frame(),
                 file = "Results/Extended Data Table 2. Beta diversity analysis.xlsx", 
                 sheetName = "Table2",
                 row.names = FALSE,
                 col.names = TRUE, 
                 append = FALSE)








