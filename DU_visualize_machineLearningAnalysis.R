

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load the packages
library("tidyverse")
library("ggthemes")


# Abundance data
countData_raw <- readRDS("RData/Interim/mOTUs_PA_prevalence10p.rds") 


# Phenotype data
phenotype_data <- readRDS("RData/Interim/Data_master.rds")


# ATC names - curate some names manually for consistence
ATC_names <- readRDS("RData/Interim/DrugUsage_ATC_names.rds")


# AUC results
AUC_data <- readRDS("RData/Results/Summarized_AUCDf_currentUsage.rds")







# Visualize the results (current usage) ----
#------------------------------------------------#
#                                                #
#             VISUALIZE THE RESULTS              # 
#                                                #
#------------------------------------------------#


# Save the results
#----------------------#
AUC_plotData <- AUC_data %>% 
  dplyr::filter(model_drug %in% c("J01CR", "J01CA", "J01FA")) %>% 
  dplyr::mutate(AUC_final = ifelse(is.na(AUC_test_clean), AUC_model, AUC_test_clean)) %>% 
  dplyr::mutate(model_drug = as.character(model_drug), 
                test_drug = as.character(test_drug)) %>% 
  dplyr::left_join(ATC_names, by = c("test_drug" = "drug_code")) %>% 
  dplyr::rename("test_drug_name" = "drug_name") %>% 
  dplyr::left_join(ATC_names, by = c("model_drug" = "drug_code")) %>% 
  dplyr::rename("model_drug_name" = "drug_name", 
                "AUROC" = "AUC_final") %>% 
  dplyr::select(model_drug_name, test_drug_name, seed, AUROC)

# Active usage
xlsx::write.xlsx(x = AUC_plotData,
                 file = "Results/Extended Data Table 4. Machine learning analysis.xlsx",
                 sheetName = "Table4",
                 col.names = TRUE,
                 row.names = FALSE, 
                 append = FALSE)






# Visualize
#----------------------------------#

# Factor order
test_drug_order <- AUC_plotData %>% 
  dplyr::group_by(test_drug_name) %>% 
  dplyr::summarise(mean_AUC = mean(AUROC)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(mean_AUC)) %>% 
  dplyr::pull(test_drug_name)
  
  
# Heatmap version - average + pairwise comparison
heatmap_plotData <- data.frame()
for (i in unique(AUC_plotData$model_drug_name)){
  for (j in unique(AUC_plotData$test_drug_name)){
    
    run_df = AUC_plotData %>% 
      dplyr::filter(model_drug_name == i & test_drug_name == j)
    
    mean_AUC = mean(run_df$AUROC)
    p_val = t.test(run_df$AUROC, mu = 0.5, alternative = "greater")$p.value
    
    run_res = data.frame(model_drug = i, test_drug = j, mean_AUC, p_val, drug_name = run_df$model_drug_name[1])
    
    heatmap_plotData = dplyr::bind_rows(heatmap_plotData, run_res)
  }
}

heatmap_plotData_final <- heatmap_plotData %>% 
  dplyr::filter(model_drug != test_drug) %>%
  dplyr::mutate(drug_name = case_when(model_drug == "Penicillins With Extended Spectrum (J01CA)" ~ "Penicillins with \nExtended Spectrum \n(J01CA)", 
                                      model_drug == "Penicillins in Combination (J01CR)" ~ "Penicillins in \nCombination \n(J01CR)", 
                                      model_drug == "Macrolides (J01FA)" ~ "Macrolides (J01FA)", 
                                      TRUE ~ drug_name), 
                mean_AUC_class = cut(mean_AUC, breaks = c(-Inf, 0.5, 0.55, 0.6, 0.65, 0.7, Inf)))


# Visualize
ML_generalizability_heatmap <- ggplot(heatmap_plotData_final, 
                                      aes(y = factor(test_drug, levels = test_drug_order), 
                                                                  x = drug_name, fill = mean_AUC_class)) +
  geom_tile(color = "white") + 
  geom_text(aes(label = round(mean_AUC, 2)), size = 5) + 
  scale_fill_manual(values = c("white", "#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c82"), name = "",) + 
  
  xlab("") + 
  ylab("") + 
  xlab("Model trained \nto detect") + 
  ylab("Model tested to detect") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12), 
        title = element_text(size = 18), 
        legend.position = "none",
        axis.title = element_text(size = 20), 
        axis.text.y = element_text(size = 12)) + 
  coord_flip()
ML_generalizability_heatmap

ggsave(ML_generalizability_heatmap, filename = "Figures/drugML_boxplot.png", height = 6, width = 18)
ggsave(ML_generalizability_heatmap, filename = "Figures/drugML_boxplot.pdf", height = 6, width = 18)


