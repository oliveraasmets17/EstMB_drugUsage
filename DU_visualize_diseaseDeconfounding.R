

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
library("ggsci")



# Reead phenotype data
phenotype_data <- readRDS("RData/Interim/Data_master.rds") 
diseases <- c(readRDS("RData/Interim/Disease_factors_analyzed.rds"))


# Read result files
disease_associations_CLR <- readRDS("RData/Interim/Naive_disease_associations_CLR.rds") 


# Deconfounding results
Disease_postHoc_assessment_CLR_base <- readRDS("RData/Results/deconfounding/Disease_postHoc_assessment_CLR_base.rds") %>% 
  dplyr::select(disease, MB_feature, confounding_status) %>% 
  dplyr::rename("confounding_status_base"= "confounding_status")

Disease_postHoc_assessment_CLR_drugs <- readRDS("RData/Results/deconfounding/Disease_postHoc_assessment_CLR_drugs.rds") %>% 
  dplyr::select(disease, MB_feature, confounding_status) %>% 
  dplyr::rename("confounding_status_drugs"= "confounding_status")

Disease_postHoc_assessment_CLR_cumulative <- readRDS("RData/Results/deconfounding/Disease_postHoc_assessment_CLR_cumulative.rds") %>% 
  dplyr::select(disease, MB_feature, confounding_status) %>% 
  dplyr::rename("confounding_status_cumulative"= "confounding_status")


# Disease names
ICD10_names <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_EstMiBiom/RData/ICD10_name_data.rds")









# Visualize the confounding associations ----
#------------------------------------------------#
#                                                #
#                VISUALIZE RESULTS               # 
#                                                #
#------------------------------------------------#

# Merge the data
disease_df_final <- disease_associations_CLR %>% 
  dplyr::left_join(Disease_postHoc_assessment_CLR_base, by = c("disease", "MB_feature")) %>% 
  dplyr::left_join(Disease_postHoc_assessment_CLR_drugs, by = c("disease", "MB_feature")) %>% 
  dplyr::left_join(Disease_postHoc_assessment_CLR_cumulative, by = c("disease", "MB_feature")) %>% 
  dplyr::left_join(ICD10_names, by = c("disease" = "ICD10_category"))


# Summarize the number of hits
deconf_aggregated <- disease_df_final %>% 
  dplyr::mutate(disease_name = paste(ICD10_name, " (", disease, ")", sep = "")) %>% 
  dplyr::select(disease_name, MB_feature, p.value, FDR, confounding_status_base, confounding_status_drugs, confounding_status_cumulative) %>% 
  dplyr::mutate(all = as.character(NA)) %>% 
  dplyr::group_by(disease_name) %>% 
  dplyr::summarise(n_nom_all = sum(p.value <= 0.05), 
                   n_FDR_all = sum(FDR <= 0.1), 
                   n_nom_base = sum(p.value <= 0.05 & confounding_status_base == "Confidently deconfounded"), 
                   n_nom_drugs = sum(p.value <= 0.05 & confounding_status_drugs == "Confidently deconfounded"), 
                   n_nom_cumulative = sum(p.value <= 0.05 & confounding_status_cumulative == "Confidently deconfounded"), 
                   n_FDR_base = sum(FDR <= 0.1 & confounding_status_base == "Confidently deconfounded"), 
                   n_FDR_drugs = sum(FDR <= 0.1 & confounding_status_drugs == "Confidently deconfounded"), 
                   n_FDR_cumulative = sum(FDR <= 0.1 & confounding_status_cumulative == "Confidently deconfounded")) %>% 
  dplyr::ungroup() %>% 
  tidyr::gather(key, value, -disease_name) %>% 
  dplyr::mutate(adjustment = ifelse(str_detect(key, "FDR"), "FDR", "nominal"), 
                confounders = substring(key, 7))


# Define disease order
disease_order <- deconf_aggregated %>% 
  dplyr::filter(confounders == "base" & adjustment != "FDR") %>% 
  dplyr::arrange(desc(value)) %>% 
  dplyr::filter(value >= 20) %>% 
  dplyr::pull(disease_name)


# Visualize nominal associations
disease_deconf_plot <- ggplot(deconf_aggregated %>% dplyr::filter(disease_name %in% disease_order & adjustment != "FDR" & confounders != "all"), 
                              aes(x = factor(disease_name, levels = disease_order),  
                                  y = value,
                                  fill = factor(confounders, levels = c("all", "base", "drugs", "cumulative"),
                                                labels = c("Nominal hits", 
                                                           "Deconfounded by lifestyle/diet/other diseases", 
                                                           "Deconfounded by lifestyle/diet/other diseases + active drug usage", 
                                                           "Deconfounded by lifestyle/diet/other diseases + active drug usage + past drug usage")))) + 
  geom_bar(stat = "identity", position = "dodge") + 
  xlab("") + 
  ylab("Number of univariate hits (p-value <= 0.05)") +
  scale_fill_manual(name = "", values = c("#509aca", "#c7d7e2","orange")) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.position = c(0.7, 0.8)) 
disease_deconf_plot
  
ggsave(disease_deconf_plot, filename = "Figures/disease_deconf_plot_nominal.png", width = 14, height = 9)
ggsave(disease_deconf_plot, filename = "Figures/disease_deconf_plot_nominal.pdf", width = 14, height = 9)




