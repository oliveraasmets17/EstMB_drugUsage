

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load the packages
library("pheatmap")
library("ggsci")
library("ggtree")
library("ggnewscale")
library("ape")
library("ggpubr")
library("tidyverse")
library("readr")
library("xlsx")
library("ggthemes")
library("ggExtra")


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
factor_names <- readRDS("RData/Interim/Factor_groups_names.rds") %>% 
  dplyr::mutate(factor_name = ifelse(factor_name == "Stool consistenct", "Stool consistency", factor_name))







# Visualize results of variance partitioning analysis ----
#------------------------------------------------#
#                                                #
#              VARIANCE PARTITIONING             # 
#                                                #
#------------------------------------------------#

# Define distance for visualization
# -------------------------------
varPart_preprocessing <- "CLR"


# Analysis results
# -------------------------------
varPart_res <- readRDS(paste("RData/Results/VariancePartitioning_", varPart_preprocessing, "_forward_SET02.rds", sep = "")) %>% 
  dplyr::select(-var.formula)


# Define factor groups
# -------------------------------
factors_diseases <- c(readRDS("RData/Interim/Disease_factors_analyzed.rds"), "gumDiseaseDiagnosed", "seasonableAllergy", "chickenpox")
factors_drugs <- c(currentlyUsed_ATC3_considered, currentlyUsed_ATC4_considered, currentlyUsed_ATC5_considered)
factors_cumulativeDrugs <- paste("Cumulative_", factors_drugs, sep = "")
factors_procedures <- c(readRDS("RData/Interim/Procedure_factors_analyzed.rds"), "cecumRemoved", "tonsilsRemoved")
factors_dietary <- c(readRDS("RData/Interim/Dietary_factors_analyzed.rds"), "eatingHabit_name")
factors_lifestyle <- c("has_smoked", "alcohol_has_used", "doesPhysicalExercise")
factors_antophometric <- readRDS("RData/Interim/Intrinsic_factors_analyzed.rds")
factors_stool <- c("usualStoolType_category", "frequencyGutEmpting_name") 



# Prepare the data
# -------------------------------
varPart_data <- varPart_res %>% 
  dplyr::filter(element == "unique") %>% 
  dplyr::mutate(var = str_replace_all(var, "`", ""), 
                var = ifelse(var == "\n    Cumulative_B01AF", "Cumulative_B01AF", var)) %>% 
  dplyr::mutate(factor_group = case_when(var %in% factors_diseases ~ "Diseases",
                                         var %in% factors_drugs ~ "Active drug \nusage",
                                         var %in% factors_procedures ~ "Medical \nprocedures",
                                         var %in% factors_dietary ~ "Diet",
                                         var %in% factors_lifestyle ~ "Lifestyle",
                                         var %in% factors_antophometric ~ "Antophometric",
                                         var %in% factors_stool ~ "Stool \ncharacteristics",
                                         var %in% factors_cumulativeDrugs ~ "Past drug \nusage",
                                         TRUE ~ "Other")) 



# Save the results
# -------------------------------
varPart_data_output <- varPart_res %>% 
  dplyr::filter(element == "unique") %>% 
  dplyr::mutate(var = str_replace_all(var, "`", ""), 
                var = ifelse(var == "\n    Cumulative_B01AF", "Cumulative_B01AF", var), 
                var_help = case_when(str_detect(var, "Cumulative") ~ substring(var, 12), 
                                     var %in% factors_drugs ~ var, 
                                     TRUE ~ as.character(NA))) %>% 
  dplyr::mutate(factor_group = case_when(var %in% factors_diseases ~ "Diseases",
                                         var %in% factors_drugs ~ "Active drug usage",
                                         var %in% factors_procedures ~ "Medical procedures",
                                         var %in% factors_dietary ~ "Diet",
                                         var %in% factors_lifestyle ~ "Lifestyle",
                                         var %in% factors_antophometric ~ "Antophometric",
                                         var %in% factors_stool ~ "Stool characteristics",
                                         var %in% factors_cumulativeDrugs ~ "Past drug usage",
                                         TRUE ~ "Other")) %>% 
  dplyr::left_join(factor_names[ ,c("factor", "factor_name")], by = c("var" = "factor")) %>% 
  dplyr::left_join(ATC_names, by = c("var_help" = "drug_code")) %>% 
  dplyr::mutate(factor_name = case_when(is.na(factor_name) == FALSE ~ factor_name, 
                                        is.na(var_help) == F & str_detect(var, "Cumulative") ~ paste("Cumulative usage of ", drug_name, sep = ""), 
                                        is.na(var_help) == F ~ drug_name, 
                                        var == "education" ~ "Highest education level", 
                                        TRUE ~ "Something else"), 
                r.squared = r.squared*100) %>% 
  dplyr::rename("R2" = "r.squared") %>%
  dplyr::select(factor_group, factor_name, R2) %>% 
  dplyr::arrange(factor_group, desc(R2))
  
xlsx::write.xlsx(x = varPart_data_output,
                 file = "Results/Extended Data Table 9. Multivariate variance partitioning analysis.xlsx",
                 sheetName = "Table9",
                 col.names = TRUE,
                 row.names = FALSE, 
                 append = FALSE)





# Visualize
# -------------------------------

# Define factor order
varPart_factorOrder <- varPart_data %>% 
  dplyr::group_by(factor_group) %>% 
  dplyr::summarize(sum_var = sum(r.squared)) %>% 
  dplyr::arrange(desc(sum_var)) %>% 
  dplyr::pull(factor_group)

varPart_plotData <- varPart_data %>% 
  dplyr::mutate(factor_group = factor(factor_group, levels = varPart_factorOrder))


# Plot
plot_varPart <- ggplot(varPart_plotData, 
                       aes(x = factor(factor_group, levels = varPart_factorOrder), 
                           y = r.squared*100, 
                           fill = factor_group)) + 
  geom_bar(stat = "identity", color = "black") + 
  theme_classic() + 
  xlab("") + 
  ylab("R2 (%)") + 
  scale_fill_futurama() +
  scale_fill_manual(values = c("#C71000FF", "#008EA0FF","#FF6F00FF", "#8A4198FF", "#5A9599FF", "#509aca", "#84D7E1FF", "#FF95A8FF", "#3D3B25FF")) + 
  theme(axis.title = element_text(size = 16), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1), 
        legend.position = "none",
        panel.grid.major.y = element_line())

ggsave(plot_varPart, filename = paste("Figures/Past_usage/VariancePartitioning_", varPart_preprocessing, ".png", sep = ""), height = 6, width = 8)
ggsave(plot_varPart, filename = paste("Figures/Past_usage/VariancePartitioning_", varPart_preprocessing, ".pdf", sep = ""), , height = 6, width = 8)
