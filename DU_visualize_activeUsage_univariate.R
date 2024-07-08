

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


# Currently used drugs
currentlyUsed_ATC3_considered <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Medications_currentlyUsed_n20_ATC3.rds")
currentlyUsed_ATC4_considered <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Medications_currentlyUsed_n20_ATC4.rds")
currentlyUsed_ATC5_considered <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Medications_currentlyUsed_n20_ATC5.rds")


# Number of drug users
n_currentlyUsed_ATC3 <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Data_medications_currentlyUsed_ATC3.rds") %>% 
  tidyr::gather(drug, value, -skood) %>% 
  dplyr::group_by(drug) %>% 
  dplyr::summarise(n_users = sum(value))

n_currentlyUsed_ATC4 <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Data_medications_currentlyUsed_ATC4.rds") %>% 
  tidyr::gather(drug, value, -skood) %>% 
  dplyr::group_by(drug) %>% 
  dplyr::summarise(n_users = sum(value))

n_currentlyUsed_ATC5 <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Data_medications_currentlyUsed_ATC5.rds") %>% 
  tidyr::gather(drug, value, -skood) %>% 
  dplyr::group_by(drug) %>% 
  dplyr::summarise(n_users = sum(value))

n_currentlyUsed <- dplyr::bind_rows(n_currentlyUsed_ATC3, n_currentlyUsed_ATC4, n_currentlyUsed_ATC5)








# Summarize the number of univariate hits ----
#------------------------------------------------#
#                                                #
#         SUMMARIZE UNIVARIATE ANALYSIS          #  
#                                                #
#------------------------------------------------#


# Read the data
#----------------------#

# Combined results - output from analyzeDeconfounding.R script
naive_univariate_results_PA <- readRDS("RData/Interim/Naive_univariateResults_V1_PA.rds") %>% dplyr::mutate(MB_preprocessing = "PA")
naive_univariate_results_CLR <- readRDS("RData/Interim/Naive_univariateResults_V1_CLR.rds") %>% dplyr::mutate(MB_preprocessing = "CLR")

naive_univariate_results <- dplyr::bind_rows(naive_univariate_results_PA, naive_univariate_results_CLR)

# Preprocessed post-hoc summary
univariate_postHoc_assessment_PA <- readRDS("RData/Results/deconfounding/V1_postHoc_assessment_PA.rds") %>% dplyr::mutate(MB_preprocessing = "PA")
univariate_postHoc_assessment_CLR <- readRDS("RData/Results/deconfounding/V1_postHoc_assessment_CLR.rds") %>% dplyr::mutate(MB_preprocessing = "CLR")

univariate_postHoc_assessment <- dplyr::bind_rows(univariate_postHoc_assessment_PA, univariate_postHoc_assessment_CLR)



# Summarize the number of hits - current drug usage
#----------------------#
Univariate_results_clean <- naive_univariate_results %>%
  dplyr::left_join(univariate_postHoc_assessment[ ,c("drug", "MB_feature", "timeGroup", "confounding_status", "confounders", "MB_preprocessing")], 
                   by = c("drug", "MB_feature", "timeGroup", "MB_preprocessing")) %>% 
  dplyr::mutate(MB_factor_group = case_when(MB_feature %in% c("observed", "shannon", "PB_ratio") ~ "Diversity", 
                                            TRUE ~ "Univariate"),
                ATC_group = case_when(nchar(as.character(drug)) == 4 ~ "ATC3",
                                      nchar(as.character(drug)) == 5 ~ "ATC4",
                                      nchar(as.character(drug)) == 7 ~ "ATC5",
                                      TRUE ~ "Something else"),
                timeGroup = case_when(as.character(timeGroup) == "Current_user" ~ "Active usage",
                                      as.character(timeGroup) == "after_1y" ~ "Last used more than 1 year ago",
                                      as.character(timeGroup) == "after_2y" ~ "Last used more than 2 years ago",
                                      as.character(timeGroup) == "after_3y" ~ "Last used more than 3 years ago",
                                      as.character(timeGroup) == "after_4y" ~ "Last used more than 4 years ago",
                                      as.character(timeGroup) == "after_5y" ~ "Last used more than 5 years ago")) %>% 
  dplyr::left_join(ATC_names, by = c("drug" = "drug_code")) 

table(Univariate_results_clean$timeGroup, Univariate_results_clean$confounding_status)
table(Univariate_results_clean$p.value <= 0.05, Univariate_results_clean$confounding_status)
table(Univariate_results_clean$FDR <= 0.1, Univariate_results_clean$confounding_status)




# Different effects for AB vs host-targeted drugs
#----------------------#
Univariate_results_clean %>% 
  dplyr::filter(timeGroup == "Active usage" & MB_preprocessing == "PA") %>%
  dplyr::filter(FDR <= 0.1 & confounding_status == "Confidently deconfounded") %>% 
  dplyr::filter(nchar(drug) == 5) %>% 
  dplyr::mutate(factor_group = ifelse(substr(drug, 1, 3) == "J01", "AB", "notAB")) %>% 
  dplyr::group_by(MB_feature, factor_group) %>% 
  dplyr::summarise(n = n()) %>% 
  tidyr::spread(factor_group, n) %>% 
  dplyr::arrange(desc(AB))


# Similarity of the profiles within a higher ATC
#----------------------#
drugSignal_hierarchical_similarity = Univariate_results_clean %>% 
  dplyr::filter(FDR <= 0.1 & confounding_status == "Confidently deconfounded" & timeGroup == "Active usage" ) %>% 
  dplyr::mutate(ATC_high = substr(drug, 1, 4)) %>% 
  dplyr::group_by(ATC_high, drug) %>% 
  dplyr::summarise(n_associations = n()) %>% 
  dplyr::left_join(n_currentlyUsed, by = "drug") %>% 
  dplyr::mutate(drug_ATC_group = case_when(nchar(drug) == 4 ~ "ATC3",
                                           nchar(drug) == 5 ~ "ATC4",
                                           nchar(drug) == 7 ~ "ATC5")) %>% 
  dplyr::group_by(ATC_high, drug_ATC_group) %>% 
  dplyr::mutate(mean_sd = mean(n_associations)/sd(n_associations)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(mean_sd)

head(drugSignal_hierarchical_similarity, 20)



# Output the univariate results
#----------------------#
CS_hits_CLR <- Univariate_results_clean %>%
  dplyr::filter(MB_preprocessing == "CLR") %>% 
  dplyr::filter(timeGroup != "<1y") %>% 
  dplyr::filter(p.value <= 0.05) %>%
  dplyr::mutate(MB_feature_group = ifelse(MB_feature %in% c("observed", "shannon"), "Diversity", "Centered log-ratio")) %>% 
  dplyr::select(ATC_group, drug_name, MB_feature_group, MB_feature, timeGroup, effect_pearson, effect_cliff, p.value, FDR, confounding_status, confounders) %>%
  dplyr::arrange(ATC_group, drug_name, MB_feature_group, MB_feature, timeGroup) %>% 
  as.data.frame()

CS_hits_PA <- Univariate_results_clean %>%
  dplyr::filter(MB_preprocessing == "PA") %>% 
  dplyr::filter(!(MB_feature %in% c("observed", "shannon", "PB_ratio"))) %>% 
  dplyr::filter(timeGroup != "<1y") %>% 
  dplyr::filter(p.value <= 0.05) %>%
  dplyr::mutate(MB_feature_group = "Presence-absence") %>% 
  dplyr::select(ATC_group, drug_name, MB_feature_group, MB_feature, timeGroup, effect_pearson, effect_cliff, p.value, FDR, confounding_status, confounders) %>%
  dplyr::arrange(ATC_group, drug_name, MB_feature_group, MB_feature, timeGroup) %>% 
  as.data.frame()

CS_hits_output <- CS_hits_CLR %>%
  dplyr::bind_rows(CS_hits_PA)





# Save the results
#----------------------#

# Active usage
xlsx::write.xlsx(x = CS_hits_output %>% dplyr::filter(timeGroup == "Active usage"),
                 file = "Results/Extended Data Table 3. Univariate analysis, active usage effects.xlsx",
                 sheetName = "Table3",
                 col.names = TRUE,
                 row.names = FALSE, 
                 append = FALSE)


# Past usage
xlsx::write.xlsx(x = CS_hits_output %>% dplyr::filter(timeGroup != "Active usage"),
                 file = "Results/Extended Data Table 7. Univariate analysis, carryover effects.xlsx",
                 sheetName = "Table7",
                 col.names = TRUE,
                 row.names = FALSE, 
                 append = FALSE)







# Visualize effect sizes for current usage of drugs ----
#------------------------------------------------#
#                                                #
#      ASSOCIATIONS WITH CURRENT DRUG USAGE      # 
#                                                #
#------------------------------------------------#

# Prepare data for visualization
#----------------------#
drugs_displayed <- Univariate_results_clean  %>% 
  dplyr::filter(ATC_group %in% c("ATC4")) %>% 
  dplyr::filter(MB_preprocessing %in% c("PA", "CLR")) %>% 
  dplyr::filter(!(MB_feature %in% c("observed", "shannon", "PB_ratio"))) %>% 
  dplyr::filter(timeGroup == "Active usage") %>% 
  dplyr::group_by(MB_preprocessing, drug_name) %>% 
  dplyr::summarise(n = sum(p.value <= 0.05 & confounding_status == "Confidently deconfounded", na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(MB_preprocessing, n) %>% 
  dplyr::mutate(total = CLR + PA) %>% 
  dplyr::arrange(desc(total)) %>% 
  tibble::rowid_to_column(var = "rank") %>% 
  dplyr::filter(rank <= 20) %>% 
  dplyr::pull(drug_name)


# Microbial factors displayed - at least 5 associations within the drugs considered for plotting
bugs_displayed <- Univariate_results_clean  %>% 
  dplyr::filter(ATC_group %in% c("ATC4")) %>% 
  dplyr::filter(MB_preprocessing %in% c("PA", "CLR")) %>% 
  dplyr::filter(!(MB_feature %in% c("observed", "shannon", "PB_ratio"))) %>% 
  dplyr::filter(timeGroup == "Active usage") %>% 
  dplyr::group_by(MB_feature, MB_preprocessing) %>% 
  dplyr::summarise(n = sum(p.value <= 0.05 & confounding_status == "Confidently deconfounded", na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(MB_preprocessing, n) %>% 
  dplyr::filter(PA + CLR >= 10) %>% 
  dplyr::pull(MB_feature)


# Final plot data
univariate_plotData <- Univariate_results_clean %>% 
  dplyr::filter(ATC_group %in% c("ATC4")) %>% 
  dplyr::filter(MB_preprocessing %in% c("PA", "CLR")) %>% 
  dplyr::filter(timeGroup == "Active usage") %>% 
  dplyr::mutate(effect_pearson = ifelse(p.value <= 0.05 & confounding_status == "Confidently deconfounded", effect_pearson, 0)) %>% 
  dplyr::filter(drug_name %in% drugs_displayed) %>% 
  dplyr::filter(MB_feature %in% bugs_displayed) %>% 
  dplyr::mutate(drug_name = factor(drug_name, levels = rev(drugs_displayed))) %>% 
  dplyr::mutate(MB_feature = factor(MB_feature, levels = c(taxonomy_raw$taxa_raw, "observed", "shannon", "PB_ratio"), labels = c(taxonomy_raw$taxa_raw, "Observed richness", "Shannon diversity", "Prevotella-Bacteroides ratio"))) %>% 
  dplyr::mutate(Significance_label = case_when(FDR <= 0.1 & confounding_status == "Confidently deconfounded" ~ "*", 
                                               p.value <= 0.05 & confounding_status == "Confidently deconfounded" ~ "", 
                                               TRUE ~ ""))



# Association plot
#--------------------------#
p_currentUsage <- ggplot(univariate_plotData, aes(x = MB_feature, y = drug_name, fill = effect_pearson)) + 
  geom_tile(color = "gray90") + 
  geom_text(aes(label = Significance_label)) + 
  scale_fill_gradient2(name = "Partial Pearson correlation", 
                       low = "cyan3", mid = "white", high = "deeppink2", midpoint = 0, limits = c(min(univariate_plotData$effect_pearson), 
                                                                                                  max(univariate_plotData$effect_pearson))) + 
  xlab("") + 
  ylab("") + 
  theme_bw() +
  facet_wrap(vars(factor(MB_preprocessing, levels = c("PA", "CLR"), labels = c("Presence-absence", "Abundance (centered log-ratio)"))), ncol = 1) + 
  theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), 
    axis.text.x = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "white"), 
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 16),
    axis.ticks = element_blank())
p_currentUsage

ggsave(plot = p_currentUsage, filename = "Figures/Active_usage/CurrentUsage_pheatmap_ATC4.png", width = 14, height = 8)
ggsave(plot = p_currentUsage, filename = "Figures/Active_usage/CurrentUsage_pheatmap_ATC4.pdf", width = 14, height = 8)



# Annotation plot
#--------------------------#
annotation_df <- univariate_plotData %>% 
  dplyr::distinct(MB_feature) %>%
  dplyr::left_join(taxonomy_raw, by = c("MB_feature" = "taxa_raw")) %>% 
  dplyr::select(MB_feature, p, f) %>% 
  dplyr::mutate(p = ifelse(is.na(p) == TRUE, "x_Composition", p), 
                f = ifelse(is.na(f) == TRUE, MB_feature, f)) %>% 
  tidyr::gather(key, value, -MB_feature) %>% 
  dplyr::mutate(key = ifelse(key == "p", "Phylum", "Family")) %>% 
  dplyr::mutate(MB_feature = factor(MB_feature, levels = c(taxonomy_raw$taxa_raw, "Observed richness", "Shannon diversity", "Prevotella-Bacteroides ratio")))
  

# Generate colour palette
n_fam <- univariate_plotData %>% 
  dplyr::distinct(MB_feature) %>%
  dplyr::left_join(taxonomy_raw, by = c("MB_feature" = "taxa_raw")) %>%
  dplyr::mutate(p = ifelse(is.na(p) == TRUE, "x_Composition", p), 
                f = ifelse(is.na(f) == TRUE, MB_feature, f)) %>% 
  dplyr::distinct(p, f) %>% 
  dplyr::group_by(p) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(p)

color_palette = c()
main_colors = c("#084081", "#EF3B2C", "#54278f", "#41ab5d", "#fed976", "black", "#dd3497")
color_list = list(c("#2b8cbe", "#7bccc4"),
                  c("#fc4e2a", "#fd8d3c"),
                  c("#6a51a3", "#9e9ac8"), 
                  c("#74c476", "#a1d99b", "#c7e9c0", "#78c679", "#addd8e", "#66c2a4", "#99d8c9"), 
                  c("#ffeda0", "#ffffcc"), 
                  c("#bdbdbd", "#737373"), 
                  c("#fa9fb5", "#f768a1"))
set.seed(1)
for (i in 1:nrow(n_fam)){
  color_palette = c(color_palette, main_colors[i])
  
  run_colors = colorRampPalette(color_list[[i]])(n = n_fam[i, ] %>% dplyr::pull(n))
  color_palette = c(color_palette, sample(run_colors))
}

# Order for the taxonomy for coloring
taxa_order_df <- univariate_plotData %>% 
  dplyr::distinct(MB_feature) %>%
  dplyr::left_join(taxonomy_raw, by = c("MB_feature" = "taxa_raw")) %>% 
  dplyr::mutate(p = ifelse(is.na(p) == TRUE, "x_Composition", p), 
                f = ifelse(is.na(f) == TRUE, MB_feature, f)) %>%
  dplyr::distinct(p, f) %>% 
  dplyr::arrange(p, f) %>% 
  dplyr::mutate(p_lag = lag(p),
                p_lag = ifelse(is.na(p_lag), "", p_lag)) 

taxa_order = c()
for (i in 1:nrow(taxa_order_df)){
  
  if (taxa_order_df$p[i] != taxa_order_df$p_lag[i]){
    taxa_order = c(taxa_order, taxa_order_df[i, ] %>% dplyr::pull(p))
  }
  taxa_order = c(taxa_order, taxa_order_df[i, ] %>% dplyr::pull(f))
}

# Visualize annotation
p_annotation <- ggplot(annotation_df, aes(x = MB_feature, y = key, fill = factor(value, levels = taxa_order))) + 
  geom_tile(color = "gray90") + 
  scale_fill_manual(name = "", values = color_palette) +  
  theme_classic() + 
  xlab("") + 
  ylab("") + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
        axis.text.y = element_text(size = 16),
        axis.ticks = element_blank(),
        axis.line = element_blank())

ggsave(plot = p_annotation, filename = "Figures/Active_usage/CurrentUsage_pheatmap_annotation_legend.png", width = 20, height = 7)
ggsave(plot = p_annotation, filename = "Figures/Active_usage/CurrentUsage_pheatmap_annotation_legend.pdf", width = 20, height = 7)







# Diversity plot
#--------------------------#

# PERMANOVA results
PERMANOVA_drugOrder <- readRDS("RData/Results/PERMANOVA_results_CLR.rds") %>% 
  dplyr::filter(analysis_factor %in% readRDS("RData/Interim/Medications_currentlyUsed_n20_ATC4.rds")) %>% 
  dplyr::left_join(ATC_names, by = c("analysis_factor" = "drug_code")) %>% 
  dplyr::rename("p.value" = "Pr(>F)") %>% 
  dplyr::mutate(FDR = p.adjust(p.value, method = "BH")) %>% 
  dplyr::select( analysis_factor, drug_name, everything()) %>% 
  dplyr::arrange(desc(R2)) %>% 
  dplyr::filter(FDR <= 0.1) %>% 
  dplyr::pull(drug_name)
  

# Alpha diversity plot data
diversity_plotData <- Univariate_results_clean %>% 
  dplyr::filter(ATC_group %in% c("ATC4")) %>% 
  dplyr::filter(MB_preprocessing %in% c("CLR")) %>% 
  dplyr::filter(timeGroup == "Active usage") %>% 
  dplyr::mutate(effect_pearson = ifelse(p.value <= 0.05, effect_pearson, 0)) %>% 
  dplyr::filter(MB_feature %in% c("observed", "shannon", "PB_ratio")) %>% 
  dplyr::filter(drug_name %in% PERMANOVA_drugOrder) %>% 
  dplyr::mutate(drug_name = factor(drug_name, levels = PERMANOVA_drugOrder)) %>% 
  dplyr::mutate(MB_feature = factor(MB_feature, levels = rev(c("observed", "shannon", "PB_ratio")), labels = rev(c("Observed richness", "Shannon diversity", "Prevotella-Bacteroides ratio")))) %>% 
  dplyr::mutate(Significance_label = case_when(FDR <= 0.1 & confounding_status == "Confidently deconfounded" ~ "*", 
                                               #FDR <= 0.1 ~ "*",
                                               p.value <= 0.05  & confounding_status == "Confidently deconfounded"~ "", 
                                               TRUE ~ ""))

# Visualize
diversity_heatmap <- ggplot(diversity_plotData, aes(x = drug_name, y = MB_feature, fill = effect_pearson)) + 
  geom_tile(color = "gray70") + 
  geom_text(aes(label = Significance_label)) + 
  scale_fill_gradient2(name = "Partial Pearson correlation", 
                       low = "cyan3", mid = "white", high = "deeppink2", midpoint = 0, limits = c(min(univariate_plotData$effect_pearson), 
                                                                                                  max(univariate_plotData$effect_pearson))) + 
  theme_bw() + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 12), 
        axis.text.y = element_text(size = 14), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 16),
        legend.position = "bottom")
diversity_heatmap


ggsave(plot = diversity_heatmap, filename = "Figures/Active_usage/Alpha_diversity_heatmap.png", width = 20, height = 4.5)
ggsave(plot = diversity_heatmap, filename = "Figures/Active_usage/Alpha_diversity_heatmap.pdf", width = 20, height = 4.5)
