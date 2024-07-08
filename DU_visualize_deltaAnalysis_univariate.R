

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
library("ggsankey")


# ATC names 
ATC_names <- readRDS("RData/Interim/DrugUsage_ATC_names.rds")


# Read MBAD analysis results
analysis1_output <- readRDS("RData/Results/drugUsage/MBAD_analysis1_output.rds") %>% 
  dplyr::select(drug, MB_feature, Delta_p.value, Delta_effect_pearson) %>% 
  dplyr::rename("Delta1_pearson" = "Delta_effect_pearson",
                "Delta1_pval" = "Delta_p.value")

analysis2_output <- readRDS("RData/Results/drugUsage/MBAD_analysis2_output.rds") %>% 
  dplyr::select(drug, MB_feature, Delta_p.value, Delta_effect_pearson) %>% 
  dplyr::rename("Delta2_pearson" = "Delta_effect_pearson",
                "Delta2_pval" = "Delta_p.value")

analysis3_output <- readRDS("RData/Results/drugUsage/MBAD_analysis3_output.rds") %>% 
  dplyr::select(drug, MB_feature, Delta_p.value, Delta_effect_pearson) %>% 
  dplyr::rename("Delta3_pearson" = "Delta_effect_pearson", 
                "Delta3_pval" = "Delta_p.value")


# Cross-sectional hits
naive_univariate_results <- readRDS("RData/Interim/Naive_univariateResults_V1_CLR.rds")
univariate_postHoc_assessment <- readRDS("RData/Results/deconfounding/V1_postHoc_assessment_CLR.rds") 


# Drug starting period
MBAD_drugInitiation_time_ATC4 <- readRDS("RData/Interim/MBAD_drugUsage_ATC4_q2.rds")

MBAD_drugInitiation_summary_ATC3 <- readRDS("RData/Interim/MBAD_drugUsage_ATC3_q3.rds") %>% dplyr::rename("ATC" = "ATC3") %>% dplyr::mutate(ATC_level = "ATC3")
MBAD_drugInitiation_summary_ATC4 <- readRDS("RData/Interim/MBAD_drugUsage_ATC4_q3.rds") %>% dplyr::rename("ATC" = "ATC4") %>% dplyr::mutate(ATC_level = "ATC4")
MBAD_drugInitiation_summary_ATC5 <- readRDS("RData/Interim/MBAD_drugUsage_ATC5_q3.rds") %>% dplyr::rename("ATC" = "ATC5") %>% dplyr::mutate(ATC_level = "ATC5")


# Count data
EstMB_countdata_raw <- readRDS("RData/Interim/mOTUs_CLR_prevalence10p.rds")
MBAD_countdata_raw <- readRDS("RData/Interim/mOTUs_TP2_CLR_prevalence10p.rds")


# ATC names
ATC_names <- readRDS("RData/Interim/DrugUsage_ATC_names.rds")


# Taxonomy information 
taxonomy_raw <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/mOTUs_taxonomy_used.rds") %>% 
  tidyr::separate(taxonomy, into = c("k", "p", "c", "o", "f", "g", "s"), sep = "\\|", remove = FALSE) %>% 
  dplyr::arrange(k, p, c, o, f, g, s)








# Save drug number of drug initiators ----
#------------------------------------------------#
#                                                #
#            SUMMARIZE DRUG USAGE DATA           # 
#                                                #
#------------------------------------------------#


# Define phenotype as non-users at all time points vs users at MBAD timepoint
scodes_all <- unique(MBAD_countdata_raw$skood)
drugs_all <- c(MBAD_drugInitiation_summary_ATC3$ATC, MBAD_drugInitiation_summary_ATC4$ATC, MBAD_drugInitiation_summary_ATC5$ATC)

shell <- data.frame(skood = rep(scodes_all, each = length(drugs_all)), 
                    ATC = rep(drugs_all, length(scodes_all)))

MBAD_drugInitiation_summary <- dplyr::bind_rows(MBAD_drugInitiation_summary_ATC3, 
                                                MBAD_drugInitiation_summary_ATC4, 
                                                MBAD_drugInitiation_summary_ATC5) %>% 
  dplyr::mutate(`Non-user-clean` = 328 - (`Continuous user` + `Initiator, active usage` + `Initiator, past usage` + `Non-user` +
                                            `Something else` + `Stopper, active usage` + `Stopper, past active`),
                `Non-users` = `Non-user-clean` + `Non-user`) %>% 
  dplyr::rename("Discontinued drug usage" = "Stopper, active usage") %>% 
  dplyr::ungroup() %>% 
  dplyr::select(ATC_level, drug_name, `Non-users`, `Initiator, active usage`, `Initiator, past usage`, `Discontinued drug usage`) %>% 
  dplyr::filter(`Initiator, active usage` >= 10 | `Initiator, past usage` >= 10 | `Discontinued drug usage` >= 10) %>% 
  as.data.frame()

xlsx::write.xlsx(x = MBAD_drugInitiation_summary,
                 file = "Results/Extended Data Table 10. Number of drug initiators.xlsx",
                 sheetName = "Table10",
                 col.names = TRUE,
                 row.names = FALSE,
                 append = FALSE)






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








# Combine and save analysis results ----
#------------------------------------------------#
#                                                #
#                COMPARE RESULTS                 # 
#                                                #
#------------------------------------------------#

# CS results
#----------------------#
Univariate_results_CLR <- naive_univariate_results %>%
  dplyr::left_join(univariate_postHoc_assessment[ ,c("drug", "MB_feature", "timeGroup", "confounding_status", "confounders")], 
                   by = c("drug", "MB_feature", "timeGroup")) %>% 
  dplyr::mutate(MB_factor_group = case_when(MB_feature %in% c("observed", "shannon", "PB_ratio") ~ as.character(MB_feature), 
                                            TRUE ~ "Univariate"),
                ATC_group = case_when(nchar(as.character(drug)) == 4 ~ "ATC3",
                                      nchar(as.character(drug)) == 5 ~ "ATC4",
                                      nchar(as.character(drug)) == 7 ~ "ATC5",
                                      TRUE ~ "Something else")) %>% 
  dplyr::left_join(ATC_names, by = c("drug" = "drug_code")) %>% 
  dplyr::filter(timeGroup == "Current_user")


# Repeated measures analysis results
#----------------------#
analysis1 <- readRDS("RData/Results/drugUsage/MBAD_analysis1_output.rds") %>% 
  dplyr::mutate(case_group = "Initiator, active usage")
  
analysis2<- readRDS("RData/Results/drugUsage/MBAD_analysis2_output.rds") %>% 
  dplyr::mutate(case_group = "Discontinued drug usage")

analysis3 <- readRDS("RData/Results/drugUsage/MBAD_analysis3_output.rds") %>% 
  dplyr::mutate(case_group = "Initiator, past usage") 


# Combine data
#----------------------#
MBAD_results_output <- bind_rows(analysis1, analysis2, analysis3) %>% 
  dplyr::left_join(Univariate_results_CLR, by = c("drug", "MB_feature")) %>% 
  dplyr::mutate(Q1_hit = ifelse(FDR <= 0.1 & confounding_status == "Confidently deconfounded", "Yes", "")) %>% 
  dplyr::select(case_group, ATC_group, drug_name, MB_feature, Delta_p.value, Delta_effect_pearson, Q1_hit) %>% 
  dplyr::group_by(case_group, ATC_group) %>% 
  dplyr::mutate(Delta_FDR = p.adjust(Delta_p.value, method = "BH")) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(case_group, ATC_group, drug_name, MB_feature, Delta_effect_pearson, Delta_p.value,  Delta_FDR, Q1_hit) %>%
  dplyr::arrange(ATC_group, drug_name, MB_feature) 

table(MBAD_results_output$Delta_p.value <= 0.05, MBAD_results_output$Q1_hit %>% is.na(), MBAD_results_output$case_group)


# Save output
xlsx::write.xlsx(x = MBAD_results_output %>% as.data.frame(),
                 file = "Results/Extended Data Table 11. Univariate analysis, drug initiation.xlsx",
                 sheetName = "Table11",
                 col.names = TRUE,
                 row.names = FALSE,
                 append = FALSE)








# Compare with univariate CLR analysis ----
#------------------------------------------------#
#                                                #
#                COMPARE RESULTS                 # 
#                                                #
#------------------------------------------------#

# Compare repeated-measures results with the univariate results
#----------------------#
comparison_df <- Univariate_results_CLR %>% 
  dplyr::left_join(analysis1_output, by = c("drug", "MB_feature")) %>% 
  dplyr::left_join(analysis2_output, by = c("drug", "MB_feature")) %>% 
  dplyr::left_join(analysis3_output, by = c("drug", "MB_feature")) %>% 
  dplyr::mutate(plot_group = ifelse(FDR <= 0.1 & confounding_status == "Confidently deconfounded", "Cross-sectional hit (FDR <= 0.1)", "No hit"))


# Analyze overlapping hits
#----------------------#
overlapping_hits <- comparison_df %>% 
  dplyr::mutate(CS_hit = ifelse(FDR <= 0.1 & confounding_status == "Confidently deconfounded", "Yes", "-")) %>%
  dplyr::mutate(consistent_estimate = sign(Delta1_pearson*Delta3_pearson))

prop.table(table(overlapping_hits$CS_hit, overlapping_hits$Delta1_pval <= 0.05), margin = 1)


overlapping_hits_CS <- overlapping_hits %>% 
  dplyr::filter(CS_hit == "Yes")

table(overlapping_hits_CS$drug)
table(overlapping_hits_CS$drug, overlapping_hits_CS$Delta1_pval <= 0.05)
table(overlapping_hits_CS$drug, overlapping_hits_CS$Delta2_pval <= 0.05)
table(overlapping_hits_CS$drug, overlapping_hits_CS$Delta3_pval <= 0.05)






# Long term effects ----
#------------------------------------------------#
#                                                #
#               LONG TERM EFFECTS                # 
#                                                #
#------------------------------------------------#

# Merge
analysis_output_merged <- analysis1_output %>% 
  dplyr::left_join(analysis2_output, by = c("drug", "MB_feature")) %>% 
  dplyr::left_join(analysis3_output, by = c("drug", "MB_feature")) %>% 
  dplyr::left_join(Univariate_results_CLR[ ,c("drug", "MB_feature", "p.value", "FDR", "confounding_status")], by = c("drug", "MB_feature")) %>% 
  dplyr::filter(FDR <= 0.1 & confounding_status == "Confidently deconfounded") 

table(analysis_output_merged$Delta1_pval <= 0.05, analysis_output_merged$Delta3_pval <= 0.05)
chisq.test(table(analysis_output_merged$Delta1_pval <= 0.05, analysis_output_merged$Delta3_pval <= 0.05))


# Hits in both analysis
both_hits <- analysis_output_merged %>% 
  dplyr::filter(Delta1_pval <= 0.05 & Delta3_pval <= 0.05) %>% 
  dplyr::filter(p.value <= 0.05)

table(both_hits$drug) %>% sort()





# Visualize effect similarity
#----------------------#
both_hits_effect <- analysis_output_merged %>% 
  dplyr::filter(Delta1_pval <= 0.05) %>% 
  dplyr::select(drug, MB_feature, Delta1_pearson, Delta2_pearson, Delta3_pearson) %>% 
  tidyr::gather(key, effect, -c("drug", "MB_feature")) %>% 
  dplyr::mutate(key = case_when(key == "Delta1_pearson" ~ "analysis1", 
                                key == "Delta2_pearson" ~ "analysis2",
                                key == "Delta3_pearson" ~ "analysis3"),
                name = paste(drug, " ", MB_feature, sep = "")) 

both_hits_p <- analysis_output_merged %>% 
  dplyr::filter(Delta1_pval <= 0.05) %>% 
  dplyr::select(drug, MB_feature, Delta1_pval, Delta2_pval, Delta3_pval) %>% 
  tidyr::gather(key, p.value, -c("drug", "MB_feature")) %>% 
  dplyr::mutate(key = case_when(key == "Delta1_pval" ~ "analysis1", 
                                key == "Delta2_pval" ~ "analysis2",
                                key == "Delta3_pval" ~ "analysis3"),
                name = paste(drug, " ", MB_feature, sep = "")) 

both_hits_plotDf <- both_hits_effect %>% 
  dplyr::left_join(both_hits_p[ ,c("name", "key", "p.value")], by = c("name", "key")) %>% 
  dplyr::left_join(ATC_names, by = c("drug" = "drug_code")) %>% 
  dplyr::mutate(bug_help = ifelse(str_detect(MB_feature, "meta"), 
                                  substr(MB_feature, 1, nchar(MB_feature) - 22),
                                  substr(MB_feature, 1, nchar(MB_feature) - 21)),
                bug_clean = paste(bug_help, " [", substr(MB_feature, nchar(MB_feature) - 5, nchar(MB_feature)), sep = ""), 
                drug_name = ifelse(drug == "A02BC", "PPI (A02BC)", drug_name))

# Naming
both_hits_plotOrder <- both_hits_plotDf %>% 
  dplyr::filter(key == "analysis1") %>% 
  dplyr::arrange(drug, desc(effect)) %>% 
  dplyr::pull(name)

both_hits_plotOrder_names <- data.frame(MB_feature = both_hits_plotOrder) %>% 
  dplyr::mutate(MB_feature = substring(MB_feature, 7), 
                bug_help = ifelse(str_detect(MB_feature, "meta"), 
                                  substr(MB_feature, 1, nchar(MB_feature) - 22),
                                  substr(MB_feature, 1, nchar(MB_feature) - 21)),
                bug_clean = paste(bug_help, " [", substr(MB_feature, nchar(MB_feature) - 5, nchar(MB_feature)), sep = "")) %>% 
  dplyr::pull(bug_clean)


# Visualize
p_effect_similarity <- ggplot() + 
  geom_path(data = both_hits_plotDf %>% dplyr::filter(drug %in% c("J01CR", "A02BC", "J01FA")) %>% dplyr::filter(key %in% c("analysis1", "analysis3")), 
            aes(x = factor(name, levels = both_hits_plotOrder, labels = both_hits_plotOrder_names), 
                y = effect), color = "gray70") +
  geom_path(data = both_hits_plotDf %>% dplyr::filter(drug %in% c("J01CR", "A02BC", "J01FA")) %>% dplyr::filter(key %in% c("analysis2", "analysis3")), 
            aes(x = factor(name, levels = both_hits_plotOrder, labels = both_hits_plotOrder_names), 
                y = effect), color = "gray70", linetype = 3) +
  geom_point(data = both_hits_plotDf %>% dplyr::filter(drug %in% c("J01CR", "A02BC", "J01FA")), 
             aes(x = factor(name, levels = both_hits_plotOrder, labels = both_hits_plotOrder_names), 
                 y = effect, 
                 color = factor(key, levels = c("analysis1", "analysis3", "analysis2"),
                                labels = c("Active usage at T2", "Last usage >1y ago from T2", "Dicontinued drug usage")),
                 shape = ifelse(p.value <= 0.05, "Statistically significant", "NS")), 
             size = 5) + 
  geom_hline(yintercept = 0, linetype = 3, size = 1, color = "gray70") + 
  scale_color_manual(name = "", values = c("#509aca", "#c7d7e2", "yellowgreen")) + 
  scale_fill_manual(name = "", values = c("#509aca", "#c7d7e2", "yellowgreen")) + 
  scale_shape_manual(values = c(21, 19)) + 
  theme_bw() + 
  facet_grid(. ~ drug_name, scales = "free", space = "free") + 
  xlab("") + 
  ylab("Effect size") + 
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1), 
        legend.position = "top",
        legend.background = element_rect(color = "black"), 
        legend.title = element_blank(),
        strip.text = element_text(size = 18), 
        strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 12))
p_effect_similarity

ggsave(p_effect_similarity, file = "Figures/MBAD_shortVsLong_effectSize_comparison.png", width = 18, height = 8)
ggsave(p_effect_similarity, file = "Figures/MBAD_shortVsLong_effectSize_comparison.pdf", width = 18, height = 8)











# Visualize some examples of change ----
#------------------------------------------------#
#                                                #
#               VISUALIZE EXAMPLES               # 
#                                                #
#------------------------------------------------#

# Define a function for visualization
#----------------------#

visualize_mbChange <- function(plot_drug, plot_bug){
  
  plot_bug_clean = paste(ifelse(str_detect(plot_bug, "ref"), substr(plot_bug, 1, nchar(plot_bug)-21), substr(plot_bug, 1, nchar(plot_bug)-22)), 
                         " [", substr(plot_bug, nchar(plot_bug) - 5, nchar(plot_bug) - 1), "]", sep = "")
  
  EstMB_drugUsage <- readRDS("RData/Interim/Data_medications_currentlyUsed_ATC4.rds") %>% 
    tidyr::gather(ATC4, drugUsage, -skood)
  
  drugGroup <- EstMB_drugUsage %>% 
    dplyr::left_join(MBAD_drugInitiation_time_ATC4, by = c("skood", "ATC4")) %>% 
    dplyr::filter(ATC4 == plot_drug) %>% 
    dplyr::filter(skood %in% MBAD_countdata$skood)
  
  table(drugGroup$starter_group)
  
  if (sum(drugGroup$starter_group == "Stopper, active usage", na.rm = T) >= 10){
    plotdata <- EstMB_countdata %>% 
      dplyr::left_join(MBAD_countdata, by = c("skood", "taxa")) %>% 
      dplyr::filter(taxa == plot_bug) %>% 
      dplyr::mutate(delta = MBAD_value - EstMB_value) %>% 
      dplyr::left_join(drugGroup, by = c("skood")) %>%
      dplyr::filter(is.na(delta) == FALSE) %>% 
      dplyr::mutate(starter_group = ifelse(is.na(starter_group) == TRUE, "Non-user", starter_group)) %>% 
      dplyr::filter(starter_group %in% c("Non-user", "Initiator, past usage", "Initiator, active usage", "Stopper, active usage"))
  } else {
    plotdata <- EstMB_countdata %>% 
      dplyr::left_join(MBAD_countdata, by = c("skood", "taxa")) %>% 
      dplyr::filter(taxa == plot_bug) %>% 
      dplyr::mutate(delta = MBAD_value - EstMB_value) %>% 
      dplyr::left_join(drugGroup, by = c("skood")) %>%
      dplyr::filter(is.na(delta) == FALSE) %>% 
      dplyr::mutate(starter_group = ifelse(is.na(starter_group) == TRUE, "Non-user", starter_group)) %>% 
      dplyr::filter(starter_group %in% c("Non-user", "Initiator, past usage", "Initiator, active usage"))
  }
  
  p2 = ggplot(plotdata, 
              aes(x = factor(starter_group, levels = c("Non-user", "Initiator, past usage", "Initiator, active usage", "Stopper, active usage"),
                             labels = c("Non-user", "Initiator\npast usage", "Initiator\nactive usage", "Dicontinued\ndrug usage")),
                  y = delta, 
                  fill = factor(starter_group, levels = c("Non-user", "Initiator, past usage", "Initiator, active usage", "Stopper, active usage"), 
                                labels = c("Non-user", "Initiator\npast usage", "Initiator\nactive usage", "Dicontinued\ndrug usage")))) + 
    geom_hline(yintercept = 0, size = 1.5, linetype = 3, color = "gray50") + 
    geom_boxplot() +
    scale_fill_manual(values = c("gray70", "#c7d7e2", "#509aca", "yellowgreen"), guide = "none") + 
    geom_jitter(color = "gray50") + 
    xlab("") + 
    ylab("Change in abundance (T2-T1)") + 
    ggtitle(plot_bug_clean) + 
    theme_classic() + 
    theme(axis.text.y = element_text(size = 14), 
          legend.position = "bottom", 
          axis.text.x = element_text(size = 16), 
          axis.title = element_text(size = 18), 
          strip.text = element_text(size = 20), 
          title = element_text(size = 18))
  
  plot_name = paste("Figures/Drug_initiation/", plot_drug, "_", plot_bug, ".png", sep = "")
  plot_name_pdf = paste("Figures/Drug_initiation/", plot_drug, "_", plot_bug, ".pdf", sep = "")
  
  if (sum(drugGroup$starter_group == "Stopper, active usage", na.rm = T) >= 10){
    ggsave(p2, file = plot_name, width = 8, height = 7)
    ggsave(p2, file = plot_name_pdf, width = 8, height = 7)
  }  else {
    ggsave(p2, file = plot_name, width = 6, height = 7)
    ggsave(p2, file = plot_name_pdf, width = 6, height = 7)
  }
}


# Apply the function for all overlaping hits
#----------------------#
example_hits <- analysis_output_merged %>% 
  dplyr::filter(nchar(drug) == 5) %>% # ATC4
  dplyr::filter(Delta1_pval <= 0.05) 

for (i in 1:nrow(example_hits)){
  run_drug = example_hits[i, ] %>% dplyr::pull(drug)
  run_bug = example_hits[i, ] %>% dplyr::pull(MB_feature)
  
  visualize_mbChange(plot_drug = run_drug, plot_bug = run_bug)
}

