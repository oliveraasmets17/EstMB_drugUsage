

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load the packages
library("tidyverse")
library("ggthemes")


# Diversity indexes based on raw data
diversity_data <- readRDS("RData/Interim/mOTUs_alphaDiversity_unfiltered.rds")


# PB ratio
PB_ratio <- readRDS("RData/Interim/PB_ratio_mOTUs.rds") %>% 
  dplyr::select(skood, PB_ratio)


# Phenotype data
phenotype_data <- readRDS("RData/Interim/Data_master.rds")


# Raw drug data
MB_medications_data <- readRDS("RData/Interim/MB_medications_data_HIF.rds") %>% 
  dplyr::filter(Time_from_MBsample <= 1826.25)


# Cell counts
cellcount_data <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_Others/Kuhn/RData/Predicted_cellcounts.rds")









# Number of different drug used ----
#------------------------------------------------#
#                                                #
#        EFFECT OF NUMBER OF DRUGS USED          # 
#                                                #
#------------------------------------------------#

# Read the data
n_drugs_used_currently <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Data_medications_currentlyUsed_ATC5.rds") %>% 
  tidyr::gather(key, value, -skood) %>% 
  dplyr::group_by(skood) %>% 
  dplyr::summarise(n_drugs = sum(value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(PB_ratio, by = "skood") %>% 
  dplyr::left_join(diversity_data, by = "skood") %>% 
  dplyr::left_join(cellcount_data, by = "skood") %>%  
  dplyr::left_join(phenotype_data[ ,c("skood", "Age_at_MBsample")], by = c("skood")) %>% 
  tidyr::gather(key, value, -c("skood", "n_drugs", "Age_at_MBsample")) %>% 
  dplyr::mutate(age_group = cut(Age_at_MBsample, breaks = c(-Inf, 40, 60, Inf), labels = c("<40", "40-60", ">60")),
                n_drugs_group = ifelse(n_drugs <= 10, as.character(n_drugs), ">10"))


median <- n_drugs_used_currently %>% 
  dplyr::filter(key %in% c("observed", "shannon")) %>% 
  dplyr::filter(n_drugs_group == 0) %>% 
  dplyr::group_by(key) %>% 
  dplyr::summarise(median = median(value, na.rm = T)) %>% 
  dplyr::ungroup()


p_nDrugs <- ggplot(n_drugs_used_currently %>% dplyr::filter(key %in% c("observed", "shannon")), 
                   aes(x = factor(n_drugs_group, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10")), y = value, fill = ifelse(n_drugs_group == 0, "A", "B"))) + 
  geom_hline(data = median, aes(yintercept = median), linetype = 3, linewidth = 1, color = "gray50") + 
  geom_boxplot() + 
  scale_fill_manual(values = c("#FF95A8FF", "#84D7E1FF", "#008EA0FF"), name = "Usage of antibiotics") + 
  facet_wrap(vars(factor(key, levels = c("observed", "shannon", "PB_ratio", "cellcount"), 
                         labels = c("Observed richness", "Shannon diversity", "Prevotella-Bacteroides ratio", "Cell counts"))), ncol = 2, scales = "free") + 
  xlab("Number of unique drugs used at T1") + 
  ylab("") + 
  theme_classic() + 
  theme(strip.text = element_text(size = 18), 
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        legend.position = "bottom", 
        axis.text = element_text(size = 12))
p_nDrugs
ggsave(p_nDrugs, filename = "Figures/Diversity_nDrugs_scatterPlot_ATC5.png", height = 5, width = 14)
ggsave(p_nDrugs, filename = "Figures/Diversity_nDrugs_scatterPlot_ATC5.pdf", height = 5, width = 14)




p_nDrugs_observed <- ggplot(n_drugs_used_currently %>% dplyr::filter(key == "observed"), aes(x = n_drugs, y = value, color = age_group)) + 
  geom_jitter() + 
  geom_smooth(method = "lm", linewidth = 1.5) + 
  scale_color_manual(values = c("#a6cee3", "#509aca", "#69b42f"), name = "Age group") + 
  xlab("Number of unique drugs used at the time of sampling") + 
  ylab("Observed richness") + 
  theme_bw() + 
  theme(strip.text = element_text(size = 18), 
        strip.background = element_rect(fill = "white"), 
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        legend.position = "bottom", 
        axis.text = element_text(size = 12))

p_nDrugs_observed
ggsave(p_nDrugs_observed, filename = "Figures/Diversity_nDrugs_scatterPlot_ATC5_observedRichness.png", height = 7, width = 12)
ggsave(p_nDrugs_observed, filename = "Figures/Diversity_nDrugs_scatterPlot_ATC5_observedRichness.pdf", height = 7, width = 12)









# Time from last antibiotics usage ----
#------------------------------------------------#
#                                                #
#          VISUALIZE TIME FROM LAST AB           # 
#                                                #
#------------------------------------------------#

# AB last used
n_AB_used <- MB_medications_data %>% 
  dplyr::filter(Time_from_MBsample > 0) %>% 
  dplyr::filter(substr(ATC3, 1, 3) == "J01") %>% 
  dplyr::group_by(skood) %>% 
  dplyr::summarize(n_AB = n()) %>% 
  dplyr::ungroup()

AB_last_used <- MB_medications_data %>% 
  dplyr::filter(Time_from_MBsample > 0) %>% 
  dplyr::filter(substr(ATC3, 1, 3) == "J01") %>% 
  dplyr::group_by(skood) %>% 
  dplyr::mutate(min_time = min(Time_from_MBsample)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(min_time == Time_from_MBsample) %>% 
  dplyr::distinct(skood, Time_from_MBsample) %>% 
  dplyr::full_join(diversity_data, by = "skood") %>% 
  dplyr::left_join(cellcount_data, by = "skood") %>%  
  dplyr::left_join(phenotype_data[ ,c("skood", "Age_at_MBsample")], by = c("skood")) %>% 
  dplyr::left_join(n_AB_used, by = "skood") %>% 
  dplyr::mutate(age_group = cut(Age_at_MBsample, breaks = c(-Inf, 40, 60, Inf), labels = c("<40", "40-60", ">60")), 
                AB_group = cut(n_AB, breaks = c(-Inf, 1, Inf)))


# Nonusers averages 
nonusers <- AB_last_used %>% 
  dplyr::filter(is.na(Time_from_MBsample) == T) %>% 
  dplyr::summarize(n = n(), 
                   mean_observed = mean(observed),
                   median_observed = median(observed),
                   sd_observed = sd(observed),
                   int_observed = 1.96*sd(observed)/sqrt(n()),
                   mean_shannon = mean(shannon), 
                   mean_cellcount = mean(cellcount))

# Plot 
p_observed_ABusage <- ggplot(AB_last_used %>% dplyr::filter(is.na(Time_from_MBsample) == F),
                             aes(x = Time_from_MBsample, y = observed, color = AB_group)) + 
  geom_jitter(size = 2, color = "gray84") +
  geom_vline(xintercept = c(180), linetype = 3, linewidth = 1, color = "black") + 
  geom_hline(yintercept = nonusers$mean_observed, linewidth = 1.5, color = "#ce1256") + 
  annotate("label", x = 180, y = 400, label = "6 months", size = 7) + 
  annotate("text", x = 365*4 + 210, y = 272, label = "Non-users", size = 5, color = "#ce1256") + 
  geom_smooth(linewidth = 1.5, se = FALSE) + 
  scale_x_continuous(breaks = c(365, 365*2, 365*3, 365*4, 365*5), labels = c("1y", "2y", "3y", "4y", "5y")) + 
  xlab("Time since last antibiotics treatment") + 
  ylab("Observed richness") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 14))

p_observed_ABusage
ggsave(p_observed_ABusage, filename = "Figures/Time_from_lastAB_observed.png", height = 6, width = 8)
ggsave(p_observed_ABusage, filename = "Figures/Time_from_lastAB_observed.pdf", height = 6, width = 8)



# Grouped plot
AB_last_used_grouped <- AB_last_used %>% 
  dplyr::mutate(Time_from_MBsample = ifelse(is.na(Time_from_MBsample) == TRUE, -1, as.numeric(Time_from_MBsample))) %>%
  dplyr::mutate(timegroup = cut(Time_from_MBsample, breaks = c(-Inf, seq(0, 5*365 - 182.5, 182.5), Inf),
                                labels = c("Non-user", "0-6m", "6m-1y", "1y-1.5y", "1.5y-2y", "2y-2.5y", "2.5y-3y", "3y-3.5y", "3.5y-4y", "4y-4.5y", "4.5y-5y"))) %>% 
  dplyr::mutate(AB_group = case_when(is.na(n_AB) == TRUE ~ "Non-user",
                                     n_AB == 1 ~ "1 course in 5 years",
                                     n_AB > 1 ~ ">1 courses in 5 years", 
                                     TRUE ~ "Something else"), 
                AB_group = factor(AB_group, levels = c("Non-user", "1 course in 5 years", ">1 courses in 5 years")))

my_comparisons <- list(c("Non-user", "0-6m"), 
                       c("Non-user", "6m-1y"), 
                       c("Non-user", "1y-1.5y"), 
                       c("Non-user", "1.5y-2y"), 
                       c("Non-user", "2y-2.5y"),
                       c("Non-user", "2.5y-3y"),
                       c("Non-user", "3.5y-4y"),
                       c("Non-user", "4.5y-5y"))

p_observed_ABusage_boxplot <- ggplot(AB_last_used_grouped %>% mutate(AB_group = ifelse(AB_group == "Non-user", "Non-user", "User")), # %>% dplyr::filter(n_AB <= 2 | is.na(n_AB) == TRUE),
                                        aes(x = timegroup, y = observed, fill = AB_group)) + 
  geom_hline(yintercept = nonusers$median_observed, linewidth = 1, color = "gray50", linetype = 3) + 
  geom_boxplot() + 
  scale_fill_manual(name = "", values = c("#FF95A8FF", "#008EA0FF")) + 
  xlab("Time since last antibiotics treatment") + 
  ylab("Observed richness") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 14), 
        legend.title = element_blank(),
        legend.position = c(0.84, 0.8),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(size = 16)) + 
  stat_compare_means(comparisons = my_comparisons, tip.length = 0.01, step.increase = 0.07)
p_observed_ABusage_boxplot

ggsave(p_observed_ABusage_boxplot, filename = "Figures/Time_from_lastAB_observed_boxplot_2classes.png", height = 7, width = 10)
ggsave(p_observed_ABusage_boxplot, filename = "Figures/Time_from_lastAB_observed_boxplot_2classes.pdf", height = 7, width = 10)

AB_last_used_grouped %>% 
  mutate(AB_group = ifelse(AB_group == "Non-user", "Non-user", "User")) %>% 
  dplyr::group_by(AB_group, timegroup) %>% 
  dplyr::summarise(n = n())
