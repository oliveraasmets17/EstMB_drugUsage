



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


# Current usage data
#----------------------#
currentlyUsed_ATC3_data <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Data_medications_currentlyUsed_ATC3.rds")
currentlyUsed_ATC4_data <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Data_medications_currentlyUsed_ATC4.rds")
currentlyUsed_ATC5_data <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Data_medications_currentlyUsed_ATC5.rds")

currentlyUsed_ATC3_considered <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Medications_currentlyUsed_n20_ATC3.rds")
currentlyUsed_ATC4_considered <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Medications_currentlyUsed_n20_ATC4.rds")
currentlyUsed_ATC5_considered <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Medications_currentlyUsed_n20_ATC5.rds")

drug_combinations_ATC4 <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Data_drugCombinations_currentlyUsed_ATC4.rds")

previouslyUsed_ATC3_data <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Cumulative_usage_ATC3.rds") %>% dplyr::rename("ATC_used" = "ATC3")
previouslyUsed_ATC4_data <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Cumulative_usage_ATC4.rds") %>% dplyr::rename("ATC_used" = "ATC4")
previouslyUsed_ATC5_data <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Cumulative_usage_ATC5.rds") %>% dplyr::rename("ATC_used" = "ATC5")

last_used_ATC3_data <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Drug_last_used_ATC3.rds") %>% dplyr::rename("ATC_used" = "ATC3")
last_used_ATC4_data <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Drug_last_used_ATC4.rds") %>% dplyr::rename("ATC_used" = "ATC4")
last_used_ATC5_data <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Drug_last_used_ATC5.rds") %>% dplyr::rename("ATC_used" = "ATC5")
  
  
# ATC names
ATC_names <- readRDS("RData/Interim/DrugUsage_ATC_names.rds")


# Phenotype data
phenotype_data <- readRDS("C:/Users/oliver17/Desktop/Doktorantuur/Projekt_DrugUsage/RData/Interim/Data_master.rds")


# Raw drug data
MB_medications_data <- readRDS("RData/Interim/MB_medications_data_HIF.rds")







# Summarize drug usage in general ----
#------------------------------------------------#
#                                                #
#             SUMMARIZE DRUG USAGE               # 
#                                                #
#------------------------------------------------#

# Number of total drugs used in the cohort
#----------------------#
length(setdiff(colnames(currentlyUsed_ATC3_data), "skood")) # 126
length(setdiff(colnames(currentlyUsed_ATC4_data), "skood")) # 225
length(setdiff(colnames(currentlyUsed_ATC5_data), "skood")) # 433

length(unique(previouslyUsed_ATC3_data$ATC_used)) # 138
length(unique(previouslyUsed_ATC4_data$ATC_used)) # 251
length(unique(previouslyUsed_ATC5_data$ATC_used)) # 507



# Summarize medication currently usage by drug class
#----------------------#

# Prepare datasets
currentlyUsed_data <- currentlyUsed_ATC3_data %>% 
  dplyr::left_join(currentlyUsed_ATC4_data, by = "skood") %>% 
  dplyr::left_join(currentlyUsed_ATC5_data, by = "skood") %>% 
  tidyr::gather(drug, current_value, -skood)

previousUsage_data <- dplyr::bind_rows(previouslyUsed_ATC3_data, previouslyUsed_ATC4_data, previouslyUsed_ATC5_data) %>% 
  dplyr::mutate(previous_value = ifelse(n_prescription > 0, 1, 0), 
                drug = ATC_used) %>% 
  dplyr::select(skood, drug, previous_value, n_prescription)
  
# Calculate summary per person
med_currentUsage_summary <- currentlyUsed_data %>% 
  dplyr::full_join(previousUsage_data, by = c("skood", "drug")) %>% 
  dplyr::mutate(group = case_when(nchar(drug) == 3 ~ "ATC2",
                                  nchar(drug) == 4 ~ "ATC3", 
                                  nchar(drug) == 5 ~ "ATC4",
                                  TRUE ~ "ATC5")) %>% 
  dplyr::group_by(skood, group) %>% 
  dplyr::summarise(n_current = sum(current_value, na.rm = T),
                   n_previous = sum(previous_value, na.rm = T)) %>% 
  dplyr::ungroup()


# Number of different medications per subject on average
med_currentUsage_summary %>% 
  dplyr::group_by(group) %>% 
  dplyr::summarise(n0_current = sum(n_current == 0),
                   mean_current = mean(n_current), 
                   median_current = median(n_current), 
                   mean_users_current = mean(n_current[n_current != 0]),
                   median_users_current = median(n_current[n_current != 0]),
                   max_current = max(n_current),
                   # Previous usage
                   n0_previous = sum(n_previous == 0),
                   mean_previous = mean(n_previous), 
                   median_previous = median(n_previous), 
                   mean_users_previous = mean(n_previous[n_previous != 0]),
                   median_users_previous = median(n_current[n_previous != 0]),
                   max_previous = max(n_previous)) %>% 
  dplyr::ungroup()




# Most used drugs  
#----------------------#
used_drugs <- currentlyUsed_data %>% 
  dplyr::group_by(drug) %>% 
  dplyr::summarise(n = sum(current_value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(ATC_names, by = c("drug" = "drug_code")) %>% 
  dplyr::mutate(p = n/2509*100,
                N_users_p = paste(n, " (", round(p, 1), "%)", sep = ""), 
                group = case_when(nchar(drug) == 3 ~ "ATC2",
                                  nchar(drug) == 4 ~ "ATC3", 
                                  nchar(drug) == 5 ~ "ATC4",
                                  TRUE ~ "ATC5")) %>% 
  dplyr::filter(n >= 20) %>% 
  dplyr::arrange(group, desc(n))

used_drugs_historical <- previousUsage_data %>% 
  dplyr::group_by(drug) %>% 
  dplyr::summarise(n_5y = sum(previous_value), 
                   mean_n_prescriptions = round(mean(n_prescription), 1),
                   sd_n_prescriptions = round(sd(n_prescription), 2)) %>% 
  dplyr::ungroup()
  
used_drugs_summaryInfo <- used_drugs %>% 
  dplyr::left_join(used_drugs_historical, by = c("drug")) %>% 
  dplyr::mutate(p_5y = paste(n_5y, " (", round(n_5y/2509*100, 1), "%)", sep = ""),
                prescr_print = paste(mean_n_prescriptions, " (", sd_n_prescriptions, ")", sep = "")) %>% 
  dplyr::arrange(group, desc(n)) %>% 
  dplyr::select(group, drug_name, N_users_p, p_5y, prescr_print) %>% 
  dplyr::rename("Drug name (ATC code)" = "drug_name", 
                "Number of active drug users (%)" = "N_users_p", 
                "Number of drug users in 5y (%)" = "p_5y", 
                "Number of prescriptions in 5y (mean (sd))" = "prescr_print") %>% 
  as.data.frame() 

xlsx::write.xlsx(x = used_drugs_summaryInfo, 
                 file = "Results/Extended Data Table 1. Number of drug users.xlsx", 
                 sheetName = "Table1",
                 col.names = TRUE, 
                 row.names = FALSE,
                 append = FALSE)




# Visualize total drug usage in the cohort
plotData1 <- data.frame(current = c(length(setdiff(colnames(currentlyUsed_ATC3_data), "skood")), 
                                    length(setdiff(colnames(currentlyUsed_ATC4_data), "skood")), 
                                    length(setdiff(colnames(currentlyUsed_ATC5_data), "skood"))), 
                        previous = c(length(unique(previouslyUsed_ATC3_data$ATC_used)), 
                                     length(unique(previouslyUsed_ATC4_data$ATC_used)), 
                                     length(unique(previouslyUsed_ATC5_data$ATC_used))), 
                        users = c(56, 63, 67), 
                        ATC = c("ATC3", "ATC4", "ATC5")) %>% 
  tidyr::gather(key, value, -ATC)

p_totalUsers <- ggplot(plotData1, aes(x = ATC, y = value, fill = factor(key, levels = c("previous", "current", "users"), 
                                                        labels = c("# drugs used within 5y preceding T1", "# drugs used at T1", "# drugs with at least 20 users at T1")))) + 
  geom_bar(stat = "identity", position = "dodge", color = "gray90") + 
  geom_text(aes(label = value, y = value + 10), position = position_dodge(width = 0.9), size = 6) + 
  xlab("") + 
  ylab("Number of distinct drugs") + 
  scale_fill_manual(values = c("#feb24c", "#fd8d3c", "#fc4e2a"), name = "") + 
  theme_classic() + 
  theme(axis.text = element_text(size = 18), 
        panel.grid.major.x = element_blank(),
        axis.title = element_text(size = 18), 
        legend.position = c(0.3, 0.9),
        legend.text = element_text(size = 14))
p_totalUsers
                                                                        
ggsave(p_totalUsers, file = "Figures/Population_summary/p_totalUsers.png", height = 6, width = 8)
ggsave(p_totalUsers, file = "Figures/Population_summary/p_totalUsers.pdf", height = 6, width = 8)







# Visualize drug usage in general ----
#------------------------------------------------#
#                                                #
#             VISUALIZE DRUG USAGE               # 
#                                                #
#------------------------------------------------#

# Visualize the number of medication usage by subjects
#----------------------#
plotData_drugPerSubject <- currentlyUsed_data %>% 
  dplyr::mutate(group = case_when(nchar(drug) == 3 ~ "ATC2",
                                  nchar(drug) == 4 ~ "ATC3", 
                                  nchar(drug) == 5 ~ "ATC4",
                                  TRUE ~ "ATC5")) %>% 
  dplyr::group_by(skood, group) %>% 
  dplyr::summarise(n = sum(current_value)) %>% 
  dplyr::ungroup()

# Non-users
n0 = plotData_drugPerSubject %>% dplyr::filter(group == "ATC4" & n == 0) %>% nrow()

# Visualize
drug_summary_plot <- ggplot(plotData_drugPerSubject, 
                                 aes(x = n, fill = group)) + 
  geom_bar(position = "dodge") +
  xlab("Number of unique drugs taken at T1") + 
  ylab("Number of individuals") + 
  scale_fill_manual(name = "", values = c("gray30", "gray60", "gray80"))  + 
  scale_y_continuous(breaks = seq(100, 1500, 100), labels = seq(100, 1500, 100)) + 
  scale_x_continuous(breaks = seq(0, 16, 2), limits = c(-0.5, 16)) + 
  theme_classic() + 
  annotate("text", x = 0, y = n0 + 50, label = n0, size = 7) + 
  theme(legend.position = c(0.8, 0.8),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20))

ggsave(drug_summary_plot, filename = "Figures/Population_summary/TotalNr_medications_used_histogram.png", height = 6, width = 8)
ggsave(drug_summary_plot, filename = "Figures/Population_summary/TotalNr_medications_used_histogram.pdf", height = 6, width = 8)



# Visualize most used drugs and combinations

#----------------------#
# Prepare data
used_drugs_help <- used_drugs %>% 
  dplyr::filter(group == "ATC4") %>% 
  dplyr::rename("one" = drug) %>% 
  dplyr::mutate(two = one) %>% 
  dplyr::select(one, two, n)

drug_combination_df <- drug_combinations_ATC4 %>% 
  dplyr::rename("two" = "one",
                "one" = "two") %>% 
  dplyr::bind_rows(drug_combinations_ATC4, used_drugs_help) %>% 
  dplyr::left_join(ATC_names, by = c("two" = "drug_code"))

# Drug order
drug_order <- drug_combination_df %>% 
  dplyr::filter(one == two) %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::filter(n >= 50)

# Visualize the drug combinations
top_drugCombinations_plot <- ggplot(drug_combination_df %>% dplyr::filter(one %in% drug_order$one & two %in% drug_order$two), 
                                    aes(x = factor(one, levels = drug_order$one),  
                                        y = factor(drug_name, levels = drug_order$drug_name), 
                                        fill = log(n))) + 
  geom_tile(color = "white") + 
  geom_text(aes(label = n)) + 
  scale_fill_gradient2_tableau(trans = "reverse", palette = "Orange-Blue-White Diverging") + 
  theme_classic() + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        legend.position = "none")

top_drugCombinations_plot

ggsave(top_drugCombinations_plot, filename = "Figures/Population_summary/Top_drugCombinations_heatmap.png", height = 8, width = 15)
ggsave(top_drugCombinations_plot, filename = "Figures/Population_summary/Top_drugCombinations_heatmap.pdf", height = 8, width = 15)



# Visualize the time drugs were last used (ATC4) 
#----------------------#
plotData_usageHistory <- last_used_ATC4_data %>% 
  dplyr::filter(ATC_used != "J07BB") %>% 
  dplyr::group_by(ATC_used, timeGroup) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(is.na(timeGroup) == FALSE & timeGroup != "<1y") %>% 
  dplyr::left_join(ATC_names, by = c("ATC_used" = "drug_code")) %>% 
  dplyr::filter(ATC_used %in% currentlyUsed_ATC4_considered)
  
drug_order <- plotData_usageHistory %>% 
  dplyr::filter(timeGroup == "Current_user") %>% 
  dplyr::arrange(n) %>% 
  dplyr::filter(n >= 50) %>% 
  dplyr::pull(drug_name)

plot_usageHistory <- ggplot(plotData_usageHistory %>% dplyr::filter(drug_name %in% drug_order), 
                            aes(x = n, 
                                y = factor(drug_name, levels = drug_order), 
                                fill = factor(timeGroup, 
                                              levels = rev(c("Current_user", "1y-2y", "2y-3y", "3y-4y", "4y-5y")),
                                              labels = rev(c("Current user", "1y-2y ago", "2y-3y ago", "3y-4y ago", "4y-5y ago"))))) + 
  geom_bar(stat = "identity", position = "stack", color = "black") + 
  scale_fill_manual(name = "Time last used", values = c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494")) + 
  scale_x_continuous(breaks = seq(0, 700, 100)) +
  xlab("Number of drug users") + 
  ylab("") + 
  theme_classic() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 18), 
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16), 
        panel.grid.major.x = element_line())
plot_usageHistory

ggsave(plot = plot_usageHistory, file = "Figures/Population_summary/ATC4_drugLastUsed_stackedBarplot.png", width = 11, height = 10)
ggsave(plot = plot_usageHistory, file = "Figures/Population_summary/ATC4_drugLastUsed_stackedBarplot.pdf", width = 11, height = 10)



