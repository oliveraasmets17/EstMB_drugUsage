

# Command line arguments ----
#------------------------------------------------#
#                                                #
#          SET UP COMMAND LINE ARGUMENTS         # 
#                                                #
#------------------------------------------------#

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)              
# Expected input:
# 1) Input set, possible outputs SET0, SET1, ..., SET4

# NB! For testing - remove
# args = c("1", "BRAY")
cat(args, sep = "\n")


# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 

# Declare phenotype set 
input_set = args[1]


# Declare distance matrix
distance_mat_name = args[2]







# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load the packages
library("vegan")
library("stringr")
library("dplyr")

# Phenotype data
old_drugs <- readRDS("RData/Interim/Medication_factors_analyzed.rds")

phenotype_data <- readRDS("RData/Interim/Data_master.rds") %>% 
  dplyr::select(-one_of(old_drugs))


# Drug usage data
currentlyUsed_ATC4 <- readRDS("RData/Interim/Data_medications_currentlyUsed_ATC4.rds")
drugs_considered <- readRDS("RData/Interim/Medications_currentlyUsed_n20_ATC4.rds")

Cumulative_drugUsage_raw <- readRDS("RData/Interim/Cumulative_usage_ATC4.rds") %>% 
  dplyr::filter(ATC4 %in% drugs_considered) %>% 
  dplyr::select(skood, ATC4, n_prescription) %>% 
  tidyr::spread(ATC4, n_prescription)


# Proprocess cumulative usage data
Cumulative_drugUsage_cont <- phenotype_data %>% 
  dplyr::select(skood) %>% 
  dplyr::full_join(Cumulative_drugUsage_raw, by = c("skood")) %>% 
  tidyr::gather(key, value, -skood) %>% 
  dplyr::mutate(value = ifelse(is.na(value), 0, value)) %>% 
  tidyr::spread(key, value)

colnames(Cumulative_drugUsage_cont) <- c("skood", paste("Cumulative_", colnames(Cumulative_drugUsage_cont)[2:ncol(Cumulative_drugUsage_cont)], sep = ""))


Cumulative_drugUsage_binary <- phenotype_data %>% 
  dplyr::select(skood) %>% 
  dplyr::full_join(Cumulative_drugUsage_raw, by = c("skood")) %>% 
  tidyr::gather(key, value, -skood) %>% 
  dplyr::mutate(value = ifelse(is.na(value), 0, value), 
                value = ifelse(value == 0, 0, 1)) %>% 
  tidyr::spread(key, value)

colnames(Cumulative_drugUsage_binary) <- c("skood", paste("Cumulative_", colnames(Cumulative_drugUsage_binary)[2:ncol(Cumulative_drugUsage_binary)], sep = ""))




# Factor sets to analyze
diseases <- readRDS("RData/Interim/Disease_factors_analyzed.rds") %>% setdiff(c("health_status_ok", "mental_health_status_ok"))

other_factors <- setdiff(c(readRDS("RData/Interim/Other_factors_analyzed.rds"), "education"),
                         c("has_smoked_lastYear", "alcohol_frequency_category", "antidepressants_history_cont", 
                           "antibiotics_history_quartile", "antidepressants_history_quartile", "antibiotics_history_cont"))

intrinsic_factors <- readRDS("RData/Interim/Intrinsic_factors_analyzed.rds")

procedure_factors <- readRDS("RData/Interim/Procedure_factors_analyzed.rds")

dietary_factors <- readRDS("RData/Interim/Dietary_factors_analyzed.rds")

drugHistory_factors <- setdiff(colnames(Cumulative_drugUsage_binary), "skood")








# Data preprocessing ----
#------------------------------------------------#
#                                                #
#               DATA PREPROCESSING               # 
#                                                #
#------------------------------------------------#

if (distance_mat_name == "BRAY"){
  # Bray curtis
  distance_mat_raw <- readRDS("RData/Interim/Bray_dist_mOTUs_filtered.rds")
  
}else if (distance_mat_name == "CLR"){
  # Aitchison distance
  distance_mat_raw <- readRDS("RData/Interim/CLR_dist_mOTUs_filtered.rds")
  
} else if (distance_mat_name == "JACCARD"){
  # Aitchison distance
  distance_mat_raw <- readRDS("RData/Interim/Jaccard_dist_mOTUs_filtered.rds")
}

# Preprocess metadata
metadata_analysisReady_cont <- phenotype_data %>% 
  dplyr::left_join(currentlyUsed_ATC4, by = "skood") %>%
  dplyr::left_join(Cumulative_drugUsage_cont, by = "skood") %>%
  dplyr::filter(skood %in% rownames(as.matrix(distance_mat_raw))) %>% 
  dplyr::mutate_at(vars(all_of(drugHistory_factors)), ~tidyr::replace_na(., 0)) %>% # Replace NA values for diagnosis
  tibble::column_to_rownames(var = "skood")


metadata_analysisReady_binary <- phenotype_data %>% 
  dplyr::left_join(currentlyUsed_ATC4, by = "skood") %>%
  dplyr::left_join(Cumulative_drugUsage_binary, by = "skood") %>%
  dplyr::filter(skood %in% rownames(as.matrix(distance_mat_raw))) %>% 
  dplyr::mutate_at(vars(all_of(drugHistory_factors)), ~tidyr::replace_na(., 0)) %>% # Replace NA values for diagnosis
  tibble::column_to_rownames(var = "skood")


# Subset metadata to contain only variables of interest
if (input_set == "SET01"){
  
  metadata_analysisReady_subset <- metadata_analysisReady_binary %>% 
    dplyr::select(all_of(diseases), all_of(drugs_considered), all_of(other_factors), 
                  all_of(intrinsic_factors), all_of(procedure_factors), 
                  all_of(dietary_factors), all_of(drugHistory_factors)) %>% 
    dplyr::filter(complete.cases(.))
  
} else if (input_set == "SET02"){
  
  metadata_analysisReady_subset <- metadata_analysisReady_cont %>% 
    dplyr::select(all_of(diseases), all_of(drugs_considered), all_of(other_factors), 
                  all_of(intrinsic_factors), all_of(procedure_factors), 
                  all_of(dietary_factors), all_of(drugHistory_factors)) %>% 
    dplyr::filter(complete.cases(.))
  
} else if (input_set == "SET1"){
  
  metadata_analysisReady_subset <- metadata_analysisReady_cont %>% 
    dplyr::select(all_of(diseases), all_of(drugs_considered), all_of(other_factors), 
                  all_of(intrinsic_factors), all_of(procedure_factors), 
                  all_of(dietary_factors)) %>% 
    dplyr::filter(complete.cases(.))
  
} else if (input_set == "SET2"){
  
  metadata_analysisReady_subset <- metadata_analysisReady_cont %>% 
    dplyr::select(all_of(diseases), all_of(drugs_considered)) %>% 
    dplyr::filter(complete.cases(.))
  
} else if (input_set == "SET3"){
  
  metadata_analysisReady_subset <- metadata_analysisReady_cont %>% 
    dplyr::select(all_of(diseases), all_of(drugs_considered), all_of(drugHistory_factors)) %>% 
    dplyr::filter(complete.cases(.))
  
} else if (input_set == "SET41"){
  
  metadata_analysisReady_subset <- metadata_analysisReady_binary %>% 
    dplyr::select(all_of(diseases), all_of(drugs_considered), all_of(drugHistory_factors), all_of(intrinsic_factors)) %>% 
    dplyr::filter(complete.cases(.))
} else if (input_set == "SET42"){
  
  metadata_analysisReady_subset <- metadata_analysisReady_cont %>% 
    dplyr::select(all_of(diseases), all_of(drugs_considered), all_of(drugHistory_factors), all_of(intrinsic_factors)) %>% 
    dplyr::filter(complete.cases(.))
}

# Reorder and subset distance mat
distance_mat_subset <- as.matrix(distance_mat_raw)
distance_mat_subset <- distance_mat_subset[rownames(metadata_analysisReady_subset), rownames(metadata_analysisReady_subset)]

distance_mat <- as.dist(distance_mat_subset)







# Modelling  ----
#------------------------------------------------#
#                                                #
#                   MODELLING                    # 
#                                                #
#------------------------------------------------#

# --------------- Select full model --------------

# Null-model
m0 <- dbrda(distance_mat ~ 1, metadata_analysisReady_subset)

# All-variable model
analysis_variables <- colnames(metadata_analysisReady_subset)
m1 <- dbrda(formula(paste("distance_mat ~ `", paste(analysis_variables, collapse = "` + `"), "`", sep = "")), metadata_analysisReady_subset)

# Model selection
m_stepwise <- ordistep(m0, scope = formula(m1), direction = "forward", trace = F, parallel = 30)

# Save outputs
selected_variables <- as.character(formula(m_stepwise))[3]
selected_variables <- str_split(selected_variables, " \\+ ")[[1]]


if (distance_mat_name == "BRAY"){
  # Bray curtis
  saveRDS(m_stepwise, paste("RData/Results/VariancePartitioning_BRAY_forward_", input_set, "_FULL.rds", sep = ""))
  
}else if (distance_mat_name == "CLR"){
  # Aitchison distance
  saveRDS(m_stepwise, paste("RData/Results/VariancePartitioning_CLR_forward_", input_set, "_FULL.rds", sep = ""))
  
}else if (distance_mat_name == "JACCARD"){
  # Aitchison distance
  saveRDS(m_stepwise, paste("RData/Results/VariancePartitioning_JACCARD_forward_", input_set, "_FULL.rds", sep = ""))
}




# ----------- Define models to analyze ----------
# Modified from https://git.embl.de/grp-bork/vpthemall

# Add main index.to.test data
model_def_df <- data.frame(var = "full",
                           element = "full",
                           var.formula = paste(selected_variables, collapse = " + "),
                           stringsAsFactors = F)

# Make index for deconfounding all variables in the model
for (var in selected_variables){
  
  # Get other variables in the vector
  remaining_var <- selected_variables[selected_variables != var]
  
  # Get formulas for X1 | X2 and X1  (a and ab) - unique
  f_controlled <- paste0(var, " + Condition(", paste(remaining_var, collapse = " + "),")")
  
  # Add controlled to data.frame
  model_def_unique <- data.frame(var = var,
                                 element = "unique",
                                 var.formula = f_controlled,
                                 stringsAsFactors = F)

  model_def_df = dplyr::bind_rows(model_def_df, model_def_unique)
}



# ----------- Define models to analyze ----------
model_def_df$p.value <- NULL
model_def_df$r.squared <- NULL
model_def_df$adj.r.squared <- NULL

for (i in 1:nrow(model_def_df)){
  
  # Set defaults
  n_permutations = 999
  
  # Make formula 
  f_run <- as.formula(paste("distance_mat ~", model_def_df$var.formula[i]))
  
  # Build dbRDA model
  m_run <- dbrda(f_run, metadata_analysisReady_subset)
  
  # Perform ANOVA type III
  anova_run <- anova(m_run, by = "margin", permutations = n_permutations)
  
  if(model_def_df$element[i] == "full"){
    p.value <- NA
    f.value <- NA
  } else{
    p.value <- anova_run$`Pr(>F)`[1]
    f.value <-  anova_run$F[1]
  }
  
  # Get adjusted R squared
  r2 <- RsquareAdj(m_run)
  
  # Add results to results dataset
  model_def_df$p.value[i] <- p.value
  model_def_df$adj.r.squared[i] <- r2$adj.r.squared
  model_def_df$r.squared[i] <- r2$r.squared
}



if (distance_mat_name == "BRAY"){
  # Bray curtis
  saveRDS(model_def_df, paste("RData/Results/VariancePartitioning_BRAY_forward_", input_set, ".rds", sep = ""))
  
}else if (distance_mat_name == "CLR"){
  # Aitchison distance
  saveRDS(model_def_df, paste("RData/Results/VariancePartitioning_CLR_forward_", input_set, ".rds", sep = ""))
  
}else if (distance_mat_name == "JACCARD"){
  # Aitchison distance
  saveRDS(model_def_df, paste("RData/Results/VariancePartitioning_JACCARD_forward_", input_set, ".rds", sep = ""))
}
