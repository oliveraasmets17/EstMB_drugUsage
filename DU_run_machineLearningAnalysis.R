

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load the packages
library("dplyr")
library("tidyr")
library("compositions")
library("tidymodels")

# Abundance data
countData_raw <- readRDS("RData/Interim/mOTUs_CLR_prevalence10p.rds") 


# Phenotype data
phenotype_data <- readRDS("RData/Interim/Data_master.rds")


# Actively used ATC4 level drugs
currentlyUsed_ATC4 <- readRDS("RData/Interim/Data_medications_currentlyUsed_ATC4.rds") %>% 
  tidyr::gather(drug, current_usage, -skood)







# ML modelling ----
#------------------------------------------------#
#                                                #
#               TIDYMODELS FUNCTION              # 
#                                                #
#------------------------------------------------#

# Preprocess medication data 
# ------------------------------#
recentUsage_data <- readRDS("RData/Interim/Drug_last_used_ATC4.rds")

recentUsage_analysis_data <- phenotype_data %>% 
  # Cover all subject - ATC pairs
  dplyr::select(skood) %>% 
  dplyr::full_join(recentUsage_data, by = "skood") %>% 
  dplyr::select(skood, ATC4) %>% 
  tidyr::expand(skood, ATC4) %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::left_join(recentUsage_data, by = c("skood", "ATC4")) %>% 
  # Define target variable
  dplyr::mutate(timeGroup = ifelse(is.na(timeGroup), "non_user", as.character(timeGroup))) %>% 
  dplyr::select(skood, ATC4, timeGroup) %>% 
  dplyr::mutate(help = 1) %>% 
  tidyr::spread(timeGroup, help)  %>% 
  tidyr::gather(key, value, -c("skood", "ATC4", "non_user")) %>% 
  dplyr::filter(key == "Current_user") %>% 
  dplyr::mutate(target = case_when(non_user == 1 ~ 0, 
                                   value == 1 ~ 1, 
                                   TRUE ~ as.numeric(NA)), 
                target = factor(target)) %>% 
  dplyr::select(skood, ATC4, target) %>% 
  dplyr::filter(skood %in% countData_raw$skood) %>% 
  dplyr::left_join(countData_raw, by = "skood")



# Define the function to run ML
# ------------------------------#
run_drug_ml = function(drug, p_seed){
  
  # ML data 
  ml_data <- recentUsage_analysis_data %>%
    dplyr::filter(ATC4 == drug) %>% 
    dplyr::select(-one_of("ATC4"))
  
  ml_data_unclean <- ml_data %>% 
    dplyr::filter(is.na(target) == T)
  
  ml_data_clean <- ml_data %>% 
    dplyr::filter(is.na(target) == F)
  
  
  # Initial data split ----
  #------------------------------------------------#
  
  # Create random data split
  set.seed(p_seed)
  data_split <- rsample::initial_split(ml_data_clean, prop = .75, strata = target) # stratified sampling
  
  # Split the data
  data_train_raw <- rsample::training(data_split)
  data_test_raw <- rsample::testing(data_split) %>% 
    dplyr::bind_rows(ml_data_unclean)
  
  data_train <- data_train_raw %>% dplyr::select(-one_of("skood"))
  data_test <- data_test_raw %>% dplyr::select(-one_of("skood"))
  
  # Create CV sets
  set.seed(p_seed)
  data_cv <- rsample::vfold_cv(data_train, v = 5, repeats = 4, strata = target) # stratified sampling
  
  
  # Preprocessing steps for the data ----
  #------------------------------------------------#
  
  # Data Preprocessing
  data_preprocessing <- 
    # Specify the outcome variable and predictor variables
    recipes::recipe(target ~ ., data = data_train) %>%
    recipes::step_dummy(all_nominal(), -all_outcomes()) %>%
    # Apply preprocessing steps to the data (either one by one or for numeric/character etc)
    recipes::step_normalize(all_predictors()) %>%
    recipes::step_nzv(all_predictors())
  
  
  # Define the model used ----
  #------------------------------------------------#
  
  # Define the model - no training done yet
  model <- parsnip::logistic_reg() %>%
    parsnip::set_args(
      penalty = tune::tune(),  
      mixture = tune::tune()) %>%
    set_engine("glmnet") %>% 
    set_mode("classification")
  
  # Construct a workflow that combines your recipe and your model
  ml_workflow <- workflows::workflow() %>%
    workflows::add_recipe(data_preprocessing) %>%
    workflows::add_model(model)
  
  
  # Start tuning the model ----
  #------------------------------------------------#
  
  # Tune the models
  set.seed(p_seed)
  tuned_model <- ml_workflow %>%
    tune::tune_grid(resamples = data_cv,
                    grid = 50, 
                    metrics = yardstick::metric_set(yardstick::roc_auc))
  
  
  # Return best parameters
  best_params <- tuned_model %>%
    tune::select_best(metric = "roc_auc")
  
  
  # Evaluate model performance ----
  #------------------------------------------------#
  
  # Describe the best models according to the best parameters gathered earlier
  best_workflow <- ml_workflow %>%
    tune::finalize_workflow(best_params) 
  
  
  # Fit the models
  model_final <- best_workflow %>%
    parsnip::fit(data_train)
  

  # Output ----
  #------------------------------------------------#
  
  # Model predictors
  predictors <- model_final %>%
    pull_workflow_fit() %>%
    tidy() %>%
    dplyr::filter(estimate != 0)
  
  return(list(model_final, data_test_raw, predictors))
  
}









# Run function for several drugs ----
#------------------------------------------------#
#                                                #
#                  RUN FUNCTION                  # 
#                                                #
#------------------------------------------------#

# antibiotics with >= 50 cases
ML_drugs <- currentlyUsed_ATC4 %>% 
  dplyr::group_by(drug) %>% 
  dplyr::summarise(n = sum(current_usage)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(n >= 50) %>% 
  dplyr::filter(substr(drug, 1, 3) == "J01") %>% 
  dplyr::pull(drug)

for (i in ML_drugs){
  for (j in 1:5){
    ML_run = run_drug_ml(drug = i, p_seed = j)
    
    saveRDS(ML_run, paste("RData/Results/drugGeneralizability/ML_currentUsage_seed", j, "_", i, ".rds", sep = "")) 
  }
}

