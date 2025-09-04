#-------------------------------------------#
# Subgroup analyses
#-------------------------------------------#
# Libraries 
library(haven)
library(foreign)
library(foreach)
library(tidyverse)
library(dplyr)
library(tidyr)
library(openxlsx)
library(ggplot2)
library(survival)
library(DescTools)
library(conflicted)
conflicts_prefer(dplyr::recode)

# I will perform Model (2) 
# Adjusted with variables selected via stepwise regression (BMI, alcohol, grade, ER, PR, HER2, year of diagnosis, country, and stage) with country and stage in strata; 

#--------------------#
# Overall mortality 
#-------------------#
# Creation of subgroups ----

# data:
set <- readRDS("/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Data/processed_data/002_OM_models_script.rds")

# Preparing subgroup data sets to loop through later:
set_pre <- set %>% filter(menop_status_dx == "Premenopausal") # 370
set_post <- set %>% filter(menop_status_dx == "Postmenopausal")
set_pre_bl <- set %>% filter(Menop_Bld == 0)  # 413, to check results among premenopausal status at blood collection 

set$TN <- ifelse(set$Er_Adj_c == "Negative" & set$Pr_Adj_c == "Negative" & set$Her2_Adj_c == "Negative", 
                 "TN", "nonTN") # 153 TN

set_pre$TN <- ifelse(set_pre$Er_Adj_c == "Negative" & set_pre$Pr_Adj_c == "Negative" & set_pre$Her2_Adj_c == "Negative", 
                 "TN", "nonTN") # 41 TN
set_post$TN <- ifelse(set_post$Er_Adj_c == "Negative" & set_post$Pr_Adj_c == "Negative" & set_post$Her2_Adj_c == "Negative", 
                     "TN", "nonTN") # 110 TN
set_pre_TN <- set_pre %>% filter(set_pre$TN == "TN") # 41
set_pre_ERp <- set_pre %>% filter(set_pre$Er_Adj_c == "Positive") # 282
set_pre_ERn <- set_pre %>% filter(set_pre$Er_Adj_c == "Negative") # 88

set_TN <- set %>% filter(set$TN == "TN")
set_nonTN <- set %>% filter(set$TN == "nonTN")

set_ERp <- set %>% filter(Er_Adj_c == "Positive")
set_ERn <- set %>% filter(Er_Adj_c == "Negative")

set_PRp <- set %>% filter(Pr_Adj_c == "Positive")
set_PRn <- set %>% filter(Pr_Adj_c == "Negative")

set_her2p <- set %>% filter(Her2_Adj_c == "Positive")
set_her2n <- set %>% filter(Her2_Adj_c == "Negative")
  
set_fast <- set %>% filter(fasting == "Yes")
set_nonfast <- set %>% filter(fasting == "No" | fasting == "In between")

set_Horm <- set %>% filter(Phrt_Bld_c == "Yes") # Nusers of pill/hrt-ert at blood collection
set_nonHorm <- set %>% filter(Phrt_Bld_c == "No") # nonsers

set$PA_2 <- ifelse(set$Pa_Index_c == "Inactive" | set$Pa_Index_c == "Mod.Inactive", "Inactive",
              ifelse(set$Pa_Index_c == "Mod.Active" | set$Pa_Index_c == "Active", "Active", NA))
set$PA_2 <- as.factor(set$PA_2)

set_active <- set %>% filter(PA_2 == "Active")
set_inactive <- set %>% filter(PA_2 == "Inactive") 

table(set$Smoke_Stat_cc)
set$smoke2 <- ifelse(set$Smoke_Stat_cc == "Current" | set$Smoke_Stat_cc == "Former", "Smoker",
                   ifelse(set$Smoke_Stat_cc == "Never","NonSmoker", NA))

set_smoke <- set %>% filter(set$smoke2 == "Smoker")
set_nonsmok <- set %>% filter(set$smoke2 == "NonSmoker")

table(set$l_school2)
set$l_sch2 <- ifelse(set$l_school2 == "None/Primary", "LowEduc",
                    ifelse(set$l_school2 == "Longer educat", "HighEduc", NA))
set$l_school2 <- as.factor(set$l_school2)

set_highEduc <- set %>% filter(set$l_sch2 == "LowEduc")
set_lowEduc <- set %>% filter(set$l_sch2 == "HighEduc")

set$stage2.2 <- ifelse(set$stage2 == "Localised", "Localised",
                    ifelse(set$stage2 == "Metastatic", "Metastatic", NA))
set$stage2.2 <- as.factor(set$stage2.2)

set_metas <- set %>% filter(stage3 == "Metastatic") 
set_nonmetas <- set %>% filter(stage3 == "Non-Metastatic")

set$stage3.1 <- ifelse(set$stage3 == "Non-Metastatic", "Non-Metastatic",
                       ifelse(set$stage3 == "Metastatic", "Metastatic", NA))
set$stage3.1 <- as.factor(set$stage3.1)

set$bmi2 <- ifelse(set$Bmi_C <25, "<25", 
                   ifelse(set$Bmi_C >=25, ">=25",NA))
set$bmi2 <- as.factor(set$bmi2)

set_above25 <- set %>% filter(bmi2 == "<25") 
set_below25 <- set %>% filter(bmi2 == ">=25")

summary(set$age_dx)
set$age_dx_2 <- ifelse(set$age_dx <55, "<55", 
                   ifelse(set$age_dx >=55, ">=55",NA)) # Based on median age
table(set$age_dx_2)
set$age_dx_2 <- as.factor(set$age_dx_2)
set_age_young55 <- set %>% filter(age_dx_2 == "<55") 
set_age_old55 <- set %>% filter(age_dx_2 == ">=55") 

#-----------------------------------------------------#
# Subgroups: fasting, hormones users, smoking, PA ----
#-----------------------------------------------------#

# Define exposures and subgroups
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     

subgroups <- list(
  set_fast = set_fast, set_nonfast = set_nonfast,
  set_Horm = set_Horm, set_nonHorm = set_nonHorm,
  set_active = set_active, set_inactive = set_inactive,
  set_smoke = set_smoke, set_nonsmok = set_nonsmok
)

# Initialize an empty results list
results_list <- list()

for (subgroup_name in names(subgroups)) {
  subgroup_data <- subgroups[[subgroup_name]]
  results_list[[subgroup_name]] <- list()  # Ensure each subgroup has its own list
  
  for (exposure in exposures) {
    # Create the formula for the Cox regression model
    formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", exposure,
                                "+ bmi_oms + Alc_Re_c + year_diagnosis + l_school2 +
                                 Er_Adj_c + Pr_Adj_c + Her2_Adj_c + gradb + Phrt_Bld_c +
                                 strata(COUNTRY, stage3, menop_status_dx)"))
    
    # Fit the Cox model for the current subgroup
    model <- coxph(formula, data = subgroup_data, id = Idepic)

    # Store model results
    results_list[[subgroup_name]][[exposure]] <- list(
      "model" = model,            
      "conf.int" = summary(model)$conf.int[1, c(1, 3:4)],
      "n" = summary(model)$n,
      "nevent" = summary(model)$nevent
    )
  }
}

# Initialize a list to store rows for the final table
result_tables <- list()

# Loop through subgroups and exposures to create result tables
for (subgroup_name in names(subgroups)) {
  for (exposure in exposures) {
    cat("Results for subgroup:", subgroup_name, "and exposure:", exposure, "\n")
    
    # Extract stored results for the current exposure
    model_results <- results_list[[subgroup_name]][[exposure]]
    
    if (!is.null(model_results)) {  # Ensure model results exist
      # Extract results into a data frame
      table <- data.frame(
        HR = round(model_results$conf.int[1], 2), 
        lowCI = round(model_results$conf.int[2], 2),
        upCI = round(model_results$conf.int[3], 2),
        n = model_results$n,
        nevent = model_results$nevent
      )
      
      # Store results table with subgroup and exposure as row names
      rownames(table) <- paste(subgroup_name, exposure)
      result_tables[[paste(subgroup_name, exposure, sep = "_")]] <- table  # Save in list
    } else {
      cat("Warning: No model results for", subgroup_name, exposure, "\n")
    }
  }
}

# Desired exposure display names & order
exposure_order <- c("TSH", "fT3", "fT4", "fT3fT4 ratio", "Anti TPO")

# Subgroup labels & order
subgroup_labels <- c(
  set_Horm   = "Hormones users",
  set_nonHorm = "Non-hormone users",
  set_active = "Physically active",
  set_inactive = "Physically inactive",
  set_fast   = "Fasting",
  set_nonfast = "Non-fasting",
  set_smoke  = "Smokers",
  set_nonsmok = "Non-smokers"
)

final_results_1 <- bind_rows(result_tables, .id = "Subgroup_Exposure") %>%
  # robustly split "set_Horm_log_TSH_R_sd" -> Subgroup = "set_Horm", Exposure = "log_TSH_R_sd"
  tidyr::extract(Subgroup_Exposure, into = c("Subgroup", "Exposure"),
                 regex = "^(set_[^_]+)_(.+)$", remove = FALSE) %>%
  # map exposure to the simplified name (use pattern matching to avoid exact-name failures)
  mutate(
    Exposure = case_when(
      grepl("TSH", Exposure, ignore.case = TRUE)                          ~ "TSH",
      grepl("fT3fT4", Exposure, ignore.case = TRUE) |
        grepl("fT3fT4r", Exposure, ignore.case = TRUE)                     ~ "fT3fT4 ratio",
      grepl("fT3", Exposure, ignore.case = TRUE)                         ~ "fT3",
      grepl("fT4", Exposure, ignore.case = TRUE)                         ~ "fT4",
      grepl("Anti[_ ]?TPO", Exposure, ignore.case = TRUE)                ~ "Anti TPO",
      TRUE                                                               ~ as.character(Exposure)  # fallback
    ),
    # map subgroup codes to readable labels
    Subgroup = recode(Subgroup, !!!subgroup_labels)
  ) %>%
  # enforce ordering
  mutate(
    Exposure = factor(Exposure, levels = exposure_order),
    Subgroup = factor(Subgroup, levels = unname(subgroup_labels))
  ) %>%
  arrange(Subgroup, Exposure) %>%
  select(Subgroup, Exposure, HR, lowCI, upCI, n, nevent)

final_results_1

# Subgroups: menopause, ER, PR, HER2, TN, stage, BMI, age ----
#-------------------------------------------------------------#

#------------------------------------------------#
# Model for pre- vs postmenopause ----
#------------------------------------------------#
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     
subgroups <- list(set_pre = set_pre,
                  set_post = set_post)

results_list <- list()

# Loop through subgroups and exposures to fit models
for (subgroup_name in names(subgroups)) {
  subgroup_data <- subgroups[[subgroup_name]]
  results_list[[subgroup_name]] <- list()  # Ensure each subgroup has its own list
  
  for (exposure in exposures) {
    # Create the formula for the Cox regression model
    formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", exposure,
                                "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb + l_school2 + Er_Adj_c +
                                 Pr_Adj_c + Her2_Adj_c + strata(COUNTRY, stage3)"))
    
    # Fit the Cox model for the current subgroup
    model <- coxph(formula, data = subgroup_data, id = Idepic)
    
    # Store model results
    results_list[[subgroup_name]][[exposure]] <- list(
      "model" = model,            
      "conf.int" = summary(model)$conf.int[1, c(1, 3:4)],
      "n" = summary(model)$n,
      "nevent" = summary(model)$nevent
    )
  }
}

# Initialize a list to store rows for the final table
result_tables <- list()

# Loop through subgroups and exposures to create result tables
for (subgroup_name in names(subgroups)) {
  for (exposure in exposures) {
    cat("Results for subgroup:", subgroup_name, "and exposure:", exposure, "\n")
    
    # Extract stored results for the current exposure
    model_results <- results_list[[subgroup_name]][[exposure]]
    
    if (!is.null(model_results)) {  # Ensure model results exist
      # Extract results into a data frame
      table <- data.frame(
        HR = round(model_results$conf.int[1], 2), 
        lowCI = round(model_results$conf.int[2], 2),
        upCI = round(model_results$conf.int[3], 2),
        n = model_results$n,
        nevent = model_results$nevent
      )
      
      # Store results table with subgroup and exposure as row names
      rownames(table) <- paste(subgroup_name, exposure)
      result_tables[[paste(subgroup_name, exposure, sep = "_")]] <- table  # Save in list
    } else {
      cat("Warning: No model results for", subgroup_name, exposure, "\n")
    }
  }
}

# View the results
result_tables

# Desired exposure display names & order
exposure_order <- c("TSH", "fT3", "fT4", "fT3fT4 ratio", "Anti TPO")

# Subgroup labels & order
subgroup_labels <- c(
  set_pre = "Premenopausal",
  set_post = "Postmenopausal"
)

final_results_prepost <- bind_rows(result_tables, .id = "Subgroup_Exposure") %>%
  tidyr::extract(Subgroup_Exposure, into = c("Subgroup", "Exposure"),
                 regex = "^(set_[^_]+)_(.+)$", remove = FALSE) %>%
  mutate(
    Exposure = case_when(
      grepl("TSH", Exposure, ignore.case = TRUE)                         ~ "TSH",
      grepl("fT3fT4", Exposure, ignore.case = TRUE) |
        grepl("fT3fT4r", Exposure, ignore.case = TRUE)                   ~ "fT3fT4 ratio",
      grepl("fT3", Exposure, ignore.case = TRUE)                         ~ "fT3",
      grepl("fT4", Exposure, ignore.case = TRUE)                         ~ "fT4",
      grepl("Anti[_ ]?TPO", Exposure, ignore.case = TRUE)                ~ "Anti TPO",
      TRUE                                                               ~ as.character(Exposure)
    ),
    # force dplyr's recode
    Subgroup = dplyr::recode(Subgroup, !!!subgroup_labels)
  ) %>%
  mutate(
    Exposure = factor(Exposure, levels = exposure_order),
    Subgroup = factor(Subgroup, levels = unname(subgroup_labels))
  ) %>%
  arrange(Subgroup, Exposure) %>%
  select(Subgroup, Exposure, HR, lowCI, upCI, n, nevent)

final_results_prepost

#---------------------------------------------#
# Models BMI >= 25 and < 25 ----
#---------------------------------------------#
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     
subgroups <- list(BMI_below25 = set_below25, 
                  BMI_above25 = set_above25)

results_list <- list()

# Loop through subgroups and exposures to fit models
for (subgroup_name in names(subgroups)) {
  subgroup_data <- subgroups[[subgroup_name]]
  results_list[[subgroup_name]] <- list()  # Ensure each subgroup has its own list
  
  for (exposure in exposures) {
    # Create the formula for the Cox regression model
    formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", exposure,
                                "+ Alc_Re_c + year_diagnosis + gradb + l_school2 + 
                                 Er_Adj_c + Pr_Adj_c + Her2_Adj_c + strata(COUNTRY, stage3, menop_status_dx)"))
    
    # Fit the Cox model for the current subgroup
    model <- coxph(formula, data = subgroup_data, id = Idepic)
    
    # Store model results
    results_list[[subgroup_name]][[exposure]] <- list(
      "model" = model,            
      "conf.int" = summary(model)$conf.int[1, c(1, 3:4)],
      "n" = summary(model)$n,
      "nevent" = summary(model)$nevent
    )
  }
}

# Initialize a list to store rows for the final table
result_tables <- list()

# Loop through subgroups and exposures to create result tables
for (subgroup_name in names(subgroups)) {
  for (exposure in exposures) {
    cat("Results for subgroup:", subgroup_name, "and exposure:", exposure, "\n")
    
    # Extract stored results for the current exposure
    model_results <- results_list[[subgroup_name]][[exposure]]
    
    if (!is.null(model_results)) {  # Ensure model results exist
      # Extract results into a data frame
      table <- data.frame(
        HR = round(model_results$conf.int[1], 2), 
        lowCI = round(model_results$conf.int[2], 2),
        upCI = round(model_results$conf.int[3], 2),
        n = model_results$n,
        nevent = model_results$nevent
      )
      
      # Store results table with subgroup and exposure as row names
      rownames(table) <- paste(subgroup_name, exposure)
      result_tables[[paste(subgroup_name, exposure, sep = "_")]] <- table  # Save in list
    } else {
      cat("Warning: No model results for", subgroup_name, exposure, "\n")
    }
  }
}

# Desired exposure display names & order
exposure_order <- c("TSH", "fT3", "fT4", "fT3fT4 ratio", "Anti TPO")

# Subgroup labels & order
subgroup_labels <- c(
  BMI_below25 = "BMI <25",
  BMI_above25 = "BMI +25"
)

final_results_BMI <- bind_rows(result_tables, .id = "Subgroup_Exposure") %>%
  tidyr::extract(Subgroup_Exposure, into = c("Subgroup", "Exposure_raw"),
                 regex = "^([^_]+_[^_]+)_(.+)$") %>%
  mutate(
    Exposure = case_when(
      grepl("TSH", Exposure_raw, ignore.case = TRUE)           ~ "TSH",
      grepl("fT3fT4", Exposure_raw, ignore.case = TRUE)        ~ "fT3fT4 ratio",
      grepl("fT3", Exposure_raw, ignore.case = TRUE)           ~ "fT3",
      grepl("fT4", Exposure_raw, ignore.case = TRUE)           ~ "fT4",
      grepl("Anti[_ ]?TPO", Exposure_raw, ignore.case = TRUE)  ~ "Anti TPO",
      TRUE                                                     ~ Exposure_raw
    ),
    Subgroup = case_when(
      Subgroup == "BMI_below25" ~ "BMI <25",
      Subgroup == "BMI_above25" ~ "BMI ≥25",
      TRUE                       ~ Subgroup
    )
  ) %>%
  mutate(
    Exposure = factor(Exposure, levels = exposure_order),
    Subgroup = factor(Subgroup, levels = c("BMI <25", "BMI ≥25"))
  ) %>%
  arrange(Subgroup, Exposure) %>%
  select(Subgroup, Exposure, HR, lowCI, upCI, n, nevent)

final_results_BMI

#------------------------------------------#
# Model for ER(+) vs ER(-) ----
#------------------------------------------#
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     
subgroups <- list(set_ERp = set_ERp, 
                  set_ERn = set_ERn)

results_list <- list()

# Loop through subgroups and exposures to fit models
for (subgroup_name in names(subgroups)) {
  subgroup_data <- subgroups[[subgroup_name]]
  results_list[[subgroup_name]] <- list()  # Ensure each subgroup has its own list
  
  for (exposure in exposures) {
    # Create the formula for the Cox regression model
    formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", exposure,
                                "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb + l_school2 + 
                                 Pr_Adj_c + Her2_Adj_c + strata(COUNTRY, stage3,menop_status_dx)"))
    
    # Fit the Cox model for the current subgroup
    model <- coxph(formula, data = subgroup_data, id = Idepic)
    
    # Store model results
    results_list[[subgroup_name]][[exposure]] <- list(
      "model" = model,            
      "conf.int" = summary(model)$conf.int[1, c(1, 3:4)],
      "n" = summary(model)$n,
      "nevent" = summary(model)$nevent
    )
  }
}

# Initialize a list to store rows for the final table
result_tables <- list()

# Loop through subgroups and exposures to create result tables
for (subgroup_name in names(subgroups)) {
  for (exposure in exposures) {
    cat("Results for subgroup:", subgroup_name, "and exposure:", exposure, "\n")
    
    # Extract stored results for the current exposure
    model_results <- results_list[[subgroup_name]][[exposure]]
    
    if (!is.null(model_results)) {  # Ensure model results exist
      # Extract results into a data frame
      table <- data.frame(
        HR = round(model_results$conf.int[1], 2), 
        lowCI = round(model_results$conf.int[2], 2),
        upCI = round(model_results$conf.int[3], 2),
        n = model_results$n,
        nevent = model_results$nevent
      )
      
      # Store results table with subgroup and exposure as row names
      rownames(table) <- paste(subgroup_name, exposure)
      result_tables[[paste(subgroup_name, exposure, sep = "_")]] <- table  # Save in list
    } else {
      cat("Warning: No model results for", subgroup_name, exposure, "\n")
    }
  }
}

# Desired exposure display names & order
exposure_order <- c("TSH", "fT3", "fT4", "fT3fT4 ratio", "Anti TPO")

# Subgroup labels & order
subgroup_labels <- c(
  set_ERp = "ER+",
  set_ERn = "ER-"
)

final_results_ER <- bind_rows(result_tables, .id = "Subgroup_Exposure") %>%
  tidyr::extract(Subgroup_Exposure, into = c("Subgroup", "Exposure_raw"),
                 regex = "^([^_]+_[^_]+)_(.+)$") %>%
  mutate(
    Exposure = case_when(
      grepl("TSH", Exposure_raw, ignore.case = TRUE)           ~ "TSH",
      grepl("fT3fT4", Exposure_raw, ignore.case = TRUE)        ~ "fT3fT4 ratio",
      grepl("fT3", Exposure_raw, ignore.case = TRUE)           ~ "fT3",
      grepl("fT4", Exposure_raw, ignore.case = TRUE)           ~ "fT4",
      grepl("Anti[_ ]?TPO", Exposure_raw, ignore.case = TRUE)  ~ "Anti TPO",
      TRUE                                                     ~ Exposure_raw
    ),
    Subgroup = recode(Subgroup,
                      set_ERp = "ER+",
                      set_ERn = "ER-")
  ) %>%
  mutate(
    Exposure = factor(Exposure, levels = exposure_order),
    Subgroup = factor(Subgroup, levels = c("ER+", "ER-"))
  ) %>%
  arrange(Subgroup, Exposure) %>%
  select(Subgroup, Exposure, HR, lowCI, upCI, n, nevent)

final_results_ER

#------------------------------------------#
# Model for Metastatic vs. Localised ----
#------------------------------------------#
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     
subgroups <- list(set_metas = set_metas, 
                  set_nonmetas = set_nonmetas)

results_list <- list()

# Loop through subgroups and exposures to fit models
for (subgroup_name in names(subgroups)) {
  subgroup_data <- subgroups[[subgroup_name]]
  results_list[[subgroup_name]] <- list()  # Ensure each subgroup has its own list
  
  for (exposure in exposures) {
    # Create the formula for the Cox regression model
    formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", exposure,
                                "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb + l_school2 +
                                 Er_Adj_c + Pr_Adj_c + Her2_Adj_c + strata(COUNTRY, menop_status_dx)"))
    
    # Fit the Cox model for the current subgroup
    model <- coxph(formula, data = subgroup_data, id = Idepic)
    
    # Store model results
    results_list[[subgroup_name]][[exposure]] <- list(
      "model" = model,            
      "conf.int" = summary(model)$conf.int[1, c(1, 3:4)],
      "n" = summary(model)$n,
      "nevent" = summary(model)$nevent
    )
  }
}

# Initialize a list to store rows for the final table
result_tables <- list()

# Loop through subgroups and exposures to create result tables
for (subgroup_name in names(subgroups)) {
  for (exposure in exposures) {
    cat("Results for subgroup:", subgroup_name, "and exposure:", exposure, "\n")
    
    # Extract stored results for the current exposure
    model_results <- results_list[[subgroup_name]][[exposure]]
    
    if (!is.null(model_results)) {  # Ensure model results exist
      # Extract results into a data frame
      table <- data.frame(
        HR = round(model_results$conf.int[1], 2), 
        lowCI = round(model_results$conf.int[2], 2),
        upCI = round(model_results$conf.int[3], 2),
        n = model_results$n,
        nevent = model_results$nevent
      )
      
      # Store results table with subgroup and exposure as row names
      rownames(table) <- paste(subgroup_name, exposure)
      result_tables[[paste(subgroup_name, exposure, sep = "_")]] <- table  # Save in list
    } else {
      cat("Warning: No model results for", subgroup_name, exposure, "\n")
    }
  }
}

# Desired exposure display names & order
exposure_order <- c("TSH", "fT3", "fT4", "fT3fT4 ratio", "Anti TPO")

# Subgroup labels & order
subgroup_labels <- c(
  set_metas = "Metastatic",
  set_nonmetas = "Localised"
)

final_results_stage <- bind_rows(result_tables, .id = "Subgroup_Exposure") %>%
  tidyr::extract(Subgroup_Exposure, into = c("Subgroup", "Exposure_raw"),
                 regex = "^([^_]+_[^_]+)_(.+)$") %>%
  mutate(
    Exposure = case_when(
      grepl("TSH", Exposure_raw, ignore.case = TRUE)           ~ "TSH",
      grepl("fT3fT4", Exposure_raw, ignore.case = TRUE)        ~ "fT3fT4 ratio",
      grepl("fT3", Exposure_raw, ignore.case = TRUE)           ~ "fT3",
      grepl("fT4", Exposure_raw, ignore.case = TRUE)           ~ "fT4",
      grepl("Anti[_ ]?TPO", Exposure_raw, ignore.case = TRUE)  ~ "Anti TPO",
      TRUE                                                     ~ Exposure_raw
    ),
    Subgroup = recode(Subgroup,
                      set_metas = "Metastatic",
                      set_nonmetas = "Localised")
  ) %>%
  mutate(
    Exposure = factor(Exposure, levels = exposure_order),
    Subgroup = factor(Subgroup, levels = c("Metastatic", "Localised"))
  ) %>%
  arrange(Subgroup, Exposure) %>%
  select(Subgroup, Exposure, HR, lowCI, upCI, n, nevent)

final_results_stage

# All tables for paper ----
final_table <- rbind(final_results_1,
                     final_results_prepost,
                     final_results_BMI,
                     final_results_ER,
                     final_results_stage)

write.table(final_table, file = paste0(output.dir,"004_subgroup_all_cause.txt"), 
            col.names = T, row.names = T, sep = "\t", quote = F)

#----------------------------------------------------------#
# Other subgroups: age, education, TN, PR, HER2, pre_ER ----
#----------------------------------------------------------#
# Model for age by median 2 groups ----
#------------------------------------------#
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     
subgroups <- list()

subgroups <- list(set_age_old55 = set_age_old55,
                  set_age_young55 = set_age_young55)

results_list <- list()

# Loop through subgroups and exposures to fit models
for (subgroup_name in names(subgroups)) {
  subgroup_data <- subgroups[[subgroup_name]]
  results_list[[subgroup_name]] <- list()  # Ensure each subgroup has its own list
  
  for (exposure in exposures) {
    # Create the formula for the Cox regression model
    formula <- as.formula(paste("Surv(time_surv, Death_Status =='Dead') ~ ", exposure,
                                "+ Alc_Re_c + year_diagnosis + gradb + l_school2 +
                           Er_Adj_c + Pr_Adj_c + Her2_Adj_c + 
                          strata(COUNTRY, stage3, menop_status_dx)"))
    
    # Fit the Cox model for the current subgroup
    model <- coxph(formula, data = subgroup_data, id = Idepic)
    
    # Store model results
    results_list[[subgroup_name]][[exposure]] <- list(
      "model" = model,            
      "conf.int" = summary(model)$conf.int[1, c(1, 3:4)],
      "n" = summary(model)$n,
      "nevent" = summary(model)$nevent
    )
  }
}

# Initialize a list to store rows for the final table
result_tables <- list()

# Loop through subgroups and exposures to create result tables
for (subgroup_name in names(subgroups)) {
  for (exposure in exposures) {
    cat("Results for subgroup:", subgroup_name, "and exposure:", exposure, "\n")
    
    # Extract stored results for the current exposure
    model_results <- results_list[[subgroup_name]][[exposure]]
    
    if (!is.null(model_results)) {  # Ensure model results exist
      # Extract results into a data frame
      table <- data.frame(
        HR = round(model_results$conf.int[1], 2), 
        lowCI = round(model_results$conf.int[2], 2),
        upCI = round(model_results$conf.int[3], 2),
        n = model_results$n,
        nevent = model_results$nevent
      )
      
      # Store results table with subgroup and exposure as row names
      rownames(table) <- paste(subgroup_name, exposure)
      result_tables[[paste(subgroup_name, exposure, sep = "_")]] <- table  # Save in list
    } else {
      cat("Warning: No model results for", subgroup_name, exposure, "\n")
    }
  }
}

# View the results
result_tables

final_results_age <- bind_rows(result_tables, .id = "Subgroup_Exposure") %>%
  separate(Subgroup_Exposure, into = c("Subgroup", "Exposure"), sep = "_", extra = "merge") %>%
  arrange(Subgroup, Exposure)  # Sort the table

final_results_age

#------------------------------------------------#
# Model for educational level ----
#------------------------------------------------#
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     
subgroups <- list(set_lowEduc = set_lowEduc,
                  set_highEduc = set_highEduc)

results_list <- list()

# Loop through subgroups and exposures to fit models
for (subgroup_name in names(subgroups)) {
  subgroup_data <- subgroups[[subgroup_name]]
  results_list[[subgroup_name]] <- list()  # Ensure each subgroup has its own list
  
  for (exposure in exposures) {
    # Create the formula for the Cox regression model
    formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", exposure,
                                "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb + 
                                 Pr_Adj_c + Her2_Adj_c + strata(COUNTRY, stage3,menop_status_dx)"))
    
    # Fit the Cox model for the current subgroup
    model <- coxph(formula, data = subgroup_data, id = Idepic)
    
    # Store model results
    results_list[[subgroup_name]][[exposure]] <- list(
      "model" = model,            
      "conf.int" = summary(model)$conf.int[1, c(1, 3:4)],
      "n" = summary(model)$n,
      "nevent" = summary(model)$nevent
    )
  }
}

# Initialize a list to store rows for the final table
result_tables <- list()

# Loop through subgroups and exposures to create result tables
for (subgroup_name in names(subgroups)) {
  for (exposure in exposures) {
    cat("Results for subgroup:", subgroup_name, "and exposure:", exposure, "\n")
    
    # Extract stored results for the current exposure
    model_results <- results_list[[subgroup_name]][[exposure]]
    
    if (!is.null(model_results)) {  # Ensure model results exist
      # Extract results into a data frame
      table <- data.frame(
        HR = round(model_results$conf.int[1], 2), 
        lowCI = round(model_results$conf.int[2], 2),
        upCI = round(model_results$conf.int[3], 2),
        n = model_results$n,
        nevent = model_results$nevent
      )
      
      # Store results table with subgroup and exposure as row names
      rownames(table) <- paste(subgroup_name, exposure)
      result_tables[[paste(subgroup_name, exposure, sep = "_")]] <- table  # Save in list
    } else {
      cat("Warning: No model results for", subgroup_name, exposure, "\n")
    }
  }
}

# View the results
result_tables

final_results_educ <- bind_rows(result_tables, .id = "Subgroup_Exposure") %>%
  separate(Subgroup_Exposure, into = c("Subgroup", "Exposure"), sep = "_", extra = "merge") %>%
  arrange(Subgroup, Exposure)  # Sort the table

final_results_educ

#------------------------------------------#
# Model for Premenopausal, ER+ and ER- ----
#------------------------------------------#

exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     
subgroups <- list(set_pre_ERp = set_pre_ERp,
                  set_pre_ERn = set_pre_ERn)

results_list <- list()

# Loop through subgroups and exposures to fit models
for (subgroup_name in names(subgroups)) {
  subgroup_data <- subgroups[[subgroup_name]]
  results_list[[subgroup_name]] <- list()  # Ensure each subgroup has its own list
  
  for (exposure in exposures) {
    # Create the formula for the Cox regression model
    formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", exposure,
                                "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb+ l_school2 +  
                                 Pr_Adj_c + Her2_Adj_c + strata(COUNTRY, stage3)"))
    
    # Fit the Cox model for the current subgroup
    model <- coxph(formula, data = subgroup_data, id = Idepic)
    
    # Store model results
    results_list[[subgroup_name]][[exposure]] <- list(
      "model" = model,            
      "conf.int" = summary(model)$conf.int[1, c(1, 3:4)],
      "n" = summary(model)$n,
      "nevent" = summary(model)$nevent
    )
  }
}

# Initialize a list to store rows for the final table
result_tables <- list()

# Loop through subgroups and exposures to create result tables
for (subgroup_name in names(subgroups)) {
  for (exposure in exposures) {
    cat("Results for subgroup:", subgroup_name, "and exposure:", exposure, "\n")
    
    # Extract stored results for the current exposure
    model_results <- results_list[[subgroup_name]][[exposure]]
    
    if (!is.null(model_results)) {  # Ensure model results exist
      # Extract results into a data frame
      table <- data.frame(
        HR = round(model_results$conf.int[1], 2), 
        lowCI = round(model_results$conf.int[2], 2),
        upCI = round(model_results$conf.int[3], 2),
        n = model_results$n,
        nevent = model_results$nevent
      )
      
      # Store results table with subgroup and exposure as row names
      rownames(table) <- paste(subgroup_name, exposure)
      result_tables[[paste(subgroup_name, exposure, sep = "_")]] <- table  # Save in list
    } else {
      cat("Warning: No model results for", subgroup_name, exposure, "\n")
    }
  }
}

# View the results
result_tables

final_results_pre_ER <- bind_rows(result_tables, .id = "Subgroup_Exposure") %>%
  separate(Subgroup_Exposure, into = c("Subgroup", "Exposure"), sep = "_", extra = "merge") %>%
  arrange(Subgroup, Exposure)  # Sort the table

final_results_pre_ER
#------------------------------------------#
# Model for TNBC vs Non-TNBC ----
#------------------------------------------#
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     
subgroups <- list(set_TN = set_TN, 
                  set_nonTN = set_nonTN)

results_list <- list()

# Loop through subgroups and exposures to fit models
for (subgroup_name in names(subgroups)) {
  subgroup_data <- subgroups[[subgroup_name]]
  results_list[[subgroup_name]] <- list()  # Ensure each subgroup has its own list
  
  for (exposure in exposures) {
    # Create the formula for the Cox regression model
    formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", exposure,
                                "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb + l_school2 + 
                                 + strata(COUNTRY, stage3, menop_status_dx)"))
    
    # Fit the Cox model for the current subgroup
    model <- coxph(formula, data = subgroup_data, id = Idepic)
    
    # Store model results
    results_list[[subgroup_name]][[exposure]] <- list(
      "model" = model,            
      "conf.int" = summary(model)$conf.int[1, c(1, 3:4)],
      "n" = summary(model)$n,
      "nevent" = summary(model)$nevent
    )
  }
}

# Initialize a list to store rows for the final table
result_tables <- list()

# Loop through subgroups and exposures to create result tables
for (subgroup_name in names(subgroups)) {
  for (exposure in exposures) {
    cat("Results for subgroup:", subgroup_name, "and exposure:", exposure, "\n")
    
    # Extract stored results for the current exposure
    model_results <- results_list[[subgroup_name]][[exposure]]
    
    if (!is.null(model_results)) {  # Ensure model results exist
      # Extract results into a data frame
      table <- data.frame(
        HR = round(model_results$conf.int[1], 2), 
        lowCI = round(model_results$conf.int[2], 2),
        upCI = round(model_results$conf.int[3], 2),
        n = model_results$n,
        nevent = model_results$nevent
      )
      
      # Store results table with subgroup and exposure as row names
      rownames(table) <- paste(subgroup_name, exposure)
      result_tables[[paste(subgroup_name, exposure, sep = "_")]] <- table  # Save in list
    } else {
      cat("Warning: No model results for", subgroup_name, exposure, "\n")
    }
  }
}

# View the results
result_tables

final_results_TN <- bind_rows(result_tables, .id = "Subgroup_Exposure") %>%
  separate(Subgroup_Exposure, into = c("Subgroup", "Exposure"), sep = "_", extra = "merge") %>%
  arrange(Subgroup, Exposure)  # Sort the table

final_results_TN

#------------------------------------------#
# Model for PR(+) vs PR(-) ----
#------------------------------------------#
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     
subgroups <- list(set_PRp = set_PRp, 
                  set_PRn = set_PRn)

results_list <- list()

# Loop through subgroups and exposures to fit models
for (subgroup_name in names(subgroups)) {
  subgroup_data <- subgroups[[subgroup_name]]
  results_list[[subgroup_name]] <- list()  # Ensure each subgroup has its own list
  
  for (exposure in exposures) {
    # Create the formula for the Cox regression model
    formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", exposure,
                                "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb + l_school2 +
                                 Er_Adj_c + Her2_Adj_c + strata(COUNTRY, stage3, menop_status_dx)"))
    
    # Fit the Cox model for the current subgroup
    model <- coxph(formula, data = subgroup_data, id = Idepic)
    
    # Store model results
    results_list[[subgroup_name]][[exposure]] <- list(
      "model" = model,            
      "conf.int" = summary(model)$conf.int[1, c(1, 3:4)],
      "n" = summary(model)$n,
      "nevent" = summary(model)$nevent
    )
  }
}

# Initialize a list to store rows for the final table
result_tables <- list()

# Loop through subgroups and exposures to create result tables
for (subgroup_name in names(subgroups)) {
  for (exposure in exposures) {
    cat("Results for subgroup:", subgroup_name, "and exposure:", exposure, "\n")
    
    # Extract stored results for the current exposure
    model_results <- results_list[[subgroup_name]][[exposure]]
    
    if (!is.null(model_results)) {  # Ensure model results exist
      # Extract results into a data frame
      table <- data.frame(
        HR = round(model_results$conf.int[1], 2), 
        lowCI = round(model_results$conf.int[2], 2),
        upCI = round(model_results$conf.int[3], 2),
        n = model_results$n,
        nevent = model_results$nevent
      )
      
      # Store results table with subgroup and exposure as row names
      rownames(table) <- paste(subgroup_name, exposure)
      result_tables[[paste(subgroup_name, exposure, sep = "_")]] <- table  # Save in list
    } else {
      cat("Warning: No model results for", subgroup_name, exposure, "\n")
    }
  }
}

# View the results
result_tables

final_results_PR <- bind_rows(result_tables, .id = "Subgroup_Exposure") %>%
  separate(Subgroup_Exposure, into = c("Subgroup", "Exposure"), sep = "_", extra = "merge") %>%
  arrange(Subgroup, Exposure)  # Sort the table

final_results_PR

#------------------------------------------#
# Model for HER2(+) vs HER2(-) ----
#------------------------------------------#
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     
subgroups <- list(set_her2p = set_her2p, 
                  set_her2n = set_her2n)
results_list <- list()

# Loop through subgroups and exposures to fit models
for (subgroup_name in names(subgroups)) {
  subgroup_data <- subgroups[[subgroup_name]]
  results_list[[subgroup_name]] <- list()  # Ensure each subgroup has its own list
  
  for (exposure in exposures) {
    # Create the formula for the Cox regression model
    formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", exposure,
                                "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb + l_school2 +
                                 Er_Adj_c + Pr_Adj_c + strata(COUNTRY, stage3,menop_status_dx)"))
    
    # Fit the Cox model for the current subgroup
    model <- coxph(formula, data = subgroup_data, id = Idepic)
    
    # Store model results
    results_list[[subgroup_name]][[exposure]] <- list(
      "model" = model,            
      "conf.int" = summary(model)$conf.int[1, c(1, 3:4)],
      "n" = summary(model)$n,
      "nevent" = summary(model)$nevent
    )
  }
}

# Initialize a list to store rows for the final table
result_tables <- list()

# Loop through subgroups and exposures to create result tables
for (subgroup_name in names(subgroups)) {
  for (exposure in exposures) {
    cat("Results for subgroup:", subgroup_name, "and exposure:", exposure, "\n")
    
    # Extract stored results for the current exposure
    model_results <- results_list[[subgroup_name]][[exposure]]
    
    if (!is.null(model_results)) {  # Ensure model results exist
      # Extract results into a data frame
      table <- data.frame(
        HR = round(model_results$conf.int[1], 2), 
        lowCI = round(model_results$conf.int[2], 2),
        upCI = round(model_results$conf.int[3], 2),
        n = model_results$n,
        nevent = model_results$nevent
      )
      
      # Store results table with subgroup and exposure as row names
      rownames(table) <- paste(subgroup_name, exposure)
      result_tables[[paste(subgroup_name, exposure, sep = "_")]] <- table  # Save in list
    } else {
      cat("Warning: No model results for", subgroup_name, exposure, "\n")
    }
  }
}

# View the results
result_tables

final_results_her2 <- bind_rows(result_tables, .id = "Subgroup_Exposure") %>%
  separate(Subgroup_Exposure, into = c("Subgroup", "Exposure"), sep = "_", extra = "merge") %>%
  arrange(Subgroup, Exposure)  # Sort the table

final_results_her2

# Join all tables 
final_table <- rbind(final_results_1,
                     final_results_prepost,
                     final_results_BMI,
                     final_results_ER,
                     final_results_stage,
                     final_results_pre_ER,
                     final_results_PR,
                     final_results_TN,
                     final_results_age)
head(final_table)

write.table(final_table, file = paste0(output.dir,"004_Subgroup.txt"), 
            col.names = T, row.names = T, sep = "\t", quote = F)

# Heterogeneity ----
# for models on: BMI, ER, PA, stage, menopausal status at dx, users hormones, education

exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     

# Het Horm Users ----
set_Horm <- set %>% filter(Phrt_Bld_c == "Yes") # Users of pill/hrt-ert at blood collection
set_nonHorm <- set %>% filter(Phrt_Bld_c == "No") # Non-users

het_list <- list()

for (i in exposures) {
  
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", i,
                              "+ Phrt_Bld_c + bmi_oms + Alc_Re_c + year_diagnosis + gradb + 
                               Er_Adj_c + Pr_Adj_c + Her2_Adj_c + l_school2 +
                               strata(COUNTRY,stage3, menop_status_dx)"))
  
  formula_2 <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", 
                                i,"+ Phrt_Bld_c +",i," *Phrt_Bld_c +
                                bmi_oms + Alc_Re_c + year_diagnosis + gradb + l_school2 +
                                Er_Adj_c + Pr_Adj_c + Her2_Adj_c +
                                strata(COUNTRY,stage3, menop_status_dx)"))
  
  model <- coxph(formula, data = set)    
  model2 <- coxph(formula_2, data = set)
  anov <- anova(model, model2, test = "Chisq")
  
  results_df <- data.frame(Group = "Horm_use", 
                           Biomarker =i,
                           het_list <- list("pvalue"= round(anov$`Pr(>|Chi|)`[2],3))) 
  
  if (i == exposures[1]) {
    final <- results_df 
  } else {
    final <- rbind(final, results_df)
  }
}
final_horm_use <- final
final_horm_use

# Het Menopause ----
het_list <- list()

for (i in exposures) {
  
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", i,
                              "+ Er_Adj_c + Alc_Re_c + year_diagnosis + gradb + l_school2 +
                               bmi_oms + menop_status_dx +
                               Er_Adj_c + Pr_Adj_c + Her2_Adj_c +
                               strata(COUNTRY,stage3)"))
  
  formula_2 <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", 
                                i,"+ menop_status_dx +",i,"* menop_status_dx + 
                                bmi_oms + Alc_Re_c + year_diagnosis + gradb + 
                                Er_Adj_c + Pr_Adj_c + Her2_Adj_c + l_school2 +
                                strata(COUNTRY,stage3)"))
  
  model <- coxph(formula, data = set)    
  model2 <- coxph(formula_2, data = set)
  anov <- anova(model, model2, test = "Chisq")
  
  results_df <- data.frame(Group = "Menop", 
                           Biomarker = i,
                           het_list <- list("pvalue" = round(anov$`Pr(>|Chi|)`[2],3))) 
  
  if (i == exposures[1]) {
    final <- results_df 
  } else {
    final <- rbind(final, results_df)
  }
}
final_Menop <- final
final_Menop

# Het BMI ----
het_list <- list()

for (exposure in exposures) {

  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", exposure,
                               "+ Alc_Re_c + year_diagnosis + gradb + l_school2 +
                           Er_Adj_c + Pr_Adj_c + Her2_Adj_c + bmi2 + 
                          strata(COUNTRY, stage3, menop_status_dx)"))
  
  formula_2 <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", 
                                exposure,"+ bmi2 +", exposure,"* bmi2 + 
                                Alc_Re_c + year_diagnosis + gradb + l_school2 +
                           Er_Adj_c + Pr_Adj_c + Her2_Adj_c + 
                          strata(COUNTRY, stage3, menop_status_dx)"))
  
  model <- coxph(formula, data = set)    
  model2 <- coxph(formula_2, data = set)
  anov <- anova(model, model2, test = "Chisq")
  
  results_df <- data.frame(Group = "BMI", 
                           Biomarker = exposure,
                           het_list <- list("pvalue"= round(anov$`Pr(>|Chi|)`[2],3))) 
  
  if (exposure == exposures[1]) {
    final <- results_df 
  } else {
    final <- rbind(final, results_df)
  }
}
final_BMI <- final
final_BMI

# Het ER ----
het_list <- list()

for (i in exposures) {
  
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", i,
                              "+ Er_Adj_c + bmi_oms + Alc_Re_c + year_diagnosis + gradb + 
                               Pr_Adj_c + Her2_Adj_c + l_school2 +
                               strata(COUNTRY,stage3, menop_status_dx)"))
  
  formula_2 <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", 
                                i,"+ Er_Adj_c +",i,"* Er_Adj_c + 
                                bmi_oms + Alc_Re_c + year_diagnosis + gradb + l_school2 +
                                Pr_Adj_c + Her2_Adj_c +
                                strata(COUNTRY,stage3, menop_status_dx)"))
  
  model <- coxph(formula, data = set)    
  model2 <- coxph(formula_2, data = set)
  anov <- anova(model, model2, test = "Chisq")
  
  results_df <- data.frame(Group = "ER", 
                           Biomarker =i,
                           het_list <- list("pvalue"= round(anov$`Pr(>|Chi|)`[2],3))) 
  
  if (i == exposures[1]) {
    final <- results_df 
  } else {
    final <- rbind(final, results_df)
  }
}
final_ER <- final
final_ER

# Het Stage ----
set_stage <- set %>% filter(stage3 %in% c("Metastatic", "Non-Metastatic"))
dim(set_stage)

het_list <- list()

for (i in exposures) {
  
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", i,
                              "+ stage3.1 + Alc_Re_c + year_diagnosis + gradb + 
                               bmi_oms + Er_Adj_c + Pr_Adj_c + Her2_Adj_c + l_school2 +
                               strata(COUNTRY,menop_status_dx)"))
  
  formula_2 <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", 
                                i,"+ stage3.1 +",i," *stage3.1 + 
                                bmi_oms + Alc_Re_c + year_diagnosis + gradb + l_school2 +
                                Er_Adj_c + Pr_Adj_c + Her2_Adj_c +
                                strata(COUNTRY,menop_status_dx)"))
  
  model <- coxph(formula, data = set_stage)    
  model2 <- coxph(formula_2, data = set_stage)
  anov <- anova(model, model2, test = "Chisq")
  
  results_df <- data.frame(Group = "stage", 
                           Biomarker =i,
                           het_list <- list("pvalue"= round(anov$`Pr(>|Chi|)`[2],3))) 
  
  if (i == exposures[1]) {
    final <- results_df 
  } else {
    final <- rbind(final, results_df)
  }
}
final_stage <- final
final_stage

# Final table Heterogeneity - paper ----
final_table_het <- rbind(final_horm_use,
                         final_Menop,
                         final_BMI,
                         final_ER,
                         final_stage
                         )
final_table_het

write.table(final_table_het, file = paste0(output.dir,"004_Subgroup_het.txt"), 
            col.names = T, row.names = T, sep = "\t", quote = F)

# Other heterogeneity tests (not included in paper):
# Het Premenopausal, ER ----
het_list <- list()

for (i in exposures) {
  
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", i,
                              "+ Er_Adj_c + bmi_oms + Alc_Re_c + year_diagnosis + gradb + 
                               Pr_Adj_c + Her2_Adj_c + strata(COUNTRY,stage3)"))
  
  formula_2 <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", 
                                i,"+ Er_Adj_c +",i,"* Er_Adj_c + 
                                bmi_oms + Alc_Re_c + year_diagnosis + gradb + 
                               Pr_Adj_c + Her2_Adj_c + strata(COUNTRY,stage3)"))
  
  model <- coxph(formula, data = set_pre)    
  model2 <- coxph(formula_2, data = set_pre)
  anov <- anova(model, model2, test = "Chisq")
  
  results_df <- data.frame(Group = "pre_ER", 
                           Biomarker =i,
                           het_list <- list("pvalue"= round(anov$`Pr(>|Chi|)`[2],3))) 
  
  if (i == exposures[1]) {
    final <- results_df 
  } else {
    final <- rbind(final, results_df)
  }
}
final_pre_ER <- final
final_pre_ER

# Het PA ----

set_pa <- set %>% filter(PA_2 %in% c("Active", "Inactive"))
dim(set_pa)

het_list <- list()

for (i in exposures) {
  
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", i,
                              "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb + 
                               PA_2 +
                               Er_Adj_c + Pr_Adj_c + Her2_Adj_c +
                               strata(COUNTRY,stage3,menop_status_dx)"))
  
  formula_2 <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", 
                                i,"+ PA_2 +",i,"* PA_2 + 
                                bmi_oms + Alc_Re_c + year_diagnosis + gradb + 
                                Er_Adj_c + Pr_Adj_c + Her2_Adj_c +
                                strata(COUNTRY,stage3,menop_status_dx)"))
  
  model <- coxph(formula, data = set_pa)    
  model2 <- coxph(formula_2, data = set_pa)
  anov <- anova(model, model2, test = "Chisq")
  
  results_df <- data.frame(Group = "PA", 
                           Biomarker =i,
                           het_list <- list("pvalue"= round(anov$`Pr(>|Chi|)`[2],3))) 
  
  if (i == exposures[1]) {
    final <- results_df 
  } else {
    final <- rbind(final, results_df)
  }
}
final_PA <- final
final_PA

# Het Education level----

set_highEduc <- set %>% filter(l_sch2== "LowEduc")
set_lowEduc <- set %>% filter(l_sch2 == "HighEduc")

set_school <- set %>% filter(l_sch2 %in% c("LowEduc", "HighEduc"))
dim(set_school)

het_list <- list()

for (i in exposures) {
  
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", i,
                              "+ l_sch2 + bmi_oms + Alc_Re_c + year_diagnosis + gradb + 
                               Er_Adj_c + Pr_Adj_c + Her2_Adj_c +
                               strata(COUNTRY,stage3,menop_status_dx)"))
  
  formula_2 <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_Status =='Dead') ~ ", 
                                i,"+ l_sch2 +",i,"* l_sch2 + 
                                bmi_oms + Alc_Re_c + year_diagnosis + gradb + 
                                Er_Adj_c + Pr_Adj_c + Her2_Adj_c +
                                strata(COUNTRY,stage3,menop_status_dx)"))
  
  model <- coxph(formula, data = set_school)    
  model2 <- coxph(formula_2, data = set_school)
  anov <- anova(model, model2, test = "Chisq")
  
  results_df <- data.frame(Group = "Lschool", 
                           Biomarker =i,
                           het_list <- list("pvalue"= round(anov$`Pr(>|Chi|)`[2],3))) 
  
  if (i == exposures[1]) {
    final <- results_df 
  } else {
    final <- rbind(final, results_df)
  }
}
final_sch <- final
final_sch

# Final table Heterogeneity ----
final_table_het <- rbind(final_Menop,
                     final_ER,
                     final_pre_ER,
                     final_horm_use,
                     final_BMI,
                     final_PA,
                     final_sch,
                     final_stage)
final_table_het

write.table(final_table, file = paste0(output.dir,"004_Subgroup_het.txt"), 
            col.names = T, row.names = T, sep = "\t", quote = F)

saveRDS(set, "/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Data/processed_data/004_subgroup.rds")
