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

set <- readRDS("/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Data/processed_data/002_OM_models_script.rds")

# P10 and P90 percentiles groups of hormones ----

# 10th<TSH<90th & 10th<fT4<90th (Ref)

# TSH ≤ 10th & fT4 ≥ 90th (Proxy of hyperthyroidism)
# TSH ≥ 90th & fT4 ≤ 10th 
# 10th < TSH < 90th & fT4 ≤ 10th 
# TSH ≤ 10th & 10th < fT4 < 90th 
# 10th < TSH < 90th & fT4 ≥ 90th 
# TSH ≥ 90th & 10th < fT4 < 90th 
# TSH ≤ 10th and/or fT4 ≥ 90th (Hyperthyroid like)
# TSH ≥ 90th and/or fT4 ≤ 10th (Hypothyroid like

set <- set %>%
  mutate(
    TSH_percentile = ntile(log_TSH_R_sd, 10),  # Splits into 10 quantile groups
    fT4_percentile = ntile(log_fT4_pmol_L_R_sd, 10)
  ) %>%
  mutate(
    TSH_fT4_Group = case_when(
      TSH_percentile <= 1 & fT4_percentile >= 9  ~ "TSH≤10th & fT4≥90th (Proxy of hypert.)",
      TSH_percentile >= 9 & fT4_percentile <= 1  ~ "TSH≥90th & fT4≤10th (Proxy of hypo.)",
      TSH_percentile > 1 & TSH_percentile < 9 & fT4_percentile <= 1  ~ "10th<TSH<90th & fT4≤10th",
      TSH_percentile <= 1 & fT4_percentile > 1 & fT4_percentile < 9  ~ "TSH≤ 10th & 10th<fT4<90th",
      TSH_percentile > 1 & TSH_percentile < 9 & fT4_percentile >= 9  ~ "10th<TSH<90th & fT4≥90th",
      TSH_percentile >= 9 & fT4_percentile > 1 & fT4_percentile < 9  ~ "TSH≥90th & 10th<fT4<90th",
      TSH_percentile > 1 & TSH_percentile < 9 & fT4_percentile > 1 & fT4_percentile < 9  ~ "10th<TSH<90th & 10th<fT4<90th (Ref.)",  # Reference group
      TRUE ~ NA_character_  # Assigns NA if no condition is met
    ),
    TSH_fT4_Group = factor(TSH_fT4_Group)  # Converts to categorical (factor) variable
  )

table(set$TSH_fT4_Group, useNA="always") # 364 NAs
levels(set$TSH_fT4_Group)

# Model for TSH_fT4 6 conditions----

# Define exposures and subgroups
exposures <- c("TSH_fT4_Group")     

results_list <- list()

for (exposure in exposures) {
  
  # Count n and events per category
  counts <- set %>%
    group_by(!!sym(exposure)) %>%
    summarise(
      n = n(),
      nevent = sum(event == levels(Death_Status)[2], na.rm = TRUE)
    ) %>%
    rename(Category = !!sym(exposure)) %>%
    as.data.frame()
  
  # Fit Cox model
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ", exposure,
                              "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb + l_school2 +
                                 Pr_Adj_c + Her2_Adj_c + strata(COUNTRY, stage3, menop_status_dx)"))
  
  model <- coxph(formula, data = set, id = Idepic)
  aic_value <- AIC(model)
  
  # Store model results
  conf <- summary(model)$conf.int[1:nrow(counts), c(1, 3, 4)]
  rownames(conf) <- levels(factor(set[[exposure]]))[2:(nrow(counts)+1)]  # Skip reference
  
  results_list[[exposure]] <- list(
    "model" = model,
    "conf.int" = conf,
    "aic" = aic_value,
    "counts" = counts
  )
}

# Crear tabla final por exposure
result_tables <- list()

for (exposure in exposures) {
  conf <- results_list[[exposure]]$conf.int
  counts <- results_list[[exposure]]$counts
  aic <- results_list[[exposure]]$aic
  
  table <- data.frame(
    Category = counts$Category,
    HR = round(conf[, 1], 2),
    lowCI = round(conf[, 2], 2),
    upCI = round(conf[, 3], 2),
    n = counts$n,
    nevent = counts$nevent,
    AIC = round(aic, 2)
  )
  
  rownames(table) <- NULL
  result_tables[[exposure]] <- table
}


# Combine all results into a final table
fi_table <- do.call(rbind, result_tables)
fi_table <- round(fi_table, 2)  # Round the entire table to 2 decimal places
print(fi_table)

# Now with another classification:
# 10th < TSH< 90th & 10th < fT4 < 90th  (Ref)
# TSH ≤ 10th and/or fT4 ≥ 90th (Hypothyroid like)
# TSH ≥ 90th and/or fT4 ≤ 10th (Hyperthyroid like)

set <- set %>%
  mutate(
    TSH_percentile = ntile(log_TSH_R_sd, 10),  # Splits into 10 quantile groups
    fT4_percentile = ntile(log_fT4_pmol_L_R_sd, 10)
  ) %>%
  mutate(
    TSH_fT4_Group_3 = case_when(
      TSH_percentile <= 1 | fT4_percentile >= 9  ~ "TSH<=10th &/or fT4>=90th hyperthyroid like",
      TSH_percentile >= 9 | fT4_percentile <= 1  ~ "10th>=TSH<90th &/or fT4<=10th hypothyroid like",
      TSH_percentile > 1 & TSH_percentile < 9 & fT4_percentile > 1 & fT4_percentile < 9  ~ "10th<TSH<90th & 10th<fT4<90th (Ref.)",  # Reference group
      TRUE ~ NA_character_  # Assigns NA if no condition is met
    ),
    TSH_fT4_Group_3 = factor(TSH_fT4_Group_3)  # Converts to categorical (factor) variable
  )

table(set$TSH_fT4_Group_3, useNA="always") # 240 NAs
levels(set$TSH_fT4_Group)

# Model for TSH_fT4 3 conditions----

# Define exposures and subgroups
# Define exposures and subgroups
exposures <- c("TSH_fT4_Group_3")     

results_list <- list()

for (exposure in exposures) {
  
  # Count n and events per category
  counts <- set %>%
    group_by(!!sym(exposure)) %>%
    summarise(
      n = n(),
      nevent = sum(event == levels(Death_Status)[2], na.rm = TRUE)
    ) %>%
    rename(Category = !!sym(exposure)) %>%
    as.data.frame()
  
  # Fit Cox model
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ", exposure,
                              "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb + l_school2 +
                                 Pr_Adj_c + Her2_Adj_c + strata(COUNTRY, stage3, menop_status_dx)"))
  
  model <- coxph(formula, data = set, id = Idepic)
  aic_value <- AIC(model)
  
  # Store model results
  conf <- summary(model)$conf.int[1:nrow(counts), c(1, 3, 4)]
  rownames(conf) <- levels(factor(set[[exposure]]))[2:(nrow(counts)+1)]  # Skip reference
  
  results_list[[exposure]] <- list(
    "model" = model,
    "conf.int" = conf,
    "aic" = aic_value,
    "counts" = counts
  )
}

# Crear tabla final por exposure
result_tables <- list()

for (exposure in exposures) {
  conf <- results_list[[exposure]]$conf.int
  counts <- results_list[[exposure]]$counts
  aic <- results_list[[exposure]]$aic
  
  table <- data.frame(
    Category = counts$Category,
    HR = round(conf[, 1], 2),
    lowCI = round(conf[, 2], 2),
    upCI = round(conf[, 3], 2),
    n = counts$n,
    nevent = counts$nevent,
    AIC = round(aic, 2)
  )
  
  rownames(table) <- NULL
  result_tables[[exposure]] <- table
}


# Combine all results into a final table
fi_table_fT4_bc <- do.call(rbind, result_tables)
fi_table_fT4_bc <- round(fi_table_fT4_bc,2)
fi_table_fT4_bc

fitable <- rbind(fi_table_TSH_bc, fi_table_fT4_bc)
fitable
