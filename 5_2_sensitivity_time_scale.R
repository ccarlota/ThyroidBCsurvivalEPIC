# Analyses using time since diagnosis to end instead of age - Sensitivity
#------------------------------------------------------------------------#
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
library(survival)

# Directories
data.dir <- "/data/Epic/subprojects/Breast_Cancer/files/Maestra/"
output.dir <- "/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Analysis/Tables/"

# Charge processed data:
set <- readRDS("/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Data/processed_data/001_Setup_script.rds")
colnames(set)
set$COUNTRY <- as.factor(set$COUNTRY)

# Hormones per 1-SD increase 
horm <- colnames(set)[grep("_R$", colnames(set))]
horm # "log_TSH_R" "log_fT3_pmol_L_R" "log_fT4_pmol_L_R" "log_fT3fT4r_R" (R stands for Residuals from Batch effect!)

for (i in horm) {
  horm_sd <- paste0(i, "_sd")
  x <- data.frame(sd = set[, i]/sd(na.omit(set[, i])))
  colnames(x) <- horm_sd
  set <- cbind(set, x)
}
colnames(set)

# Model OM ----
# age at diagnostic, BMI, alcohol, year of diagnosis, ER, PR, HER2, grade, education and stratified by menopause, stage, country

exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     

results_list <- list()

for (exposure in exposures) {
  formula <- as.formula(paste("Surv(time_surv, event == levels(Death_Status)[2]) ~ ", exposure,
                              "+ age_dx + bmi_oms + Alc_Re_c + year_diagnosis + 
                               Er_Adj_c + Pr_Adj_c + Her2_Adj_c + gradb + 
                               + l_school2 + 
                               strata(COUNTRY, stage3, menop_status_dx)"))
  model <- coxph(formula, data=set, id=Idepic)
  aic_value <- AIC(model)  # Compute AIC
  
  # Store model results, including AIC
  results_list[[exposure]] <- list(
    "model" = model,            
    "conf.int" = summary(model)$conf.int[1, c(1, 3:4)],
    "n" = summary(model)$n,
    "nevent" = summary(model)$nevent,
    "AIC" = aic_value  # Store AIC inside results_list
  )
}

# Initialize a list to store rows for final table
result_tables <- list()

# Loop through exposures and extract stored results
for (exposure in exposures) {
  cat("Results for exposure: ", exposure, "\n")
  print(results_list[[exposure]]$conf.int) # Print confidence intervals
  
  # Extract results into a data frame
  table <- data.frame(
    HR = round(results_list[[exposure]]$conf.int[1], 2), 
    lowCI = round(results_list[[exposure]]$conf.int[2], 2),
    upCI = round(results_list[[exposure]]$conf.int[3], 2),
    n = results_list[[exposure]]$n,
    nevent = results_list[[exposure]]$nevent,
    AIC = round(results_list[[exposure]]$AIC, 2)  # Retrieve the stored AIC
  )
  
  rownames(table) <- exposure
  result_tables[[exposure]] <- table  # Store table for binding
}

# Combine all results into a final table
fi_table <- do.call(rbind, result_tables)
fi_table <- round(fi_table, 2)
fi_table_1 <- fi_table
fi_table_1

# Model BC-specific ----
# age at diagnostic, BMI, alcohol, year of diagnosis, ER, PR, HER2, grade, education and stratified by menopause, stage, country

exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     

results_list <- list()

for (exposure in exposures) {
  formula <- as.formula(paste("Surv(time_surv, Death_specific =='Dead') ~ ", exposure,
                              "+ age_dx + bmi_oms + Alc_Re_c + year_diagnosis + 
                               Er_Adj_c + Pr_Adj_c + Her2_Adj_c + gradb +  
                               + l_school2 + 
                               strata(COUNTRY, stage3, menop_status_dx)"))
  model <- coxph(formula, data=set, id=Idepic)
  aic_value <- AIC(model)  # Compute AIC
  
  # Store model results, including AIC
  results_list[[exposure]] <- list(
    "model" = model,            
    "conf.int" = summary(model)$conf.int[1, c(1, 3:4)],
    "n" = summary(model)$n,
    "nevent" = summary(model)$nevent,
    "AIC" = aic_value  # Store AIC inside results_list
  )
}

# Initialize a list to store rows for final table
result_tables <- list()

# Loop through exposures and extract stored results
for (exposure in exposures) {
  cat("Results for exposure: ", exposure, "\n")
  print(results_list[[exposure]]$conf.int) # Print confidence intervals
  
  # Extract results into a data frame
  table <- data.frame(
    HR = round(results_list[[exposure]]$conf.int[1], 2), 
    lowCI = round(results_list[[exposure]]$conf.int[2], 2),
    upCI = round(results_list[[exposure]]$conf.int[3], 2),
    n = results_list[[exposure]]$n,
    nevent = results_list[[exposure]]$nevent,
    AIC = round(results_list[[exposure]]$AIC, 2)  # Retrieve the stored AIC
  )
  
  rownames(table) <- exposure
  result_tables[[exposure]] <- table  # Store table for binding
}

# Combine all results into a final table
fi_table <- do.call(rbind, result_tables)
fi_table <- round(fi_table, 2)
fi_table_2 <- fi_table
fi_table_2

# Both OM and BCspecific tables:
final_table <- rbind(fi_table_1, fi_table_2)
final_table
