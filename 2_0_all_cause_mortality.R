#-------------------------#
# Survival models 
# 21/02/2025
#-------------------------#
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

# Hormones per 1-SD increase ----
horm <- colnames(set)[grep("_R$", colnames(set))]
horm # "log_TSH_R" "log_fT3_pmol_L_R" "log_fT4_pmol_L_R" "log_fT3fT4r_R" (R stands for Residuals from Batch effect!)

for (i in horm) {
  horm_sd <- paste0(i, "_sd")
  x <- data.frame(sd = set[, i]/sd(na.omit(set[, i])))
  colnames(x) <- horm_sd
  set <- cbind(set, x)
}
colnames(set)

#--------------------------------#
# All-cause mortality models ----
#--------------------------------#

# Model (1): Unadjusted including country in strata
# Model (2) Adjusted with variables selected via stepwise regression (BMI, alcohol, grade, ER, PR, HER2, year of diagnosis, country, and stage) with country and stage in strata; 
# Model (3) further adjusted by adding menopausal status at diagnosis, education, physical activity, and smoking.

# (1) Unadjusted model (already adjusted for batch) ----
# List of exposures to loop through
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     

results_list <- list()

# Fit Cox models and store results
for (exposure in exposures) {
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ", exposure,
                              "+ strata(COUNTRY, stage3, menop_status_dx)"))
  
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
fi_table <- round(fi_table,2)
fi_table
fi_table_1 <- fi_table

# (2) Model with variables selected via stepwise regression ----

# Exposures to test
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     

# Initialize a list to store retained variables per exposure
stepwise_results <- list()

for (exp in exposures) {
  # Subset to rows with complete data for this exposure
  set_sub <- set[complete.cases(set[[exp]]), ]
  
  # Full model including the exposure + covariates
  formula <- as.formula(paste(
    "Surv(age_dx, Age_Exit_Vs, Death_Status=='Dead') ~",
    exp,
    "+ l_school2 + bmi_oms + Pa_Index_c + Alc_Re_c + gradb + year_diagnosis + Phrt_Bld_c +",
    "Smoke_Stat_cc + Er_Adj_c + Pr_Adj_c + Her2_Adj_c + COUNTRY + menop_status_dx + stage3"
  ))
  
  model <- coxph(formula, data = set_sub, id = Idepic)
  
  # Stepwise selection
  stepwise_model <- step(model, direction = "both", trace = 0, k = 2)
  
  # Extract retained predictors (excluding the outcome)
  retained <- attr(terms(stepwise_model), "term.labels")
  
  # Store in list
  stepwise_results[[exp]] <- retained
}

# Convert the list to a tidy table
stepwise_summary <- tibble(
  Exposure = names(stepwise_results),
  Retained_Covariates = sapply(stepwise_results, function(x) paste(x, collapse = ", "))
)

# View the summary
stepwise_summary

  

# (2) Model incorporating stepwise retained variables + menop status at dx and level education:
results_list <- list()

# Fit Cox models and store results
for (exposure in exposures) {
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ", exposure,
                              "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb + l_school2 +
                              Er_Adj_c + Pr_Adj_c + Her2_Adj_c + strata(COUNTRY, stage3, menop_status_dx)"))
  
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
fi_table <- round(fi_table,2)
fi_table
fi_table_2 <- fi_table

# (3) Further adjustments ----
# Education, menopausal status at dx, use of exogenous hormones at blood (Phrt_Bld_c)
#------------------------------------------------------------------------------------------------------#
results_list <- list()

# Fit Cox models and store results
for (exposure in exposures) {
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ", exposure,
                              "+ bmi_oms + Alc_Re_c + year_diagnosis + 
                               Er_Adj_c + Pr_Adj_c + Her2_Adj_c + gradb + 
                               + l_school2 + strata(COUNTRY, stage3, menop_status_dx)"))
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
fi_table_3 <- fi_table
fi_table_3

## Save this dataset to start in 002 script:
saveRDS(set, "/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Data/processed_data/002_OM_models_script.rds")

final_table <- rbind(fi_table_1, fi_table_2, fi_table_3)
final_table

# save:
write.table(final_table, file = paste0(output.dir,"002_OM.txt"), 
            col.names = T, row.names = T, sep = "\t", quote = F)
