# Test of Cox proportional hazards assumption
#---------------------------------------------#

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
library(survminer)

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
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     

results_list <- list()

for (exposure in exposures) {
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ", exposure,
                              "+ bmi_oms + Alc_Re_c + gradb + year_diagnosis + 
                               Er_Adj_c + Pr_Adj_c + Her2_Adj_c + 
                               + l_school2 + COUNTRY + stage3 + menop_status_dx"))
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

# PH test results for OM ----
ph_test_results <- list()

# Create a folder to save plots (optional)
dir.create("PH_plots", showWarnings = FALSE)

for (exposure in exposures) {
  model <- results_list[[exposure]]$model
  
  # Test proportional hazards assumption
  ph_test <- cox.zph(model, transform = "rank")  # 'rank' es lo más común
  ph_test_results[[exposure]] <- ph_test
  
  # Print summary in console
  cat("\nProportional hazards test for exposure:", exposure, "\n")
  print(ph_test)
  
  # Plot PH test (Schoenfeld residuals)
  png(filename = paste0("PH_plots/PH_", exposure, ".png"), width = 800, height = 600)
  plot(ph_test, main = paste("PH test plots for", exposure))
  dev.off()
}

# GLOBAL and individual variables: >0.05: no violation of proportional hazards assumption

# Model BC-specific ----
# Specific-cause deaths
table(set$Death_Status_Csm, useNA="always") #217 

set %>%
  filter(Death_Status_Csm == "Dead") %>%
  select(Idepic, starts_with("C_Death_")) %>%
  filter(!apply(., 1, function(x) any(startsWith(as.character(x), "C5")))) 

# Death_specific: 0 = Alive, 1 = Specific death (BC-related), 2 = Other death
set$Death_specific <- if_else(set$Death_Status_Csm == "Alive", "Alive", 
                              if_else(set$Death_Status_Csm == "Dead" & 
                                        set$C_Death_U %in% c('C500','C501','C502','C503','C504','C505','C506','C508','C509'), 
                                      "Dead", "Othercauses")
)

table(set$Death_specific, useNA="always")
#   Alive Dead  Othercauses  <NA> 
#   1304  161   56            0 
set$Death_specific <- as.factor(set$Death_specific)

# Exposures
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     

results_list <- list()

for (exposure in exposures) {
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Csm, Death_specific =='Dead') ~ ", exposure,
                              "+ bmi_oms + Alc_Re_c + gradb + year_diagnosis + 
                               Er_Adj_c + Pr_Adj_c + Her2_Adj_c + 
                               + l_school2 + COUNTRY + stage3 + menop_status_dx"))
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

# PH test results for BCspecific ----
ph_test_results <- list()

# Create a folder to save plots (now for bc-specific)
dir.create("PH_plots_specific", showWarnings = FALSE)

for (exposure in exposures) {
  model <- results_list[[exposure]]$model
  
  # Test proportional hazards assumption
  ph_test <- cox.zph(model, transform = "rank")  # 'rank' es lo más común
  ph_test_results[[exposure]] <- ph_test
  
  # Print summary in console
  cat("\nProportional hazards test for exposure:", exposure, "\n")
  print(ph_test)
  
  # Plot PH test (Schoenfeld residuals)
  png(filename = paste0("PH_plots_specific/PH_", exposure, ".png"),
      width = 800, height = 600)
  plot(ph_test, main = paste("PH test plots for", exposure))
  dev.off()
}

# test
test_model <- results_list[["Anti_TPO_p"]]$model
cox.zph(test_model)  

library(car)

# Build a simple linear model with the same predictors
lm_check <- lm(as.numeric(event) ~ log_TSH_R_sd + bmi_oms + Alc_Re_c + gradb +
                 year_diagnosis + Er_Adj_c + Pr_Adj_c + Her2_Adj_c +
                 l_school2 + COUNTRY + stage3 + menop_status_dx,
               data = set)

vif(lm_check)
