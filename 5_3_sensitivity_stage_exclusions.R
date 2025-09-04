# Excluding "unknown" stage from main analyses
# ---------------------------------------------#

set_stage <- set %>% filter(stage2 %in% c("Metastatic", "Localised"))
dim(set_stage)

set_nonmetas <- set %>% filter(stage3 %in% c("Non-Metastatic"))
dim(set_nonmetas)

results_list <- list()

# Set: Metastatic + Localised (set_stage) using stage2.2 as covariate ----

# Fit Cox models and store results
for (exposure in exposures) {
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ", exposure,
                              "+ bmi_oms + Alc_Re_c + year_diagnosis + 
                               Er_Adj_c + Pr_Adj_c + Her2_Adj_c + gradb + 
                               + l_school2 + strata(COUNTRY, stage2.2, menop_status_dx)"))
  model <- coxph(formula, data=set_stage, id=Idepic)
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

# Set: Metastatic + Non-metastatic (set_stage) using stage2.2 as covariate ----
set$stage3.1 <- ifelse(set$stage3 == "Non-Metastatic", "Non-Metastatic",
                       ifelse(set$stage3 == "Metastatic", "Metastatic", NA))
set$stage3.1 <- as.factor(set$stage3.1)

set_metas <- set %>% filter(stage3 == "Metastatic") 
set_nonmetas <- set %>% filter(stage3 == "Non-Metastatic")

set_stage_met_nonmet <- set %>% filter(stage3 %in% c("Metastatic", "Non-Metastatic"))
dim(set_stage_met_nonmet)

results_list <- list()

# Fit Cox models and store results
for (exposure in exposures) {
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ", exposure,
                              "+ bmi_oms + Alc_Re_c + year_diagnosis + 
                               Er_Adj_c + Pr_Adj_c + Her2_Adj_c + gradb + 
                               + l_school2 + strata(COUNTRY, stage3.1, menop_status_dx)"))
  model <- coxph(formula, data=set_stage_met_nonmet, id=Idepic)
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

