exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     

results_list <- list()

# Fit Cox models and store results
for (exposure in exposures) {
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_specific =='Dead') ~ ", exposure,
                              "+ bmi_oms + Alc_Re_c + gradb + year_diagnosis + log_fT4_pmol_L_R +
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
fi_table_3 <- fi_table
fi_table_3
