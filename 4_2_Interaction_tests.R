library(survival)

exposures <- c("log_TSH_R_sd", "log_fT3_pmol_L_R_sd", "log_fT4_pmol_L_R_sd", "log_fT3fT4r_R_sd", "Anti_TPO_p")

results_list <- list()
lrt_results <- data.frame(Exposure=character(), LRT_pval=numeric(), stringsAsFactors=FALSE)

# Fit Cox models and store results
for (exposure in exposures) {
  
  # Model without interaction
  formula_no_int <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ", exposure,
                                     "+ age_dx_2 + bmi_oms + Alc_Re_c + year_diagnosis + gradb + 
                                     Er_Adj_c + Pr_Adj_c + Her2_Adj_c + strata(COUNTRY,stage3)"))
  
  # Model with interaction
  formula_int <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ", exposure,
                                  "*age_dx_2 + age_dx_2 + bmi_oms + Alc_Re_c + year_diagnosis + gradb + 
                                   Er_Adj_c + Pr_Adj_c + Her2_Adj_c + strata(COUNTRY,stage3)"))
  
  model_no_int <- coxph(formula_no_int, data=set)
  model_int <- coxph(formula_int, data=set)
  
  # Compute LRT
  lrt_test <- anova(model_no_int, model_int, test="LRT")
  lrt_pval <- lrt_test$`Pr(>|Chi|)`[2]  # Extract p-value
  
  # Extract AIC
  aic_value <- AIC(model_int)
  
  # Extract HR and CIs for the exposure, ensuring correct indexing
  conf_int <- NA  # Default to NA in case the exposure is not found
  conf_table <- summary(model_int)$conf.int
  
  if (exposure %in% rownames(conf_table)) {
    conf_int <- conf_table[exposure, c(1, 3:4)]
  }
  
  # Store results
  results_list[[exposure]] <- list(
    "model" = model_int,            
    "conf.int" = conf_int,  
    "n" = summary(model_int)$n,
    "nevent" = summary(model_int)$nevent,
    "AIC" = aic_value
  )
  
  # Store LRT results
  lrt_results <- rbind(lrt_results, data.frame(Exposure=exposure, LRT_pval=round(lrt_pval, 4)))
}

# Print LRT results
print(lrt_results)

rownames(table) <- exposure
result_tables[[exposure]] <- table  # Store table for binding


# Combine all results into a final table
fi_table <- do.call(rbind, result_tables)
fi_table <- round(fi_table, 2)

# Print results
print(fi_table)
print("Likelihood Ratio Test results:")
print(lrt_results)
