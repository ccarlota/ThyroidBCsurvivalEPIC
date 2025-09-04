
#-------------------------------------------#
# Survival models (specific cause of death)
# 21/02/2025
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

# Charge processed data:
set <- readRDS("/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Data/processed_data/002_OM_models_script.rds")

#-------------------------------------#
# Specific-cause deaths ----
#-------------------------------------#
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

#-------------------------------------#
# Models for BC-specific mortality ----
#-------------------------------------#

# List of exposures to loop through
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","Anti_TPO_p")     

# Model (1): Unadjusted including country in strata
# Model (2) Adjusted with variables selected via stepwise regression (BMI, alcohol, grade, ER, PR, HER2, year of diagnosis, country, and stage) with country and stage in strata; 

results_list <- list()

# Fit Cox models and store results
for (exposure in exposures) {
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Csm, Death_specific =='Dead') ~ ", exposure,
                              "+ strata(COUNTRY,stage3,menop_status_dx)"))
  
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
results_list <- list()

# Fit Cox models and store results
for (exposure in exposures) {
  formula <- as.formula(paste("Surv(age_dx, Age_Exit_Csm, Death_specific =='Dead') ~ ", exposure,
                              "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb +
                                  Er_Adj_c + Pr_Adj_c + Her2_Adj_c + l_school2 +
                                  strata(COUNTRY,stage3,menop_status_dx)"))
  
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

final_table <- rbind(fi_table_1, fi_table_2)
final_table

# save:
write.table(final_table, file = paste0(output.dir,"003_BCM.txt"), 
            col.names = T, row.names = T, sep = "\t", quote = F)


