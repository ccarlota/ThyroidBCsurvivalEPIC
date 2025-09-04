#-----------------------------------------------------------------------------------#
# Analyses of OM and BC-specific mortality by tertiles/quartiles of thyroid hormones
#-----------------------------------------------------------------------------------#
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
library(survival)

# Directories
data.dir <- "/data/Epic/subprojects/Breast_Cancer/files/Maestra/"
output.dir <- "/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Analysis/Tables/"

# Charge processed data:
set <- readRDS("/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Data/processed_data/002_OM_models_script.rds")

horm <- c("log_TSH_R","log_fT3_pmol_L_R","log_fT4_pmol_L_R","log_fT3fT4r_R")     

for (i in horm) {
  
  temp_data <- set[!is.na(set[[i]]), ]
  
  T_breaks <- unique(quantile(temp_data[[i]], probs = seq(0, 1, 1/3), na.rm = TRUE))
  # Q_breaks <- unique(quantile(temp_data[[i]], probs = seq(0, 1, 1/4), na.rm = TRUE))
  
  if (length(T_breaks) > 2) {
    set[, paste0(i, "_T")] <- cut(set[, i], T_breaks, include.lowest = TRUE, labels = c("T1", "T2", "T3"))
  } else {
    set[, paste0(i, "_T")] <- NA
  }
}

#-------------------------------------------------------#
# All-cause mortality by Tertiles UNADJUSTED model ----
#------------------------------------------------------#
exposures <- c("log_TSH_R_T","log_fT3_pmol_L_R_T","log_fT4_pmol_L_R_T","log_fT3fT4r_R_T")     

results_list <- list()

# Create a copy of your dataset (if needed) or work directly with 'set'
for (exposure in exposures) {
  
  # First, fit your original model by tertiles (the categorical variable)
  formula_cat <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ",
                                  exposure, "+ strata(COUNTRY,stage3,menop_status_dx)"))
  model_cat <- coxph(formula_cat, data = set, id = Idepic)
  
  results_list[[exposure]] <- list(
    model = model_cat,
    conf.int = summary(model_cat)$conf.int[1:2, c(1, 3:4)],
    n = summary(model_cat)$n,
    nevent = summary(model_cat)$nevent
  )
  
  # Now, create an ordinal (numeric) trend variable from your tertile factor.
  trend_var_name <- paste0(exposure, "_trend")
  set[[trend_var_name]] <- with(set, ifelse(get(exposure) == "T1", 1,
                                            ifelse(get(exposure) == "T2", 2,
                                                   ifelse(get(exposure) == "T3", 3, NA))))
  
  # Fit a trend Cox model using the numeric trend variable.
  formula_trend <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ",
                                    trend_var_name, "+ strata(COUNTRY,stage3,menop_status_dx)"))
  model_trend <- coxph(formula_trend, data = set, id = Idepic)
  
  # Extract the p-value for the trend from the model summary.
  # The p-value for the continuous predictor (trend variable) is in the coefficients table.
  p_trend <- summary(model_trend)$coefficients[1, "Pr(>|z|)"]
  
  # Store the p_trend in your results list
  results_list[[exposure]]$p_trend <- p_trend
}

# Initialize a list to store rows for final table with p_trend
result_tables <- list()

for (exposure in exposures) {
  cat("Results for exposure: ", exposure, "\n")
  print(results_list[[exposure]]$conf.int)  # Print confidence intervals
  
  # Extract results into a data frame and include p_trend
  table <- data.frame(
    HR_T2 = round(results_list[[exposure]]$conf.int[1, 1], 2), 
    lowCI_T2 = round(results_list[[exposure]]$conf.int[1, 2], 2),
    upCI_T2 = round(results_list[[exposure]]$conf.int[1, 3], 2),
    HR_T3 = round(results_list[[exposure]]$conf.int[2, 1], 2), 
    lowCI_T3 = round(results_list[[exposure]]$conf.int[2, 2], 2),
    upCI_T3 = round(results_list[[exposure]]$conf.int[2, 3], 2),
    p_trend = round(results_list[[exposure]]$p_trend, 3),
    n = results_list[[exposure]]$n,
    nevent = results_list[[exposure]]$nevent
  )
  
  rownames(table) <- exposure
  result_tables[[exposure]] <- table  # Store table for binding
}

# Combine all results into a final table
fi_table <- do.call(rbind, result_tables)
fi_table <- round(fi_table, 2)
fi_table

write.table(fi_table, file = paste0(output.dir,"003_BC_tertiles.txt"), 
            col.names = T, row.names = T, sep = "\t", quote = F)

#---------------------------------------------------#
# All-cause mortality by Tertiles Model Adjusted ----
#-----------------------------------------------------#
exposures <- c("log_TSH_R_T","log_fT3_pmol_L_R_T","log_fT4_pmol_L_R_T","log_fT3fT4r_R_T")     

results_list <- list()

# Create a copy of your dataset (if needed) or work directly with 'set'
for (exposure in exposures) {
  
  # First, fit your original model by tertiles (the categorical variable)
  formula_cat <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ",
                                  exposure, "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb +
                                  Er_Adj_c + Pr_Adj_c + Her2_Adj_c + l_school2 +
                                  strata(COUNTRY,stage3,menop_status_dx)"))
  model_cat <- coxph(formula_cat, data = set, id = Idepic)
  
  results_list[[exposure]] <- list(
    model = model_cat,
    conf.int = summary(model_cat)$conf.int[1:2, c(1, 3:4)],
    n = summary(model_cat)$n,
    nevent = summary(model_cat)$nevent
  )
  
  # Now, create an ordinal (numeric) trend variable from your tertile factor.
  trend_var_name <- paste0(exposure, "_trend")
  set[[trend_var_name]] <- with(set, ifelse(get(exposure) == "T1", 1,
                                            ifelse(get(exposure) == "T2", 2,
                                                   ifelse(get(exposure) == "T3", 3, NA))))
  
  # Fit a trend Cox model using the numeric trend variable.
  formula_trend <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ",
                                    trend_var_name, "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb +
                                  Er_Adj_c + Pr_Adj_c + Her2_Adj_c + l_school2 +
                                  strata(COUNTRY,stage3,menop_status_dx)"))
  
  model_trend <- coxph(formula_trend, data = set, id = Idepic)
  
  # Extract the p-value for the trend from the model summary.
  # The p-value for the continuous predictor (trend variable) is in the coefficients table.
  p_trend <- summary(model_trend)$coefficients[1, "Pr(>|z|)"]
  
  # Store the p_trend in your results list
  results_list[[exposure]]$p_trend <- p_trend
}

# Initialize a list to store rows for final table with p_trend
result_tables <- list()

for (exposure in exposures) {
  cat("Results for exposure: ", exposure, "\n")
  print(results_list[[exposure]]$conf.int)  # Print confidence intervals
  
  # Extract results into a data frame and include p_trend
  table <- data.frame(
    HR_T2 = round(results_list[[exposure]]$conf.int[1, 1], 2), 
    lowCI_T2 = round(results_list[[exposure]]$conf.int[1, 2], 2),
    upCI_T2 = round(results_list[[exposure]]$conf.int[1, 3], 2),
    HR_T3 = round(results_list[[exposure]]$conf.int[2, 1], 2), 
    lowCI_T3 = round(results_list[[exposure]]$conf.int[2, 2], 2),
    upCI_T3 = round(results_list[[exposure]]$conf.int[2, 3], 2),
    p_trend = round(results_list[[exposure]]$p_trend, 3),
    n = results_list[[exposure]]$n,
    nevent = results_list[[exposure]]$nevent
  )
  
  rownames(table) <- exposure
  result_tables[[exposure]] <- table  # Store table for binding
}

# Combine all results into a final table
fi_table <- do.call(rbind, result_tables)
fi_table <- round(fi_table, 2)
fi_table

write.table(fi_table, file = paste0(output.dir,"003_BC_tertiles.txt"), 
            col.names = T, row.names = T, sep = "\t", quote = F)

#---------------------------------------#
# BC-specific mortality by Tertiles UNADJUSTED model ----
#---------------------------------------#

# Fit Cox models and store results
results_list <- list()

# Create a copy of your dataset (if needed) or work directly with 'set'
for (exposure in exposures) {
  
  # First, fit your original model by tertiles (the categorical variable)
  formula_cat <- as.formula(paste("Surv(age_dx, Age_Exit_Csm, Death_specific =='Dead') ~ ",
                                  exposure, "+ strata(COUNTRY,stage3,menop_status_dx)"))
  model_cat <- coxph(formula_cat, data = set, id = Idepic)
  
  results_list[[exposure]] <- list(
    model = model_cat,
    conf.int = summary(model_cat)$conf.int[1:2, c(1, 3:4)],
    n = summary(model_cat)$n,
    nevent = summary(model_cat)$nevent
  )
  
  # Now, create an ordinal (numeric) trend variable from your tertile factor.
  trend_var_name <- paste0(exposure, "_trend")
  set[[trend_var_name]] <- with(set, ifelse(get(exposure) == "T1", 1,
                                            ifelse(get(exposure) == "T2", 2,
                                                   ifelse(get(exposure) == "T3", 3, NA))))
  
  # Fit a trend Cox model using the numeric trend variable.
  formula_trend <- as.formula(paste("Surv(age_dx, Age_Exit_Csm, Death_specific =='Dead') ~ ",
                                    trend_var_name, "+ strata(COUNTRY,stage3,menop_status_dx)"))
  model_trend <- coxph(formula_trend, data = set, id = Idepic)
  
  # Extract the p-value for the trend from the model summary.
  # The p-value for the continuous predictor (trend variable) is in the coefficients table.
  p_trend <- summary(model_trend)$coefficients[1, "Pr(>|z|)"]
  
  # Store the p_trend in your results list
  results_list[[exposure]]$p_trend <- p_trend
}

# Initialize a list to store rows for final table with p_trend
result_tables <- list()

for (exposure in exposures) {
  cat("Results for exposure: ", exposure, "\n")
  print(results_list[[exposure]]$conf.int)  # Print confidence intervals
  
  # Extract results into a data frame and include p_trend
  table <- data.frame(
    HR_T2 = round(results_list[[exposure]]$conf.int[1, 1], 2), 
    lowCI_T2 = round(results_list[[exposure]]$conf.int[1, 2], 2),
    upCI_T2 = round(results_list[[exposure]]$conf.int[1, 3], 2),
    HR_T3 = round(results_list[[exposure]]$conf.int[2, 1], 2), 
    lowCI_T3 = round(results_list[[exposure]]$conf.int[2, 2], 2),
    upCI_T3 = round(results_list[[exposure]]$conf.int[2, 3], 2),
    p_trend = round(results_list[[exposure]]$p_trend, 3),
    n = results_list[[exposure]]$n,
    nevent = results_list[[exposure]]$nevent
  )
  
  rownames(table) <- exposure
  result_tables[[exposure]] <- table  # Store table for binding
}

# Combine all results into a final table
fi_table <- do.call(rbind, result_tables)
fi_table <- round(fi_table, 2)
fi_table

write.table(fi_table, file = paste0(output.dir,"003_BC_tertiles.txt"), 
            col.names = T, row.names = T, sep = "\t", quote = F)

# Crear una lista vacía para guardar las filas en formato largo
long_format_table <- list()

for (exposure in exposures) {
  res <- results_list[[exposure]]
  
  n_event <- paste0(res$n, " (", res$nevent, ")")
  
  # Fila para T2
  row_T2 <- data.frame(
    Exposure = exposure,
    Tertile = "T2",
    HR = round(res$conf.int[1, 1], 2),
    CI_low = round(res$conf.int[1, 2], 2),
    CI_up = round(res$conf.int[1, 3], 2),
    p_trend = NA,
    N_event = n_event
  )
  
  # Fila para T3
  row_T3 <- data.frame(
    Exposure = exposure,
    Tertile = "T3",
    HR = round(res$conf.int[2, 1], 2),
    CI_low = round(res$conf.int[2, 2], 2),
    CI_up = round(res$conf.int[2, 3], 2),
    p_trend = NA,
    N_event = n_event
  )
  
  # Fila para p-trend
  row_trend <- data.frame(
    Exposure = exposure,
    Tertile = "p-trend",
    HR = NA,
    CI_low = NA,
    CI_up = NA,
    p_trend = round(res$p_trend, 3),
    N_event = n_event
  )
  
  # Añadir las tres filas
  long_format_table[[exposure]] <- rbind(row_T2, row_T3, row_trend)
}

# Combinar todo en una tabla
final_table_long <- do.call(rbind, long_format_table)
rownames(final_table_long) <- NULL

# Mostrar la tabla
final_table_long
write.xlsx(final_table_long, file = "cox_results_by_tertile.xlsx", rowNames = FALSE)

#-----------------------------------------------------#
# BC-specific mortality by Tertiles Model Adjusted ----
#-----------------------------------------------------#
results_list <- list()

# Create a copy of your dataset (if needed) or work directly with 'set'
for (exposure in exposures) {
  
  # First, fit your original model by tertiles (the categorical variable)
  formula_cat <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_specific =='Dead') ~ ",
                                  exposure, "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb +
                                  Er_Adj_c + Pr_Adj_c + Her2_Adj_c + l_school2 +
                                  strata(COUNTRY,stage3,menop_status_dx)"))
  model_cat <- coxph(formula_cat, data = set, id = Idepic)
  
  results_list[[exposure]] <- list(
    model = model_cat,
    conf.int = summary(model_cat)$conf.int[1:2, c(1, 3:4)],
    n = summary(model_cat)$n,
    nevent = summary(model_cat)$nevent
  )
  
  # Now, create an ordinal (numeric) trend variable from your tertile factor.
  trend_var_name <- paste0(exposure, "_trend")
  set[[trend_var_name]] <- with(set, ifelse(get(exposure) == "T1", 1,
                                            ifelse(get(exposure) == "T2", 2,
                                                   ifelse(get(exposure) == "T3", 3, NA))))
  
  # Fit a trend Cox model using the numeric trend variable.
  formula_trend <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, Death_specific =='Dead') ~ ",
                                    trend_var_name, "+ bmi_oms + Alc_Re_c + year_diagnosis + gradb +
                                  Er_Adj_c + Pr_Adj_c + Her2_Adj_c + l_school2 +
                                  strata(COUNTRY,stage3,menop_status_dx)"))
  
  model_trend <- coxph(formula_trend, data = set, id = Idepic)
  
  # Extract the p-value for the trend from the model summary.
  # The p-value for the continuous predictor (trend variable) is in the coefficients table.
  p_trend <- summary(model_trend)$coefficients[1, "Pr(>|z|)"]
  
  # Store the p_trend in your results list
  results_list[[exposure]]$p_trend <- p_trend
}

# Initialize a list to store rows for final table with p_trend
result_tables <- list()

for (exposure in exposures) {
  cat("Results for exposure: ", exposure, "\n")
  print(results_list[[exposure]]$conf.int)  # Print confidence intervals
  
  # Extract results into a data frame and include p_trend
  table <- data.frame(
    HR_T2 = round(results_list[[exposure]]$conf.int[1, 1], 2), 
    lowCI_T2 = round(results_list[[exposure]]$conf.int[1, 2], 2),
    upCI_T2 = round(results_list[[exposure]]$conf.int[1, 3], 2),
    HR_T3 = round(results_list[[exposure]]$conf.int[2, 1], 2), 
    lowCI_T3 = round(results_list[[exposure]]$conf.int[2, 2], 2),
    upCI_T3 = round(results_list[[exposure]]$conf.int[2, 3], 2),
    p_trend = round(results_list[[exposure]]$p_trend, 3),
    n = results_list[[exposure]]$n,
    nevent = results_list[[exposure]]$nevent
  )
  
  rownames(table) <- exposure
  result_tables[[exposure]] <- table  # Store table for binding
}

# Combine all results into a final table
fi_table <- do.call(rbind, result_tables)
fi_table <- round(fi_table, 2)
fi_table

long_format_table <- list()

for (exposure in exposures) {
  res <- results_list[[exposure]]
  
  n_event <- paste0(res$n, " (", res$nevent, ")")
  
  # Fila para T2
  row_T2 <- data.frame(
    Exposure = exposure,
    Tertile = "T2",
    HR = round(res$conf.int[1, 1], 2),
    CI_low = round(res$conf.int[1, 2], 2),
    CI_up = round(res$conf.int[1, 3], 2),
    p_trend = NA,
    N_event = n_event
  )
  
  # Fila para T3
  row_T3 <- data.frame(
    Exposure = exposure,
    Tertile = "T3",
    HR = round(res$conf.int[2, 1], 2),
    CI_low = round(res$conf.int[2, 2], 2),
    CI_up = round(res$conf.int[2, 3], 2),
    p_trend = NA,
    N_event = n_event
  )
  
  # Fila para p-trend
  row_trend <- data.frame(
    Exposure = exposure,
    Tertile = "p-trend",
    HR = NA,
    CI_low = NA,
    CI_up = NA,
    p_trend = round(res$p_trend, 3),
    N_event = n_event
  )
  
  # Añadir las tres filas
  long_format_table[[exposure]] <- rbind(row_T2, row_T3, row_trend)
}

# Combinar todo en una tabla
final_table_long <- do.call(rbind, long_format_table)
rownames(final_table_long) <- NULL

# Mostrar la tabla
final_table_long
write.xlsx(final_table_long, file = "cox_results_by_tertile.xlsx", rowNames = FALSE)
