#-----------------------------------------------#
# Check the linearity assumption using splines
#-----------------------------------------------#

set <- readRDS("/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Data/processed_data/002_OM_models_script.rds")

# Libraries:
library(splines)
library(rms)
library(lmtest)
library(haven)
library(foreign)
library(foreach)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(DescTools)
library(devtools)

# List of exposures
exposures <- c("log_TSH_R_sd","log_fT3_pmol_L_R_sd","log_fT4_pmol_L_R_sd","log_fT3fT4r_R_sd","log_TSH_R_sd")
               # "Anti_TPO_p")     

# Initialize an empty list to store results
results_list <- list()

# Loop through each exposure and fit Cox regression models with linear and spline terms
for (exposure in exposures) {

  formula_linear <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ", exposure,
                                     "+ bmi_oms + Alc_Re_c + gradb + year_diagnosis +
                               Er_Adj_c + Pr_Adj_c + Her2_Adj_c + 
                               strata(COUNTRY, stage3)"))
  
  # Create the formula for the Cox regression model with spline term
  formula_spline <- as.formula(paste("Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ rcs(", exposure, ", df = 4)",
                                     "+ bmi_oms + Alc_Re_c + gradb + year_diagnosis +
                               Er_Adj_c + Pr_Adj_c + Her2_Adj_c + 
                               strata(COUNTRY, stage3)"))
  
  # Fit the Cox regression models
  model_linear <- coxph(formula_linear, data = set, id = Idepic)
  model_spline <- coxph(formula_spline, data = set, id = Idepic)
  
  # Perform ANOVA between the models
  anova_result <- lrtest(model_linear, model_spline)
  
  # Store the results in the list
  results_list[[exposure]] <- list("formula_linear" = formula_linear,
                                   "formula_spline" = formula_spline,
                                   "model_linear" = model_linear,
                                   "model_spline" = model_spline,
                                   "anova_result" = anova_result)
}

# Display ANOVA results for each exposure
for (exposure in exposures) {
  cat("Exposure:", exposure, "\n")
  print(results_list[[exposure]]$anova_result)
  cat("\n")
}


#-------------------------#
# Plot: fT3 splines & OM ----
#------------------------#
# Define the exposure for fT3
exposure <- "log_fT3_pmol_L_R_sd"  # Exposure variable

# Create the linear and spline formulas for the model (for fT3)
formula_linear <- Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ 
  get(exposure) + bmi_oms + Alc_Re_c + gradb + year_diagnosis +
  Er_Adj_c + Pr_Adj_c + Her2_Adj_c + 
  strata(COUNTRY, stage3)

formula_spline <- Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ 
  rcs(get(exposure), df = 4) + bmi_oms + Alc_Re_c + gradb + year_diagnosis +
  Er_Adj_c + Pr_Adj_c + Her2_Adj_c + 
  strata(COUNTRY, stage3)

# Fit the Cox models for both linear and spline (for fT3)
model_linear <- coxph(formula_linear, data = set, id = Idepic)
model_spline <- coxph(formula_spline, data = set, id = Idepic)

# Perform ANOVA between the models for fT3
anova_result <- anova(model_linear, model_spline)

# Extract p-value for the comparison between linear and spline models (for fT3)
p_value_fT3 <- anova_result$`Pr(>|Chi|)`[2]  # Extract p-value for fT3

# Print p-value for fT3
print(p_value_fT3)

# Create a dataset with complete cases for the variable log_fT3_pmol_L_R_sd
set_complete <- set[complete.cases(set[[exposure]]), ]

# Ensure the model is fitted with complete data
pred3 <- predict(model_spline, newdata = set_complete, type = "terms", se = TRUE)

# Extract the fitted values and standard errors
hfit <- pred3$fit[, 1]  # Fitted values
hse <- pred3$se.fit[, 1]  # Standard errors

# Create confidence intervals (CI)
hmat_om <- cbind(hfit, hfit + 2 * hse, hfit - 2 * hse)  # CIs

# Order the rows based on log_fT3_pmol_L_R_sd
o <- order(set_complete[[exposure]], na.last = TRUE)

# Dynamically adjust position for p-value
x_pos <- min(set_complete[[exposure]]) + 0.2 * (max(set_complete[[exposure]]) - min(set_complete[[exposure]]))  # 20% to the right of min x
y_pos <- max(hmat_om) * 0.85  # 85% of the maximum y value for positioning

# Create a new dataframe for ggplot
plot_data <- data.frame(exposure = set_complete[[exposure]][o], 
                        fit = hmat_om[o, 1], 
                        upper = hmat_om[o, 2], 
                        lower = hmat_om[o, 3])

# ggplot version with less visible grid lines and p-value annotation for fT3
plot_fT3 <- ggplot(plot_data, aes(x = exposure, y = fit)) +
  geom_line(color = "lightblue", size = 1.2) +  # Line for the spline fit
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightgrey", alpha = 0.5) +  # Confidence interval shaded region
  labs(x = "fT3", y = "Log HR", title = "fT3 and Overall Mortality") +
  theme_minimal(base_size = 15) +  # Minimal theme with base size for text
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Title settings
    axis.title = element_text(size = 14),  # Axis labels size
    axis.text = element_text(size = 12),  # Axis text size
    panel.grid.major = element_line(color = "grey80", size = 0.5),  # Less visible major grid lines
    panel.grid.minor = element_line(color = "grey90", size = 0.25),  # Less visible minor grid lines
    plot.background = element_rect(fill = "white", color = "white"),  # Set the plot background to white
    panel.background = element_rect(fill = "white", color = "white")  # Set the panel background to white
  ) +
  # Annotate the p-value extracted from anova_result for fT3
  annotate("text", x = x_pos, y = y_pos,  # Position the p-value annotation
           label = paste("p-value: ", round(p_value_fT3, 5)),  # Use extracted p-value for fT3
           color = "darkblue", size = 6, fontface = "bold", hjust = 0)  # p-value in dark blue, bold font

# Print the plot for fT3
print(plot_fT3)

#-------------------------#
# Plot: fT4 splines & OM ----
#------------------------#

# Define the exposure for fT4
exposure <- "log_fT4_pmol_L_R_sd"  # Exposure variable for fT4

# Create the linear and spline formulas for the model (for fT4)
formula_linear <- Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ 
  get(exposure) + bmi_oms + Alc_Re_c + gradb + year_diagnosis +
  Er_Adj_c + Pr_Adj_c + Her2_Adj_c + 
  strata(COUNTRY, stage3)

formula_spline <- Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ 
  rcs(get(exposure), df = 4) + bmi_oms + Alc_Re_c + gradb + year_diagnosis +
  Er_Adj_c + Pr_Adj_c + Her2_Adj_c + 
  strata(COUNTRY, stage3)

# Fit the Cox models for both linear and spline (for fT4)
model_linear <- coxph(formula_linear, data = set, id = Idepic)
model_spline <- coxph(formula_spline, data = set, id = Idepic)

# Perform ANOVA between the models for fT4
anova_result <- anova(model_linear, model_spline)

# Extract p-value for the comparison between linear and spline models (for fT4)
p_value_fT4 <- anova_result$`Pr(>|Chi|)`[2]  # Extract p-value for fT4

# Print p-value for fT4
print(p_value_fT4)

# Create a dataset with complete cases for the variable log_fT4_pmol_L_R_sd
set_complete <- set[complete.cases(set[[exposure]]), ]

# Ensure the model is fitted with complete data
pred3 <- predict(model_spline, newdata = set_complete, type = "terms", se = TRUE)

# Extract the fitted values and standard errors
hfit <- pred3$fit[, 1]  # Fitted values
hse <- pred3$se.fit[, 1]  # Standard errors

# Create confidence intervals (CI)
hmat_om <- cbind(hfit, hfit + 2 * hse, hfit - 2 * hse)  # CIs

# Order the rows based on log_fT4_pmol_L_R_sd
o <- order(set_complete[[exposure]], na.last = TRUE)

# Dynamically adjust position for p-value
x_pos <- min(set_complete[[exposure]]) + 0.2 * (max(set_complete[[exposure]]) - min(set_complete[[exposure]]))  # 20% to the right of min x
y_pos <- max(hmat_om) * 0.85  # 85% of the maximum y value for positioning

# Create a new dataframe for ggplot
plot_data <- data.frame(exposure = set_complete[[exposure]][o], 
                        fit = hmat_om[o, 1], 
                        upper = hmat_om[o, 2], 
                        lower = hmat_om[o, 3])

# ggplot version with less visible grid lines and p-value annotation for fT4
plot_fT4 <- ggplot(plot_data, aes(x = exposure, y = fit)) +
  geom_line(color = "lightblue", size = 1.2) +  # Line for the spline fit
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightgrey", alpha = 0.5) +  # Confidence interval shaded region
  labs(x = "fT4", y = "Log HR", title = "fT4 and Overall Mortality") +
  theme_minimal(base_size = 15) +  # Minimal theme with base size for text
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Title settings
    axis.title = element_text(size = 14),  # Axis labels size
    axis.text = element_text(size = 12),  # Axis text size
    panel.grid.major = element_line(color = "grey80", size = 0.5),  # Less visible major grid lines
    panel.grid.minor = element_line(color = "grey90", size = 0.25),  # Less visible minor grid lines
    plot.background = element_rect(fill = "white", color = "white"),  # Set the plot background to white
    panel.background = element_rect(fill = "white", color = "white")  # Set the panel background to white
  ) +
  # Annotate the p-value extracted from anova_result for fT4
  annotate("text", x = x_pos, y = y_pos,  # Position the p-value annotation
           label = paste("p-value: ", round(p_value_fT4, 5)),  # Use extracted p-value for fT4
           color = "darkblue", size = 6, fontface = "bold", hjust = 0)  # p-value in dark blue, bold font

# Print the plot for fT4
print(plot_fT4)

#--------------------------------------#
# Plot fT4/fT3 ration splines & OM ----
#--------------------------------------#

# Define the exposure for fT3
exposure <- "log_fT3fT4r_R_sd"  # Exposure variable for fT3 and fT4 combined

# Create the linear and spline formulas for the model (for fT3 and fT4 combined)
formula_linear <- Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ 
  get(exposure) + bmi_oms + Alc_Re_c + gradb + year_diagnosis +
  Er_Adj_c + Pr_Adj_c + Her2_Adj_c + 
  strata(COUNTRY, stage3)

formula_spline <- Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ 
  rcs(get(exposure), df = 4) + bmi_oms + Alc_Re_c + gradb + year_diagnosis +
  Er_Adj_c + Pr_Adj_c + Her2_Adj_c + 
  strata(COUNTRY, stage3)

# Fit the Cox models for both linear and spline (for fT3 and fT4 combined)
model_linear <- coxph(formula_linear, data = set, id = Idepic)
model_spline <- coxph(formula_spline, data = set, id = Idepic)

# Perform ANOVA between the models for fT3 and fT4 combined
anova_result <- anova(model_linear, model_spline)

# Extract p-value for the comparison between linear and spline models (for fT3 and fT4 combined)
p_value_fT3fT4 <- anova_result$`Pr(>|Chi|)`[2]  # Extract p-value for fT3 and fT4 combined

# Print p-value for fT3 and fT4 combined
print(p_value_fT3fT4)

# Create a dataset with complete cases for the variable log_fT3fT4r_R_sd
set_complete <- set[complete.cases(set[[exposure]]), ]

# Ensure the model is fitted with complete data
pred3 <- predict(model_spline, newdata = set_complete, type = "terms", se = TRUE)

# Extract the fitted values and standard errors
hfit <- pred3$fit[, 1]  # Fitted values
hse <- pred3$se.fit[, 1]  # Standard errors

# Create confidence intervals (CI)
hmat_om <- cbind(hfit, hfit + 2 * hse, hfit - 2 * hse)  # CIs

# Order the rows based on log_fT3fT4r_R_sd
o <- order(set_complete[[exposure]], na.last = TRUE)

# Dynamically adjust position for p-value
x_pos <- min(set_complete[[exposure]]) + 0.2 * (max(set_complete[[exposure]]) - min(set_complete[[exposure]]))  # 20% to the right of min x
y_pos <- max(hmat_om) * 0.85  # 85% of the maximum y value for positioning

# Create a new dataframe for ggplot
plot_data <- data.frame(exposure = set_complete[[exposure]][o], 
                        fit = hmat_om[o, 1], 
                        upper = hmat_om[o, 2], 
                        lower = hmat_om[o, 3])

# ggplot version with less visible grid lines and p-value annotation for fT3 and fT4 combined
plot_fT3fT4 <- ggplot(plot_data, aes(x = exposure, y = fit)) +
  geom_line(color = "lightblue", size = 1.2) +  # Line for the spline fit
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightgrey", alpha = 0.5) +  # Confidence interval shaded region
  labs(x = "fT3/fT4 ratio", y = "Log HR", title = "fT3/fT4 ratio and Overall Mortality") +
  theme_minimal(base_size = 15) +  # Minimal theme with base size for text
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Title settings
    axis.title = element_text(size = 14),  # Axis labels size
    axis.text = element_text(size = 12),  # Axis text size
    panel.grid.major = element_line(color = "grey80", size = 0.5),  # Less visible major grid lines
    panel.grid.minor = element_line(color = "grey90", size = 0.25),  # Less visible minor grid lines
    plot.background = element_rect(fill = "white", color = "white"),  # Set the plot background to white
    panel.background = element_rect(fill = "white", color = "white")  # Set the panel background to white
  ) +
  # Annotate the p-value extracted from anova_result for fT3 and fT4 combined
  annotate("text", x = x_pos, y = y_pos,  # Position the p-value annotation
           label = paste("p-value: ", round(p_value_fT3fT4, 5)),  # Use extracted p-value for fT3 and fT4 combined
           color = "darkblue", size = 6, fontface = "bold", hjust = 0)  # p-value in dark blue, bold font

# Print the plot for fT3 and fT4 combined
print(plot_fT3fT4)
#---------------------------#
# Plot TSH splines & OM ----
#---------------------------#

# Define the exposure for TSH
exposure <- "log_TSH_R_sd"  # Exposure variable for TSH

# Create the linear and spline formulas for the model (for TSH)
formula_linear <- Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ 
  get(exposure) + bmi_oms + Alc_Re_c + gradb + year_diagnosis +
  Er_Adj_c + Pr_Adj_c + Her2_Adj_c + 
  strata(COUNTRY, stage3)

formula_spline <- Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ 
  rcs(get(exposure), df = 4) + bmi_oms + Alc_Re_c + gradb + year_diagnosis +
  Er_Adj_c + Pr_Adj_c + Her2_Adj_c + 
  strata(COUNTRY, stage3)

# Fit the Cox models for both linear and spline (for TSH)
model_linear <- coxph(formula_linear, data = set, id = Idepic)
model_spline <- coxph(formula_spline, data = set, id = Idepic)

# Perform ANOVA between the models for TSH
anova_result_tsh <- anova(model_linear, model_spline)

# Extract p-value for the comparison between linear and spline models (for TSH)
p_value_tsh <- anova_result_tsh$`Pr(>|Chi|)`[2]  # Extract p-value for TSH

# Print p-value for TSH
print(p_value_tsh)

# Create a dataset with complete cases for the variable log_TSH_R_sd
set_complete_tsh <- set[complete.cases(set[[exposure]]), ]

# Ensure the model is fitted with complete data
pred3_tsh <- predict(model_spline, newdata = set_complete_tsh, type = "terms", se = TRUE)

# Extract the fitted values and standard errors
hfit_tsh <- pred3_tsh$fit[, 1]  # Fitted values
hse_tsh <- pred3_tsh$se.fit[, 1]  # Standard errors

# Create confidence intervals (CI) for TSH
hmat_om_tsh <- cbind(hfit_tsh, hfit_tsh + 2 * hse_tsh, hfit_tsh - 2 * hse_tsh)  # CIs

# Order the rows based on log_TSH_R_sd
o_tsh <- order(set_complete_tsh[[exposure]], na.last = TRUE)

# Dynamically adjust position for p-value
x_pos_tsh <- min(set_complete_tsh[[exposure]]) + 0.2 * (max(set_complete_tsh[[exposure]]) - min(set_complete_tsh[[exposure]]))  # 20% to the right of min x
y_pos_tsh <- max(hmat_om_tsh) * 0.85  # 85% of the maximum y value for positioning

# Create a new dataframe for ggplot (for TSH)
plot_data_tsh <- data.frame(exposure = set_complete_tsh[[exposure]][o_tsh], 
                            fit = hmat_om_tsh[o_tsh, 1], 
                            upper = hmat_om_tsh[o_tsh, 2], 
                            lower = hmat_om_tsh[o_tsh, 3])

# ggplot version with less visible grid lines and p-value annotation for TSH
plot_tsh <- ggplot(plot_data_tsh, aes(x = exposure, y = fit)) +
  geom_line(color = "lightblue", size = 1.2) +  # Line for the spline fit
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightgrey", alpha = 0.5) +  # Confidence interval shaded region
  labs(x = "TSH", y = "Log HR", title = "TSH and Overall Mortality") +
  theme_minimal(base_size = 15) +  # Minimal theme with base size for text
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Title settings
    axis.title = element_text(size = 14),  # Axis labels size
    axis.text = element_text(size = 12),  # Axis text size
    panel.grid.major = element_line(color = "grey80", size = 0.5),  # Less visible major grid lines
    panel.grid.minor = element_line(color = "grey90", size = 0.25),  # Less visible minor grid lines
    plot.background = element_rect(fill = "white", color = "white"),  # Set the plot background to white
    panel.background = element_rect(fill = "white", color = "white")  # Set the panel background to white
  ) +
  # Annotate the p-value extracted from anova_result for TSH
  annotate("text", x = x_pos_tsh, y = y_pos_tsh,  # Position the p-value annotation
           label = paste("p-value: ", round(p_value_tsh, 5)),  # Use extracted p-value for TSH
           color = "darkblue", size = 6, fontface = "bold", hjust = 0)  # p-value in dark blue, bold font

# Print the plot for TSH
print(plot_tsh)

# All plots in one grid ----
grid_plot <- grid.arrange(plot_fT4, plot_ft3, plot_fT3fT4, plot_tsh, ncol = 2)

# save
ggsave("Splines_Plots.png", plot = grid_plot, width = 12, height = 10, dpi = 300)

