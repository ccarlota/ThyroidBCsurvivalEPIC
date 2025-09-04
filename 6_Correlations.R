#------------------------------------#
# Correlations among thyroid hormones
#------------------------------------#

library(ggplot2)
library(GGally)
library(ggcorrplot)
library(dplyr)
library(ppcor)
library(tidyr)

complete_data <- set[complete.cases(set[c("log_TSH", "log_fT3_pmol_L", "log_fT4_pmol_L", "log_fT3fT4r", "log_Anti_TPO")]), ]
dim(complete_data)

# I will use residuals of hormones regressed on age at blood, batch and country to perform partial correlations.

# Define the variables
hormones <- c("log_fT3_pmol_L", "log_fT4_pmol_L", "log_fT3fT4r")
adjust_vars <- c("Batch_Ft", "COUNTRY", "Age_Blood")  # Variables to adjust for
hormones1 <- c("log_TSH", "log_Anti_TPO")
adjust_vars1 <- c("Batch_T", "COUNTRY", "Age_Blood")  # Different variables to adjust for

# Ensure that Batch_Ft is treated as a factor
complete_data$Batch_Ft_f <- as.factor(complete_data$Batch_Ft)
complete_data$Batch_T_f <- as.factor(complete_data$Batch_T)

# Loop over each hormone in the first set ('hormones') to perform regression and get residuals
for (hormone in hormones) {
  # Create the formula for the regression
  formula <- as.formula(paste(hormone, "~", paste(adjust_vars, collapse = "+")))
  
  # Fit the linear model
  model <- lm(formula, data = complete_data)
  
  # Extract the residuals (adjusted values) and add them to the dataframe
  complete_data[[paste(hormone, "_adj", sep = "")]] <- residuals(model)
}

# Loop over each hormone in the second set ('hormones1') to perform regression and get residuals
for (hormone in hormones1) {
  # Create the formula for the regression
  formula <- as.formula(paste(hormone, "~", paste(adjust_vars1, collapse = "+")))
  
  # Fit the linear model
  model <- lm(formula, data = complete_data)
  
  # Extract the residuals (adjusted values) and add them to the dataframe
  complete_data[[paste(hormone, "_adj", sep = "")]] <- residuals(model)
}

#--------------------------------------------------------#
# Correlation matrices for all and pre/postmenopausal ----
#--------------------------------------------------------#

hormones <- c("log_fT3_pmol_L_adj", "log_fT4_pmol_L_adj", "log_fT3fT4r_adj", "log_TSH_adj", "log_Anti_TPO_adj")
cor_all <- cor(complete_data[, hormones], use = "complete.obs")
cor_all

complete_data_premenopausal <- subset(complete_data, menop_status_dx == "Premenopausal")
complete_data_postmenopausal <- subset(complete_data, menop_status_dx == "Postmenopausal")

# Compute the correlation matrix for premenopausal group
cor_premenopausal <- cor(complete_data_premenopausal[, hormones], use = "complete.obs")

# Compute the correlation matrix for postmenopausal group
cor_postmenopausal <- cor(complete_data_postmenopausal[, hormones], use = "complete.obs")

print("Correlation Matrix for Premenopausal Group:")
print(cor_premenopausal)

print("Correlation Matrix for Postmenopausal Group:")
print(cor_postmenopausal)

# Now to calculate pvalues ----
#------------------------------#

# Define the hormones you are interested in
hormones <- c("log_fT3_pmol_L", "log_fT4_pmol_L", "log_fT3fT4r", "log_TSH", "log_Anti_TPO")

# Function to calculate correlation and p-value matrix
get_correlation_matrix_with_pvalues <- function(data, hormones) {
  # Initialize matrices to store correlations and p-values
  cor_matrix <- matrix(NA, nrow = length(hormones), ncol = length(hormones))
  pvalue_matrix <- matrix(NA, nrow = length(hormones), ncol = length(hormones))
  
  # Loop over each pair of hormones to calculate correlation and p-value
  for (i in 1:length(hormones)) {
    for (j in 1:length(hormones)) {
      if (i != j) {
        # Perform the correlation test
        test_result <- cor.test(data[[hormones[i]]], data[[hormones[j]]], use = "complete.obs")
        
        # Store the correlation coefficient and p-value in the matrices
        cor_matrix[i, j] <- test_result$estimate
        pvalue_matrix[i, j] <- test_result$p.value
      } else {
        # Correlation of a variable with itself is 1
        cor_matrix[i, j] <- 1
        pvalue_matrix[i, j] <- NA  # No p-value for self-correlation
      }
    }
  }
  
  # Set row and column names
  rownames(cor_matrix) <- colnames(cor_matrix) <- hormones
  rownames(pvalue_matrix) <- colnames(pvalue_matrix) <- hormones
  
  return(list(cor_matrix = cor_matrix, pvalue_matrix = pvalue_matrix))
}

# Compute for all cases
all_results <- get_correlation_matrix_with_pvalues(complete_data, hormones)
cor_all <- all_results$cor_matrix
pvalues_all <- all_results$pvalue_matrix

# Subset the data for premenopausal and postmenopausal groups
complete_data_premenopausal <- subset(complete_data, menop_status_dx == "Premenopausal")
complete_data_postmenopausal <- subset(complete_data, menop_status_dx == "Postmenopausal")

# Compute for premenopausal group
premenopausal_results <- get_correlation_matrix_with_pvalues(complete_data_premenopausal, hormones)
cor_premenopausal <- premenopausal_results$cor_matrix
pvalues_premenopausal <- premenopausal_results$pvalue_matrix

# Compute for postmenopausal group
postmenopausal_results <- get_correlation_matrix_with_pvalues(complete_data_postmenopausal, hormones)
cor_postmenopausal <- postmenopausal_results$cor_matrix
pvalues_postmenopausal <- postmenopausal_results$pvalue_matrix

# View the correlation matrices and p-values
print("Correlation Matrix for All Cases:")
print(cor_all)
print("P-value Matrix for All Cases:")
print(pvalues_all)

print("Correlation Matrix for Premenopausal Group:")
print(cor_premenopausal)
print("P-value Matrix for Premenopausal Group:")
print(pvalues_premenopausal)

print("Correlation Matrix for Postmenopausal Group:")
print(cor_postmenopausal)
print("P-value Matrix for Postmenopausal Group:")
print(pvalues_postmenopausal)

# To visualize the correlation matrices with heatmaps 
#----------------------------------------------------#

# Install and load corrplot package if you haven't already
# install.packages("corrplot")
library(corrplot)

# Create heatmap for the correlation matrix for all cases
par(mfrow = c(1,1))  
png("correlation_all_cases.png", width = 800, height = 800, res = 150)

corrplot(cor_all, 
         method = "color",    # Type of plot (other options: "square", "number", etc.)
         type = "upper",       # Display only the upper triangle
         tl.col = "black",     # Color of the text labels
         tl.cex = 1.1,
         p.mat = pvalues_all,  # p-values matrix to show significance levels
         addCoef.col = "white",
         sig.level = 0.05,     # Show only significant correlations (p-value <= 0.05)
         insig = "blank") # Option to hide insignificant correlations
         # title = "Correlation Heatmap for All Cases")

dev.off()

# For premenopausal group
png("correlation_all_cases.png", width = 800, height = 800, res = 150)
corrplot(cor_premenopausal, 
         method = "color",
         type = "upper",       # Display only the upper triangle
         tl.col = "black",     # Color of the text labels
         tl.cex = 1.1,
         p.mat = pvalues_premenopausal,  # p-values matrix to show significance levels
         addCoef.col = "white",
         sig.level = 0.05,     # Show only significant correlations (p-value <= 0.05)
         insig = "blank")
dev.off()

# For postmenopausal group
png("correlation_all_cases.png", width = 800, height = 800, res = 150)
corrplot(cor_postmenopausal, 
         method = "color",
         type = "upper",       # Display only the upper triangle
         tl.col = "black",     # Color of the text labels
         tl.cex = 1.1,
         p.mat = pvalues_postmenopausal,  # p-values matrix to show significance levels
         addCoef.col = "white",
         sig.level = 0.05,     # Show only significant correlations (p-value <= 0.05)
         insig = "blank")
dev.off()
