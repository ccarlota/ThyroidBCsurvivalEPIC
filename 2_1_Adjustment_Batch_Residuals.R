#-----------------------------------------#
# Adjustment of hormones for batch effect: 
#-----------------------------------------#

library(ggplot2)

set$set$Batch_Ft_f <- as.factor(set$Batch_Ft)

table(set$Batch_Ft, set$COUNTRY)
chisq.test(table(set$Batch_T, set$COUNTRY))

# distribution of Batch_T by COUNTRY:
ggplot(set, aes(x = Batch_T, fill = COUNTRY)) +
  geom_bar(position = "dodge") +
  labs(x = "Batch_T", y = "Count", title = "Distribution of Batch_T by COUNTRY") +
  theme_minimal()

library(nlme)

set <- readRDS("/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Data/processed_data/001_Setup_script.rds")

# Compute residuals on log-transformed data:
#-------------------------------------------#
# Concentrations of biomarkers were log-transformed to approximate 
# normal distribution and we used residuals of log-transformed variables 
# regressed on analytical batch 

# Log-transformation of variables ----
data <- set[,c("fT3_pmol_L","fT4_pmol_L","fT3fT4r", "TSH")] 

for (i in names(data)) {
  horm <- data[, i]
  set2 <- data.frame(horm = log(horm))
  colnames(set2) <- paste0("log_",i)
  set <- cbind(set,set2)
}

#  Biomarker data -->
# Log_biomarker = COUNTRY + Batch -->
# Batch -> random effect --> Residuals with random effect-->

Log_horm_Ft <- c("log_fT3_pmol_L","log_fT4_pmol_L","log_fT3fT4r") # with Batch_Ft
# "log_TSH" with Batch_T

Residuals <- list()

for (i in Log_horm_Ft) {
  
  Log_horm_R <- paste0(i, "_R")
  Residuals[[Log_horm_R]] <- as.numeric(residuals(lme(as.formula(paste(i, "~ COUNTRY")), 
                                          random = ~ 1 | Batch_Ft, 
                                            data = set, 
                                              na.action = "na.exclude"), 
                                                  type = "pearson"))
  set[[Log_horm_R]] <- Residuals[[Log_horm_R]]
}
summary(set$log_fT3_pmol_L_R)

set$log_TSH_R = as.numeric(residuals(lme(log_TSH ~ COUNTRY, 
                              random = ~ 1 | Batch_T, data = set, 
                                na.action = "na.exclude"), 
                                  type = "pearson"))
summary(set$log_TSH_R)

# Compare results now using lmer and not adjusting for COUNTRY:
library(lme4)

set_clean <- na.omit(set)
model <- lmer(log_TSH ~ 1 + (1 | Batch_T), data = set, 
                                          na.action = na.omit,
                                          control=lmerControl(optimizer="bobyqa", 
                                          optCtrl=list(maxfun=2e5)))
                                       
log_TSH_R_2 <- residuals(model)
summary(set$log_TSH_R)

set$log_TSH_R_2[!is.na(set$log_TSH)] <- log_TSH_R_2

summary(set$log_TSH_R) # lme
summary(set$log_TSH_R_2) # lmer

