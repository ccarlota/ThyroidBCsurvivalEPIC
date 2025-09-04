#--------------------#
# Checking Missings 
#--------------------#

# Libraries 
library(mice)
library(VIM)

library(tidyverse)
library(dplyr)
library(tidyr)
library(openxlsx)
library(survival)
library(DescTools)

# Directories
data.dir <- "/data/Epic/subprojects/Breast_Cancer/files/Maestra/"
output.dir <- "/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Analysis/Tables/"

# Charge processed data:
set <- readRDS("/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Data/processed_data/001_Setup_script.rds")

# Select columns in data
set2 <- set %>%
  select(Idepic, age_dx, age_at_event, Death_Status, Death_Status_Csm, 
         L_School, Bmi_Adj, Pa_Index, Alc_Re, Fasting_C, Smoke_Stat_cc, 
         Phrt_Bld, Er_Adj, Pr_Adj, Her2_Adj, COUNTRY, Menop_Bld)
dim(set2)

set2 <- set %>%
  select(
    id = Idepic, 
    agedx = age_dx, 
    ageEv = age_at_event, 
    DStat = Death_Status, 
    DStatCsm = Death_Status_Csm, 
    school = L_School, 
    bmi = Bmi_Adj, 
    pa = Pa_Index, 
    alc = Alc_Re, 
    fasting = Fasting_C, 
    smoke = Smoke_Stat_cc, 
    bld = Phrt_Bld, 
    #grad = Gradbrea, 
    er = Er_Adj, 
    pr = Pr_Adj, 
    her2 = Her2_Adj, 
    country = COUNTRY, 
    meno = Menopause
  )
# Checking missing data patterns/misstable patterns*/
md.pattern(set2[, c("agedx", "alc","bmi","ageEv")]) # few variables only
md.pattern(set2)

# Shrink the pattern
aggr(set2, numbers = TRUE, cex.axis = 0.7)

# Calculates percentage of complete cases (i.e., rows with no missing values) 
complete_cases <- sum(complete.cases(set2)) / nrow(set2) * 100
cat("Percentage of complete cases before imputation:", round(complete_cases, 2), "%\n")
# 97.37 % (2.6% missings)

missing_counts <- colSums(is.na(set2))
missing_vars <- names(missing_counts[missing_counts > 0])  # Select columns with missing values
print(missing_vars) # "school"  "bmi"     "alc"     "fasting"

# Replace missing values with the median of complete cases for each variable
set2 <- set2 %>%
  mutate(across(everything(), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))

# Check now grade of tumour and stage of cancer:
set2 <- set %>%
  select(Idepic,Gradbrea)
md.pattern(set2)
complete_cases <- sum(complete.cases(set2)) / nrow(set2) * 100
cat("Percentage", round(complete_cases, 2), "%\n")
# 80.7% (19.3% Missings)

set2 <- set %>%
  select(Idepic,Stagbrea)
md.pattern(set2)
complete_cases <- sum(complete.cases(set2)) / nrow(set2) * 100
cat("Percentage", round(complete_cases, 2), "%\n")
# 78.78 % (21.2% Missings)
