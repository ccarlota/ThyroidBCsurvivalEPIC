#------------------------------------------------------#
# Formatting variables and descriptes
# 10/02/2025
#------------------------------------------------------#

# Libraries 
library(nlme)
library(haven)
library(foreign)
library(foreach)
library(tidyverse)
library(dplyr)
library(tidyr)
library(openxlsx)
# library(xlsx)
library(ggplot2)

library(mice)
library(VIM)

# Directories
data.dir <- "/data/Epic/subprojects/Breast_Cancer/files/Maestra/"
output.dir <- "/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Analysis/Tables/"

# Data
dat <- read_sas(paste0(data.dir, "surv_brea.sas7bdat")) # 1551
table(dat$Cncr_Caco_Brea, useNA = "always") # 1524 BC cases
set <- dat %>% filter(Cncr_Caco_Brea == 1) # 1524 BC cases
table(set$Cncr_Caco_Brea, useNA = "always") 

# Missing missing ER, PR, HER2 checked and there're no missings.
table(set$Er_Adj, useNA = "always") 
table(set$Pr_Adj, useNA = "always")
table(set$Her2_Adj, useNA = "always") 

# Check unique Idepic:
anyDuplicated(set$Idepic) > 0  # ok, no duplicates

#-----------------------#
# Checking missings ----
#-----------------------#

# Select columns in set that I will use and can contain missings:
set2 <- set %>%
  select(Idepic, Age_Blood, Death_Status, Death_Status_Csm, 
         L_School, Bmi_C, Pa_Index, Alc_Re, Fasting_C, Smoke_Stat, 
         Phrt_Bld, Er_Adj, Pr_Adj, Her2_Adj, COUNTRY, Menop_Bld)

# Checking missing data patterns/misstable patterns*/
md.pattern(set2)

# Calculates percentage of complete cases (i.e., rows with no missing values) 
complete_cases <- sum(complete.cases(set2)) / nrow(set2) * 100
cat("Percentage of complete cases before imputation:", round(complete_cases, 2), "%\n")
# 97.64 % (2.4% missings)

missing_counts <- colSums(is.na(set2))
missing_vars <- names(missing_counts[missing_counts > 0])  # Columns with missing values
missing_vars # "L_School"  "Alc_Re"    "Fasting_C"

# Replace missing values with the median of complete cases for each variable
vars_to_impute <- c("L_School","Alc_Re","Fasting_C")

for (var in vars_to_impute) {
  set2[[var]][is.na(set2[[var]])] <- median(set2[[var]], na.rm = TRUE)
}

md.pattern(set2) # no missings now.

# Replace these complete variables in set2 with imputed values into the original set:
set[vars_to_impute] <- set2[vars_to_impute]

#----------------------------#
# Duration time ----
#----------------------------#
# (1) Date diagnosis to Date end of study ----

summary(is.na(set$D_Dgbrea)) # Date of diagnosis of the tumour
summary(is.na(set$D_Exit_Vs)) # (for all-cause mortality) Date end of follow-up for vital status, i.e. when the last vital status is known
summary(is.na(set$D_Exit_Csm)) # (for specific-cause mortalitY)

# To calculate this time period:
set$date_dx <- as.Date(set$D_Dgbrea, format= "%Y-%m-%d")
set$year_diagnosis <- as.numeric(format(set$date_dx, "%Y"))
set$date_exit_vit_status <- as.Date(set$D_Exit_Vs, format = "%Y-%m-%d") 

set$time_surv <- set$date_exit_vit_status - set$date_dx 
set$time_surv <- as.numeric(set$time_surv)
summary(is.na(set$time_surv)) 

summary(set$time_surv)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -4305    2021    2543    2657    3410    5999 

# Remove variables with inconsistent period (e.g. date at diagnosis after date of death. therefore negative 'time_surv')
table(set$time_surv <= 0) # 3 cases have time_surv < 0 and 1 cases = 0.
test <- set %>% filter(time_surv <= 0)
check <- test[,c("Idepic","D_Dgbrea","date_exit_vit_status")]
check

set <- subset(set, !(Idepic %in% check$Idepic))
dim(set) # 1521

table(set$time_surv < 0) # 0
summary(set$time_surv)

as.data.frame(set %>% filter(time_surv == 0) %>% select(Idepic, D_Dgbrea, date_exit_vit_status)) # 1 case

summary(set$time_surv/365.25) # IN YEARS
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.03559  5.54689  6.97057  7.29723  9.34702 16.42437 

#---------------------------#
# (2) Age at diagnosis ----
#---------------------------#
# 1) difference between the date at blood collection and the date at diagnosis
# 2) this period (1) in years will be summed to the age at blood collection to get the age at diagnosis.

# To calculate time period (in years) between date entry and date of diagnosis:
set$date_dx <- as.Date(format(set$D_Dgbrea, format= "%Y-%m-%d")) # Date at diagnosis of breast cancer
set$date_recr <- as.Date(format(set$D_Recrui, format = "%Y-%m-%d")) # Data at recruitment
set$date_blo <- as.Date(format(set$D_Bld_Coll, format = "%Y-%m-%d")) # Date at blood collection

summary(is.na(set$date_dx)) # 0 NAs
summary(is.na(set$date_recr)) # 0 NAs
summary(is.na(set$date_blo)) # 0 NAs

set$blo_dx <- (set$date_dx - set$date_blo)/365.25 # years between blood collection and diagnosis
head(set$blo_dx)
class(set$blo_dx)
set$blo_dx <- as.numeric(set$blo_dx)
summary(set$blo_dx) # mean time between blood collection to diagnosis: 8.3 years (median 8.6y) SD 2.7y

# Age at diagnosis will be age at blood collection + duration time (years) from blood to diagnosis
set$age_dx <- set$Age_Blood + set$blo_dx 
summary(set$age_dx) # median age at diagnosis: 61 years
summary(is.na(set$age_dx)) # 0 NAs

# Age at diagnosis (as a proxy of menopausal status at diagnosis) ≤ 55 years and > 55 years:
set$age_dx_2 <- if_else(set$age_dx <55, "<55y", ">=55y",NA)
table(set$age_dx_2, useNA="always") 

# Age in 5 years categories: 
set$age_dx_5y <- cut(set$age_dx,
                     breaks = c(35, seq(40, max(set$age_dx, na.rm = TRUE), by = 5), Inf),  # Add Inf for ages over 80
                     labels = c(paste(seq(35, max(set$age_dx, na.rm = TRUE) - 5, by = 5), "-", 
                                      seq(39, max(set$age_dx, na.rm = TRUE), by = 5)), 
                                "> 80"),  # Label for ages over 80
                     right = FALSE,  # Left-closed intervals
                     include.lowest = TRUE)  # Include the lowest value (35) in the first interval

table(set$age_dx_5y, useNA = "always")
set$age_dx_5y <- as.factor(set$age_dx_5y)

#--------------------#
# Age at event ----
#--------------------#
# Age as time-scale in cox models:
# Create a new variable 'age_at_event' which represents the age at the time of event

set$age_at_event <- as.numeric(difftime(set$date_exit_vit_status, set$date_dx, units = "days") / 365.25) + set$age_dx
summary(set$age_at_event)

summary(set$age_dx)

set$age_time_surv <- set$age_at_event - set$age_dx
summary(set$age_time_surv)

#------------------------------------------------------------------------#
# Time from blood collection to diagnosis and age at blood collection ----
#------------------------------------------------------------------------#
# 1) difference between the date at recruitment and date at blood collection
# 2) this period (1) in years is summed to age_recruitment to get the age_blood collection.

# variable date at blood collection: 'D_Bld_Coll'
# variable date of recruitment: 'D_Recrui'

# To calculate time period (in years) between date blood collection and recruitment:
set$date_blo <- as.Date(format(set$D_Bld_Coll, format= "%Y-%m-%d"))
set$date_recr <- as.Date(format(set$D_Recrui, format = "%Y-%m-%d"))
set$date_dx <- as.Date(format(set$D_Dgbrea, format= "%Y-%m-%d"))

summary(is.na(set$date_blo)) # 0 NAs
summary(is.na(set$date_recr)) # 0 NAs
summary(is.na(set$date_dx)) # 0 NAs

# Time in years between blood collection & recruitment:
set$blo_recr <- (set$date_blo - set$date_recr)/365.25
head(set$blo_recr)
class(set$blo_recr)
set$blo_recr <- as.numeric(set$blo_recr)
summary(set$blo_recr, useNA=="always") # median between blood collection and recruitment is 0 years

summary(is.na(set$blo_recr)) # ok 0 NAs
set %>% group_by(COUNTRY) %>% summarise(min = min(blo_recr),
                                        median = median(blo_recr),
                                        mean= mean(blo_recr)) # France is the country where the period is longer (median 3.9 y from blood to recruit)
set %>%
  filter(blo_recr > 3) %>%
  summarise(n_casos = n())

# Check breast malignant cancer:
table(set$Cncr_Mal_Brea, useNA = "always") # ok, all cases are malignant tumours

# Check deaths:
table(set$Death_Status, useNA = "always")
#   0    1   <NA> 
# 1297  224    0 

# Check specific deaths:
table(set$Death_Status_Csm, useNA = "always")
#   0    1   <NA> 
# 1304  217    0 

#--------------------------------------#
# Death and specific death status ----
#-------------------------------------#
# For mortality models: 0 Alive - 1 Dead 
table(set$Death_Status,useNA = "always")
set$Death_Status <- as.factor(set$Death_Status) 
table(set$Death_Status)
set$Death_Status <- factor(set$Death_Status, labels = c("Alive","Dead"))
table(set$Death_Status, useNA="always") 
set$event <- set$Death_Status
table(set$event,useNA = "always")

# For specific cause mortality models: 
table(set$Death_Status_Csm,useNA = "always")
set$Death_Status_Csm <- as.factor(set$Death_Status_Csm) 
table(set$Death_Status_Csm)
set$Death_Status_Csm <- factor(set$Death_Status_Csm, labels = c("Alive","Dead"))
table(set$Death_Status_Csm, useNA="always") 
set$event_Csm <- set$Death_Status_Csm
table(set$event_Csm,useNA = "always")

# Descriptive of deaths by country ----
#-------------------------------------------------#
country_name <- c("7"="Germany",
                  "2"="Italy",
                  "9"="Denmark",
                  "3"="Spain",
                  "5"="The Netherlands",
                  "4"="United Kingdom",
                  "1"="France")

table_death_country <- set %>%
  filter(Death_Status == 1) %>% 
  group_by(COUNTRY) %>% 
  summarise(
    death_count = n(),
    total_count = nrow(set %>% filter(Death_Status == 1)),
    death_percent = round(death_count/total_count*100,1),
  ) %>% arrange(desc(death_count)) %>%
  mutate(
    cumulative_freq = cumsum(death_count),
    cumulative_percent = cumsum(round(death_percent,1)),
    country_names = country_name[as.character(COUNTRY)]
  ) %>%
  select(COUNTRY, country_names, 
         death_count,cumulative_freq,
         death_percent,cumulative_percent)

table_death_country

# save table ofsummary statistics for thyroid horm,their ratio and TSH:
write.table(table_death_country, file = paste0(output.dir,"001_table_death_country.txt"), 
            col.names = T, row.names = F, sep = "\t", quote = F)

# Descriptive of specific deaths by country ----
#----------------------------------------------------------#
country_name <- c("7"="Germany",
                  "2"="Italy",
                  "9"="Denmark",
                  "3"="Spain",
                  "5"="The Netherlands",
                  "4"="United Kingdom",
                  "1"="France")

table_death_country <- set %>%
  filter(Death_Status_Csm == 1) %>% 
  group_by(COUNTRY) %>% 
  summarise(
    death_count = n(),
    total_count = nrow(set %>% filter(Death_Status_Csm == 1)),
    death_percent = round(death_count/total_count*100,1),
  ) %>% arrange(desc(death_count)) %>%
  mutate(
    cumulative_freq = cumsum(death_count),
    cumulative_percent = cumsum(round(death_percent,1)),
    country_names = country_name[as.character(COUNTRY)]
  ) %>%
  select(COUNTRY, country_names, 
         death_count,cumulative_freq,
         death_percent,cumulative_percent)

table_death_country

# save table ofsummary statistics for thyroid horm,their ratio and TSH:
write.table(table_death_country, file = paste0(output.dir,"001_table_s.death_country.txt"), 
            col.names = T, row.names = F, sep = "\t", quote = F)

# Follow-up time since diagnosis ----
summary(set$Length_Vs)
5986/365.25 # mean 16.4 days

# Thyroid hormones units ----
# fT3, fT4, TSH, fT3fT4r
#----------------------------#
# fT3 pg/mL 
# fT4 ng/dL 
# Mathilde had them in pmol/L. 
# To convert fT3 values from pg/mL to pmol/L: fT3(pmol/L)=fT3(pg/mL)× 1.536
# To convert fT4 values from ng/dL to pmol/L: fT4(pmol/L)=fT4(ng/dL)× 12.87

set$fT3_pmol_L <- set$fT3 * 1.536
set$fT4_pmol_L <-set$fT4 * 12.87

set <- set %>%
  mutate(fT3fT4r = fT3_pmol_L / fT4_pmol_L)


# Table of horm concentrations by country ----
#---------------------------------------------#
table_horm_country <- set %>%
  group_by(COUNTRY) %>% 
  summarise(across(c(fT3_pmol_L, fT4_pmol_L, TSH, fT3fT4r, Anti_TPO), list(
    n = \(x) sum(!is.na(x)),  
    na = \(x) sum(is.na(x)),  
    min = \(x) ifelse(all(is.na(x)), NA, min(x, na.rm = TRUE)),  
    mean = \(x) ifelse(all(is.na(x)), NA, mean(x, na.rm = TRUE)),  
    geomean = \(x) exp(mean(log(ifelse(x > 0, x, NA)), na.rm = TRUE)),  
    median = \(x) ifelse(all(is.na(x)), NA, median(x, na.rm = TRUE)),  
    sd = \(x) ifelse(all(is.na(x)), NA, sd(x, na.rm = TRUE)),  
    max = \(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE))),  
    .names = "{.col}_{.fn}")) %>%
  
  # Ungroup to avoid selection issues
  ungroup() %>%
  
  # Check if COUNTRY exists
  relocate(COUNTRY, .before = everything()) %>%
  
  pivot_longer(-COUNTRY, names_to = c("Hormone", "Stat"), names_pattern = "^(.*)_(.*)$") %>%
  
  # Pivot wider to create separate columns for each country
  pivot_wider(names_from = COUNTRY, values_from = value) %>%
  
  # Arrange hormones in order
  arrange(Hormone) %>%
  
  mutate(across(where(is.numeric), ~ round(., 2)))

table_horm_country

# save 
write.table(table_horm_country, file = paste0(output.dir,"001_table_horm_country.txt"), 
            col.names = T, row.names = F, sep = "\t", quote = F)

# TPO-Ab status, % Positive ----
#---------------------------#
# His M, et al. "For TPO-Ab, values were dichotomized according to cut-off recommended by the vendor (40 IU/mL), 
# and sample showing TPO-Ab concentration above this cut-off were considered TPO-Ab-positive"

summary(set$Anti_TPO)
set$Anti_TPO_p <- if_else(set$Anti_TPO <= 40, "Negative",
                          if_else(set$Anti_TPO >40, "Positive",NA))
table(set$Anti_TPO_p, useNA = "always")
# set$Anti_TPO_p <- ifelse(is.na(set$Anti_TPO_p),"Missing", set$Anti_TPO_p)

counts <- table(set$Anti_TPO_p, useNA = "always")
percentages <- round(prop.table(counts),3) * 100
percentages
# Negative  Positive   Unknown
#   87.9     11.6      0.5 
set$Anti_TPO_p <- as.factor(set$Anti_TPO_p)

# Different cut-offs to check if any differences
#- - - - - - - - - - - - - - - - - - - - - - - - -#
# 30 IU/mL:
summary(set$Anti_TPO)
set$Anti_TPO_30 <- if_else(set$Anti_TPO <= 30, "Negative",
                          if_else(set$Anti_TPO >30, "Positive",NA))
table(set$Anti_TPO_30, useNA = "always")
# set$Anti_TPO_p <- ifelse(is.na(set$Anti_TPO_p),"Missing", set$Anti_TPO_p)

counts <- table(set$Anti_TPO_30, useNA = "always")
percentages <- round(prop.table(counts),3) * 100
percentages
# Negative  Positive   Unknown
#   86.4    13.1      0.5 
set$Anti_TPO_30 <- as.factor(set$Anti_TPO_30)

# 60 IU/mL:
summary(set$Anti_TPO)
set$Anti_TPO_60 <- if_else(set$Anti_TPO <= 60, "Negative",
                           if_else(set$Anti_TPO >60, "Positive",NA))
table(set$Anti_TPO_60, useNA = "always")
# set$Anti_TPO_p <- ifelse(is.na(set$Anti_TPO_p),"Missing", set$Anti_TPO_p)

counts <- table(set$Anti_TPO_60, useNA = "always")
percentages <- round(prop.table(counts),3) * 100
percentages
# Negative  Positive   Unknown
#   87.9     11.6      0.5 
set$Anti_TPO_60 <- as.factor(set$Anti_TPO_60)

#---------------------------------------------------------#
# Descriptive of sociodemographic-reproductive variables:
#---------------------------------------------------------#

# Descriptive Continuous variable ----
#-------------------------------------#
table_desc <- set %>% 
  select(Age_Blood, Age_Recr, age_dx, blo_dx, Height_C, Bmi_C) %>% 
  summarise(across(everything(), list(mean = ~round(mean(.), 2),
                                      median = ~round(median(.), 2),
                                      sd = ~round(sd(.), 2)))) %>%
  pivot_longer(cols = everything(), names_to = c("Variable", "Statistic"), 
               names_pattern = "(.+)_(mean|median|sd)") %>%
  pivot_wider(names_from = Statistic, values_from = value)

print(table_desc)

# save table:
write.xlsx(table_desc, file = paste0(output.dir,"table_desc_continuous.csv"))

# Descriptive Categorical variables:
#-----------------------------------#
# Labeling name categories:

# Estrogen receptor status 
set$Er_Adj_c <- as.factor(set$Er_Adj)
levels(set$Er_Adj_c) <- c("Negative","Positive")
table(set$Er_Adj_c, useNA = "always")

# Progesterone receptor status 
set$Pr_Adj_c <- as.factor(set$Pr_Adj)
levels(set$Pr_Adj_c) <- c("Negative","Positive")
table(set$Pr_Adj_c, useNA = "always")

# Human epidermal receptor 2 status 
set$Her2_Adj_c <- as.factor(set$Her2_Adj)
levels(set$Her2_Adj_c) <- c("Negative","Positive")
table(set$Her2_Adj_c, useNA = "always")

# Fasting ----
set$fasting <- if_else(set$Fasting_C==0,"No",
                       if_else(set$Fasting_C==1,"In between",
                               if_else(set$Fasting_C==2,"Yes",NA)))
set$fasting <- factor(set$fasting, levels = c("No","In between","Yes"))
table(set$fasting, useNA="always")

# Batch Residual method ----

# 1) Concentrations of hormones were log-transformed to approximate 
# normal distribution and we used residuals of log-transformed variables 
# regressed on analytical batch 

# Log-transformation of variables ----
data <- set[,c("fT3_pmol_L","fT4_pmol_L","fT3fT4r", "TSH", "Anti_TPO")] 

for (i in names(data)) {
  horm <- data[, i]
  set2 <- data.frame(horm = log(horm))
  colnames(set2) <- paste0("log_",i)
  set <- cbind(set,set2)
}

# Horm. data --> Log_horm = Batch + Intercept --> random effect --> Residuals with random effect
Log_horm_Ft <- c("log_fT3_pmol_L","log_fT4_pmol_L","log_fT3fT4r") # with Batch_Ft 

# fT3,Ft4,ft3ft4r:Batch residuals
Residuals <- list()
for (i in Log_horm_Ft) {
  
  Log_horm_R <- paste0(i, "_R")
  Residuals[[Log_horm_R]] <- as.numeric(residuals(lme(as.formula(paste(i, "~ 1")), 
                                                      random = ~ 1 | Batch_Ft, 
                                                      data = set, 
                                                      na.action = "na.exclude"), 
                                                  type = "pearson"))
  set[[Log_horm_R]] <- Residuals[[Log_horm_R]]
}
summary(set$log_fT3_pmol_L_R)

# TSH:Batch_T residuals
set$log_TSH_R = as.numeric(residuals(lme(log_TSH ~ 1, 
                                         random = ~ 1 | Batch_T, data = set, 
                                         na.action = "na.exclude"), 
                                     type = "pearson"))
summary(set$log_TSH_R)

# Anti_TPO:Batch_T
set$log_Anti_TPO_R = as.numeric(residuals(lme(log_Anti_TPO ~ 1, 
                                         random = ~ 1 | Batch_T, data = set, 
                                         na.action = "na.exclude"), 
                                     type = "pearson"))
summary(set$log_Anti_TPO_R)

#----------------------#
# Country ----
#----------------------#
set$COUNTRY_c <- if_else(set$COUNTRY == 1, "France", 
                       if_else(set$COUNTRY == 2, "Italy",
                               if_else(set$COUNTRY == 3, "Spain",
                                       if_else(set$COUNTRY == 4, "United Kingdom",
                                               if_else(set$COUNTRY == 5, "The Netherlands",
                                                       if_else(set$COUNTRY == 7, "Germany",
                                                               if_else(set$COUNTRY == 9, "Denmark",NA)))))))
table(set$COUNTRY_c, useNA="always")
set$COUNTRY_c <- factor(set$COUNTRY_c, levels = c("France","Italy","Spain","United Kingdom", "The Netherlands", "Germany", "Denmark"))

# Menopause at blood collection:
table(set$Menop_Bld, useNA="always") # 0 premenopausal; 1 postmenopausal; 2 perimenopausal; 3 surgical postmenopausal
set$menop <- if_else(set$Menop_Bld == 1 | set$Menop_Bld == 3, 1, set$Menop_Bld) # post + surgical in 1 category
table(set$menop, useNA="always") # 0 premenopausal; 1 postmenopausal; 2 perimenopausal

set$menop <- if_else(set$menop == 0, "Premenopausal",
                            if_else(set$menop == 1, "Postmenopausal",
                                    "Perimenopausal"))
set$menop <- factor(set$menop, levels = c("Premenopausal","Perimenopausal","Postmenopausal"))
table(set$menop, useNA = "always")

# Menopause at diagnosis (created) ----
# Women who were postmenopausal at blood collection or were 55 years or older at diagnosis (regardless of their menopausal status at blood collection) 
# will be classified as postmenopausal at diagnosis. All others will be classified as premenopausal.
set$menop_status_dx <- if_else(set$menop == "Postmenopausal", "Postmenopausal",
                               if_else((set$menop == "Premenopausal" | set$menop == "Perimenopausal") &
                                         set$age_dx >= 55, "Postmenopausal","Premenopausal"))
table(set$menop_status_dx, useNA="always")
set$menop_status_dx <- factor(set$menop_status_dx, levels = c("Premenopausal","Postmenopausal"))

length(set$age_dx[set$age_dx >=55 & set$menop == "Premenopausal"]) # 123 were premenopausal at recruitment with =>55 years at dx, so we estimate that at diagnosis they will be postmenopausal
length(set$age_dx[set$age_dx >=55 & set$menop == "Perimenopausal"]) # 210 were perimenopausal at recruitment with =>55 years, so we estimate that at diagnosis they will be postmenopausal
123+210 # 333 were pre or peri at recruitment with >=55y so at dx we count them as postmenopausal.

length(set$age_dx[set$menop == "Postmenopausal"]) # 818
819 + 333 # 1154 are categorized as postmenopausal at diagnosis. 

# Use of hormones at blood collection ----
# Carine has added these variables:

# Horm_Bld     Num       8    YES_NO.      Use of hormones for menopause at blood collection
# Pill_Bld     Num       8    YES_NO.      Use of pill at blood collection
# Hormon       Num       8    YES_NO.      Have you taken hormone
# Phrt_Bld     Num       8    YES_NO.      Use of pill/hrt-ert at blood collection

table(set$Phrt_Bld, useNA = "always") # we will use this in our survival models
set$Phrt_Bld_c <- as.factor(set$Phrt_Bld)
levels(set$Phrt_Bld_c) <- c("No", "Yes")
table(set$Phrt_Bld_c, useNA = "always")

# BMI 'bmi_oms' ----
summary(is.na(set$Bmi_C)) # 0 NAs

set$bmi_oms <- if_else(set$Bmi_C<18.5, "UnderWeight",
                       if_else(set$Bmi_C >=18.5 & set$Bmi_C <25, "NormalWeight",
                               if_else(set$Bmi_C >=25 & set$Bmi_C <30, "Overweight",
                                       if_else(set$Bmi_C >=30, "Obesity", ""))))
table(set$bmi_oms, useNA = "always")  
set$bmi_oms <- factor(set$bmi_oms, levels = c("NormalWeight","Overweight","Obesity","UnderWeight"))
table(set$bmi_oms, useNA = "always")  

# Education level ----
table(set$L_School, useNA = "always")
set$l_school2 <- if_else(set$L_School == 0 | set$L_School == 1, "None/Primary", 
                         if_else(set$L_School == 2, "Tech/Prof",
                                 if_else(set$L_School == 3, "Secondary",
                                         if_else(set$L_School == 4, "Longer educat", 
                                                 if_else(set$L_School==5, "Unknown","")))))   

table(set$l_school2, useNA = "always")
set$l_school2 <- if_else(is.na(set$l_school2), "Unknown", set$l_school2)
table(set$l_school2, useNA = "always")
set$l_school2 <- factor(set$l_school2, levels = c("None/Primary", "Tech/Prof", "Secondary", "Longer educat", "Unknown"))
table(set$l_school2, useNA = "always")

# Physical activity ----
table(set$Pa_Index, useNA = "always")
class(set$Pa_Index)
set$Pa_Index_c <- factor(set$Pa_Index)
levels(set$Pa_Index_c) <- c("Inactive","Mod.Inactive","Mod.Active","Active","Unknown")
table(set$Pa_Index_c, useNA = "always")

# Smoking status ----
table(set$Smoke_Intensity, useNA = "always")
set$Smoke_Intensity_cc <- if_else(set$Smoke_Intensity==1, "Never",
                                  if_else(set$Smoke_Intensity==2, "Current, 1-15 cig/d",
                                          if_else(set$Smoke_Intensity==3, "Current, 16-25 cig/d",
                                                  if_else(set$Smoke_Intensity==4, "Current, 26+ cig/d",
                                                          if_else(set$Smoke_Intensity==5, "Former, quit <= 10 years",
                                                                  if_else(set$Smoke_Intensity==6, "Former, quit 11-20 years",
                                                                          if_else(set$Smoke_Intensity==7, "Former, quit 20+ years",
                                                                                  if_else(set$Smoke_Intensity==8| set$Smoke_Intensity==9, "Miscellaneous","Unknown"))))))))
table(set$Smoke_Intensity_cc, useNA = "always")

set$Smoke_Intensity_cc <- factor(set$Smoke_Intensity_cc, levels = c("Never","Current, 1-15 cig/d","Current, 16-25 cig/d",
                                                                    "Current, 26+ cig/d","Former, quit <= 10 years",
                                                                    "Former, quit 11-20 years","Former, quit 20+ years",
                                                                    "Miscellaneous","Unknown"))
table(set$Smoke_Intensity_cc, useNA = "always")

set$Smoke_Stat_cc <- as.factor(set$Smoke_Stat); levels(set$Smoke_Stat_cc) <- c("Never","Former","Current","Unknown")
table(set$Smoke_Stat_cc, useNA="always")

# Alcohol ----
summary(is.na(set$Alc_Re))
class(set$Alc_Re)
set$Alc_Re_c <- if_else(set$Alc_Re == 0, "Non drinker",
                        if_else(set$Alc_Re >0 & set$Alc_Re <3, ">0-3",
                                if_else(set$Alc_Re >= 3 & set$Alc_Re <12, ">3-12",
                                        if_else(set$Alc_Re >= 12 & set$Alc_Re <24, ">12-24",
                                                if_else(set$Alc_Re >=24, ">24","")))))
set$Alc_Re_c <- ifelse(is.na(set$Alc_Re_c),"Unknown", set$Alc_Re_c)
table(set$Alc_Re_c, useNA = "always")

set$Alc_Re_c <- factor(set$Alc_Re_c, levels = c("Non drinker",">0-3",">3-12",">12-24",">24","Unknown"))

#-------------------------------#
# Grade of tumour ----
#-------------------------------#
table(set$Gradbrea, useNA = "always") # 294 NAs

set$gradb <- if_else(set$Gradbrea == 1, "Well.diff", 
                    if_else(set$Gradbrea == 2, "Mod.diff",
                          if_else(set$Gradbrea == 3 | set$Gradbrea == 4, "Poorly diff/Undiff", 
                                  if_else(set$Gradbrea == 9, "Not determined", NA))))

set$gradb <- if_else(is.na(set$gradb), "Not determined", set$gradb)
table(set$gradb, useNA = "always")

set$gradb <- factor(set$gradb, levels = c("Well.diff","Mod.diff","Poorly diff/Undiff","Not determined"))
table(set$gradb, useNA = "always")

# Stage of tumour ----
#---------------------#
table(set$Stagbrea)
set$Stagbrea <- as.factor(set$Stagbrea)
set$stage2 <- if_else(set$Stagbrea == 2, "Localised",
                      if_else(set$Stagbrea == 3 | set$Stagbrea == 4 | set$Stagbrea == 5, "Metastatic",
                              if_else(set$Stagbrea == 9, "Unknown", set$Stagbrea)))

table(set$stage2, useNA = "always")
set$stage2 <- if_else(is.na(set$stage2), "Missing", set$stage2)
table(set$stage2, useNA = "always")

# When M0 is no metastasic:
set$stage3 <- case_when(
  set$stage2 == "Metastatic" ~ "Metastatic",  # Keep "Metastatic" from stage2
  grepl("M1", set$Tnmbrea) ~ "Metastatic",  # Ensure any M1 in Tnmbrea is Metastatic
  grepl("M0", set$Tnmbrea) ~ "Non-Metastatic", 
  is.na(set$Tnmbrea) | set$Tnmbrea == "" | !grepl("M[01]", set$Tnmbrea) ~ "Unknown"
)

table(set$stage3, useNA = "always")                      
set$stage3 <- factor(set$stage3, levels = c("Non-Metastatic","Metastatic","Unknown"))

# To see the TNM's within the "Unknown" category of stage3:
set %>% filter(stage3 == "Unknown") %>% count(Tnmbrea, TRUE)

# check that all metastatic are indeed M1 or stage2=="Metastatic" and non-metastatic are M0 or stage2 "non-metastatic"               
set %>% 
  filter(stage3 == "Metastatic") %>% 
  count(grepl("M1", Tnmbrea) | stage2 == "Metastatic") # ok

set %>% 
  filter(stage3 == "Non-Metastatic") %>% 
  count(grepl("M0", Tnmbrea) | stage2 == "Metastatic") # ok


# Table with categorical variables ----
#--------------------------------------#
table_categ <- set %>%
  select(menop_status_dx, Phrt_Bld_c, Er_Adj_c, Pr_Adj_c, Her2_Adj_c,
         fasting, Pa_Index_c, Smoke_Stat_cc, Alc_Re_c, l_school2) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "categ") %>%
  group_by(Variable, categ) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Variable) %>%
  mutate(Total = sum(N),  
         percent = (N / Total) * 100) %>%
  ungroup() %>%
  arrange(Variable, categ) 

print(table_categ, n = Inf)

# save table:
write.xlsx(table_categ, file = paste0(output.dir,"table_desc_catego.csv"))

## Save this dataset to start in 002 script:
saveRDS(set, "/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Data/processed_data/001_Setup_script.rds")


