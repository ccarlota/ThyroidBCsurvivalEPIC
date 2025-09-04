# Distribution of thyroid hormones among AbTPO+ and AbTPO-
#---------------------------------------------------------#
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
set <- readRDS("/data/Epic/subprojects/Breast_Cancer/work/Thyroid_Horm_BC_survival/Data/processed_data/001_Setup_script.rds")
colnames(set)
set$COUNTRY <- as.factor(set$COUNTRY)

# Descriptive of thyroid hormones across all BC cases ----
#--------------------------------------------------------#

table_thy <- set %>% 
  dplyr::select(fT3_pmol_L, fT4_pmol_L, TSH, fT3fT4r) %>%
  summarise(across(everything(), list(
    n = \(x) sum(!is.na(x)),  # Count non-NA values
    na = \(x) sum(is.na(x)),  # Count NA values
    min = \(x) min(x, na.rm = TRUE),
    P5 = \(x) quantile(x, 0.05, na.rm = TRUE),  
    P10 = \(x) quantile(x, 0.10, na.rm = TRUE),
    P25 = \(x) quantile(x, 0.25, na.rm = TRUE), 
    mean = \(x) mean(x, na.rm = TRUE), 
    geomean = \(x) exp(mean(log(x), na.rm = TRUE)),
    median = \(x) median(x, na.rm = TRUE),  
    sd = \(x) sd(x, na.rm = TRUE),
    P75 = \(x) quantile(x, 0.75, na.rm = TRUE),
    P90 = \(x) quantile(x, 0.90, na.rm = TRUE),
    P95 = \(x) quantile(x, 0.95, na.rm = TRUE),
    max = \(x) max(x, na.rm = TRUE)), 
    .names = "{.col} ({.fn})")) %>%  
  
  # Pivot longer to get a tidy format
  pivot_longer(everything(),
               names_to = c("variable", "stat"),
               names_pattern = "^(.*) \\((.*)\\)$") %>%
  
  pivot_wider(names_from = stat, values_from = value) %>%
  
  mutate(across(where(is.numeric), ~ round(., 2)))

table_thy

# save table of summary statistics for thyroid horm, their ratio and TSH:
write.table(table_thy, file = paste0(output.dir,"001.1_table_horm_concentration.txt"), 
            col.names = T, row.names = F, sep = "\t", quote = F)

# Descriptive of Thyroid variables across antiTPO (+) and antiTPO (-) ----
#------------------------------------------------------------------------#
# Continuous variables
thy_vars <- c("fT3_pmol_L", "fT4_pmol_L", "TSH", "fT3fT4r")

table_thy_strat <- set %>%
  dplyr::filter(Anti_TPO_p %in% c("Positive", "Negative")) %>%
  dplyr::group_by(Anti_TPO_p) %>%
  summarise(across(all_of(thy_vars), list(
    n = \(x) sum(!is.na(x)),
    na = \(x) sum(is.na(x)),
    min = \(x) min(x, na.rm = TRUE),
    P5 = \(x) quantile(x, 0.05, na.rm = TRUE),
    P10 = \(x) quantile(x, 0.10, na.rm = TRUE),
    P25 = \(x) quantile(x, 0.25, na.rm = TRUE),
    mean = \(x) mean(x, na.rm = TRUE),
    geomean = \(x) exp(mean(log(x), na.rm = TRUE)),
    median = \(x) median(x, na.rm = TRUE),
    sd = \(x) sd(x, na.rm = TRUE),
    P75 = \(x) quantile(x, 0.75, na.rm = TRUE),
    P90 = \(x) quantile(x, 0.90, na.rm = TRUE),
    P95 = \(x) quantile(x, 0.95, na.rm = TRUE),
    max = \(x) max(x, na.rm = TRUE)
  ), .names = "{.col} ({.fn})")) %>%
  dplyr::ungroup()


# Optional: pivot to long format for tidy table
# Desired order
# Desired order of variables
var_order <- c("TSH", "fT3_pmol_L", "fT4_pmol_L", "fT3fT4r")

table_thy_long <- table_thy_long %>%
  mutate(variable = factor(variable, levels = var_order)) %>%
  arrange(Anti_TPO_p, variable)  # first by Anti_TPO_p, then by variable order

table_thy_long

write.table(table_thy_long, file = paste0(output.dir,"001.1_table_horm_across_TPO.txt"), 
            col.names = T, row.names = F, sep = "\t", quote = F)



