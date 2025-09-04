# Check whether use as a confounder: covariate use of exogenous horm at blood collection (Phrt_Bld_c)
# --------------------------------------------------------------------------------------------------#

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

# Hormones per 1-SD increase ----
horm <- colnames(set)[grep("_R$", colnames(set))]
horm # "log_TSH_R" "log_fT3_pmol_L_R" "log_fT4_pmol_L_R" "log_fT3fT4r_R" (R stands for Residuals from Batch effect!)

for (i in horm) {
  horm_sd <- paste0(i, "_sd")
  x <- data.frame(sd = set[, i]/sd(na.omit(set[, i])))
  colnames(x) <- horm_sd
  set <- cbind(set, x)
}
colnames(set)

#--------------------------------#
# All-cause mortality models ----
#--------------------------------#

# (1) exposure ↔ covariate association ----
exposure_cov_assoc <- lapply(exposures, function(x){
  message("Running exposure: ", x)
  res <- list(exposure = x)
  
  xv <- set[[x]]
  
  if (is.numeric(xv)) {
    fit <- lm(reformulate("Phrt_Bld_c", response = x), data=set)
    coefs <- coef(summary(fit))
    print(coefs)  # show all coefficients
    
    # extract all Phrt_Bld_c-related rows
    idx <- grep("^Phrt_Bld_c", rownames(coefs))
    est <- coefs[idx, "Estimate"]
    p <- coefs[idx, "Pr(>|t|)"]
    
    res$type <- "numeric-lm"
    res$estimate <- paste(round(est,3), collapse=", ")
    res$p <- paste(signif(p,3), collapse=", ")
    
  } else {
    xv <- as.factor(xv)
    if (nlevels(xv) == 2) {
      fit <- glm(reformulate("Phrt_Bld_c", response = x), data=set, family=binomial())
      coefs <- coef(summary(fit))
      print(coefs)
      
      idx <- grep("^Phrt_Bld_c", rownames(coefs))
      or <- exp(coef(fit)[idx])
      p <- coefs[idx, "Pr(>|z|)"]
      
      res$type <- "binary-glm"
      res$estimate <- paste(round(or,3), collapse=", ")
      res$p <- paste(signif(p,3), collapse=", ")
      
    } else {
      # fallback: chi-square (crude)
      tab <- table(set[[x]], set[["Phrt_Bld_c"]], useNA="no")
      ch <- suppressWarnings(chisq.test(tab))
      res$type <- "categorical-chi2"
      res$estimate <- NA
      res$p <- ch$p.value
    }
  }
  
  return(res)
})

exposure_cov_table <- do.call(rbind, lapply(exposure_cov_assoc, as.data.frame))
row.names(exposure_cov_table) <- NULL

# (2) covariate ↔ outcome association (Cox) ----
get_hr_ci_multi <- function(fit, term_prefix){
  ci <- as.data.frame(summary(fit)$conf.int)  # force data frame
  rn <- rownames(ci)
  idx <- grep(paste0("^", term_prefix), rn)
  if(length(idx) == 0) return(NULL)
  
  out <- ci[idx, c("exp(coef)", "lower .95", "upper .95"), drop=FALSE]
  out$term <- rn[idx]
  colnames(out)[1:3] <- c("HR","lowCI","upCI")
  rownames(out) <- NULL
  return(out)
}

uni <- get_hr_ci_multi(cox_cov_uni, "Phrt_Bld_c")
adj <- get_hr_ci_multi(cox_cov_adj, "Phrt_Bld_c")

covariate_outcome_table <- rbind(
  cbind(model="Univariate", uni,
        n=summary(cox_cov_uni)$n,
        nevent=summary(cox_cov_uni)$nevent,
        AIC=AIC(cox_cov_uni)),
  cbind(model="Adjusted", adj,
        n=summary(cox_cov_adj)$n,
        nevent=summary(cox_cov_adj)$nevent,
        AIC=AIC(cox_cov_adj))
)

covariate_outcome_table

# (3) change-in-estimate + model fit: Cox without vs with Phrt_Bld_c, flagged at ≥10% HR change ----

# Helper: get HR and CI for a coefficient (or first dummy if factor)
get_hr_ci_multi <- function(fit, term_prefix){
  ci <- as.data.frame(summary(fit)$conf.int)
  rn <- rownames(ci)
  idx <- grep(paste0("^", term_prefix), rn)
  if(length(idx) == 0) return(NULL)
  
  out <- ci[idx, c("exp(coef)", "lower .95", "upper .95"), drop=FALSE]
  out$term <- rn[idx]
  colnames(out)[1:3] <- c("HR","lowCI","upCI")
  rownames(out) <- NULL
  return(out)
}

# Helper: always return a 1-row data.frame
safe_first_hr <- function(hr_table){
  if(is.null(hr_table) || nrow(hr_table) == 0){
    return(data.frame(HR=NA, lowCI=NA, upCI=NA))
  } 
  out <- hr_table[1, c("HR","lowCI","upCI"), drop=FALSE]
  rownames(out) <- NULL
  return(out)
}

results_list <- list()
result_tables <- list()

for (exposure in exposures) {
  # Base model (without Phrt_Bld_c)
  form_base <- as.formula(paste0(
    "Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ",
    exposure, " + bmi_oms + Alc_Re_c + year_diagnosis + ",
    "Er_Adj_c + Pr_Adj_c + Her2_Adj_c + stage3 + menop_status_dx + COUNTRY + gradb + l_school2"
  ))
  mod_base <- coxph(form_base, data=set, id=Idepic)
  
  # Adjusted model (with Phrt_Bld_c)
  form_adj <- update(form_base, . ~ . + Phrt_Bld_c)
  mod_adj <- coxph(form_adj, data=set, id=Idepic)
  
  # Extract HRs safely
  hr_base <- safe_first_hr(get_hr_ci_multi(mod_base, exposure))
  hr_adj  <- safe_first_hr(get_hr_ci_multi(mod_adj, exposure))
  
  # Percent change in HR
  pct_change <- if(!is.na(hr_base$HR) && !is.na(hr_adj$HR)){
    abs(hr_adj$HR - hr_base$HR) / hr_base$HR
  } else NA
  
  # Likelihood ratio test (safe)
  lrt <- tryCatch(anova(mod_base, mod_adj, test="LRT"), error=function(e) NULL)
  lrt_p <- NA
  
  if(!is.null(lrt) && nrow(lrt) >= 2){
    # Look for possible LRT columns
    # Often it is named "Pr(>Chi)", but can also be "P(>|Chi|)" depending on R version
    pcol <- grep("P.*Chi", colnames(lrt), value=TRUE)
    
    if(length(pcol) == 0){
      # If no p-value column exists, compute manually
      chisq <- lrt$Chisq[2]        # second row, chi-square statistic
      df    <- lrt$Df[2]           # degrees of freedom
      if(!is.na(chisq) && !is.na(df) && df > 0){
        lrt_p <- 1 - pchisq(chisq, df)
      }
    } else {
      lrt_p <- lrt[2, pcol[1]]
    }
  }
  lrt_p
  
  # Compile table row
  tbl <- data.frame(
    exposure = exposure,
    HR_base = hr_base$HR, lowCI_base = hr_base$lowCI, upCI_base = hr_base$upCI,
    HR_adj  = hr_adj$HR,  lowCI_adj  = hr_adj$lowCI,  upCI_adj  = hr_adj$upCI,
    pct_change = pct_change,
    AIC_base = AIC(mod_base), AIC_adj = AIC(mod_adj),
    LRT_p = lrt_p,
    n = summary(mod_adj)$n,
    nevent = summary(mod_adj)$nevent,
    confounder_flag_10pct = ifelse(!is.na(pct_change) & pct_change >= 0.10, TRUE, FALSE)
  )
  
  # Store results
  result_tables[[exposure]] <- tbl
  results_list[[exposure]] <- list(base=mod_base, adj=mod_adj)
}

# Combine all exposures into a single table
change_in_estimate_table <- do.call(rbind, result_tables)

# Optional: round numeric columns
num_cols <- c("HR_base","lowCI_base","upCI_base","HR_adj","lowCI_adj","upCI_adj",
              "pct_change","AIC_base","AIC_adj","LRT_p")
change_in_estimate_table[num_cols] <- 
  lapply(change_in_estimate_table[num_cols], function(x) round(as.numeric(x), 3))

# --- OUTPUT ---
change_in_estimate_table
results_list  # all fitted Cox models

#-------------------------------------------------#
# Change-in-estimate: exposure-only vs +Phrt_Bld_c ----
#-------------------------------------------------#

results_list_simple <- list()
result_tables_simple <- list()

for (exposure in exposures) {
  # Base model (exposure only)
  form_base <- as.formula(paste0(
    "Surv(age_dx, Age_Exit_Vs, event == levels(Death_Status)[2]) ~ ",
    exposure
  ))
  mod_base <- coxph(form_base, data=set, id=Idepic)
  
  # Adjusted model (exposure + Phrt_Bld_c only)
  form_adj <- update(form_base, . ~ . + Phrt_Bld_c)
  mod_adj <- coxph(form_adj, data=set, id=Idepic)
  
  # Extract HRs safely
  hr_base <- safe_first_hr(get_hr_ci_multi(mod_base, exposure))
  hr_adj  <- safe_first_hr(get_hr_ci_multi(mod_adj, exposure))
  
  # Percent change in HR
  pct_change <- if(!is.na(hr_base$HR) && !is.na(hr_adj$HR)){
    abs(hr_adj$HR - hr_base$HR) / hr_base$HR
  } else NA
  
  # Likelihood ratio test
  lrt <- tryCatch(anova(mod_base, mod_adj, test="LRT"), error=function(e) NULL)
  lrt_p <- NA
  if(!is.null(lrt) && nrow(lrt) >= 2){
    pcol <- grep("P.*Chi", colnames(lrt), value=TRUE)
    if(length(pcol) == 0){
      chisq <- lrt$Chisq[2]
      df    <- lrt$Df[2]
      if(!is.na(chisq) && !is.na(df) && df > 0){
        lrt_p <- 1 - pchisq(chisq, df)
      }
    } else {
      lrt_p <- lrt[2, pcol[1]]
    }
  }
  
  # Compile table row
  tbl <- data.frame(
    exposure = exposure,
    HR_base = hr_base$HR, lowCI_base = hr_base$lowCI, upCI_base = hr_base$upCI,
    HR_adj  = hr_adj$HR,  lowCI_adj  = hr_adj$lowCI,  upCI_adj  = hr_adj$upCI,
    pct_change = pct_change,
    AIC_base = AIC(mod_base), AIC_adj = AIC(mod_adj),
    LRT_p = lrt_p,
    n = summary(mod_adj)$n,
    nevent = summary(mod_adj)$nevent,
    confounder_flag_10pct = ifelse(!is.na(pct_change) & pct_change >= 0.10, TRUE, FALSE)
  )
  
  result_tables_simple[[exposure]] <- tbl
  results_list_simple[[exposure]] <- list(base=mod_base, adj=mod_adj)
}

# Combine into one table
change_in_estimate_table_simple <- do.call(rbind, result_tables_simple)

# Round numeric columns
num_cols <- c("HR_base","lowCI_base","upCI_base","HR_adj","lowCI_adj","upCI_adj",
              "pct_change","AIC_base","AIC_adj","LRT_p")
change_in_estimate_table_simple[num_cols] <- 
  lapply(change_in_estimate_table_simple[num_cols], function(x) round(as.numeric(x), 3))

# --- OUTPUT ---
change_in_estimate_table_simple

