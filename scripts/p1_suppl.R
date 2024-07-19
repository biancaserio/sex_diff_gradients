# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Project: Sex differences in brain organization

# Supplementary Analyses

# Content: 
# - GSP
# - Linear Model analyses for Age effects in functional gradients and local CT
# - Local SA

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  




# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SET UP 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load packages
require(rstudioapi) # gets path of script  
require(tidyverse)
require(plyr)  # ddply
require(lme4)  # for lmer (linear mixed effects model)
require(lmerTest)  # to obtain p-values for lmer - this actually overrides lme4's lmer() and prints the p-values for the fixed effects (which aren't present in lme4's lmer())
# note that no p-values are printed for the random effects because p values cannot be estimated for random effects because these are latent variables without standard deviations


# Clear environment
rm(list = ls())

# set up directories
codedir = dirname(getActiveDocumentContext()$path)  # get path to current script
resdir_gsp = '/data/p_02667/sex_diff_gradients/results/GSP/'
resdir_hcp = '/data/p_02667/sex_diff_gradients/results/HCP/'
resdir_suppl_gsp = '/data/p_02667/sex_diff_gradients/results/GSP/suppl/'
resdir_suppl_hcp = '/data/p_02667/sex_diff_gradients/results/HCP/suppl/'

# set directory to path of current script
setwd(codedir) 



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DEFINING FUNCTIONS 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# note: these functions are hardcoded for given datasets (variable names, variables included in the regression)


### Linear Regression GSP Functional Connectivity

lm.gsp_fc_sex_contrast <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for SEX effects, including sex, age and ICV in the model as covariates (relevant to functional connectivity)
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for SEX contrast
  '
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val_sex = vector(mode = "double", length = ncol(df_dv))  
  p_val_sex = vector(mode = "double", length = ncol(df_dv))
  
  # Degrees of Freedom = N (subjects) - number of IV (3 *** HARD-CODED FOR THIS SPECIFIC LM MODEL ***) - 1 (mean)
  #DoF = nrow(df_dv) - 3 -1
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    # Fit a linear model: lm = Gradient_Eigenvalues ~ Sex + Age + ICV
    lm_fit = lm(df_dv[[i]] ~ df_iv$Sex + df_iv$Age_Bin + df_iv$ICV)
    
    # Extract from summary of lm_fit the t- and p-values
    # summary(lm_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 Age, 4 ICV; columns = 1 Estimate, 2 Std. Error, 3 t-value, 4 p-value
    t_val_sex[[i]] = summary(lm_fit)$coefficients[2,3]
    p_val_sex[[i]] = summary(lm_fit)$coefficients[2,4]  # if want to calculate p value by hand: p_val_sex[[i]] = 2*pt(abs(t_val_sex[[i]]), DoF, lower.tail = F)
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val_sex = p.adjust(p_val_sex, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val_sex, p_val_sex, q_val_sex)
  
  return(output_df)
}



#### Linear Regression GSP Functional Connectivity - AGE effect

lm.gsp_fc_age_contrast <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for AGE effects, including sex, age and ICV in the model as covariates (relevant to functional connectivity)
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for AGE contrast
  '
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val_age = vector(mode = "double", length = ncol(df_dv))  
  p_val_age = vector(mode = "double", length = ncol(df_dv))
  
  # Degrees of Freedom = N (subjects) - number of IV (3 *** HARD-CODED FOR THIS SPECIFIC LM MODEL ***) - 1 (mean)
  #DoF = nrow(df_dv) - 3 -1
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    # Fit a linear model: lm = Gradient_Eigenvalues ~ Sex + Age + ICV
    lm_fit = lm(df_dv[[i]] ~ df_iv$Sex + df_iv$Age_Bin + df_iv$ICV)
    
    # Extract from summary of lm_fit the t- and p-values
    # summary(lm_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 Age, 4 ICV; columns = 1 Estimate, 2 Std. Error, 3 t-value, 4 p-value
    t_val_age[[i]] = summary(lm_fit)$coefficients[3,3]
    p_val_age[[i]] = summary(lm_fit)$coefficients[3,4]  # if want to calculate p value by hand: p_val_sex[[i]] = 2*pt(abs(t_val_sex[[i]]), DoF, lower.tail = F)
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val_age = p.adjust(p_val_age, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val_age, p_val_age, q_val_age)
  
  return(output_df)
}


### Linear Regression HCP Functional Connectivity - AGE effect

lmer.hcp_fc_age_contrast <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for AGE effects, including sex, age and ICV, as well as random effects family id, twin status and family id * twin status, in the model as covariates (relevant to functional connectivity)
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for AGE contrast
  '
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val_age = vector(mode = "double", length = ncol(df_dv))  
  p_val_age = vector(mode = "double", length = ncol(df_dv))
  beta_val_age = vector(mode = "double", length = ncol(df_dv))
  
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    family_id = df_iv$Family_ID
    twin_status = df_iv$TwinStatus
    
    # Fit a linear mixed effects model: lm = Gradient_Eigenvalues ~ Sex + Age + ICV + random effect(family relatedness) + random effect (family relatedness * twin status)
    # error message: https://stackoverflow.com/questions/60028673/lme4-error-boundary-singular-fit-see-issingular -> Your model did fit, but it generated that warning because your random effects are very small
    # https://stats.stackexchange.com/questions/96600/interactions-between-random-effects
    # https://stackoverflow.com/questions/71340764/interaction-between-two-factors-as-random-effects-in-mixed-model-in-r
    
    # Model including only one of the "single" random effects included in the interaction random effect (i.e., family ID, not including twin status because would group together unrelated subjects who are e.g., also twins) -> I think this would be the correct model to use
    #lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol + (1 | family_id) + (1 | family_id:twin_status), REML = FALSE)
    
    # Model including both "single" random effects included in the interaction random effect (Sofie's decision)
    lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol + (1 | family_id) + (1 | twin_status) + (1 | family_id:twin_status), REML = FALSE)
    
    # Extract from summary of lmer_fit the t- and p-values
    # summary(lmer_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 Age, 4 ICV; columns = 1 Estimate, 2 Std. Error, 3 df, 4 t-value, 5 p-value
    t_val_age[[i]] = summary(lmer_fit)$coefficients[3,4]
    p_val_age[[i]] = summary(lmer_fit)$coefficients[3,5]
    beta_val_age[[i]] = summary(lmer_fit)$coefficients[3,1]
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val_age = p.adjust(p_val_age, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val_age, p_val_age, q_val_age, beta_val_age)
  
  return(output_df)
}



### Linear Regression GSP CT (both local and gradient)

lm.gsp_ct_sex_contrast <- function(df_dv, df_iv) {
  
  '
  - fits and runs linear model to test for SEX effects, including sex, age and global CT in the model as covariates (relevant to CT)
  - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
  - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for SEX contrast
'
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val_sex = vector(mode = "double", length = ncol(df_dv))  
  p_val_sex = vector(mode = "double", length = ncol(df_dv))
  
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    # Fit a linear model: lm = Gradient_Eigenvalues ~ Sex + Age + Global CT
    lm_fit = lm(df_dv[[i]] ~ df_iv$Sex + df_iv$Age_Bin + df_iv$global_ct)
    
    # Extract from summary of lm_fit the t- and p-values
    # summary(lm_fit)$coefficients[row, column]; row = 1 intercept, 2 Sex, 3 Age, 4 Global CT; columns = 1 Estimate, 2 Std. Error, 3 t-value, 4 p-value
    t_val_sex[[i]] = summary(lm_fit)$coefficients[2,3]
    p_val_sex[[i]] = summary(lm_fit)$coefficients[2,4]  # if want to calculate p value by hand: p_val_sex[[i]] = 2*pt(abs(t_val_sex[[i]]), DoF, lower.tail = F)
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val_sex = p.adjust(p_val_sex, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val_sex, p_val_sex, q_val_sex)
  
  return(output_df)
}


### Linear Regression GSP CT (both local and gradient) - AGE effect

lm.gsp_ct_age_contrast <- function(df_dv, df_iv) {
  
  '
  - fits and runs linear model to test for AGE effects, including sex, age and global CT in the model as covariates (relevant to CT)
  - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
  - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for AGE contrast
'
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val_age = vector(mode = "double", length = ncol(df_dv))  
  p_val_age = vector(mode = "double", length = ncol(df_dv))
  
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    # Fit a linear model: lm = Gradient_Eigenvalues ~ Sex + Age + Global CT
    lm_fit = lm(df_dv[[i]] ~ df_iv$Sex + df_iv$Age_Bin + df_iv$global_ct)
    
    # Extract from summary of lm_fit the t- and p-values
    # summary(lm_fit)$coefficients[row, column]; row = 1 intercept, 2 Sex, 3 Age, 4 Global CT; columns = 1 Estimate, 2 Std. Error, 3 t-value, 4 p-value
    t_val_age[[i]] = summary(lm_fit)$coefficients[3,3]
    p_val_age[[i]] = summary(lm_fit)$coefficients[3,4]  # if want to calculate p value by hand: p_val_sex[[i]] = 2*pt(abs(t_val_sex[[i]]), DoF, lower.tail = F)
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val_age = p.adjust(p_val_age, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val_age, p_val_age, q_val_age)
  
  return(output_df)
}





### Linear Regression HCP CT (both local and gradient)

lmer.hcp_ct_sex_contrast <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for SEX effects, including sex, age and global CT, as well as random effects family id, twin status and family id * twin status, in the model as covariates (relevant to CT)
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for SEX contrast
  '
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val_sex = vector(mode = "double", length = ncol(df_dv))  
  p_val_sex = vector(mode = "double", length = ncol(df_dv))
  beta_val_sex = vector(mode = "double", length = ncol(df_dv))
  
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    family_id = df_iv$Family_ID
    twin_status = df_iv$TwinStatus
    
    # Fit a linear mixed effects model: lm = Gradient_Eigenvalues ~ Sex + Age + Global CT + random effect(family relatedness) + random effect(twin status) + random effect (family relatedness * twin status)
    # error message: https://stackoverflow.com/questions/60028673/lme4-error-boundary-singular-fit-see-issingular -> Your model did fit, but it generated that warning because your random effects are very small
    # https://stats.stackexchange.com/questions/96600/interactions-between-random-effects
    # https://stackoverflow.com/questions/71340764/interaction-between-two-factors-as-random-effects-in-mixed-model-in-r
    
    # Model including only one of the "single" random effects included in the interaction random effect (i.e., family ID, not including twin status because would group together unrelated subjects who are e.g., also twins) -> I think this would be the correct model to use
    #lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$global_ct + (1 | family_id) + (1 | family_id:twin_status), REML = FALSE)
    
    # Model including both "single" random effects included in the interaction random effect (Sofie's decision)
    lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$global_ct + (1 | family_id) + (1 | twin_status) + (1 | family_id:twin_status), REML = FALSE)
    
    # Extract from summary of lmer_fit the t- and p-values
    # summary(lmer_fit)$coefficients[row, column]; row = 1 intercept, 2 Sex, 3 Age, 4 Global CT; columns = 1 Estimate, 2 Std. Error, 3 df, 4 t-value, 5 p-value
    t_val_sex[[i]] = summary(lmer_fit)$coefficients[2,4]
    p_val_sex[[i]] = summary(lmer_fit)$coefficients[2,5]
    beta_val_sex[[i]] = summary(lmer_fit)$coefficients[2,1]
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val_sex = p.adjust(p_val_sex, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val_sex, p_val_sex, q_val_sex, beta_val_sex)
  
  return(output_df)
}




### Linear Regression HCP CT (both local and gradient) - AGE effect

lmer.hcp_ct_age_contrast <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for AGE effects, including sex, age and global CT, as well as random effects family id, twin status and family id * twin status, in the model as covariates (relevant to CT)
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for AGE contrast
  '
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val_age = vector(mode = "double", length = ncol(df_dv))  
  p_val_age = vector(mode = "double", length = ncol(df_dv))
  beta_val_age = vector(mode = "double", length = ncol(df_dv))
  
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    family_id = df_iv$Family_ID
    twin_status = df_iv$TwinStatus
    
    # Fit a linear mixed effects model: lm = Gradient_Eigenvalues ~ Sex + Age + Global CT + random effect(family relatedness) + random effect(twin status) + random effect (family relatedness * twin status)
    # error message: https://stackoverflow.com/questions/60028673/lme4-error-boundary-singular-fit-see-issingular -> Your model did fit, but it generated that warning because your random effects are very small
    # https://stats.stackexchange.com/questions/96600/interactions-between-random-effects
    # https://stackoverflow.com/questions/71340764/interaction-between-two-factors-as-random-effects-in-mixed-model-in-r
    
    # Model including only one of the "single" random effects included in the interaction random effect (i.e., family ID, not including twin status because would group together unrelated subjects who are e.g., also twins) -> I think this would be the correct model to use
    #lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$global_ct + (1 | family_id) + (1 | family_id:twin_status), REML = FALSE)
    
    # Model including both "single" random effects included in the interaction random effect (Sofie's decision)
    lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$global_ct + (1 | family_id) + (1 | twin_status) + (1 | family_id:twin_status), REML = FALSE)
    
    # Extract from summary of lmer_fit the t- and p-values
    # summary(lmer_fit)$coefficients[row, column]; row = 1 intercept, 2 Sex, 3 Age, 4 Global CT; columns = 1 Estimate, 2 Std. Error, 3 df, 4 t-value, 5 p-value
    t_val_age[[i]] = summary(lmer_fit)$coefficients[3,4]
    p_val_age[[i]] = summary(lmer_fit)$coefficients[3,5]
    beta_val_age[[i]] = summary(lmer_fit)$coefficients[3,1]
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val_age = p.adjust(p_val_age, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val_age, p_val_age, q_val_age, beta_val_age)
  
  return(output_df)
}




# Sex contrast
lmer.hcp_sex_contrast <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for SEX effects, including sex, age and total surface area, as well as random nested effect(family relatedness/twin status in the model as covariates
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for SEX contrast
  '
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val = vector(mode = "double", length = ncol(df_dv))  
  p_val = vector(mode = "double", length = ncol(df_dv))
  beta_val = vector(mode = "double", length = ncol(df_dv))
  
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    family_id = df_iv$Family_ID
    twin_status = df_iv$TwinStatus
    
    # Fit a linear mixed effects model: DV ~ Sex + Age + tot_SA + random nested effect(family relatedness/twin status)
    lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$tot_SA + (1 | family_id/twin_status), REML = FALSE)
    
    # Extract from summary of lmer_fit the t- and p-values
    # summary(lmer_fit)$coefficients[row, column]; row = 1 intercept, 2 Sex, 3 Age, 4 tot_SA; columns = 1 Estimate, 2 Std. Error, 3 df, 4 t-value, 5 p-value
    t_val[[i]] = summary(lmer_fit)$coefficients[2,4]
    p_val[[i]] = summary(lmer_fit)$coefficients[2,5]
    beta_val[[i]] = summary(lmer_fit)$coefficients[2,1]
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val = p.adjust(p_val, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val, p_val, q_val, beta_val)
  
  return(output_df)
}



### Linear Regression HCP Connectivity lmer = Gradient_Eigenvalues ~ Sex + Age + tot_SA + local_SA + random nested effect(family relatedness/twin status) 

# Local SA contrast
lmer.hcp_local_SA_contrast <- function(df_dv, df_iv, df_local_SA) {
  
  '
    - fits and runs linear model to test for LOCAL SA effects, including sex, age, tot_SA and local SA, as well as random nested effect(family relatedness/twin status in the model as covariates
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables), df_geo (dataframe containing geodesic distance)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for local SA contrast
  '
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val = vector(mode = "double", length = ncol(df_dv))  
  p_val = vector(mode = "double", length = ncol(df_dv))
  beta_val = vector(mode = "double", length = ncol(df_dv))
  
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    family_id = df_iv$Family_ID
    twin_status = df_iv$TwinStatus
    
    # Fit a linear mixed effects model: DV ~ Sex + Age + tot_SA + random nested effect(family relatedness/twin status)
    
    # Model including both "single" random effects included in the interaction random effect (Sofie's decision)
    lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$tot_SA + df_local_SA[[i]] + (1 | family_id/twin_status), REML = FALSE)
    
    # Extract from summary of lmer_fit the t- and p-values
    # summary(lmer_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 age, 4 tot_SA, 5 local_SA; columns = 1 Estimate, 2 Std. Error, 3 df, 4 t-value, 5 p-value
    t_val[[i]] = summary(lmer_fit)$coefficients[5,4]
    p_val[[i]] = summary(lmer_fit)$coefficients[5,5]
    beta_val[[i]] = summary(lmer_fit)$coefficients[5,1]
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val = p.adjust(p_val, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val, p_val, q_val, beta_val)
  
  return(output_df)
}




# ICV contrast in lmer = Gradient_Eigenvalues ~ Sex + Age + ICV + random nested effect(family relatedness/twin status) 
lmer.hcp_icv_contrast <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for ICV effects, including sex, age and ICV, as well as random nested effect(family relatedness/twin status in the model as covariates
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for ICV contrast
  '
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val = vector(mode = "double", length = ncol(df_dv))  
  p_val = vector(mode = "double", length = ncol(df_dv))
  beta_val = vector(mode = "double", length = ncol(df_dv))
  
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    family_id = df_iv$Family_ID
    twin_status = df_iv$TwinStatus
    
    # Fit a linear mixed effects model: DV ~ Sex + Age + ICV + random nested effect(family relatedness/twin status)
    lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol + (1 | family_id/twin_status), REML = FALSE)
    
    # Extract from summary of lmer_fit the t- and p-values
    # summary(lmer_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 age, 4 ICV; columns = 1 Estimate, 2 Std. Error, 3 df, 4 t-value, 5 p-value
    t_val[[i]] = summary(lmer_fit)$coefficients[4,4]
    p_val[[i]] = summary(lmer_fit)$coefficients[4,5]
    beta_val[[i]] = summary(lmer_fit)$coefficients[4,1]
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val = p.adjust(p_val, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val, p_val, q_val, beta_val)
  
  return(output_df)
}





# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PREPARE DATA 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Local CT data
GSP_ct_schaefer400 = read.csv(paste(resdir_gsp, 'ct_schaefer400.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_ct_schaefer400 = read.csv(paste(resdir_hcp, 'ct_schaefer400.csv', sep = ''), fileEncoding = 'UTF-8-BOM')


# Aligned CT gradient loadings
GSP_array_aligned_ct_G1 = read.csv(paste(resdir_gsp, 'array_aligned_ct_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
GSP_array_aligned_ct_G2 = read.csv(paste(resdir_gsp, 'array_aligned_ct_G2.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

HCP_array_aligned_ct_G1 = read.csv(paste(resdir_hcp, 'array_aligned_ct_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_array_aligned_ct_G2 = read.csv(paste(resdir_hcp, 'array_aligned_ct_G2.csv', sep = ''), fileEncoding = 'UTF-8-BOM')



# Aligned GSP functional gradient loadings 
GSP_array_aligned_fc_G1 = read.csv(paste(resdir_gsp, 'array_aligned_fc_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
GSP_array_aligned_fc_G2 = read.csv(paste(resdir_gsp, 'array_aligned_fc_G2.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
GSP_array_aligned_fc_G3 = read.csv(paste(resdir_gsp, 'array_aligned_fc_G3.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

# Local CT data
GSP_ct_schaefer400 = read.csv(paste(resdir_gsp, 'ct_schaefer400.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_ct_schaefer400 = read.csv(paste(resdir_hcp, 'ct_schaefer400.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

# Local SA
local_SA = read.csv(paste(resdir_hcp, 'local_SA.csv', sep = ''), fileEncoding = 'UTF-8-BOM')


# class(array_aligned_G1)
# typeof(array_aligned_G1)
# dim(GSP_array_aligned_G1)


# Descriptives 
GSP_demographics_cleaned = read.csv(paste(resdir_gsp, 'demographics_cleaned.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_demographics_cleaned_final = read.csv(paste(resdir_hcp, 'demographics_cleaned_final.csv', sep = ''), fileEncoding = 'UTF-8-BOM')  # this doesn't include mean CT (so would need to run that script again if I want to run CT analyses)

#str(HCP_demographics_cleaned, list.len=ncol(HCP_demographics_cleaned))


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LINEAR REGRESSIONS
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#### LOCAL CT

### GSP: SEX effects in model = Local CT ~ Sex + Age + Global CT 

# run model
GSP_lm_ct_local_sex_contrast_res = lm.gsp_ct_sex_contrast(df_dv = GSP_ct_schaefer400, df_iv = GSP_demographics_cleaned)

# number of significant parcels
sum(GSP_lm_ct_local_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(GSP_lm_G1_res$q_val_sex < 0.05))


### HCP: SEX effects in model = Local CT ~ Sex + Age + Global CT + random effect (family relatedness) + random effect (twin status) + random effect (family relatedness * twin status) 

# run model
HCP_lmer_ct_local_sex_contrast_res = lmer.hcp_ct_sex_contrast(df_dv = HCP_ct_schaefer400, df_iv = HCP_demographics_cleaned)

# number of significant parcels
sum(HCP_lmer_ct_local_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(GSP_lm_G1_res$q_val_sex < 0.05))



#### CT GRADIENTS

### GSP: SEX effects in model = Gradient Loadings ~ Sex + Age + Global CT 

# run model

GSP_lm_ct_G1_sex_contrast_res = lm.gsp_ct_sex_contrast(df_dv = GSP_array_aligned_ct_G1, df_iv = GSP_demographics_cleaned)
GSP_lm_ct_G2_sex_contrast_res = lm.gsp_ct_sex_contrast(df_dv = GSP_array_aligned_ct_G2, df_iv = GSP_demographics_cleaned)

# number of significant parcels
sum(GSP_lm_ct_G1_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(GSP_lm_G1_res$q_val_sex < 0.05))
sum(GSP_lm_ct_G2_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)


### HCP: SEX effects in model = Gradient Loadings ~ Sex + Age + Global CT + random effect (family relatedness) + random effect (twin status) + random effect (family relatedness * twin status) 

# run model
HCP_lmer_ct_G1_sex_contrast_res = lmer.hcp_ct_sex_contrast(df_dv = HCP_array_aligned_ct_G1, df_iv = HCP_demographics_cleaned)
HCP_lmer_ct_G2_sex_contrast_res = lmer.hcp_ct_sex_contrast(df_dv = HCP_array_aligned_ct_G2, df_iv = HCP_demographics_cleaned)

# number of significant parcels
sum(HCP_lmer_ct_G1_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(GSP_lm_G1_res$q_val_sex < 0.05))
sum(HCP_lmer_ct_G2_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE) 




##### FUNCTIONAL 

### GSP: SEX effects in model = Gradient_Eigenvalues ~ Sex + Age + ICV 

# run model
GSP_lm_fc_G1_sex_contrast_res = lm.gsp_fc_sex_contrast(df_dv = GSP_array_aligned_fc_G1, df_iv = GSP_demographics_cleaned)
GSP_lm_fc_G2_sex_contrast_res = lm.gsp_fc_sex_contrast(df_dv = GSP_array_aligned_fc_G2, df_iv = GSP_demographics_cleaned)
GSP_lm_fc_G3_sex_contrast_res = lm.gsp_fc_sex_contrast(df_dv = GSP_array_aligned_fc_G3, df_iv = GSP_demographics_cleaned)

# number of significant parcels
sum(GSP_lm_fc_G1_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(GSP_lm_G1_res$q_val_sex < 0.05))
sum(GSP_lm_fc_G2_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  
sum(GSP_lm_fc_G3_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  



### GSP: AGE effects in model = Gradient_Eigenvalues ~ Sex + Age + ICV 

# run model
GSP_lm_fc_G1_age_contrast_res = lm.gsp_fc_age_contrast(df_dv = GSP_array_aligned_fc_G1, df_iv = GSP_demographics_cleaned)
GSP_lm_fc_G2_age_contrast_res = lm.gsp_fc_age_contrast(df_dv = GSP_array_aligned_fc_G2, df_iv = GSP_demographics_cleaned)
GSP_lm_fc_G3_age_contrast_res = lm.gsp_fc_age_contrast(df_dv = GSP_array_aligned_fc_G3, df_iv = GSP_demographics_cleaned)

# number of significant parcels
sum(GSP_lm_fc_G1_age_contrast_res$q_val_age < 0.05, na.rm=TRUE)  # other way: length(which(GSP_lm_G1_res$q_val_age < 0.05))
sum(GSP_lm_fc_G2_age_contrast_res$q_val_age < 0.05, na.rm=TRUE)  
sum(GSP_lm_fc_G3_age_contrast_res$q_val_age < 0.05, na.rm=TRUE)  



### HCP: AGE effects in model = Gradient_Eigenvalues ~ Sex + Age + ICV + random effect (family relatedness) + random effect (twin status) + random effect (family relatedness * twin status) 

# run model
HCP_lmer_fc_G1_age_contrast_res = lmer.hcp_fc_age_contrast(df_dv = HCP_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned)
HCP_lmer_fc_G2_age_contrast_res = lmer.hcp_fc_age_contrast(df_dv = HCP_array_aligned_fc_G2, df_iv = HCP_demographics_cleaned)
HCP_lmer_fc_G3_age_contrast_res = lmer.hcp_fc_age_contrast(df_dv = HCP_array_aligned_fc_G3, df_iv = HCP_demographics_cleaned)

# number of significant parcels
sum(HCP_lmer_fc_G1_age_contrast_res$q_val_age < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_age < 0.05))
sum(HCP_lmer_fc_G2_age_contrast_res$q_val_age < 0.05, na.rm=TRUE)  
sum(HCP_lmer_fc_G3_age_contrast_res$q_val_age < 0.05, na.rm=TRUE)  



##### STRUCTURAL

#### LOCAL CT

### GSP: AGE effects in model = Local CT ~ Sex + Age + Global CT 

# run model
GSP_lm_ct_local_age_contrast_res = lm.gsp_ct_age_contrast(df_dv = GSP_ct_schaefer400, df_iv = GSP_demographics_cleaned)

# number of significant parcels
sum(GSP_lm_ct_local_age_contrast_res$q_val_age < 0.05, na.rm=TRUE)  # other way: length(which(GSP_lm_G1_res$q_val_age < 0.05))


### HCP: AGE effects in model = Local CT ~ Sex + Age + Global CT + random effect (family relatedness) + random effect (twin status) + random effect (family relatedness * twin status) 

# run model
HCP_lmer_ct_local_age_contrast_res = lmer.hcp_ct_age_contrast(df_dv = HCP_ct_schaefer400, df_iv = HCP_demographics_cleaned)

# number of significant parcels
sum(HCP_lmer_ct_local_age_contrast_res$q_val_age < 0.05, na.rm=TRUE)  # other way: length(which(GSP_lm_G1_res$q_val_age < 0.05))





#### LOCAL SA

## sex effects in local SA

# run model
HCP_lmer_local_SA_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = local_SA, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_local_SA_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  



## Local SA contrast in lmer = Gradient_Eigenvalues ~ Sex + Age + tot_SA + local_SA + random nested effect(family relatedness/twin status) 

# run model
HCP_lmer_hemi_fc_G1_local_SA_contrast_res = lmer.hcp_local_SA_contrast(df_dv = HCP_hemi_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final, df_local_SA = local_SA)
HCP_lmer_fc_G1_local_SA_contrast_res = lmer.hcp_local_SA_contrast(df_dv = HCP_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final, df_local_SA = local_SA)

# number of significant parcels
sum(HCP_lmer_hemi_fc_G1_local_SA_contrast_res$q_val < 0.05, na.rm=TRUE)  
sum(HCP_lmer_fc_G1_local_SA_contrast_res$q_val < 0.05, na.rm=TRUE)  




############ not used, from p1_main.R -> might need at some point? ##############

# run model effects of ICV on mean geodesic distance of top10% connections computed at the individual level
HCP_lmer_geo_icv_contrast_res = lmer.hcp_icv_contrast(df_dv = HCP_mean_geodesic_distances, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_geo_icv_contrast_res$q_val < 0.05, na.rm=TRUE) 





# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# EXPORT RESULTS 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Local CT SEX effects

write.csv(GSP_lm_ct_local_sex_contrast_res, paste(resdir_gsp, 'R_lm_ct_local_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_ct_local_sex_contrast_res, paste(resdir_hcp, 'R_lmer_ct_local_sex_contrast_res.csv', sep = ''), row.names = FALSE)


# Gradients CT SEX effects

write.csv(GSP_lm_ct_G1_sex_contrast_res, paste(resdir_gsp, 'R_lm_ct_G1_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(GSP_lm_ct_G2_sex_contrast_res, paste(resdir_gsp, 'R_lm_ct_G2_sex_contrast_res.csv', sep = ''), row.names = FALSE)

write.csv(HCP_lmer_ct_G1_sex_contrast_res, paste(resdir_hcp, 'R_lmer_ct_G1_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_ct_G2_sex_contrast_res, paste(resdir_hcp, 'R_lmer_ct_G2_sex_contrast_res.csv', sep = ''), row.names = FALSE)



### Functional (Gradients) AGE effects

write.csv(GSP_lm_fc_G1_sex_contrast_res, paste(resdir_gsp, 'R_lm_fc_G1_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(GSP_lm_fc_G2_sex_contrast_res, paste(resdir_gsp, 'R_lm_fc_G2_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(GSP_lm_fc_G3_sex_contrast_res, paste(resdir_gsp, 'R_lm_fc_G3_sex_contrast_res.csv', sep = ''), row.names = FALSE)

write.csv(GSP_lm_fc_G1_age_contrast_res, paste(resdir_suppl_gsp, 'R_lm_fc_G1_age_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(GSP_lm_fc_G2_age_contrast_res, paste(resdir_suppl_gsp, 'R_lm_fc_G2_age_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(GSP_lm_fc_G3_age_contrast_res, paste(resdir_suppl_gsp, 'R_lm_fc_G3_age_contrast_res.csv', sep = ''), row.names = FALSE)

write.csv(HCP_lmer_fc_G1_age_contrast_res, paste(resdir_suppl_hcp, 'R_lmer_fc_G1_age_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_fc_G2_age_contrast_res, paste(resdir_suppl_hcp, 'R_lmer_fc_G2_age_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_fc_G3_age_contrast_res, paste(resdir_suppl_hcp, 'R_lmer_fc_G3_age_contrast_res.csv', sep = ''), row.names = FALSE)



### Structural

# Local CT AGE effects

write.csv(GSP_lm_ct_local_age_contrast_res, paste(resdir_suppl_gsp, 'R_lm_ct_local_age_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_ct_local_age_contrast_res, paste(resdir_suppl_hcp, 'R_lmer_ct_local_age_contrast_res.csv', sep = ''), row.names = FALSE)



### Local SA

# sex effects in local SA
write.csv(HCP_lmer_local_SA_sex_contrast_res, paste(resdir_hcp, 'R_lmer_local_SA_sex_contrast_res.csv', sep = ''), row.names = FALSE)

# Local SA contrast in lmer = Gradient_Eigenvalues ~ Sex + Age + tot_SA + local_SA + random nested effect(family relatedness/twin status) 
write.csv(HCP_lmer_hemi_fc_G1_local_SA_contrast_res, paste(resdir_hcp, 'R_lmer_hemi_fc_G1_local_SA_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_fc_G1_local_SA_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G1_local_SA_contrast_res.csv', sep = ''), row.names = FALSE)




############ not used, from p1_main.R -> might need at some point? ##############
write.csv(HCP_lmer_geo_icv_contrast_res, paste(resdir_hcp, 'R_lmer_geo_icv_contrast_res.csv', sep = ''), row.names = FALSE)






### power analysis
library(pwrss)

power.z.test(ncp = 0.05, alpha = 0.05, 
             alternative = "not equal", plot = TRUE)
