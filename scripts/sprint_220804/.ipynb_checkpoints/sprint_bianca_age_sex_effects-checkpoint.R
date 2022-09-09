---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  
  # Lab Sprint: Multiscale organization of social brain networks
  # 2-4 August 2022
  # Lina H. Schaare, Şeyma Bayrak, Yun-Shuang Fan, Benjamin Hänisch, Meike D. Hettwer, Alexandra John, Svenja Küchenhoff, Neville Magielse, Katerina Manoli, Amin Saberi, Bianca Serio, Bin Wan, Matthias Schurz, Sofie L. Valk
  
  # Content: Linear regression analyses for sex and age contrast in gradient eigenvalues
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  
  
  
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SET UP 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
# Load packages
require(rstudioapi) # gets path of script  
require(tidyverse)
require(plyr)  # ddply


# Clear environment
rm(list = ls())

# set up directories
codedir = dirname(getActiveDocumentContext()$path)  # get path to current script
datadir = '/data/pt_02542/social_networks/data/' 
resdir = '/data/pt_02542/social_networks/results/age_sex/'

# set directory to path of current script
setwd(codedir) 



----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DEFINING FUNCTIONS 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Descriptives

summarise.descriptives <- function(df, variable){
  
  '
  - summarizes main descriptives: mean, sd, min, max, sem, IQR
  - to supply: dataframe and variable of interest
  - was not able to make summarise by Gender work
  '
  
  summarise(df,
            mean = mean(variable),
            sd = sd(variable),
            min = min(variable),
            max = max(variable),
            sem = sd(variable)/sqrt(length(variable)),
            IQR = IQR(variable))
}




### Linear Regression Models : Function definitions 

lm.SEX_age_icv <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for SEX effects, including sex, age and ICV in the model as covariates
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for SEX contrast
  '
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val_sex = vector(mode = "double", length = ncol(df_dv))  
  p_val_sex = vector(mode = "double", length = ncol(df_dv))
  
  # Degrees of Freedom = N (subjects) - number of IV (3 *** HARD-CODED FOR THIS SPECIFIC LM MODEL ***) - 1 (mean)
  DoF = nrow(df_dv) - 3 -1
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    # Fit a linear model: lm = Gradient_Eigenvalues ~ Sex + Age + ICV
    lm_fit = lm(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol)
    
    # Extract from summary of lm_fit the t- and p-values
    # summary(lm_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 Age, 4 ICV; columns = 1 Estimate, 2 Std. Error, 3 t-value, 4 p-value
    t_val_sex[[i]] = summary(lm_fit)$coefficients[2,3]
    p_val_sex[[i]] = 2*pt(abs(t_val_sex[[i]]), DoF, lower.tail = F)  # calculating p value by hand but could directly obtain it from summary(lm_fit)$coefficients[2,4]
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val_sex = p.adjust(p_val_sex, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val_sex, p_val_sex, q_val_sex)
  
  return(output_df)
}


lm.sex_AGE_icv <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for AGE effects, including sex, age and ICV in the model as covariates
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for AGE contrast
  '
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val_sex = vector(mode = "double", length = ncol(df_dv))  
  p_val_sex = vector(mode = "double", length = ncol(df_dv))
  
  # Degrees of Freedom = N (subjects) - number of IV (3 *** HARD-CODED FOR THIS SPECIFIC LM MODEL ***) - 1 (mean)
  DoF = nrow(df_dv) - 3 -1
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    # Fit a linear model: lm = Gradient_Eigenvalues ~ Sex + Age + ICV
    lm_fit = lm(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol)
    
    # Extract from summary of lm_fit the t- and p-values
    # summary(lm_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 Age, 4 ICV; columns = 1 Estimate, 2 Std. Error, 3 t-value, 4 p-value
    t_val_sex[[i]] = summary(lm_fit)$coefficients[3,3]
    p_val_sex[[i]] = 2*pt(abs(t_val_sex[[i]]), DoF, lower.tail = F)  # calculating p value by hand but could directly obtain it from summary(lm_fit)$coefficients[2,4]
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val_sex = p.adjust(p_val_sex, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val_sex, p_val_sex, q_val_sex)
  
  return(output_df)
}



# functions below are specifcially for MPC analyses (different because need to inlcude mean MPC as covariate in model)

lm.SEX_age_icv_meanmpc <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for SEX effects *MPC analyses*, including sex, age, ICV and mean MPC in the model as covariates
       - mean MPC (mean across all parcels/vertices (covariate) per subject) is included as a covaiate in model because standard procedure in MPC analyses to control for it (because of parcel artifacts for each specific subject (intensity of t1/t2 could vary by subject))
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for SEX contrast
  '
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val_sex = vector(mode = "double", length = ncol(df_dv))  
  p_val_sex = vector(mode = "double", length = ncol(df_dv))
  
  # Degrees of Freedom = N (subjects) - number of IV (4 *** HARD-CODED FOR THIS SPECIFIC LM MODEL ***) - 1 (mean)
  DoF = nrow(df_dv) - 4 -1
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    # Fit a linear model: lm = Gradient_Eigenvalues ~ Sex + Age + ICV + mean MPC
    lm_fit = lm(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol + df_iv$mpc_mean)
    
    # Extract from summary of lm_fit the t- and p-values
    # summary(lm_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 Age, 4 ICV, 5 mean MPC; columns = 1 Estimate, 2 Std. Error, 3 t-value, 4 p-value
    t_val_sex[[i]] = summary(lm_fit)$coefficients[2,3]
    p_val_sex[[i]] = 2*pt(abs(t_val_sex[[i]]), DoF, lower.tail = F)  # calculating p value by hand but could directly obtain it from summary(lm_fit)$coefficients[2,4]
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val_sex = p.adjust(p_val_sex, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val_sex, p_val_sex, q_val_sex)
  
  return(output_df)
}



lm.sex_AGE_icv_meanmpc <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for AGE effects in *MPC analyses*, including sex, age, ICV and mean MPC in the model as covariates
       - mean MPC (mean across all parcels/vertices (covariate) per subject) is included as a covaiate in model because standard procedure in MPC analyses to control for it (because of parcel artifacts for each specific subject (intensity of t1/t2 could vary by subject))
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for AGE contrast
  '
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val_sex = vector(mode = "double", length = ncol(df_dv))  
  p_val_sex = vector(mode = "double", length = ncol(df_dv))
  
  # Degrees of Freedom = N (subjects) - number of IV (4 *** HARD-CODED FOR THIS SPECIFIC LM MODEL ***) - 1 (mean)
  DoF = nrow(df_dv) - 4 -1
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    # Fit a linear model: lm = Gradient_Eigenvalues ~ Sex + Age + ICV + mean MPC
    lm_fit = lm(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol + df_iv$mpc_mean)
    
    # Extract from summary of lm_fit the t- and p-values
    # summary(lm_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 Age, 4 ICV, 5 mean MPC; columns = 1 Estimate, 2 Std. Error, 3 t-value, 4 p-value
    t_val_sex[[i]] = summary(lm_fit)$coefficients[3,3]
    p_val_sex[[i]] = 2*pt(abs(t_val_sex[[i]]), DoF, lower.tail = F)  # calculating p value by hand but could directly obtain it from summary(lm_fit)$coefficients[2,4]
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val_sex = p.adjust(p_val_sex, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val_sex, p_val_sex, q_val_sex)
  
  return(output_df)
}



----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PREPARE DATA 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Functional gradient eigenvalues for sjh parcellation (462 parcels)
array_G1_fc_between = read.csv(paste(resdir, 'array_G1_fc_between.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
array_G1_fc_within = read.csv(paste(resdir, 'array_G1_fc_within.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

# MPC pairwise gradient eigenvalues for sjh parcellation (462 parcels)
array_G1_mpc_pairwise_between = read.csv(paste(resdir, 'array_G1_mpc_pairwise_between.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
array_G1_mpc_pairwise_within = read.csv(paste(resdir, 'array_G1_mpc_pairwise_within.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

# class(array_aligned_G1)
# typeof(array_aligned_G1)
dim(array_G1_fc_within)
dim(array_G1_mpc_pairwise_within)

length(array_G1_mpc_pairwise_within$X1)  # parcel 1 out of 462 (length = n subjects = 426)

means_mpc_9197v = read.csv(paste(resdir, 'means_mpc_9197v.csv', sep = ''), fileEncoding = 'UTF-8-BOM')$X..mpc_means


# Descriptives 
merged_demographics_cleaned = read.csv(paste(datadir, 'HCP_426unrelatedSub_phenotype.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
merged_demographics_cleaned$mpc_mean <- means_mpc_9197v


----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DESCRIPTIVES 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
str(merged_demographics_cleaned, list.len=ncol(merged_demographics_cleaned))  # untruncated output
dim(merged_demographics_cleaned)


# Sample size by sex
ftable(merged_demographics_cleaned$Gender)


# Age
summarise.descriptives(merged_demographics_cleaned, merged_demographics_cleaned$Age_in_Yrs)

# Age by sex
ddply(merged_demographics_cleaned, 'Gender', summarise,
      mean = mean(Age_in_Yrs),
      sd = sd(Age_in_Yrs),
      min = min(Age_in_Yrs),
      max = max(Age_in_Yrs),
      sem = sd(Age_in_Yrs)/sqrt(length(Age_in_Yrs)),
      IQR = IQR(Age_in_Yrs))


# ICV
summarise.descriptives(merged_demographics_cleaned, merged_demographics_cleaned$FS_IntraCranial_Vol)

# ICV by sex
ddply(merged_demographics_cleaned, 'Gender', summarise,
      mean = mean(FS_IntraCranial_Vol),
      sd = sd(FS_IntraCranial_Vol),
      min = min(FS_IntraCranial_Vol),
      max = max(FS_IntraCranial_Vol),
      sem = sd(FS_IntraCranial_Vol)/sqrt(length(FS_IntraCranial_Vol)),
      IQR = IQR(FS_IntraCranial_Vol))



----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LINEAR REGRESSION 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


### Run models

## fc G1 -> model = G1 ~ Sex + Age + ICV 
  
# sex effects
lm_G1_fc_between_SEX_age_icv_res = lm.SEX_age_icv(df_dv = array_G1_fc_between, df_iv = merged_demographics_cleaned)
lm_G1_fc_within_SEX_age_icv_res = lm.SEX_age_icv(df_dv = array_G1_fc_within, df_iv = merged_demographics_cleaned)

# age effects
lm_G1_fc_between_sex_AGE_icv_res = lm.sex_AGE_icv(df_dv = array_G1_fc_between, df_iv = merged_demographics_cleaned)
lm_G1_fc_within_sex_AGE_icv_res = lm.sex_AGE_icv(df_dv = array_G1_fc_within, df_iv = merged_demographics_cleaned)


## mpc G1 -> model = G1 ~ Sex + Age + ICV + mean(MPC)

# sex effects
lm_G1_mpc_pairwise_between_SEX_age_icv_meanmpc_res = lm.SEX_age_icv_meanmpc(df_dv = array_G1_mpc_pairwise_between, df_iv = merged_demographics_cleaned)
lm_G1_mpc_pairwise_within_SEX_age_icv_meanmpc_res = lm.SEX_age_icv_meanmpc(df_dv = array_G1_mpc_pairwise_within, df_iv = merged_demographics_cleaned)


# age effects
lm_G1_mpc_pairwise_between_sex_AGE_icv_meanmpc_res = lm.sex_AGE_icv_meanmpc(df_dv = array_G1_mpc_pairwise_between, df_iv = merged_demographics_cleaned)
lm_G1_mpc_pairwise_within_sex_AGE_icv_meanmpc_res = lm.sex_AGE_icv_meanmpc(df_dv = array_G1_mpc_pairwise_within, df_iv = merged_demographics_cleaned)



### number of significant parcels

## fc G1 -> model = G1 ~ Sex + Age + ICV 

# sex effects
sum(lm_G1_fc_between_SEX_age_icv_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))
sum(lm_G1_fc_within_SEX_age_icv_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))

# age effects
sum(lm_G1_fc_between_sex_AGE_icv_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))
sum(lm_G1_fc_within_sex_AGE_icv_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))


## mpc G1 -> model = G1 ~ Sex + Age + ICV + mean(MPC)

# sex effects
sum(lm_G1_mpc_pairwise_between_SEX_age_icv_meanmpc_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))
sum(lm_G1_mpc_pairwise_within_SEX_age_icv_meanmpc_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))


# age effects
sum(lm_G1_mpc_pairwise_between_sex_AGE_icv_meanmpc_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))
sum(lm_G1_mpc_pairwise_within_sex_AGE_icv_meanmpc_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))


----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# EXPORT RESULTS 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

write.csv(lm_G1_fc_between_SEX_age_icv_res, paste(resdir, 'lm_G1_fc_between_SEX_age_icv_res.csv', sep = ''), row.names = FALSE)
write.csv(lm_G1_fc_within_SEX_age_icv_res, paste(resdir, 'lm_G1_fc_within_SEX_age_icv_res.csv', sep = ''), row.names = FALSE)
write.csv(lm_G1_fc_between_sex_AGE_icv_res, paste(resdir, 'lm_G1_fc_between_sex_AGE_icv_res.csv', sep = ''), row.names = FALSE)
write.csv(lm_G1_fc_within_sex_AGE_icv_res, paste(resdir, 'lm_G1_fc_within_sex_AGE_icv_res.csv', sep = ''), row.names = FALSE)

write.csv(lm_G1_mpc_pairwise_between_SEX_age_icv_meanmpc_res, paste(resdir, 'lm_G1_mpc_pairwise_between_SEX_age_icv_meanmpc_res.csv', sep = ''), row.names = FALSE)
write.csv(lm_G1_mpc_pairwise_within_SEX_age_icv_meanmpc_res, paste(resdir, 'lm_G1_mpc_pairwise_within_SEX_age_icv_meanmpc_res.csv', sep = ''), row.names = FALSE)
write.csv(lm_G1_mpc_pairwise_between_sex_AGE_icv_meanmpc_res, paste(resdir, 'lm_G1_mpc_pairwise_between_sex_AGE_icv_meanmpc_res.csv', sep = ''), row.names = FALSE)
write.csv(lm_G1_mpc_pairwise_within_sex_AGE_icv_meanmpc_res, paste(resdir, 'lm_G1_mpc_pairwise_within_sex_AGE_icv_meanmpc_res.csv', sep = ''), row.names = FALSE)
