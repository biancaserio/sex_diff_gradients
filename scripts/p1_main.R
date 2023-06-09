# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Project: Sex differences in brain organization

# Main script

# Content: Linear Mixed Effects Model for sex contrast in gradient eigenvalues controlling for sex, age, *ICV*, nested random effect of family relatedness and twin status

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
datadir = '/data/p_02667/sex_diff_gradients/data/'
resdir_hcp = '/data/p_02667/sex_diff_gradients/results/HCP/'

# set directory to path of current script
setwd(codedir) 



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DEFINING FUNCTIONS 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# note: these functions are hardcoded for given datasets (variable names, variables included in the regression)


### Linear Regression HCP Connectivity lmer = Gradient_Eigenvalues ~ Sex + Age + ICV + random nested effect(family relatedness/twin status) 

# Sex contrast
lmer.hcp_sex_contrast <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for SEX effects, including sex, age and ICV, as well as random nested effect(family relatedness/twin status in the model as covariates
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
    
    # Fit a linear mixed effects model: DV ~ Sex + Age + ICV + random nested effect(family relatedness/twin status)
    
    # Model including both "single" random effects included in the interaction random effect (Sofie's decision)
    lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol + (1 | family_id/twin_status), REML = FALSE)
    
    # Extract from summary of lmer_fit the t- and p-values
    # summary(lmer_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 Age, 4 ICV; columns = 1 Estimate, 2 Std. Error, 3 df, 4 t-value, 5 p-value
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




# ICV contrast
lmer.hcp_icv_contrast <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for ICV effects, including sex, age and ICV, as well as random nested effect(family relatedness/twin status in the model as covariates
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
    
    # Fit a linear mixed effects model: DV ~ Sex + Age + ICV + random nested effect(family relatedness/twin status)
    
    # Model including both "single" random effects included in the interaction random effect (Sofie's decision)
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



### Linear Regression HCP Connectivity lmer = Gradient_Eigenvalues ~ Sex + Age + ICV + geodesic distance + random nested effect(family relatedness/twin status) 

# Geodesic distance contrast
lmer.hcp_geo_contrast <- function(df_dv, df_iv, df_geo) {
  
  '
    - fits and runs linear model to test for GEODESIC DISTANCE effects, including sex, age, geodesic distance and ICV, as well as random nested effect(family relatedness/twin status in the model as covariates
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables), df_geo (dataframe containing geodesic distance)
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
    
    # Fit a linear mixed effects model: DV ~ Sex + Age + ICV + random nested effect(family relatedness/twin status)
    
    # Model including both "single" random effects included in the interaction random effect (Sofie's decision)
    lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol + df_geo[[i]] + (1 | family_id/twin_status), REML = FALSE)
    
    # Extract from summary of lmer_fit the t- and p-values
    # summary(lmer_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 age, 4 ICV, 5 geodesic distance; columns = 1 Estimate, 2 Std. Error, 3 df, 4 t-value, 5 p-value
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




### Linear Regression HCP Connectivity lmer = Gradient_Eigenvalues ~ Sex + Age + ICV + MPC G1 + random nested effect(family relatedness/twin status) 

# MPC G1 contrast
lmer.hcp_mpc_contrast <- function(df_dv, df_iv, df_mpc) {
  
  '
    - fits and runs linear model to test for MPC G1 effects, including sex, age, MPC G1 and ICV, as well as random nested effect(family relatedness/twin status in the model as covariates
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables), df_mpc (dataframe containing MPC G1)
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
    
    # Fit a linear mixed effects model: DV ~ Sex + Age + ICV + random nested effect(family relatedness/twin status)
    
    # Model including both "single" random effects included in the interaction random effect (Sofie's decision)
    lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol + df_mpc[[i]] + (1 | family_id/twin_status), REML = FALSE)
    
    # Extract from summary of lmer_fit the t- and p-values
    # summary(lmer_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 age, 4 ICV, 5 MPC G1; columns = 1 Estimate, 2 Std. Error, 3 df, 4 t-value, 5 p-value
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



### Linear Regression HCP Connectivity lmer = Gradient_Eigenvalues ~ Age + ICV + random nested effect(family relatedness/twin status)

# ICV contrast 
lmer.hcp_icv_contrast_within_sex <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for ICV effects, including age and ICV, as well as random nested effect(family relatedness/twin status in the model as covariates
    - not controlling for sex because running this within sex
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
    
    # Fit a linear mixed effects model: DV ~ Age + ICV + random nested effect(family relatedness/twin status)
    
    # Model including both "single" random effects included in the interaction random effect (Sofie's decision)
    lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol + (1 | family_id/twin_status), REML = FALSE)
    
    # Extract from summary of lmer_fit the t- and p-values
    # summary(lmer_fit)$coefficients[row, column]; row = 1 intercept, 2 Age, 3 ICV; columns = 1 Estimate, 2 Std. Error, 3 df, 4 t-value, 5 p-value
    t_val[[i]] = summary(lmer_fit)$coefficients[3,4]
    p_val[[i]] = summary(lmer_fit)$coefficients[3,5]
    beta_val[[i]] = summary(lmer_fit)$coefficients[3,1]
    
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

# Aligned functional gradient loadings 
HCP_array_aligned_fc_G1 = read.csv(paste(resdir_hcp, 'array_aligned_fc_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_array_aligned_fc_G2 = read.csv(paste(resdir_hcp, 'array_aligned_fc_G2.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_array_aligned_fc_G3 = read.csv(paste(resdir_hcp, 'array_aligned_fc_G3.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

HCP_array_aligned_fc_G1_rescaled = read.csv(paste(resdir_hcp, 'array_aligned_fc_G1_rescaled.csv', sep = ''), fileEncoding = 'UTF-8-BOM')  # rescaled between 0-1

HCP_hemi_array_aligned_fc_G1 = read.csv(paste(resdir_hcp, 'hemi_array_aligned_fc_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')


# Aligned MPC gradient loadings
HCP_array_aligned_mpc_G1 = read.csv(paste(resdir_hcp, 'array_aligned_mpc_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_array_aligned_mpc_G2 = read.csv(paste(resdir_hcp, 'array_aligned_mpc_G2.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_array_aligned_mpc_G3 = read.csv(paste(resdir_hcp, 'array_aligned_mpc_G3.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

# Geodesic distances (mean of top 10% of functional connections)
HCP_mean_geodesic_distances = read.csv(paste(resdir_hcp, 'mean_geodesic_distances.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

# Mean geodesic distances of mean connectivity profiles by sex (ie each subject has the mean (by seed region) of the geodesic distance of the most frequent top 10% connections for their sex)
HCP_recomp_sub_mean_geodesic_distances_most_frequent_connections_bysexprofile = read.csv(paste(resdir_hcp, 'recomp_sub_mean_geodesic_distances_most_frequent_connections_bysexprofile.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

# Mean fc strength by seed for top 10% connections computed at the individual level
HCP_mean_fc_strengths_top10 = read.csv(paste(resdir_hcp, 'mean_fc_strengths_top10.csv', sep = ''), fileEncoding = 'UTF-8-BOM')


# class(array_aligned_G1)
# typeof(array_aligned_G1)
dim(HCP_array_aligned_fc_G1)
dim(HCP_hemi_array_aligned_fc_G1)
dim(HCP_array_aligned_mpc_G1)
dim(HCP_mean_geodesic_distances)



# Descriptives 
HCP_demographics_cleaned_final = read.csv(paste(resdir_hcp, 'demographics_cleaned_final.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

dim(HCP_demographics_cleaned_final)
str(HCP_demographics_cleaned_final)


### Make separate dataframes for different sexes to run within sex ICV analyses

# FC G1
HCP_array_aligned_fc_G1_M = read.csv(paste(resdir_hcp, 'array_aligned_fc_G1_M.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_array_aligned_fc_G1_F = read.csv(paste(resdir_hcp, 'array_aligned_fc_G1_F.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

# demographics
HCP_demographics_cleaned_final_M = HCP_demographics_cleaned_final %>%
  filter(Gender == "M")

HCP_demographics_cleaned_final_F = HCP_demographics_cleaned_final %>%
  filter(Gender == "F")



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LINEAR REGRESSIONS
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### SEX effects in model = DV ~ Sex + Age + ICV + random nested effect(family relatedness/twin status)

### FUNCTIONAL -> fc

# run model
HCP_lmer_fc_G1_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final)
HCP_lmer_fc_G2_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_fc_G2, df_iv = HCP_demographics_cleaned_final)
HCP_lmer_fc_G3_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_fc_G3, df_iv = HCP_demographics_cleaned_final)

HCP_lmer_fc_G1_rescaled_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_fc_G1_rescaled, df_iv = HCP_demographics_cleaned_final)

HCP_lmer_hemi_fc_G1_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_hemi_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final)


# number of significant parcels
sum(HCP_lmer_fc_G1_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val < 0.05))
sum(HCP_lmer_fc_G2_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  
sum(HCP_lmer_fc_G3_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  

sum(HCP_lmer_fc_G1_rescaled_sex_contrast_res$q_val < 0.05, na.rm=TRUE)

sum(HCP_lmer_hemi_fc_G1_sex_contrast_res$q_val < 0.05, na.rm=TRUE)


### STRUCTURAL -> MPC

# run model
HCP_lmer_mpc_G1_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_mpc_G1, df_iv = HCP_demographics_cleaned_final)
HCP_lmer_mpc_G2_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_mpc_G2, df_iv = HCP_demographics_cleaned_final)
HCP_lmer_mpc_G3_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_mpc_G3, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_mpc_G1_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))
sum(HCP_lmer_mpc_G2_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  
sum(HCP_lmer_mpc_G3_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  



### GEODESIC DISTANCE

# run model
HCP_lmer_geo_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_mean_geodesic_distances, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_geo_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))



##### ICV effects in model = DV ~ Sex + Age + ICV + random nested effect(family relatedness/twin status)

# run model
HCP_lmer_fc_G1_icv_contrast_res = lmer.hcp_icv_contrast(df_dv = HCP_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final)
HCP_lmer_hemi_fc_G1_icv_contrast_res = lmer.hcp_icv_contrast(df_dv = HCP_hemi_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_fc_G1_icv_contrast_res$q_val < 0.05, na.rm=TRUE) 
sum(HCP_lmer_hemi_fc_G1_icv_contrast_res$q_val < 0.05, na.rm=TRUE) 


# run model effects of ICV on mean geodesic distance of top10% connections computed at the individual level
HCP_lmer_geo_icv_contrast_res = lmer.hcp_icv_contrast(df_dv = HCP_mean_geodesic_distances, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_geo_icv_contrast_res$q_val < 0.05, na.rm=TRUE) 


##### Geodesic distance effects in model = DV ~ Sex + Age + ICV + geodesic distance + random nested effect(family relatedness/twin status)

# run model
HCP_lmer_fc_G1_geo_contrast_res = lmer.hcp_geo_contrast(df_dv = HCP_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final, df_geo = HCP_mean_geodesic_distances)
HCP_lmer_hemi_fc_G1_geo_contrast_res = lmer.hcp_geo_contrast(df_dv = HCP_hemi_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final, df_geo = HCP_mean_geodesic_distances)


# number of significant parcels
sum(HCP_lmer_fc_G1_geo_contrast_res$q_val < 0.05, na.rm=TRUE) 
sum(HCP_lmer_hemi_fc_G1_geo_contrast_res$q_val < 0.05, na.rm=TRUE) 


##### MPC G1 effects in model = DV ~ Sex + Age + ICV + MPC G1 + random nested effect(family relatedness/twin status)

# run model
HCP_lmer_fc_G1_mpc_G1_contrast_res = lmer.hcp_mpc_contrast(df_dv = HCP_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final, df_mpc = HCP_array_aligned_mpc_G1)
HCP_lmer_hemi_fc_G1_mpc_G1_contrast_res = lmer.hcp_mpc_contrast(df_dv = HCP_hemi_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final, df_mpc = HCP_array_aligned_mpc_G1)


# number of significant parcels
sum(HCP_lmer_fc_G1_mpc_G1_contrast_res$q_val < 0.05, na.rm=TRUE) 
sum(HCP_lmer_hemi_fc_G1_mpc_G1_contrast_res$q_val < 0.05, na.rm=TRUE) 




##### ICV effects within sexes in model = DV ~ Age + ICV + random nested effect(family relatedness/twin status)

# run models
HCP_lmer_fc_G1_icv_contrast_within_sex_M_res = lmer.hcp_icv_contrast_within_sex(df_dv = HCP_array_aligned_fc_G1_M, df_iv = HCP_demographics_cleaned_final_M)
HCP_lmer_fc_G1_icv_contrast_within_sex_F_res = lmer.hcp_icv_contrast_within_sex(df_dv = HCP_array_aligned_fc_G1_F, df_iv = HCP_demographics_cleaned_final_F)

# number of significant parcels
sum(HCP_lmer_fc_G1_icv_contrast_within_sex_M_res$q_val < 0.05, na.rm=TRUE)  
sum(HCP_lmer_fc_G1_icv_contrast_within_sex_F_res$q_val < 0.05, na.rm=TRUE) 




##### sex effects for mean geodesic distance by seed region of top 10% connections (appearing most frequently by sex)

# run model
HCP_lmer_mean_geo_most_frequent_connections_bysexprofile_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_recomp_sub_mean_geodesic_distances_most_frequent_connections_bysexprofile, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_mean_geo_most_frequent_connections_bysexprofile_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  



##### sex effects for mean fc strength by seed for top 10% connections computed at the individual level

# run model
HCP_lmer_mean_fc_strengths_top10_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_mean_fc_strengths_top10, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_mean_fc_strengths_top10_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  





# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# EXPORT RESULTS 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### SEX effects in model = DV ~ Sex + Age + ICV + random nested effect(family relatedness/twin status)

### FUNCTIONAL -> fc

write.csv(HCP_lmer_fc_G1_sex_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G1_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_fc_G2_sex_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G2_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_fc_G3_sex_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G3_sex_contrast_res.csv', sep = ''), row.names = FALSE)

write.csv(HCP_lmer_fc_G1_rescaled_sex_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G1_rescaled_sex_contrast_res.csv', sep = ''), row.names = FALSE)

write.csv(HCP_lmer_hemi_fc_G1_sex_contrast_res, paste(resdir_hcp, 'R_lmer_hemi_fc_G1_sex_contrast_res.csv', sep = ''), row.names = FALSE)


write.csv(HCP_lmer_fc_G1_icv_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G1_icv_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_hemi_fc_G1_icv_contrast_res, paste(resdir_hcp, 'R_lmer_hemi_fc_G1_icv_contrast_res.csv', sep = ''), row.names = FALSE)

write.csv(HCP_lmer_geo_icv_contrast_res, paste(resdir_hcp, 'R_lmer_geo_icv_contrast_res.csv', sep = ''), row.names = FALSE)


write.csv(HCP_lmer_fc_G1_geo_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G1_geo_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_hemi_fc_G1_geo_contrast_res, paste(resdir_hcp, 'R_lmer_hemi_fc_G1_geo_contrast_res.csv', sep = ''), row.names = FALSE)

write.csv(HCP_lmer_fc_G1_mpc_G1_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G1_mpc_G1_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_hemi_fc_G1_mpc_G1_contrast_res, paste(resdir_hcp, 'R_lmer_hemi_fc_G1_mpc_G1_contrast_res.csv', sep = ''), row.names = FALSE)


### STRUCTURAL -> MPC
write.csv(HCP_lmer_mpc_G1_sex_contrast_res, paste(resdir_hcp, 'R_lmer_mpc_G1_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_mpc_G2_sex_contrast_res, paste(resdir_hcp, 'R_lmer_mpc_G2_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_mpc_G3_sex_contrast_res, paste(resdir_hcp, 'R_lmer_mpc_G3_sex_contrast_res.csv', sep = ''), row.names = FALSE)


### GEODESIC DISTANCE
write.csv(HCP_lmer_geo_sex_contrast_res, paste(resdir_hcp, 'R_lmer_geo_sex_contrast_res.csv', sep = ''), row.names = FALSE)



##### sex effects for mean geodesic distance by seed region of top 10% connections (appearing most frequently by sex)
write.csv(HCP_lmer_mean_geo_most_frequent_connections_bysexprofile_sex_contrast_res, paste(resdir_hcp, 'R_lmer_mean_geo_most_frequent_connections_bysexprofile_sex_contrast_res.csv', sep = ''), row.names = FALSE)

##### sex effects for mean fc strength by seed for top 10% connections computed at the individual level
write.csv(HCP_lmer_mean_fc_strengths_top10_sex_contrast_res, paste(resdir_hcp, 'R_lmer_mean_fc_strengths_top10_sex_contrast_res.csv', sep = ''), row.names = FALSE)




############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################

# extra code

############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DESCRIPTIVES 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Function displaying discriptive statistics
summarise.descriptives <- function(df, variable){
  
  '
- summarizes main descriptives: mean, sd, min, max, sem, IQR
- to supply: dataframe and variable of interest
- was not able to make summarise by Sex work
'
  
  summarise(df,
            mean = mean(variable),
            sd = sd(variable),
            min = min(variable),
            max = max(variable),
            sem = sd(variable)/sqrt(length(variable)),
            IQR = IQR(variable))
}




str(demographics_cleaned, list.len=ncol(demographics_cleaned))  # untruncated output
dim(demographics_cleaned)


# Sample size by sex
ftable(demographics_cleaned$Sex)


# Age
summarise.descriptives(demographics_cleaned, demographics_cleaned$Age_Bin)

# Age by sex
ddply(demographics_cleaned, 'Sex', summarise,
      mean = mean(Age_Bin),
      sd = sd(Age_Bin),
      min = min(Age_Bin),
      max = max(Age_Bin),
      sem = sd(Age_Bin)/sqrt(length(Age_Bin)),
      IQR = IQR(Age_Bin))


# ICV
summarise.descriptives(demographics_cleaned, demographics_cleaned$ICV)

# ICV by sex
ddply(demographics_cleaned, 'Sex', summarise,
      mean = mean(ICV),
      sd = sd(ICV),
      min = min(ICV),
      max = max(ICV),
      sem = sd(ICV)/sqrt(length(ICV)),
      IQR = IQR(ICV))




