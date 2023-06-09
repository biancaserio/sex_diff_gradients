# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Project: Sex differences in brain organization

# Main script

# Content: Linear Mixed Effects Model for sex contrast in gradient eigenvalues controlling for sex, age, total *SA*, nested random effect of family relatedness and twin status

# Analyses conducted and outputted in SA subfolders (see set up directories)

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
resdir_hcp = '/data/p_02667/sex_diff_gradients/results/HCP/SA/'

# set directory to path of current script
setwd(codedir) 



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DEFINING FUNCTIONS 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# note: these functions are hardcoded for given datasets (variable names, variables included in the regression)


### Linear Regression HCP Connectivity lmer = Gradient_Eigenvalues ~ Sex + Age + tot_SA + random nested effect(family relatedness/twin status) 

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




### Linear Regression HCP Connectivity lmer = Gradient_Eigenvalues ~ Sex + Age + tot_SA + MPC G1 + random nested effect(family relatedness/twin status) 

# MPC G1 contrast
lmer.hcp_mpc_contrast <- function(df_dv, df_iv, df_mpc) {
  
  '
    - fits and runs linear model to test for MPC G1 effects, including sex, age, MPC G1 and tot_SA, as well as random nested effect(family relatedness/twin status in the model as covariates
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables), df_mpc (dataframe containing MPC G1)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for MPC G1 contrast
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
    lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$tot_SA + df_mpc[[i]] + (1 | family_id/twin_status), REML = FALSE)
    
    # Extract from summary of lmer_fit the t- and p-values
    # summary(lmer_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 age, 4 tot_SA, 5 MPC G1; columns = 1 Estimate, 2 Std. Error, 3 df, 4 t-value, 5 p-value
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




### Linear Regression HCP Connectivity lmer = Gradient_Eigenvalues ~ Sex + Age + tot_SA + geodesic distance + random nested effect(family relatedness/twin status) 

# Geodesic distance contrast
lmer.hcp_geo_contrast <- function(df_dv, df_iv, df_geo) {

'
    - fits and runs linear model to test for GEODESIC DISTANCE effects, including sex, age, geodesic distance and tot_SA, as well as random nested effect(family relatedness/twin status in the model as covariates
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables), df_geo (dataframe containing geodesic distance)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for GEODESIC DISTANCE contrast
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
  lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$tot_SA + df_geo[[i]] + (1 | family_id/twin_status), REML = FALSE)
  
  # Extract from summary of lmer_fit the t- and p-values
  # summary(lmer_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 age, 4 tot_SA, 5 geodesic distance; columns = 1 Estimate, 2 Std. Error, 3 df, 4 t-value, 5 p-value
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





# total surface area contrast in lmer = Gradient_Eigenvalues ~ Sex + Age + tot_SA + random nested effect(family relatedness/twin status) 


lmer.hcp_totSA_contrast <- function(df_dv, df_iv) {

'
    - fits and runs linear model to test for total SA effects, including sex, age and tot SA, as well as random nested effect(family relatedness/twin status in the model as covariates
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for TOT SA contrast
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
  # summary(lmer_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 age, 4 tot_SA; columns = 1 Estimate, 2 Std. Error, 3 df, 4 t-value, 5 p-value
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


















# ICV contrast
#lmer.hcp_icv_contrast <- function(df_dv, df_iv) {
  
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

# Descriptives 
HCP_demographics_cleaned_final = read.csv(paste(resdir_hcp, 'demographics_cleaned_final.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

dim(HCP_demographics_cleaned_final)
str(HCP_demographics_cleaned_final)


# Aligned functional gradient loadings (wholebrain and hemi)
HCP_hemi_array_aligned_fc_G1 = read.csv(paste(resdir_hcp, 'hemi_array_aligned_fc_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

HCP_array_aligned_fc_G1 = read.csv(paste(resdir_hcp, 'array_aligned_fc_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')


# Aligned MPC gradient loadings
HCP_array_aligned_mpc_G1 = read.csv(paste(resdir_hcp, 'array_aligned_mpc_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')


# Geodesic distances (mean of top 10% of functional connections)
HCP_mean_geodesic_distances = read.csv(paste(resdir_hcp, 'mean_geodesic_distances.csv', sep = ''), fileEncoding = 'UTF-8-BOM')



#HCP_array_aligned_fc_G1_rescaled = read.csv(paste(resdir_hcp, 'array_aligned_fc_G1_rescaled.csv', sep = ''), fileEncoding = 'UTF-8-BOM')  # rescaled between 0-1


# Mean fc strength by seed for top 10% connections computed at the individual level
HCP_mean_fc_strengths_top10 = read.csv(paste(resdir_hcp, 'mean_fc_strengths_top10.csv', sep = ''), fileEncoding = 'UTF-8-BOM')


# Local SA
local_SA = read.csv(paste(resdir_hcp, 'local_SA.csv', sep = ''), fileEncoding = 'UTF-8-BOM')








# Mean geodesic distances of mean connectivity profiles by sex (ie each subject has the mean (by seed region) of the geodesic distance of the most frequent top 10% connections for their sex)
HCP_recomp_sub_mean_geodesic_distances_most_frequent_connections_bysexprofile = read.csv(paste(resdir_hcp, 'recomp_sub_mean_geodesic_distances_most_frequent_connections_bysexprofile.csv', sep = ''), fileEncoding = 'UTF-8-BOM')



# class(array_aligned_G1)
# typeof(array_aligned_G1)
dim(HCP_array_aligned_fc_G1)
dim(HCP_hemi_array_aligned_fc_G1)
dim(HCP_array_aligned_mpc_G1)
dim(HCP_mean_geodesic_distances)







# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LINEAR REGRESSIONS
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### SEX effects in model = DV ~ Sex + Age + tot_SA + random nested effect(family relatedness/twin status)

### FUNCTIONAL -> fc

# run models (with hemi and wholebrain computed gradients)
HCP_lmer_hemi_fc_G1_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_hemi_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final)
HCP_lmer_fc_G1_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final)


#HCP_lmer_fc_G1_rescaled_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_fc_G1_rescaled, df_iv = HCP_demographics_cleaned_final)


# number of significant parcels
sum(HCP_lmer_hemi_fc_G1_sex_contrast_res$q_val < 0.05, na.rm=TRUE)
sum(HCP_lmer_fc_G1_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val < 0.05))


#sum(HCP_lmer_fc_G1_rescaled_sex_contrast_res$q_val < 0.05, na.rm=TRUE)




### STRUCTURAL -> MPC

# run model
HCP_lmer_mpc_G1_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_mpc_G1, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_mpc_G1_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))




### GEODESIC DISTANCE

# run model
HCP_lmer_geo_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_mean_geodesic_distances, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_geo_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))




##### MPC G1 effects in model = DV ~ Sex + Age + tot_SA + MPC G1 + random nested effect(family relatedness/twin status)

# run model
HCP_lmer_fc_G1_mpc_G1_contrast_res = lmer.hcp_mpc_contrast(df_dv = HCP_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final, df_mpc = HCP_array_aligned_mpc_G1)
HCP_lmer_hemi_fc_G1_mpc_G1_contrast_res = lmer.hcp_mpc_contrast(df_dv = HCP_hemi_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final, df_mpc = HCP_array_aligned_mpc_G1)

# number of significant parcels
sum(HCP_lmer_fc_G1_mpc_G1_contrast_res$q_val < 0.05, na.rm=TRUE) 
sum(HCP_lmer_hemi_fc_G1_mpc_G1_contrast_res$q_val < 0.05, na.rm=TRUE) 




##### Geodesic distance effects in model = DV ~ Sex + Age + tot_SA + geodesic distance + random nested effect(family relatedness/twin status)

# run model
HCP_lmer_fc_G1_geo_contrast_res = lmer.hcp_geo_contrast(df_dv = HCP_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final, df_geo = HCP_mean_geodesic_distances)
HCP_lmer_hemi_fc_G1_geo_contrast_res = lmer.hcp_geo_contrast(df_dv = HCP_hemi_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final, df_geo = HCP_mean_geodesic_distances)

# number of significant parcels
sum(HCP_lmer_fc_G1_geo_contrast_res$q_val < 0.05, na.rm=TRUE) 
sum(HCP_lmer_hemi_fc_G1_geo_contrast_res$q_val < 0.05, na.rm=TRUE) 



##### Total SA contrast in model = DV ~ Sex + Age + tot_SA + geodesic distance + random nested effect(family relatedness/twin status)

# run model
HCP_lmer_hemi_fc_G1_totSA_contrast_res = lmer.hcp_totSA_contrast(df_dv = HCP_hemi_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final)
HCP_lmer_fc_G1_totSA_contrast_res = lmer.hcp_totSA_contrast(df_dv = HCP_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_hemi_fc_G1_totSA_contrast_res$q_val < 0.05, na.rm=TRUE) 
sum(HCP_lmer_fc_G1_totSA_contrast_res$q_val < 0.05, na.rm=TRUE) 


##### sex effects for mean fc strength by seed for top 10% connections computed at the individual level

# run model
HCP_lmer_mean_fc_strengths_top10_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_mean_fc_strengths_top10, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_mean_fc_strengths_top10_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  


##### sex effects in local SA

# run model
HCP_lmer_local_SA_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = local_SA, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_local_SA_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  





##### Local SA contrast in lmer = Gradient_Eigenvalues ~ Sex + Age + tot_SA + local_SA + random nested effect(family relatedness/twin status) 

# run model
HCP_lmer_hemi_fc_G1_local_SA_contrast_res = lmer.hcp_local_SA_contrast(df_dv = HCP_hemi_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final, df_local_SA = local_SA)
HCP_lmer_fc_G1_local_SA_contrast_res = lmer.hcp_local_SA_contrast(df_dv = HCP_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final, df_local_SA = local_SA)

# number of significant parcels
sum(HCP_lmer_hemi_fc_G1_local_SA_contrast_res$q_val < 0.05, na.rm=TRUE)  
sum(HCP_lmer_fc_G1_local_SA_contrast_res$q_val < 0.05, na.rm=TRUE)  














############ not used, from p1_main.R -> might need at some point? ##############


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





##### sex effects for mean geodesic distance by seed region of top 10% connections (appearing most frequently by sex)

# run model
HCP_lmer_mean_geo_most_frequent_connections_bysexprofile_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_recomp_sub_mean_geodesic_distances_most_frequent_connections_bysexprofile, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_mean_geo_most_frequent_connections_bysexprofile_sex_contrast_res$q_val < 0.05, na.rm=TRUE)  








# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# EXPORT RESULTS 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### SEX effects in model = DV ~ Sex + Age + tot_SA + random nested effect(family relatedness/twin status)

### FUNCTIONAL -> fc
write.csv(HCP_lmer_hemi_fc_G1_sex_contrast_res, paste(resdir_hcp, 'R_lmer_hemi_fc_G1_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_fc_G1_sex_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G1_sex_contrast_res.csv', sep = ''), row.names = FALSE)

#write.csv(HCP_lmer_fc_G1_rescaled_sex_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G1_rescaled_sex_contrast_res.csv', sep = ''), row.names = FALSE)



### STRUCTURAL -> MPC
write.csv(HCP_lmer_mpc_G1_sex_contrast_res, paste(resdir_hcp, 'R_lmer_mpc_G1_sex_contrast_res.csv', sep = ''), row.names = FALSE)



### GEODESIC DISTANCE
write.csv(HCP_lmer_geo_sex_contrast_res, paste(resdir_hcp, 'R_lmer_geo_sex_contrast_res.csv', sep = ''), row.names = FALSE)



# MPC G1 effects in model = DV ~ Sex + Age + tot_SA + MPC G1 + random nested effect(family relatedness/twin status)
write.csv(HCP_lmer_fc_G1_mpc_G1_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G1_mpc_G1_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_hemi_fc_G1_mpc_G1_contrast_res, paste(resdir_hcp, 'R_lmer_hemi_fc_G1_mpc_G1_contrast_res.csv', sep = ''), row.names = FALSE)


# Geodesic distance effects in model = DV ~ Sex + Age + tot_SA + geodesic distance + random nested effect(family relatedness/twin status)
write.csv(HCP_lmer_fc_G1_geo_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G1_geo_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_hemi_fc_G1_geo_contrast_res, paste(resdir_hcp, 'R_lmer_hemi_fc_G1_geo_contrast_res.csv', sep = ''), row.names = FALSE)


# Total SA contrast in model = DV ~ Sex + Age + tot_SA + geodesic distance + random nested effect(family relatedness/twin status)
write.csv(HCP_lmer_hemi_fc_G1_totSA_contrast_res, paste(resdir_hcp, 'R_lmer_hemi_fc_G1_totSA_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_fc_G1_totSA_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G1_totSA_contrast_res.csv', sep = ''), row.names = FALSE)


# sex effects for mean fc strength by seed for top 10% connections computed at the individual level
write.csv(HCP_lmer_mean_fc_strengths_top10_sex_contrast_res, paste(resdir_hcp, 'R_lmer_mean_fc_strengths_top10_sex_contrast_res.csv', sep = ''), row.names = FALSE)


# sex effects in local SA
write.csv(HCP_lmer_local_SA_sex_contrast_res, paste(resdir_hcp, 'R_lmer_local_SA_sex_contrast_res.csv', sep = ''), row.names = FALSE)


# Local SA contrast in lmer = Gradient_Eigenvalues ~ Sex + Age + tot_SA + local_SA + random nested effect(family relatedness/twin status) 
write.csv(HCP_lmer_hemi_fc_G1_local_SA_contrast_res, paste(resdir_hcp, 'R_lmer_hemi_fc_G1_local_SA_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_fc_G1_local_SA_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G1_local_SA_contrast_res.csv', sep = ''), row.names = FALSE)








############ not used, from p1_main.R -> might need at some point? ##############






write.csv(HCP_lmer_fc_G1_icv_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G1_icv_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_hemi_fc_G1_icv_contrast_res, paste(resdir_hcp, 'R_lmer_hemi_fc_G1_icv_contrast_res.csv', sep = ''), row.names = FALSE)

write.csv(HCP_lmer_geo_icv_contrast_res, paste(resdir_hcp, 'R_lmer_geo_icv_contrast_res.csv', sep = ''), row.names = FALSE)





##### sex effects for mean geodesic distance by seed region of top 10% connections (appearing most frequently by sex)
write.csv(HCP_lmer_mean_geo_most_frequent_connections_bysexprofile_sex_contrast_res, paste(resdir_hcp, 'R_lmer_mean_geo_most_frequent_connections_bysexprofile_sex_contrast_res.csv', sep = ''), row.names = FALSE)


