# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Project: Sex differences in brain organization

# Main script

# Content: Linear Model analyses for sex contrast in gradient eigenvalues in GSP, Linear Mixed Effects Model for sex contrast in gradient eigenvalues controlling for family relatedness and twin status

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


### Linear Regression HCP Connectivity lm = Gradient_Eigenvalues ~ Sex + Age + ICV + random nested effect(family relatedness/twin status)

lmer.hcp_sex_contrast <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for SEX effects, including sex, age and ICV, as well as random effects family id, twin status and family id * twin status, in the model as covariates (relevant to functional connectivity)
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
    
    # Fit a linear mixed effects model: DV ~ Sex + Age + ICV + random nested effect(family relatedness/twin status)
    
    # Model including both "single" random effects included in the interaction random effect (Sofie's decision)
    lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol + (1 | family_id/twin_status), REML = FALSE)
    
    # Extract from summary of lmer_fit the t- and p-values
    # summary(lmer_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 Age, 4 ICV; columns = 1 Estimate, 2 Std. Error, 3 df, 4 t-value, 5 p-value
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






# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PREPARE DATA 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Aligned functional gradient loadings 
HCP_array_aligned_fc_G1 = read.csv(paste(resdir_hcp, 'array_aligned_fc_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_array_aligned_fc_G2 = read.csv(paste(resdir_hcp, 'array_aligned_fc_G2.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_array_aligned_fc_G3 = read.csv(paste(resdir_hcp, 'array_aligned_fc_G3.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

# Aligned MPC gradient loadings
HCP_array_aligned_mpc_G1 = read.csv(paste(resdir_hcp, 'array_aligned_mpc_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_array_aligned_mpc_G2 = read.csv(paste(resdir_hcp, 'array_aligned_mpc_G2.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_array_aligned_mpc_G3 = read.csv(paste(resdir_hcp, 'array_aligned_mpc_G3.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

# Geodesic distances (mean of top 10% of functional connections)
HCP_mean_geodesic_distances = read.csv(paste(resdir_hcp, 'mean_geodesic_distances.csv', sep = ''), fileEncoding = 'UTF-8-BOM')


# class(array_aligned_G1)
# typeof(array_aligned_G1)
dim(HCP_array_aligned_fc_G1)
dim(HCP_array_aligned_mpc_G1)
dim(HCP_mean_geodesic_distances)



# Descriptives 
HCP_demographics_cleaned_final = read.csv(paste(resdir_hcp, 'HCP_demographics_cleaned_final.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

dim(HCP_demographics_cleaned_final)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LINEAR REGRESSIONS
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### SEX effects in model = DV ~ Sex + Age + ICV + random nested effect(family relatedness/twin status)

### FUNCTIONAL -> fc

# run model
HCP_lmer_fc_G1_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_fc_G1, df_iv = HCP_demographics_cleaned_final)
HCP_lmer_fc_G2_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_fc_G2, df_iv = HCP_demographics_cleaned_final)
HCP_lmer_fc_G3_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_fc_G3, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_fc_G1_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))
sum(HCP_lmer_fc_G2_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  
sum(HCP_lmer_fc_G3_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  



### STRUCTURAL -> MPC

# run model
HCP_lmer_mpc_G1_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_mpc_G1, df_iv = HCP_demographics_cleaned_final)
HCP_lmer_mpc_G2_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_mpc_G2, df_iv = HCP_demographics_cleaned_final)
HCP_lmer_mpc_G3_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_mpc_G3, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_mpc_G1_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))
sum(HCP_lmer_mpc_G2_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  
sum(HCP_lmer_mpc_G3_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  



### GEODESIC DISTANCE

# run model
HCP_lmer_geo_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_mean_geodesic_distances, df_iv = HCP_demographics_cleaned_final)

# number of significant parcels
sum(HCP_lmer_geo_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))




# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# EXPORT RESULTS 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### SEX effects in model = DV ~ Sex + Age + ICV + random nested effect(family relatedness/twin status)

### FUNCTIONAL -> fc

write.csv(HCP_lmer_fc_G1_sex_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G1_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_fc_G2_sex_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G2_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_fc_G3_sex_contrast_res, paste(resdir_hcp, 'R_lmer_fc_G3_sex_contrast_res.csv', sep = ''), row.names = FALSE)


### STRUCTURAL -> MPC
write.csv(HCP_lmer_mpc_G1_sex_contrast_res, paste(resdir_hcp, 'R_lmer_mpc_G1_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_mpc_G2_sex_contrast_res, paste(resdir_hcp, 'R_lmer_mpc_G2_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_mpc_G3_sex_contrast_res, paste(resdir_hcp, 'R_lmer_mpc_G3_sex_contrast_res.csv', sep = ''), row.names = FALSE)


### GEODESIC DISTANCE
write.csv(HCP_lmer_geo_sex_contrast_res, paste(resdir_hcp, 'R_lmer_geo_sex_contrast_res.csv', sep = ''), row.names = FALSE)






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




