# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Project: Sex differences in brain organization

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
resdir_gsp = '/data/p_02667/sex_diff_gradients/results/GSP/'
resdir_hcp = '/data/p_02667/sex_diff_gradients/results/HCP/'

# set directory to path of current script
setwd(codedir) 



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DEFINING FUNCTIONS 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# note: these functions are hardcoded for given datasets (variable names, variables included in the regression)

### Linear Regression GSP

lm.gsp_sex_contrast <- function(df_dv, df_iv) {

'
  - fits and runs linear model to test for SEX effects, including sex, age and ICV in the model as covariates
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




lmer.hcp_sex_contrast <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for SEX effects, including sex, age and ICV, as well as random effects family id and family id * twin status, in the model as covariates
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
    
    # Fit a linear mixed effects model: lm = Gradient_Eigenvalues ~ Sex + Age + ICV + random effect(family relatedness) + random effect (family relatedness * twin status)
    # error message: https://stackoverflow.com/questions/60028673/lme4-error-boundary-singular-fit-see-issingular -> Your model did fit, but it generated that warning because your random effects are very small
    # https://stats.stackexchange.com/questions/96600/interactions-between-random-effects
    # https://stackoverflow.com/questions/71340764/interaction-between-two-factors-as-random-effects-in-mixed-model-in-r
    
    lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol + (1 | family_id) + (1 | family_id:twin_status), REML = FALSE)
    #lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol + (1 | family_id) + (1 | twin_status) + (1 | family_id:twin_status), REML = FALSE)
    
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

# Aligned gradient values 
GSP_array_aligned_G1 = read.csv(paste(resdir_gsp, 'array_aligned_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
GSP_array_aligned_G2 = read.csv(paste(resdir_gsp, 'array_aligned_G2.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
GSP_array_aligned_G3 = read.csv(paste(resdir_gsp, 'array_aligned_G3.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

HCP_array_aligned_G1 = read.csv(paste(resdir_hcp, 'array_aligned_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_array_aligned_G2 = read.csv(paste(resdir_hcp, 'array_aligned_G2.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_array_aligned_G3 = read.csv(paste(resdir_hcp, 'array_aligned_G3.csv', sep = ''), fileEncoding = 'UTF-8-BOM')


# class(array_aligned_G1)
# typeof(array_aligned_G1)
dim(GSP_array_aligned_G1)
dim(HCP_array_aligned_G1)


# Descriptives 
GSP_demographics_cleaned = read.csv(paste(resdir_gsp, 'demographics_cleaned.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
HCP_demographics_cleaned = read.csv(paste(resdir_hcp, 'demographics_cleaned.csv', sep = ''), fileEncoding = 'UTF-8-BOM')


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LINEAR REGRESSIONS
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### GSP: model = Gradient_Eigenvalues ~ Sex + Age + ICV 

# run model
GSP_lm_G1_sex_contrast_res = lm.gsp_sex_contrast(df_dv = GSP_array_aligned_G1, df_iv = GSP_demographics_cleaned)
GSP_lm_G2_sex_contrast_res = lm.gsp_sex_contrast(df_dv = GSP_array_aligned_G2, df_iv = GSP_demographics_cleaned)
GSP_lm_G3_sex_contrast_res = lm.gsp_sex_contrast(df_dv = GSP_array_aligned_G3, df_iv = GSP_demographics_cleaned)


# number of significant parcels
sum(GSP_lm_G1_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(GSP_lm_G1_res$q_val_sex < 0.05))
sum(GSP_lm_G2_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  
sum(GSP_lm_G3_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  



### HCP: model = Gradient_Eigenvalues ~ Sex + Age + ICV + random effect family relatedness + random effect (family relatedness * twin status) 

# run model
HCP_lmer_G1_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_G1, df_iv = HCP_demographics_cleaned)
HCP_lmer_G2_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_G2, df_iv = HCP_demographics_cleaned)
HCP_lmer_G3_sex_contrast_res = lmer.hcp_sex_contrast(df_dv = HCP_array_aligned_G3, df_iv = HCP_demographics_cleaned)

# number of significant parcels
sum(HCP_lmer_G1_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))
sum(HCP_lmer_G2_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  
sum(HCP_lmer_G3_sex_contrast_res$q_val_sex < 0.05, na.rm=TRUE)  



# tests below to delete!!!!
family_id = HCP_demographics_cleaned$Family_ID
twin_status = HCP_demographics_cleaned$TwinStatus

test = lmer(HCP_array_aligned_G1[[1]] ~ HCP_demographics_cleaned$Gender + HCP_demographics_cleaned$Age_in_Yrs + HCP_demographics_cleaned$FS_IntraCranial_Vol + (1 | family_id) + (1 | family_id:twin_status), REML = FALSE)
summary(test)$coefficients

test = lmer(HCP_array_aligned_G1[[1]] ~ HCP_demographics_cleaned$Gender + HCP_demographics_cleaned$Age_in_Yrs + HCP_demographics_cleaned$FS_IntraCranial_Vol + (1 | family_id) + (1 | twin_status) + (1 | family_id:twin_status), REML = FALSE)
summary(test)$coefficients


----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# EXPORT RESULTS 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

write.csv(GSP_lm_G1_sex_contrast_res, paste(resdir_gsp, 'R_lm_G1_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(GSP_lm_G2_sex_contrast_res, paste(resdir_gsp, 'R_lm_G2_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(GSP_lm_G3_sex_contrast_res, paste(resdir_gsp, 'R_lm_G3_sex_contrast_res.csv', sep = ''), row.names = FALSE)


write.csv(HCP_lmer_G1_sex_contrast_res, paste(resdir_hcp, 'R_lmer_G1_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_G2_sex_contrast_res, paste(resdir_hcp, 'R_lmer_G2_sex_contrast_res.csv', sep = ''), row.names = FALSE)
write.csv(HCP_lmer_G3_sex_contrast_res, paste(resdir_hcp, 'R_lmer_G3_sex_contrast_res.csv', sep = ''), row.names = FALSE)








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




