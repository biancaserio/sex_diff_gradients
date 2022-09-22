---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Project: Sex differences in brain organization
  
# Content: Linear regression analyses for sex contrast in gradient eigenvalues in GSP
  
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
##require(ggplot2) #this should be included in tidyverse
#require(qwraps2)
# options(qwraps2_markup = "markdown") # define the markup language we are working in
# require(reshape2)
# require(psych)
# require(moments) # for skewness() kurtosis()
# require(ggfortify) # to plot PCA
# require(PCAmixdata)
# require(factoextra) # visualisation of PCA (fviz_pca_var)
# require(RColorBrewer)
# require(openxlsx) # to create workbook (excel)
# require(nortest) # for Anderson Darlin test (of normality) ad.test() (https://www.itl.nist.gov/div898/handbook/eda/section3/eda35e.htm H0: normal distribution)
# require(effectsize) # for rank_biserial() (effect size Mann-Whitney U)
# require(REdaS) # for Bartlett's test of sphericity (bart_spher) and Kaiser-Meyer-Olkin Measure of Sampling Adequacy (KMOS)
# require(car) # Levene's test
# require(rstatix) # Cohen's d (cohens_d) https://www.datanovia.com/en/lessons/t-test-effect-size-using-cohens-d-measure/#:~:text=To%20calculate%20an%20effect%20size,)%20%3D%20sd(x)%20.&text=%CE%BC%20is%20the%20theoretical%20mean,value%20is%20mu%20%3D%200).
# require(ppcor) # Partial correlation (pcor.test)
# require(chron) # to convert character variables in the format of "XX:XX:XX" to times variables (sing chron::times)
# # display.brewer.all() # displays all color schemes
# 
# require(lme4)  # for lmer (linear mixed effects model)
# require(lmerTest)  # to obtain p-values for lmer - this actually overrides lme4's lmer() and prints the p-values for the fixed effects (which aren't present in lme4's lmer())
#                       # note that no p-values are printed for the random effects because p values cannot be estimated for random effects because these are latent variables without standard deviations
# require(emmeans)  # for Estimated Marginal Means, aka Least-Squares Means


# Clear environment
rm(list = ls())

# set up directories
codedir = dirname(getActiveDocumentContext()$path)  # get path to current script
datadir = '/data/p_02667/sex_diff_gradients/data/'
resdir = '/data/p_02667/sex_diff_gradients/results/GSP/'

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




### Linear Regression

lm.sex_age_icv <- function(df_dv, df_iv) {
  
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





----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PREPARE DATA 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Aligned gradient values 
array_aligned_G1 = read.csv(paste(resdir, 'array_aligned_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
array_aligned_G2 = read.csv(paste(resdir, 'array_aligned_G2.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
array_aligned_G3 = read.csv(paste(resdir, 'array_aligned_G3.csv', sep = ''), fileEncoding = 'UTF-8-BOM')


# class(array_aligned_G1)
# typeof(array_aligned_G1)
dim(array_aligned_G1)
length(array_aligned_G3$X1)

# Descriptives 
demographics_cleaned = read.csv(paste(resdir, 'demographics_cleaned.csv', sep = ''), fileEncoding = 'UTF-8-BOM')


----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DESCRIPTIVES 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
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



----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LINEAR REGRESSION 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
### model = Gradient_Eigenvalues ~ Sex + Age + ICV 

# run model
lm_G1_sex_age_icv_res = lm.sex_age_icv(df_dv = array_aligned_G1, df_iv = demographics_cleaned)
lm_G2_sex_age_icv_res = lm.sex_age_icv(df_dv = array_aligned_G2, df_iv = demographics_cleaned)
lm_G3_sex_age_icv_res = lm.sex_age_icv(df_dv = array_aligned_G3, df_iv = demographics_cleaned)


# number of significant parcels
sum(lm_G1_sex_age_icv_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))
sum(lm_G2_sex_age_icv_res$q_val_sex < 0.05, na.rm=TRUE)  
sum(lm_G3_sex_age_icv_res$q_val_sex < 0.05, na.rm=TRUE)  





----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# EXPORT RESULTS 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

write.csv(lm_G1_sex_age_icv_res, paste(resdir, 'R_lm_G1_sex_age_icv_res.csv', sep = ''), row.names = FALSE)
write.csv(lm_G2_sex_age_icv_res, paste(resdir, 'R_lm_G2_sex_age_icv_res.csv', sep = ''), row.names = FALSE)
write.csv(lm_G3_sex_age_icv_res, paste(resdir, 'R_lm_G3_sex_age_icv_res.csv', sep = ''), row.names = FALSE)








----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# EXPLORATORY:  Gradients aligned to HCP mean gradient (instead of same data (GSP) mean gradient)
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
### PREPARE THE DATA
  
# HCP-Aligned gradient values (standard -> aligned to mean HCP gradient)
array_aligned_hcp_G1 = read.csv(paste(resdir, 'exploratory/array_aligned_hcp_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
array_aligned_hcp_G2 = read.csv(paste(resdir, 'exploratory/array_aligned_hcp_G2.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
array_aligned_hcp_G3 = read.csv(paste(resdir, 'exploratory/array_aligned_hcp_G3.csv', sep = ''), fileEncoding = 'UTF-8-BOM')



### LINEAR REGRESSION

## model = Gradient_Eigenvalues ~ Sex + Age + ICV 

# run model
lm_G1_sex_age_icv_res_hcp_aligned = lm.sex_age_icv(df_dv = array_aligned_hcp_G1, df_iv = demographics_cleaned)
lm_G2_sex_age_icv_res_hcp_aligned = lm.sex_age_icv(df_dv = array_aligned_hcp_G2, df_iv = demographics_cleaned)
lm_G3_sex_age_icv_res_hcp_aligned = lm.sex_age_icv(df_dv = array_aligned_hcp_G3, df_iv = demographics_cleaned)

# number of significant parcels
sum(lm_G1_sex_age_icv_res_hcp_aligned$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))
sum(lm_G2_sex_age_icv_res_hcp_aligned$q_val_sex < 0.05, na.rm=TRUE)  
sum(lm_G3_sex_age_icv_res_hcp_aligned$q_val_sex < 0.05, na.rm=TRUE)  



### EXPORT RESULTS

write.csv(lm_G1_sex_age_icv_res_hcp_aligned, paste(resdir, 'exploratory/R_lm_G1_sex_age_icv_res_hcp_aligned.csv', sep = ''), row.names = FALSE)
write.csv(lm_G2_sex_age_icv_res_hcp_aligned, paste(resdir, 'exploratory/R_lm_G2_sex_age_icv_res_hcp_aligned.csv', sep = ''), row.names = FALSE)
write.csv(lm_G3_sex_age_icv_res_hcp_aligned, paste(resdir, 'exploratory/R_lm_G3_sex_age_icv_res_hcp_aligned.csv', sep = ''), row.names = FALSE)




