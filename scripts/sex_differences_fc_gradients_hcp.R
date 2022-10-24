---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Project: Sex differences in brain organization
  
# Content: Linear regression analyses for sex contrast in gradient eigenvalues in HCP S1200
  
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
require(lme4)  # for lmer (linear mixed effects model)
require(lmerTest)  # to obtain p-values for lmer - this actually overrides lme4's lmer() and prints the p-values for the fixed effects (which aren't present in lme4's lmer())
        # note that no p-values are printed for the random effects because p values cannot be estimated for random effects because these are latent variables without standard deviations

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


# require(emmeans)  # for Estimated Marginal Means, aka Least-Squares Means


# Clear environment
rm(list = ls())

# set up directories
codedir = dirname(getActiveDocumentContext()$path)  # get path to current script
datadir = '/data/p_02667/sex_diff_gradients/data/'
resdir = '/data/p_02667/sex_diff_gradients/results/HCP/'

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
  beta_val_sex = vector(mode = "double", length = ncol(df_dv))
  
  # Degrees of Freedom = N (subjects) - number of IV (3 *** HARD-CODED FOR THIS SPECIFIC LM MODEL ***) - 1 (mean)
  #DoF = nrow(df_dv) - 3 -1
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    # Fit a linear model: lm = Gradient_Eigenvalues ~ Sex + Age + ICV
    lm_fit = lm(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol)
    
    # Extract from summary of lm_fit the t- and p-values
    # summary(lm_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 Age, 4 ICV; columns = 1 Estimate, 2 Std. Error, 3 t-value, 4 p-value
    t_val_sex[[i]] = summary(lm_fit)$coefficients[2,3]
    p_val_sex[[i]] = summary(lm_fit)$coefficients[2,4]  # if want to calculate p value by hand: p_val_sex[[i]] = 2*pt(abs(t_val_sex[[i]]), DoF, lower.tail = F)
    beta_val_sex[[i]] = summary(lm_fit)$coefficients[2,1]
    
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val_sex = p.adjust(p_val_sex, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val_sex, p_val_sex, q_val_sex, beta_val_sex)
  
  return(output_df)
}




lmer.sex_age_icv_ctrl_related <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for SEX effects, including sex, age and ICV in the model as covariates
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
    
    # Fit a linear mixed effects model: lm = Gradient_Eigenvalues ~ Sex + Age + ICV + controlling for family relatedness 
    lmer_fit <- lmer(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol + (1 | family_id), REML = FALSE)  
    #ProcSpeed_raw ~ age_months + sex_female + (1 | Pair) + (1 | M), REML = FALSE, data = ses01)
    
    # controling for twin pairs
    #lmer_fit <- lmer(ProcSpeed_raw ~ age_months + sex_female + (1 | Pair) + (1 | M), REML = FALSE, data = ses01)
    
    
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






----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PREPARE DATA 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###### Full Sample
  
# Aligned gradient values 
array_aligned_G1 = read.csv(paste(resdir, 'array_aligned_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
array_aligned_G2 = read.csv(paste(resdir, 'array_aligned_G2.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
array_aligned_G3 = read.csv(paste(resdir, 'array_aligned_G3.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

# class(array_aligned_G1)
# typeof(array_aligned_G1)
dim(array_aligned_G1)
length(array_aligned_G1$X1)

# Descriptives 
merged_demographics_cleaned = read.csv(paste(resdir, 'demographics_cleaned.csv', sep = ''), fileEncoding = 'UTF-8-BOM')



###### Unrelated Sample

# Aligned gradient values (unrelated sample)
array_aligned_G1_unrel_1 = read.csv(paste(resdir, 'unrelated_sample/array_aligned_G1_unrel_1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
array_aligned_G2_unrel_1 = read.csv(paste(resdir, 'unrelated_sample/array_aligned_G2_unrel_1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
array_aligned_G3_unrel_1 = read.csv(paste(resdir, 'unrelated_sample/array_aligned_G3_unrel_1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

# Descriptives
merged_demographics_cleaned_unrel_1 = read.csv(paste(resdir, 'unrelated_sample/merged_demographics_cleaned_unrel_1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')



###### Random n = 430 Sample

# Aligned gradient values (unrelated sample)
array_aligned_G1_rand430 = read.csv(paste(resdir, 'random_430/array_aligned_G1_rand430.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
array_aligned_G2_rand430 = read.csv(paste(resdir, 'random_430/array_aligned_G2_rand430.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
array_aligned_G3_rand430 = read.csv(paste(resdir, 'random_430/array_aligned_G3_rand430.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

# Descriptives
merged_demographics_cleaned_rand430 = read.csv(paste(resdir, 'random_430/demographics_cleaned_rand430.csv', sep = ''), fileEncoding = 'UTF-8-BOM')




----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DESCRIPTIVES 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###### Full Sample
  
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



###### Unrelated Sample

str(merged_demographics_cleaned_unrel_1, list.len=ncol(merged_demographics_cleaned_unrel_1))  # untruncated output
dim(merged_demographics_cleaned_unrel_1)


# Sample size by sex
ftable(merged_demographics_cleaned_unrel_1$Gender)


# Age
summarise.descriptives(merged_demographics_cleaned_unrel_1, merged_demographics_cleaned_unrel_1$Age_in_Yrs)

# Age by sex
ddply(merged_demographics_cleaned_unrel_1, 'Gender', summarise,
      mean = mean(Age_in_Yrs),
      sd = sd(Age_in_Yrs),
      min = min(Age_in_Yrs),
      max = max(Age_in_Yrs),
      sem = sd(Age_in_Yrs)/sqrt(length(Age_in_Yrs)),
      IQR = IQR(Age_in_Yrs))


# ICV
summarise.descriptives(merged_demographics_cleaned_unrel_1, merged_demographics_cleaned_unrel_1$FS_IntraCranial_Vol)

# ICV by sex
ddply(merged_demographics_cleaned_unrel_1, 'Gender', summarise,
      mean = mean(FS_IntraCranial_Vol),
      sd = sd(FS_IntraCranial_Vol),
      min = min(FS_IntraCranial_Vol),
      max = max(FS_IntraCranial_Vol),
      sem = sd(FS_IntraCranial_Vol)/sqrt(length(FS_IntraCranial_Vol)),
      IQR = IQR(FS_IntraCranial_Vol))



----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LINEAR REGRESSION 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###### Full Sample
  
### model = Gradient_Eigenvalues ~ Sex + Age + ICV 

# run model
lm_G1_sex_age_icv_res = lm.sex_age_icv(df_dv = array_aligned_G1, df_iv = merged_demographics_cleaned)
lm_G2_sex_age_icv_res = lm.sex_age_icv(df_dv = array_aligned_G2, df_iv = merged_demographics_cleaned)
lm_G3_sex_age_icv_res = lm.sex_age_icv(df_dv = array_aligned_G3, df_iv = merged_demographics_cleaned)





# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

##### effect size calculation: try https://stats.stackexchange.com/questions/71816/calculating-effect-size-for-variables-in-a-multiple-regression-in-r 


# to test on 1 value:

# Fit a linear model: lm = Gradient_Eigenvalues ~ Sex + Age + ICV
lm_fit = lm(array_aligned_G1[[1]] ~ merged_demographics_cleaned$Gender + merged_demographics_cleaned$Age_in_Yrs + merged_demographics_cleaned$FS_IntraCranial_Vol)
summary(lm_fit)



# function with integrated output of semi parcellation squared (error y must be numeric)

lm.sex_age_icv_TEST <- function(df_dv, df_iv) {
  
  '
    - fits and runs linear model to test for SEX effects, including sex, age and ICV in the model as covariates
    - to supply: df_dv (dataframe containing the dependent variable), df_iv (dataframe containing the independent variables)
    - outputs dataframe containing t-values, p-values, and FDR-corrected q-values for SEX contrast
  '
  
  # Create empty vectors (0s) of type "double precision" and length of len(df_dv) 
  t_val_sex = vector(mode = "double", length = ncol(df_dv))  
  p_val_sex = vector(mode = "double", length = ncol(df_dv))
  beta_val_sex = vector(mode = "double", length = ncol(df_dv))
  
  sp_sq_sex = vector(mode = "double", length = ncol(df_dv))
  
  # Degrees of Freedom = N (subjects) - number of IV (3 *** HARD-CODED FOR THIS SPECIFIC LM MODEL ***) - 1 (mean)
  #DoF = nrow(df_dv) - 3 -1
  
  # Loop over the df_dv columns (= parcels)
  for (i in seq_along(df_dv)) {
    
    # Fit a linear model: lm = Gradient_Eigenvalues ~ Sex + Age + ICV
    lm_fit = lm(df_dv[[i]] ~ df_iv$Gender + df_iv$Age_in_Yrs + df_iv$FS_IntraCranial_Vol)
    
    # Extract from summary of lm_fit the t- and p-values
    # summary(lm_fit)$coefficients[row, column]; row = 1 intercept, 2 sex, 3 Age, 4 ICV; columns = 1 Estimate, 2 Std. Error, 3 t-value, 4 p-value
    t_val_sex[[i]] = summary(lm_fit)$coefficients[2,3]
    p_val_sex[[i]] = summary(lm_fit)$coefficients[2,4]  # if want to calculate p value by hand: p_val_sex[[i]] = 2*pt(abs(t_val_sex[[i]]), DoF, lower.tail = F)
    beta_val_sex[[i]] = summary(lm_fit)$coefficients[2,1]
    
    sp_sq_sex[[i]] = (semi.r(y = df_dv[[i]], x = df_iv$Gender, given = df_iv$Age_in_Yrs))^2
 
  }
  
  # Calculate FDR-corrected q-values from p-values
  q_val_sex = p.adjust(p_val_sex, method = "fdr")
  
  # Create output dataframe containing t-values, p-values, and q-values
  output_df = data.frame(t_val_sex, p_val_sex, q_val_sex, beta_val_sex, sp_sq_sex)
  
  return(output_df)
}





# Effect size measure: semi-partial correlation squared (gives you the proportion of variance in y accounted for by x1 having controlled for x2) 

# function to compute the semi-partial r (see explanation at https://stats.stackexchange.com/questions/71816/calculating-effect-size-for-variables-in-a-multiple-regression-in-r)
semi.r = function(y, x, given){  
  ryx  = cor(y, x)
  ryg  = cor(y, given)
  rxg  = cor(x, given)
  num  = ryx - (ryg*rxg)
  dnm  = sqrt( (1-rxg^2) )
  sp.r = num/dnm
  return(sp.r)
}



# demonstration -> works

set.seed(9503)                   # this makes the example exactly reproducible
x1 = rnorm(10)                   # these variables are uncorrelated in the population
x2 = rnorm(10)                   # but not perfectly uncorrelated in this sample:
cor(x1, x2)                      # [1]  0.1265472
y  = 4 + .5*x1 - .3*x2 + rnorm(10, mean=0, sd=1)
model = lm(y~x1+x2)
summary(model)
# ...
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   4.1363     0.4127  10.022 2.11e-05 ***
# x1            0.1754     0.3800   0.461    0.658    
# x2           -0.6181     0.3604  -1.715    0.130    
# ...
sp.x1 = semi.r(y=y, x=x1, given=x2);  sp.x1  # [1]  0.1459061
sp.x1^2                                      # [1]  0.02128858
c.x2 = cor(x2, y);  c.x2                     # [1] -0.5280958
c.x2^2                                       # [1]  0.2788852
c.x2^2 + sp.x1^2                             # [1]  0.3001738
summary(model)$r.squared                     # [1]  0.3001738



# trial on my array G1 -> error y must be numeric

sp.sex = semi.r(y = array_aligned_G1[[1]], x = merged_demographics_cleaned$Gender, given = merged_demographics_cleaned$Age_in_Yrs)
class(array_aligned_G1[[1]])

lm.sex_age_icv_TEST(df_dv = array_aligned_G1, df_iv = merged_demographics_cleaned)






# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------




# number of significant parcels
sum(lm_G1_sex_age_icv_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))
sum(lm_G2_sex_age_icv_res$q_val_sex < 0.05, na.rm=TRUE)  
sum(lm_G3_sex_age_icv_res$q_val_sex < 0.05, na.rm=TRUE)  





### model = Gradient_Eigenvalues ~ Sex + Age + ICV + controlling for family relatedness

# run model

lm_G1_sex_age_icv_ctrl_related_res = lmer.sex_age_icv_ctrl_related(df_dv = array_aligned_G1, df_iv = merged_demographics_cleaned)
lm_G2_sex_age_icv_ctrl_related_res = lmer.sex_age_icv_ctrl_related(df_dv = array_aligned_G2, df_iv = merged_demographics_cleaned)
lm_G3_sex_age_icv_ctrl_related_res = lmer.sex_age_icv_ctrl_related(df_dv = array_aligned_G3, df_iv = merged_demographics_cleaned)

# number of significant parcels
sum(lm_G1_sex_age_icv_ctrl_related_res$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))
sum(lm_G2_sex_age_icv_ctrl_related_res$q_val_sex < 0.05, na.rm=TRUE)  
sum(lm_G3_sex_age_icv_ctrl_related_res$q_val_sex < 0.05, na.rm=TRUE)  


#test = lmer(array_aligned_G1[[1]] ~ merged_demographics_cleaned$Gender + merged_demographics_cleaned$Age_in_Yrs + merged_demographics_cleaned$FS_IntraCranial_Vol + (1 | merged_demographics_cleaned$Family_ID), REML = FALSE)  
#lm(array_aligned_G1[[1]] ~ merged_demographics_cleaned$Gender + merged_demographics_cleaned$Age_in_Yrs + merged_demographics_cleaned$FS_IntraCranial_Vol)






###### Unrelated Sample

### model = Gradient_Eigenvalues ~ Sex + Age + ICV 

# run model
lm_G1_sex_age_icv_res_unrel_1 = lm.sex_age_icv(df_dv = array_aligned_G1_unrel_1, df_iv = merged_demographics_cleaned_unrel_1)
lm_G2_sex_age_icv_res_unrel_1 = lm.sex_age_icv(df_dv = array_aligned_G2_unrel_1, df_iv = merged_demographics_cleaned_unrel_1)
lm_G3_sex_age_icv_res_unrel_1 = lm.sex_age_icv(df_dv = array_aligned_G3_unrel_1, df_iv = merged_demographics_cleaned_unrel_1)


# number of significant parcels
sum(lm_G1_sex_age_icv_res_unrel_1$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))
sum(lm_G2_sex_age_icv_res_unrel_1$q_val_sex < 0.05, na.rm=TRUE)  
sum(lm_G3_sex_age_icv_res_unrel_1$q_val_sex < 0.05, na.rm=TRUE)  



###### Random n=430 Sample

### model = Gradient_Eigenvalues ~ Sex + Age + ICV 

# run model
lm_G1_sex_age_icv_res_rand430 = lm.sex_age_icv(df_dv = array_aligned_G1_rand430, df_iv = merged_demographics_cleaned_rand430)
lm_G2_sex_age_icv_res_rand430 = lm.sex_age_icv(df_dv = array_aligned_G2_rand430, df_iv = merged_demographics_cleaned_rand430)
lm_G3_sex_age_icv_res_rand430 = lm.sex_age_icv(df_dv = array_aligned_G3_rand430, df_iv = merged_demographics_cleaned_rand430)


# number of significant parcels
sum(lm_G1_sex_age_icv_res_rand430$q_val_sex < 0.05, na.rm=TRUE)  # other way: length(which(G1_lm_res$q_val_sex < 0.05))
sum(lm_G2_sex_age_icv_res_rand430$q_val_sex < 0.05, na.rm=TRUE)  
sum(lm_G3_sex_age_icv_res_rand430$q_val_sex < 0.05, na.rm=TRUE)  



----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# COMPARISON RELATED AND UNRELATED SAMPLES
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
# Spearman correlation of t-values
cor.test(lm_G1_sex_age_icv_res$t_val_sex, lm_G1_sex_age_icv_res_unrel_1$t_val_sex, method="spearman")
cor.test(lm_G2_sex_age_icv_res$t_val_sex, lm_G2_sex_age_icv_res_unrel_1$t_val_sex, method="spearman")
cor.test(lm_G3_sex_age_icv_res$t_val_sex, lm_G3_sex_age_icv_res_unrel_1$t_val_sex, method="spearman")

# Spearman correlation of beta values
    # more indicative than comparing t-values because betas should be similar (regardless of the sample size), whereas t will be different (in function of sample size) because there is a different SE
    # in true null, SE will increase and will create a lot of “noise”
    # if betas are very different, it is a red flag -> betas are indeed very different
cor.test(lm_G1_sex_age_icv_res$beta_val_sex, lm_G1_sex_age_icv_res_unrel_1$beta_val_sex, method="spearman")
cor.test(lm_G2_sex_age_icv_res$beta_val_sex, lm_G2_sex_age_icv_res_unrel_1$beta_val_sex, method="spearman")
cor.test(lm_G3_sex_age_icv_res$beta_val_sex, lm_G3_sex_age_icv_res_unrel_1$beta_val_sex, method="spearman")




----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# EXPORT RESULTS 
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
###### Full Sample
write.csv(lm_G1_sex_age_icv_res, paste(resdir, 'R_lm_G1_sex_age_icv_res.csv', sep = ''), row.names = FALSE)
write.csv(lm_G2_sex_age_icv_res, paste(resdir, 'R_lm_G2_sex_age_icv_res.csv', sep = ''), row.names = FALSE)
write.csv(lm_G3_sex_age_icv_res, paste(resdir, 'R_lm_G3_sex_age_icv_res.csv', sep = ''), row.names = FALSE)


###### Full Sample - model controlling for family relatedness

write.csv(lm_G1_sex_age_icv_ctrl_related_res, paste(resdir, 'full_sample_ctrl_related/R_lm_G1_sex_age_icv_ctrl_related_res.csv', sep = ''), row.names = FALSE)
write.csv(lm_G2_sex_age_icv_ctrl_related_res, paste(resdir, 'full_sample_ctrl_related/R_lm_G2_sex_age_icv_ctrl_related_res.csv', sep = ''), row.names = FALSE)
write.csv(lm_G3_sex_age_icv_ctrl_related_res, paste(resdir, 'full_sample_ctrl_related/R_lm_G3_sex_age_icv_ctrl_related_res.csv', sep = ''), row.names = FALSE)


###### Unrelated Sample
write.csv(lm_G1_sex_age_icv_res_unrel_1, paste(resdir, 'unrelated_sample/R_lm_G1_sex_age_icv_res_unrel_1.csv', sep = ''), row.names = FALSE)
write.csv(lm_G2_sex_age_icv_res_unrel_1, paste(resdir, 'unrelated_sample/R_lm_G2_sex_age_icv_res_unrel_1.csv', sep = ''), row.names = FALSE)
write.csv(lm_G3_sex_age_icv_res_unrel_1, paste(resdir, 'unrelated_sample/R_lm_G3_sex_age_icv_res_unrel_1.csv', sep = ''), row.names = FALSE)


###### Random n=430 Sample
write.csv(lm_G1_sex_age_icv_res_rand430, paste(resdir, 'random_430/R_lm_G1_sex_age_icv_res_rand430.csv', sep = ''), row.names = FALSE)
write.csv(lm_G2_sex_age_icv_res_rand430, paste(resdir, 'random_430/R_lm_G2_sex_age_icv_res_rand430.csv', sep = ''), row.names = FALSE)
write.csv(lm_G3_sex_age_icv_res_rand430, paste(resdir, 'random_430/R_lm_G3_sex_age_icv_res_rand430.csv', sep = ''), row.names = FALSE)