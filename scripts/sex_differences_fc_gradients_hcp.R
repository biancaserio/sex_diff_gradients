


####### Load packages #######
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





###### Directories #######

# set up directories
codedir = dirname(getActiveDocumentContext()$path)  # get path to current script
datadir = '/data/p_02667/sex_diff_gradients/data/'
resdir = '/data/p_02667/sex_diff_gradients/results/hcp/'


# set directory to path of current script
setwd(codedir) 




###### Get data ######

# Aligned gradient values
array_aligned_G1 = read.csv(paste(resdir, 'array_aligned_G1.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
array_aligned_G2 = read.csv(paste(resdir, 'array_aligned_G2.csv', sep = ''), fileEncoding = 'UTF-8-BOM')
array_aligned_G3 = read.csv(paste(resdir, 'array_aligned_G3.csv', sep = ''), fileEncoding = 'UTF-8-BOM')

dim(array_aligned_G1)
length(array_aligned_G1$X1)

# # Dataframe containing all aligned gradients
# aligned_gradients = do.call(rbind, Map(data.frame, G1 = array_aligned_G1, G2 = array_aligned_G2, G3 = array_aligned_G3))
# str(aligned_gradients)
# ## PROBLEM: concatenates everything, doesnt't keep the 400 parcels separate for each subject - think about this...

# Descriptives
merged_demographics_cleaned = read.csv('merged_demographics_cleaned.csv', fileEncoding = 'UTF-8-BOM')




###### Descriptives ######

str(merged_demographics_cleaned, list.len=ncol(merged_demographics_cleaned))  # untruncated output


### Age

summarise(merged_demographics_cleaned,
          mean = mean(Age_in_Yrs),
          sd = sd(Age_in_Yrs),
          min = min(Age_in_Yrs),
          max = max(Age_in_Yrs),
          sem = sd(Age_in_Yrs)/sqrt(length(Age_in_Yrs)),
          IQR = IQR(Age_in_Yrs))


# Age by sex

ddply(merged_demographics_cleaned, 'Gender', summarise,
          mean = mean(Age_in_Yrs),
          sd = sd(Age_in_Yrs),
          min = min(Age_in_Yrs),
          max = max(Age_in_Yrs),
          sem = sd(Age_in_Yrs)/sqrt(length(Age_in_Yrs)),
          IQR = IQR(Age_in_Yrs))

ftable(merged_demographics_cleaned$Gender)


### ICV

summarise(merged_demographics_cleaned,
          mean = mean(FS_IntraCranial_Vol),
          sd = sd(FS_IntraCranial_Vol),
          min = min(FS_IntraCranial_Vol),
          max = max(FS_IntraCranial_Vol),
          sem = sd(FS_IntraCranial_Vol)/sqrt(length(FS_IntraCranial_Vol)),
          IQR = IQR(FS_IntraCranial_Vol))




###### Linear models ######

lm_fit = lm(array_aligned_G1$X1 ~ merged_demographics_cleaned$Gender + merged_demographics_cleaned$Age_in_Yrs)
summary(lm_fit)

summary(lm_fit)$terms

names(summary(lm_fit))

p_val2*pt(abs())


# ERROR: tried unlist() but then says variable lenth differs