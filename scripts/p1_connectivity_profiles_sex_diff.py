#!/usr/bin/env python3


########################################
### Load packages
########################################

print("----- Loading packages -----")

# General
import numpy as np
import pandas as pd
import os
import csv

# Computing
import scipy.io
from scipy import stats
from heapq import nlargest  # gives you the largest values of a list
from statsmodels.stats.multitest import fdrcorrection

# Visualization
import matplotlib.pyplot as plt

# Gradients
import brainspace
from brainspace.datasets import load_parcellation, load_conte69
from brainspace.plotting import plot_hemispheres
from brainspace.gradient import GradientMaps
from brainspace.utils.parcellation import map_to_labels



########################################
### Define directories
########################################

print("----- Defining directories -----")

datadir = '/data/p_02667/sex_diff_gradients/data/'
datadir_geodesic = '/data/p_02721/geodesic_HCP/mica_pipe/output/micapipe/'

resdir_hcp = '/data/p_02667/sex_diff_gradients/results/HCP/'
resdir_fig = '/data/p_02667/sex_diff_gradients/results/figures/'


path_list = os.listdir(datadir_geodesic)
path_list.sort()






################################################################################
################################################################################
################################################################################

############################        Data        ################################

################################################################################
################################################################################
################################################################################

print("----- Loading data -----")

########################################
### Load binary matrices
########################################

# loading the full sample binary matrices
binary_top_fc = scipy.io.loadmat(resdir_hcp+'binary_top_fc.mat')['binary_top_fc']

# loading the binary matrices by sex
binary_top_fc_M = scipy.io.loadmat(resdir_hcp+'binary_top_fc.mat')['binary_top_fc_M']
binary_top_fc_F = scipy.io.loadmat(resdir_hcp+'binary_top_fc.mat')['binary_top_fc_F']



################################################################################
################################################################################
################################################################################

#############    Test sex differences in connectivity profiles     #############

################################################################################
################################################################################
################################################################################

# Via:
# - Chi-square test of indepndence (H0) of variables in contingency table to test the indepdence of variables "sex" (rows) and "area is/is not in top 10% of functional connections" (columns) -> non-independence == sex differences in top 10% of connections
# - Odds Ratio (effect size)

# Explanation:
# - for every row (parcel) in binary matrix, for ever parcel in row of binary matrix -> make a contingency table 
# - for contingency table, need to a) count the number of males and females (separately) that have this area marked (with a 1) as being one of the top 10% of connections (simply take the sum across subjects for this cell in matrix given that if area is not in top 10% of connections, it is marked with 0) and b) deduce the counts of subjects for which this areas is not part of top 10% connections
# - With contingency table, can compute a) Chi-square test of independence of variables, b) (from contents of table) the Odds Ratio


print("----- Computing chisquare_pval_matrix and OR_matrix-----")




# list that will contain the matrix of p values for the Chi-square test of independence of variables in contingency table to see sex differences in top 10% of connections
chisquare_pval_matrix = []

# list that will contain the matrix of Odds Ratio from contingency table (to-be-used as effect size for Chi-square test of indepndence of variables)
OR_matrix = []

        
# transposing the binary_top_fc array in order to have it in shape 400 x 400 x N -> iterating over 400 parcels
for i in range(len(binary_top_fc.T)):
    
    print(f"-------- matrix row: {i} --------")
    
    # temporary lists that will contail the rows of the chisquare_pval_matrix and OR_matrix
    temp_row_chisquare_pval_matrix = []
    temp_row_OR_matrix = []
    
    
    # iterating over 400 parcels
    for j in range(len(binary_top_fc.T[0])):
        
        print(f"matrix column: {j}")
        
        # only run Chi-square test if i and j are not the same (otherwise we're at the intersection of the parcel's connection with itself in the matrix, and the value for it is 0 (i.e., not chosen as a top connection))
        if i != j:
        
            ### creating contingency table with counts of occurances for Chi-square test of indepndence of variables to test if there are statistically significant sex differences in top 10% of fc connections
            # C: connection (area) in top 10%; NC: connection (area) not in top 10%; m: male; f: female
            # row: male (top), female (bottom)
            # columns: C (left), NC right

            contingency_table = []
            Cm = 0
            NCm = 0
            Cf = 0
            NCf = 0
            
            
            ## male row
            Cm = sum(binary_top_fc_M.T[i][j])
            NCm = len(binary_top_fc_M) - Cm  # number of male subjects minus sum of counts in the top 10%

            # append male row to table
            contingency_table.append([Cm, NCm])


            ## female row
            Cf = sum(binary_top_fc_F.T[i][j])
            NCf = len(binary_top_fc_F) - Cf  # number of male subjects minus sum of counts in the top 10%

            # append female row to table
            contingency_table.append([Cf, NCf])

            contingency_table = np.array(contingency_table)
            
            
            # compute Odds Ratio (regardless of pvalue of Chi-square test)
            OR = (Cm/NCm)/(Cf/NCf)
            
            # in case the OR == infinity (because dividing by zero), recompute the OR: can replace the 0 at denominator by 0.5 for more interpretability (https://www.researchgate.net/post/How_to_calculate_OR_odd_ratio_if_one_of_groups_is_0_in_a_case-control_study)
            if OR == float('inf'):
                if NCm == 0:
                    NCm = 0.5
                if NCf == 0:
                    NCf = 0.5
                               
                # recompute OR 
                OR = (Cm/NCm)/(Cf/NCf)
                
                # if still == inf, it means that the whole denominator is == 0, replace it by 0.5
                if OR == float('inf'):
                    
                    # recompute OR
                    OR = (Cm/NCm)/0.5
            
            # append OR to row
            temp_row_OR_matrix.append(OR)
            
            
            # only run Chi-square test if there aren't just zeros in both columns (ie either male or female has area as one of top connections (totC > 0) or not all males and all females have area chosen as top connection (totNC > 0) - otherwise throws an error
            if (Cm + Cf > 0) and (NCm + NCf > 0):
                ### Chi-square test of independence of variables in a contingency table, [1] indexes the p val
                p_val = stats.chi2_contingency(contingency_table)[1]

                # append p value to row
                temp_row_chisquare_pval_matrix.append(p_val)
                
            else:
                # append value of 1 to matrix row because both males and females have a value of 0 for this cell in matrix therefore p > 0.05 anyway (no difference) 
                temp_row_chisquare_pval_matrix.append(1)
        
        else:
            # append value of 1 to matrix of p val row because we're at the intersection of the parcel's connection with itself in the matrix
            temp_row_chisquare_pval_matrix.append(1)
            
            # append value of 1 to matrix of OR row because we're at the intersection of the parcel's connection with itself in the matrix - fyi an odds ratio of exactly 1 means the odds of the event happening are the exact same in the exposed versus the non-exposed group.
            temp_row_OR_matrix.append(1)
                
        
    # append row to full matrix
    chisquare_pval_matrix.append(temp_row_chisquare_pval_matrix)
    OR_matrix.append(temp_row_OR_matrix)
    
    
# make into array
chisquare_pval_matrix = np.array(chisquare_pval_matrix)
OR_matrix = np.array(OR_matrix)



## compute FDR correction

# flatten() because function only takes a 1D array; [1] retrieves FDR-corrected q values
chisquare_qvals = fdrcorrection(chisquare_pval_matrix.flatten(), alpha=0.05, method='indep')[1]

# reshape the q values in the 400x400 array format
chisquare_qval_matrix = chisquare_qvals.reshape((400, 400))





### Export sex differences results matrices in arrays

print("----- Exporting results at /data/p_02667/sex_diff_gradients/results/HCP/top_fc_sex_diff_...  -----")

# arrays
np.savetxt(resdir_hcp+'top_fc_sex_diff_chisquare_pval_matrix.csv', chisquare_pval_matrix, delimiter=',', fmt = '%.16g', comments = '')
np.savetxt(resdir_hcp+'top_fc_sex_diff_chisquare_qval_matrix.csv', chisquare_qval_matrix, delimiter=',', fmt = '%.16g', comments = '')
np.savetxt(resdir_hcp+'top_fc_sex_diff_OR_matrix.csv', OR_matrix, delimiter=',', fmt = '%.16g', comments = '')
