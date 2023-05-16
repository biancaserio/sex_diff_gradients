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

resdir_hcp = '/data/p_02667/sex_diff_gradients/results/HCP/'
resdir_fig = '/data/p_02667/sex_diff_gradients/results/figures/'

dataout = '/data/p_02667/sex_diff_gradients/data/geodesic_distance/'





################################################################################
################################################################################
################################################################################

############################        Data        ################################

################################################################################
################################################################################
################################################################################

print("----- Loading data -----")

########################################
### Load fc data 
########################################

# taking the list of subjects that have fc data
HCP_sub_list_fc = scipy.io.loadmat(datadir+'fc_matrices/fc_matrices.mat')['HCP_sub_list_fc']

# taking the fc matrices from the mat file
fc_matrices_full = scipy.io.loadmat(datadir+'fc_matrices/fc_matrices.mat')['HCP_fc_matrices']



### Get subsample of fc matrices that has full data (matching HCP_sub_list_final)

HCP_sub_list_final = np.genfromtxt(datadir+'HCP_sub_list_final.csv', delimiter=',')
HCP_sub_list_final = np.array([str(int(e)) for e in HCP_sub_list_final])


# list that will contain the fc matrices in final sample
fc_matrices = []

for i in range(len(HCP_sub_list_fc)):
    if HCP_sub_list_fc[i] in HCP_sub_list_final:
        
        # append this subjects' list of fc matrices to the list of fc matrices in final sample
        fc_matrices.append(fc_matrices_full[i])
    
fc_matrices = np.array(fc_matrices)





########################################
### Load demographics data to retrieve sex of subjects
########################################

HCP_demographics_cleaned_final = pd.read_csv(resdir_hcp+'demographics_cleaned_final.csv', index_col=0)





################################################################################
################################################################################
################################################################################

###########    Compute functional connectivity strength at individual subject level    ################

################################################################################
################################################################################
################################################################################


print("----- Computing the fc strength of top 10% connections at the individual level  -----")


### lists that will contain the fc strengths for 10% connections at the individual level
# for all subjects (N x 400 x 20)
fc_strengths_top10 = []



# iterate over subjects
for i in range(len(fc_matrices)):
    print(f"----- subject = {i} ----")

    # temporary list for each subject's fc_strengths_top10 matrix
    sub_fc_strengths_top10 = []

    # iterate over each row of the subject's fc 400x400 matrix 
    for j in range(len(fc_matrices[0])):

        # identify the indices of the top 10% functional connections i.e., top 20 connections (10% of 200 (400/2) given that geodesic distance was calculated across single hemispheres only)
        # code explanation available at: https://www.educative.io/answers/how-to-get-indices-of-n-maximum-values-in-numpy-arrays
        # retrieving indices from subject i and row j of the fc matrix
        
        # if left hemisphere
        if j <= 199:
            # take top funnctional connections from left hemisphere
            indices_top_fc = np.argsort(fc_matrices[i][j][:200])[::-1][:20]
            
        else:
            # take top funnctional connections from right hemisphere (need to add 200 to the indices outputted given that we are taking them from 200-parcel right hemisphere but applying them to 400 array in fc_matrices)
            indices_top_fc = 200 + np.argsort(fc_matrices[i][j][200:])[::-1][:20]
            
            
        # temporary list for row of subject's fc_strengths_top10 matrix
        sub_row_fc_strengths_top10 = []
        
        # iterate over each parcel of subject's 400-parcel-long row
        for k in range(len(fc_matrices[0][0])):
                
            # if parcel belongs to one of the top 10% connections, assign it's value to list for row of subject's fc_strengths_most_frequent_connections matrix 
            if k in indices_top_fc:
                sub_row_fc_strengths_top10.append(fc_matrices[i][j][k])
        

        # append temporary list for row of subject's binary_top_fc matrix to subject's temporary list
        sub_fc_strengths_top10.append(sub_row_fc_strengths_top10)

    # append current subject's list of binary_top_fc to list containing all subjects
    fc_strengths_top10.append(sub_fc_strengths_top10)
            
            
fc_strengths_top10 = np.array(fc_strengths_top10)





### Export sex differences results matrices in matfile

print("----- Exporting results at /data/p_02667/sex_diff_gradients/results/HCP/fc_strengths_top10_individual_level.mat  -----")


mdict = {'fc_strengths_top10': fc_strengths_top10}

scipy.io.savemat(resdir_hcp+'fc_strengths_top10_individual_level.mat', mdict)


