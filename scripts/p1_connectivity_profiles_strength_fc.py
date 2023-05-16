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
### Load binary matrices
########################################

# loading the full sample binary matrices
binary_top_fc = scipy.io.loadmat(resdir_hcp+'binary_top_fc.mat')['binary_top_fc']

# loading the binary matrices by sex
binary_top_fc_M = scipy.io.loadmat(resdir_hcp+'binary_top_fc.mat')['binary_top_fc_M']
binary_top_fc_F = scipy.io.loadmat(resdir_hcp+'binary_top_fc.mat')['binary_top_fc_F']




########################################
### Load demographics data to retrieve sex of subjects
########################################

HCP_demographics_cleaned_final = pd.read_csv(resdir_hcp+'demographics_cleaned_final.csv', index_col=0)





################################################################################
################################################################################
################################################################################

###########    Compute functional connectivity strength for the mean connectivity profiles    ################

################################################################################
################################################################################
################################################################################


### Identify the areas that most frequently come up as top 10% connection (for each seed)

print("----- Identifying the areas that most frequently come up as top 10% connection  -----")



## adding up the binary matrices across subjects to get "counts" matrices (across all subjects, males and females)
# counting the connections marked as top 10% (yields a 400x400 matrix of the summed counts) -> across all subjects
binary_top_fc_counts = np.sum(binary_top_fc, axis = 0)

# counting the connections marked as top 10% (yields a 400x400 matrix of the summed counts) -> across males
binary_top_fc_counts_M = np.sum(binary_top_fc_M, axis = 0)

# counting the connections marked as top 10% (yields a 400x400 matrix of the summed counts) -> across females
binary_top_fc_counts_F = np.sum(binary_top_fc_F, axis = 0)



## identifying the top 10% indices (across all subjects, males and females)

# list that will contain the "mean connectivity profile" (i.e. the most recurrent top 10% connections) across all subjects -> group level
most_frequent_connections_indices = []
most_frequent_connections_indices_M = []
most_frequent_connections_indices_F = []

# iterate over rows of the 400x400 matrix
for i in range(len(binary_top_fc_counts)):

    # identify the indices of the connections that most often appeared as top 10% functional connections (with highest count) i.e., top 20 connections (10% of 200 (400/2) given that geodesic distance was calculated across single hemispheres only)
    # code explanation available at: https://www.educative.io/answers/how-to-get-indices-of-n-maximum-values-in-numpy-arrays
    # retrieving indices from row i of the  matrix
    
    # if left hemisphere
    if i <= 199:
        # take top funnctional connections from left hemisphere (across subjects)
        binary_top_fc_indices_top10 = np.argsort(binary_top_fc_counts[i][:200])[::-1][:20]
        
        # in males
        binary_top_fc_indices_top10_M = np.argsort(binary_top_fc_counts_M[i][:200])[::-1][:20]
        
        # in females
        binary_top_fc_indices_top10_F = np.argsort(binary_top_fc_counts_F[i][:200])[::-1][:20]
        

    else:
        # take top funnctional connections from right hemisphere (across subjects) (need to add 200 to the indices outputted given that we are taking them from 200-parcel right hemisphere but applying them to 400 array)
        binary_top_fc_indices_top10 = 200 + np.argsort(binary_top_fc_counts[i][200:])[::-1][:20]
        
        # in males
        binary_top_fc_indices_top10_M = 200 + np.argsort(binary_top_fc_counts_M[i][200:])[::-1][:20]
        
        # in females
        binary_top_fc_indices_top10_F = 200 + np.argsort(binary_top_fc_counts_F[i][200:])[::-1][:20]
        
        
    
    most_frequent_connections_indices.append(binary_top_fc_indices_top10)
    most_frequent_connections_indices_M.append(binary_top_fc_indices_top10_M)
    most_frequent_connections_indices_F.append(binary_top_fc_indices_top10_F)
    

# make as arrays    
most_frequent_connections_indices = np.array(most_frequent_connections_indices)
most_frequent_connections_indices_M = np.array(most_frequent_connections_indices_M)
most_frequent_connections_indices_F = np.array(most_frequent_connections_indices_F)







##### Compute the functional connectivity strengths of the mean connectivity profiles (i.e., for each seed region, the connections with areas that recur most frequently) across all subjects, across males, and across females

print("----- Computing the fc strengths of the mean connectivity profiles  -----")



### lists that will contain the fc strengths for the most frequent connections 
# for all subjects (N x 400 x 20)
fc_strengths_most_frequent_connections = []

# for males (n x 400 x 20)
fc_strengths_most_frequent_connections_M = []

# for females (n x 400 x 20)
fc_strengths_most_frequent_connections_F = []


### lists that will contain the mean fc strengths for the most frequent connections 
# for all subjects (N x 400)
mean_fc_strengths_most_frequent_connections = []

# for males (n x 400)
mean_fc_strengths_most_frequent_connections_M = []

# for females (n x 400)
mean_fc_strengths_most_frequent_connections_F = []



# iterate over subjects (N) 
for i in range(len(fc_matrices)):

    # temporary lists for each subject's fc strengths for the most frequent connections and the mean of those
    sub_fc_strengths_most_frequent_connections = []
    sub_mean_fc_strength_most_frequent_connections = []
    
    M_sub_fc_strengths_most_frequent_connections = []
    M_sub_mean_fc_strength_most_frequent_connections = []
    
    F_sub_fc_strengths_most_frequent_connections = []
    F_sub_mean_fc_strength_most_frequent_connections = []
    
    

    # iterate over each row of the subject's FC 400x400 matrix 
    for j in range(len(fc_matrices[0])):
        
        # retrieve the fc strengths of each of the most_frequent_connections_indices 
        fc_strength_most_frequent_connections = fc_matrices[i][j][most_frequent_connections_indices[j]]

        # compute the mean fc strength of the most frequent connections 
        mean_fc_strength_most_frequent_connections = np.mean(fc_strength_most_frequent_connections)

        # append to respective lists (temp subject list)
        sub_fc_strengths_most_frequent_connections.append(fc_strength_most_frequent_connections)
        sub_mean_fc_strength_most_frequent_connections.append(mean_fc_strength_most_frequent_connections)
        
        
        # do the same as above but in respective sex list
        if HCP_demographics_cleaned_final.Gender.tolist()[i] == 'M':
            
            # retrieve the fc strengths of each of the most_frequent_connections_indices for males!
            fc_strength_most_frequent_connections = fc_matrices[i][j][most_frequent_connections_indices_M[j]]

            # compute the mean fc strength of the most frequent connections 
            mean_fc_strength_most_frequent_connections = np.mean(fc_strength_most_frequent_connections)

            # append to respective lists (temp subject list)
            M_sub_fc_strengths_most_frequent_connections.append(fc_strength_most_frequent_connections)
            M_sub_mean_fc_strength_most_frequent_connections.append(mean_fc_strength_most_frequent_connections)
            
            
        else:
            
            # retrieve the fc strengths of each of the most_frequent_connections_indices for females!
            fc_strength_most_frequent_connections = fc_matrices[i][j][most_frequent_connections_indices_F[j]]

            # compute the mean fc strength of the most frequent connections 
            mean_fc_strength_most_frequent_connections = np.mean(fc_strength_most_frequent_connections)

            # append to respective lists (temp subject list)
            F_sub_fc_strengths_most_frequent_connections.append(fc_strength_most_frequent_connections)
            F_sub_mean_fc_strength_most_frequent_connections.append(mean_fc_strength_most_frequent_connections)
            

            
    # append current subject's lists to lists containing all subjects
    fc_strengths_most_frequent_connections.append(sub_fc_strengths_most_frequent_connections)
    mean_fc_strengths_most_frequent_connections.append(sub_mean_fc_strength_most_frequent_connections)
    
    # append curent subject's lists to lists containing only male or female depending on what the subject is
    if HCP_demographics_cleaned_final.Gender.tolist()[i] == 'M':
        fc_strengths_most_frequent_connections_M.append(M_sub_fc_strengths_most_frequent_connections)
        mean_fc_strengths_most_frequent_connections_M.append(M_sub_mean_fc_strength_most_frequent_connections)
        
    else:
        fc_strengths_most_frequent_connections_F.append(F_sub_fc_strengths_most_frequent_connections)
        mean_fc_strengths_most_frequent_connections_F.append(F_sub_mean_fc_strength_most_frequent_connections)
            
            
            
# make as arrays
fc_strengths_most_frequent_connections = np.array(fc_strengths_most_frequent_connections)
mean_fc_strengths_most_frequent_connections = np.array(mean_fc_strengths_most_frequent_connections)
fc_strengths_most_frequent_connections_M = np.array(fc_strengths_most_frequent_connections_M)
mean_fc_strengths_most_frequent_connections_M = np.array(mean_fc_strengths_most_frequent_connections_M)
fc_strengths_most_frequent_connections_F = np.array(fc_strengths_most_frequent_connections_F)
mean_fc_strengths_most_frequent_connections_F = np.array(mean_fc_strengths_most_frequent_connections_F)






### Export sex differences results matrices in matfile

print("----- Exporting results at /data/p_02667/sex_diff_gradients/results/HCP/connectivity_profiles_fc_strengths.mat  -----")


mdict = {'fc_strengths_most_frequent_connections': fc_strengths_most_frequent_connections, 'mean_fc_strengths_most_frequent_connections': mean_fc_strengths_most_frequent_connections, 'most_frequent_connections_indices': most_frequent_connections_indices, 'fc_strengths_most_frequent_connections_M': fc_strengths_most_frequent_connections_M, 'mean_fc_strengths_most_frequent_connections_M': mean_fc_strengths_most_frequent_connections_M, 'most_frequent_connections_indices_M': most_frequent_connections_indices_M, 'fc_strengths_most_frequent_connections_F': fc_strengths_most_frequent_connections_F, 'mean_fc_strengths_most_frequent_connections_F': mean_fc_strengths_most_frequent_connections_F, 'most_frequent_connections_indices_F': most_frequent_connections_indices_F}

scipy.io.savemat(resdir_hcp+'connectivity_profiles_fc_strengths.mat', mdict)

