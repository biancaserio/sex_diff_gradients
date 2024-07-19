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

print("----- Loading Data -----")

########################################
### Load fc data (needed to compute mean geodesic distance for top 10% of fc connections)
########################################

# load fc matrices
scipy.io.whosmat(datadir+'fc_matrices/fc_matrices.mat')

# taking the list of subjects that have fc data
HCP_sub_list_fc = scipy.io.loadmat(datadir+'fc_matrices/fc_matrices.mat')['HCP_sub_list_fc']

# taking the fc matrices from the mat file
fc_matrices_full = scipy.io.loadmat(datadir+'fc_matrices/fc_matrices.mat')['HCP_fc_matrices']



########################################
### Get subsample of fc matrices that has full data (matching HCP_sub_list_final)
########################################


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

################    Compute functional connectivity profiles    ################

################################################################################
################################################################################
################################################################################


########################################
### Binary matrices
########################################

print("----- Computing binary matrices -----")

# list that will contain the binary coded matrices indicating if area represents one of top 10% connections (1) or not (0)
binary_top_fc = []


# iterate over subjects
for i in range(len(fc_matrices)):
    print(f"----- subject = {i} ----")

    # temporary list for each subject's binary_top_fc matrix
    sub_binary_top_fc = []

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
            # take top funnctional connections from right hemisphere (need to add 200 to the indices outputted given that we are taking them from 200-parcel right hemisphere but applying them to 400 array in geodesic_distances_matrices_corrected)
            indices_top_fc = 200 + np.argsort(fc_matrices[i][j][200:])[::-1][:20]
            
            
        # temporary list for row of subject's binary_top_fc matrix
        sub_row_binary_top_fc = []
        
        # iterate over each parcel of subject's 400-parcel-long row
        for k in range(len(fc_matrices[0][0])):
                
            # if parcel belongs to one of the top 10% connections, assign a 1 to it; else assign a 0
            if k in indices_top_fc:
                sub_row_binary_top_fc.append(1)

            else:
                sub_row_binary_top_fc.append(0)
        

        # append temporary list for row of subject's binary_top_fc matrix to subject's temporary list
        sub_binary_top_fc.append(sub_row_binary_top_fc)

    # append current subject's list of binary_top_fc to list containing all subjects
    binary_top_fc.append(sub_binary_top_fc)
            
            
binary_top_fc = np.array(binary_top_fc)





### Separate the binary matrices into male and female

print("----- Separating binary matrices into male and female binary matrices -----")

binary_top_fc_M = []
binary_top_fc_F = []

# iterate over subject (matrices)
for i in range(len(binary_top_fc)):

    if HCP_demographics_cleaned_final.Gender.tolist()[i] == 'M':
        binary_top_fc_M.append(binary_top_fc[i])
        
    else:
        binary_top_fc_F.append(binary_top_fc[i])

binary_top_fc_M = np.array(binary_top_fc_M)
binary_top_fc_F = np.array(binary_top_fc_F)





### Export binary matrices in mat file

print("----- Exporting binary matrices at /data/p_02667/sex_diff_gradients/results/HCP/binary_top_fc.mat  -----")

mdict = {'binary_top_fc': binary_top_fc, 'binary_top_fc_M': binary_top_fc_M, 'binary_top_fc_F': binary_top_fc_F}

scipy.io.savemat(resdir_hcp+'binary_top_fc.mat', mdict)



