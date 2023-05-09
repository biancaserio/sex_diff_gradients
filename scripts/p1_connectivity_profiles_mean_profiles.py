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
resdir_test = '/data/p_02667/sex_diff_gradients/results/test/'

dataout = '/data/p_02667/sex_diff_gradients/data/geodesic_distance/'


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
### Load geodesic distance (subject-level, i.e., based on their individual top fc connections)
########################################

mean_geodesic_distances_full = np.genfromtxt(datadir+'geodesic_distance/mean_geodesic_distances.csv', delimiter=',')


### Get subsample of fc matrices that has full data (matching HCP_sub_list_final)

HCP_sub_list_geodesic = pd.read_csv(datadir+'geodesic_distance/HCP_sub_list_geodesic.csv', header = None)[0].to_list()
HCP_sub_list_geodesic = np.array([str(int(e)) for e in HCP_sub_list_geodesic])


# list that will contain the mean of the top geodesic distances in final sample
mean_geodesic_distances = []

for i in range(len(HCP_sub_list_geodesic)):
    if HCP_sub_list_geodesic[i] in HCP_sub_list_final:
        
        # append this subjects' list of mean top geodesic distances to the list of mean of the top geodesic distances in final sample
        mean_geodesic_distances.append(mean_geodesic_distances_full[i])
        
    
mean_geodesic_distances = np.array(mean_geodesic_distances)


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

#######    Raw geodesic distance data -> geodesic distance matrices    #########

################################################################################
################################################################################
################################################################################


### Retrieve geodesic distance data from original directory and shaping it into a readable format (from txt file to 400x400 matrix) whilst removing medial mask (original number of datapoints 402x402)

print("----- Retrieving raw geodesic distance data and obtaining geodesic distance matrices -----")


# list that will contain the matrices of the raw geodesic distance (pre-computed pro hemisphere - cross hemisphere cells are coded with 0)
geodesic_distances_matrices = []

# list that will contain the subjects IDs with geodesic distance data
HCP_sub_list_geodesic = []


# iterate over subjects in directory
for dir in path_list:
    if os.path.exists(datadir_geodesic+dir+'/ses-01/anat/surfaces/geo_dist/'+dir+'_ses-01_space-fsnative_atlas-schaefer-400_GD.txt'):
        
        print('---------------- executing: '+dir+' ----------------')

        ### save subject ID to list of subjects with geodesic distance data (save just the number (removing "sub" before "-")
        HCP_sub_list_geodesic.append(dir.split("-")[1])


        ### Get geodesic distance data

        # temporary list for each subject containing a matrix of the raw geodesic distances
        temp_matrix_geod_list = []

        # read txt file
        reader = csv.reader(open(datadir_geodesic+dir+'/ses-01/anat/surfaces/geo_dist/'+dir+'_ses-01_space-fsnative_atlas-schaefer-400_GD.txt'), delimiter="\t")

        # creating a counter for index
        count_i = 0

        # iterate over each row of the text file
        for row in reader:

            # there are 402 values instead of 400 -> indices 0 and 201 are medial wall masks so need to be removed -> only include the rest
            if count_i != 0 and count_i != 201:

                # data is separated by space so need to split by space to make each number an element of the list
                splitted_row = row[0].split(" ")

                # convert every number into a float (from string)
                float_splitted_row = [float(current_integer) for current_integer in splitted_row]

                # there are 402 values instead of 400 -> indices 0 and 201 are medial wall masks so need to be removed -> only include the rest
                float_splitted_row_corrected = float_splitted_row[1:201]+float_splitted_row[202:402]

                # append the whole raw geodesic distances that will constitute a geodesic distance matrix
                temp_matrix_geod_list.append(float_splitted_row_corrected)

            # +1 to counter
            count_i += 1

        # append current subject's matrix of raw geodesic distances to list containing all subjects
        geodesic_distances_matrices.append(temp_matrix_geod_list)
        

            
    else:
        print('---------------- '+dir+' does not have geodesic distance data ----------------') 
        
        
# transform lists to array format
geodesic_distances_matrices = np.array(geodesic_distances_matrices)




### Get subsample of geodesic distance matrices that has full data (matching HCP_sub_list_final)

# list that will contain the geodesic distance matrices in final sample
geodesic_distances_matrices_corrected = []

for i in range(len(HCP_sub_list_geodesic)):
    if HCP_sub_list_geodesic[i] in HCP_sub_list_final:
        
        # append this subjects' list of fc matrices to the list of fc matrices in final sample
        geodesic_distances_matrices_corrected.append(geodesic_distances_matrices[i])
    
geodesic_distances_matrices_corrected = np.array(geodesic_distances_matrices_corrected)






################################################################################
################################################################################
################################################################################

###########    Compute mean functional connectivity profiles    ################

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






##### Compute the geodesic distances of the mean connectivity profiles (i.e., for each seed region, the connections with areas that recur most frequently) across all subjects, across males, and across females

print("----- Computing the geodesic distances of the mean connectivity profiles  -----")



### lists that will contain the geodesic distances for the most frequent connections 
# for all subjects (N x 400 x 20)
geodesic_distances_most_frequent_connections = []

# for males (n x 400 x 20)
geodesic_distances_most_frequent_connections_M = []

# for females (n x 400 x 20)
geodesic_distances_most_frequent_connections_F = []


### lists that will contain the mean geodesic distance for the most frequent connections 
# for all subjects (N x 400)
mean_geodesic_distances_most_frequent_connections = []

# for males (n x 400)
mean_geodesic_distances_most_frequent_connections_M = []

# for females (n x 400)
mean_geodesic_distances_most_frequent_connections_F = []



# iterate over subjects (N) 
for i in range(len(geodesic_distances_matrices_corrected)):

    # temporary lists for each subject's geodesic distances for the most frequent connections and the mean of those
    sub_geodesic_distances_most_frequent_connections = []
    sub_mean_geodesic_distance_most_frequent_connections = []
    
    M_sub_geodesic_distances_most_frequent_connections = []
    M_sub_mean_geodesic_distance_most_frequent_connections = []
    
    F_sub_geodesic_distances_most_frequent_connections = []
    F_sub_mean_geodesic_distance_most_frequent_connections = []
    
    

    # iterate over each row of the subject's geodesic distance 400x400 matrix 
    for j in range(len(geodesic_distances_matrices_corrected[0])):
        
        # retrieve the geodesic distance of each of the most_frequent_connections_indices 
        geodesic_distance_most_frequent_connections = geodesic_distances_matrices_corrected[i][j][most_frequent_connections_indices[j]]

        # compute the mean geodesic distance of the most frequent connections 
        mean_geodesic_distance_most_frequent_connections = np.mean(geodesic_distance_most_frequent_connections)

        # append to respective lists (temp subject list)
        sub_geodesic_distances_most_frequent_connections.append(geodesic_distance_most_frequent_connections)
        sub_mean_geodesic_distance_most_frequent_connections.append(mean_geodesic_distance_most_frequent_connections)
        
        
        # do the same as above but in respective sex list
        if HCP_demographics_cleaned_final.Gender.tolist()[i] == 'M':
            
            # retrieve the geodesic distance of each of the most_frequent_connections_indices for males!
            geodesic_distance_most_frequent_connections = geodesic_distances_matrices_corrected[i][j][most_frequent_connections_indices_M[j]]

            # compute the mean geodesic distance of the most frequent connections 
            mean_geodesic_distance_most_frequent_connections = np.mean(geodesic_distance_most_frequent_connections)

            # append to respective lists (temp subject list)
            M_sub_geodesic_distances_most_frequent_connections.append(geodesic_distance_most_frequent_connections)
            M_sub_mean_geodesic_distance_most_frequent_connections.append(mean_geodesic_distance_most_frequent_connections)
            
            
        else:
            
            # retrieve the geodesic distance of each of the most_frequent_connections_indices for females!
            geodesic_distance_most_frequent_connections = geodesic_distances_matrices_corrected[i][j][most_frequent_connections_indices_F[j]]

            # compute the mean geodesic distance of the most frequent connections 
            mean_geodesic_distance_most_frequent_connections = np.mean(geodesic_distance_most_frequent_connections)

            # append to respective lists (temp subject list)
            F_sub_geodesic_distances_most_frequent_connections.append(geodesic_distance_most_frequent_connections)
            F_sub_mean_geodesic_distance_most_frequent_connections.append(mean_geodesic_distance_most_frequent_connections)
            

            
    # append current subject's lists to lists containing all subjects
    geodesic_distances_most_frequent_connections.append(sub_geodesic_distances_most_frequent_connections)
    mean_geodesic_distances_most_frequent_connections.append(sub_mean_geodesic_distance_most_frequent_connections)
    
    # append curent subject's lists to lists containing only male or female depending on what the subject is
    if HCP_demographics_cleaned_final.Gender.tolist()[i] == 'M':
        geodesic_distances_most_frequent_connections_M.append(M_sub_geodesic_distances_most_frequent_connections)
        mean_geodesic_distances_most_frequent_connections_M.append(M_sub_mean_geodesic_distance_most_frequent_connections)
        
    else:
        geodesic_distances_most_frequent_connections_F.append(F_sub_geodesic_distances_most_frequent_connections)
        mean_geodesic_distances_most_frequent_connections_F.append(F_sub_mean_geodesic_distance_most_frequent_connections)
            
            
            
# make as arrays
geodesic_distances_most_frequent_connections = np.array(geodesic_distances_most_frequent_connections)
mean_geodesic_distances_most_frequent_connections = np.array(mean_geodesic_distances_most_frequent_connections)
geodesic_distances_most_frequent_connections_M = np.array(geodesic_distances_most_frequent_connections_M)
mean_geodesic_distances_most_frequent_connections_M = np.array(mean_geodesic_distances_most_frequent_connections_M)
geodesic_distances_most_frequent_connections_F = np.array(geodesic_distances_most_frequent_connections_F)
mean_geodesic_distances_most_frequent_connections_F = np.array(mean_geodesic_distances_most_frequent_connections_F)






### Export sex differences results matrices in matfile

print("----- Exporting results at /data/p_02667/sex_diff_gradients/results/HCP/connectivity_profiles_geodesic_distances.mat  -----")


mdict = {'geodesic_distances_most_frequent_connections': geodesic_distances_most_frequent_connections, 'mean_geodesic_distances_most_frequent_connections': mean_geodesic_distances_most_frequent_connections, 'most_frequent_connections_indices': most_frequent_connections_indices, 'geodesic_distances_most_frequent_connections_M': geodesic_distances_most_frequent_connections_M, 'mean_geodesic_distances_most_frequent_connections_M': mean_geodesic_distances_most_frequent_connections_M, 'most_frequent_connections_indices_M': most_frequent_connections_indices_M, 'geodesic_distances_most_frequent_connections_F': geodesic_distances_most_frequent_connections_F, 'mean_geodesic_distances_most_frequent_connections_F': mean_geodesic_distances_most_frequent_connections_F, 'most_frequent_connections_indices_F': most_frequent_connections_indices_F}

scipy.io.savemat(resdir_hcp+'connectivity_profiles_geodesic_distances.mat', mdict)


mdict = {'geodesic_distances_matrices_corrected': geodesic_distances_matrices_corrected}

scipy.io.savemat(dataout+'geodesic_distances_matrices.mat', mdict)


