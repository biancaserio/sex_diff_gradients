#!/usr/bin/env python3


# Runs spin permutation with for within and between network dispersion analyses with lmer controling for ICV


### Define parameters

surface_name='fsa5'
parcellation_name='schaefer_400'
n_rot=1000


### Load required packages
import os
import sys
import numpy as np
import pandas as pd

from enigmatoolbox.permutation_testing import centroid_extraction_sphere
from enigmatoolbox.permutation_testing import rotate_parcellation
import statsmodels.regression.mixed_linear_model as sm
from enigmatoolbox.datasets import load_fsa5



### Define directories

datadir = '/data/p_02667/sex_diff_gradients/data/'
resdir_hcp = '/data/p_02667/sex_diff_gradients/results/HCP/'



### Load the required data

print(f'\nLoad the required data')


# functional G1 hemi alligned gradient array
hemi_array_aligned_fc_G1 = np.genfromtxt(resdir_hcp+'hemi_array_aligned_fc_G1.csv', delimiter=',', skip_header=1)

# demographics data
demographics_cleaned_final = pd.read_csv(resdir_hcp+'demographics_cleaned_final.csv')




# Yeo network array
# labels: 1=visual, 2=sensory motor, 3=dorsal attention, 4=ventral attention, 5=limbic, 6=fronto parietal, 7= DMN
yeo7_networks_array = np.genfromtxt(datadir+'yeo_7.csv', delimiter=',', skip_header=0)



# root_pth = os.path.dirname(__file__)  

# root_pth via os.pathdirname(__file__) did not work so me so manually directed to path to downlowaded enigmatoolbox 
# ie via 3 commands (see https://enigma-toolbox.readthedocs.io/en/latest/pages/01.install/index.html): 
## git clone https://github.com/MICA-MNI/ENIGMA.git
## cd ENIGMA
## python setup.py install
root_pth = '/data/p_02667/sex_diff_gradients/scripts/ENIGMA/enigmatoolbox/permutation_testing/'



# load fsa5 surface
sphere_lh, sphere_rh = load_fsa5(as_sphere=True)

# get sphere coordinates of parcels
annotfile_lh = os.path.join(root_pth, 'annot', surface_name + '_lh_' + parcellation_name + '.annot')
annotfile_rh = os.path.join(root_pth, 'annot', surface_name + '_rh_' + parcellation_name + '.annot')

lh_centroid = centroid_extraction_sphere(sphere_lh.Points, annotfile_lh)
rh_centroid = centroid_extraction_sphere(sphere_rh.Points, annotfile_rh)




### Generate permutation maps (rotated centroids n_rot times)

print(f'\nGenerate permutation maps')

perm_id = rotate_parcellation(lh_centroid, rh_centroid, n_rot)

nroi = perm_id.shape[0]   # number of regions
nperm = perm_id.shape[1]  # number of permutations




### Permutation of yeo7 network labels

print(f'\nPermutate Yeo network labels')

yeo7_networks_array_perm = np.empty((0, nroi))

for rr in range(nperm):
    
    print(f'---- Permuation No. {rr+1} ----')
    
    yeo7_networks_array_perm2 = []
    
    for ii in range(nroi):
        yeo7_networks_array_perm2 = np.append(yeo7_networks_array_perm2, yeo7_networks_array[int(perm_id[ii, rr])])

    yeo7_networks_array_perm = np.vstack((yeo7_networks_array_perm, yeo7_networks_array_perm2))




### Compute within and between network dispersion & sex contrast -- same code as for original, only differences: 
# WN: don't take the real yeo7_networks_array labels, but iterate over the permuted ones in yeo7_networks_array_perm)
# BN: how data is saved

print(f'\nCompute within and between network dispersion - as well as sex contrast - for {n_rot} permutations')

# lists that will contain the permutation t-values (null distribution of t-values) for WN and BN dispersion sex contrasts
WN_dispersion_perm_tval_sex_contrast = []
BN_dispersion_perm_tval_sex_contrast = []



# iterate over the i permutations of the yeo network labels
for i in range(len(yeo7_networks_array_perm)):
    
    print(f'\n---- Permuation No. {i+1} ----')
    
    
    print(f'---- Computed WN dispersion ----')
    
    yeo_cog_perm = []  # center of gravity (median) for each network (i.e., network centroid position) (7 x 1000)
    WN_dispersion_perm = []  # Within network dispersion: sum squared Euclidean distance of network nodes to the network centroid at individual level (7 x 1000)


    # gradient values
    g1 = hemi_array_aligned_fc_G1.T  # transpose to obtain shape (400 x 1000) in order to access/index the relevant network nodes


    # list that will contain the current permuation's t values (len 7 because will contain all 7 networks)
    WN_dispersion_perm_tval_sex_contrast_temp = []
    
    # iterate over the 7 Yeo networks
    for n in range(7):
        
        print(f'Network: {n+1}')

        # identify the nodes of given Yeo network
        netNodes = np.where(yeo7_networks_array_perm[i] == (n + 1))  # here we take the ith permutation of the yeo network array labels !!
        netNodes = np.squeeze(netNodes)

        # get the gradient loadings of the nodes of the given Yeo network, for each subject (shape: number of nodes in network x N)
        G1_net = g1[netNodes]


        ### identify the centroid / center of gravity (= median) of the given Yeo network for each subject (shape: N)
        yeo_cog_perm_net = np.median(G1_net, axis=0)  
        yeo_cog_perm.append(yeo_cog_perm_net)


        ### within network dispersion: 1 within network dispersion value per subject (per network)

        # compute (per subject) the Eucledian distance between each gradient loading (in Yeo network) and that network's centroid
        dist_nodes_to_cog = G1_net - yeo_cog_perm_net  # shape: number of nodes in network x N

        # take the sum of squares of this Eucledian distance 
        sum_of_squares = np.sum((dist_nodes_to_cog**2), axis = 0)  # shape: N

        # append to list
        WN_dispersion_perm.append(sum_of_squares)


        
        ### compute linear model for current network (WN dispersion)
        print(f'lmer in progress...')
        
        # make a dataframe that will contain the data for linear mixed effects model (only change is the WN dispersion metric (1-7 networks)
        df = pd.DataFrame({'WN_dispersion': sum_of_squares, 'sex': demographics_cleaned_final.Gender, 'age': demographics_cleaned_final.Age_in_Yrs, 'ICV': demographics_cleaned_final.FS_IntraCranial_Vol, 'Family_ID': demographics_cleaned_final.Family_ID, 'TwinStatus': demographics_cleaned_final.TwinStatus})

        # define model
        model = sm.MixedLM.from_formula(
            formula = "WN_dispersion ~ 1 + sex + age + ICV", 
            data = df, 
            re_formula="1",
            groups="Family_ID",
            vc_formula={"TwinStatus": "0 + C(TwinStatus)"})  # random effect of TwinStatus nested in Family_ID

        # fit model
        results = model.fit()
        #results.summary()

        # save results
        WN_dispersion_perm_tval_sex_contrast_temp.append(results.tvalues[1])
        
        
    # append to permutation results
    WN_dispersion_perm_tval_sex_contrast.append(WN_dispersion_perm_tval_sex_contrast_temp)

    
    
    
    ### Compute between network dispersion (Eucledian distance between two network centroids) only once for each pariwise network comparison, whilst keeping track what networks are being compared
    
    print(f'---- Computed BN dispersion ----')
    
    # list that will contain the current permuation's t values (len 21 because will contain all 21 comparisons of networks)
    BN_dispersion_perm_tval_sex_contrast_temp = []
    
    # to keep track the order in which the pairwise comparisons are computed
    pairwise_comparison_order = []

    # iterate over 7 Yeo networks
    for n1 in range(7):
    
        for n2 in range(7):

            current_pairwise_comparison = [n1+1, n2+1]


            if n1!=n2 and [n1+1, n2+1] not in pairwise_comparison_order and list(reversed([n1+1, n2+1])) not in pairwise_comparison_order:
                
                print(f'Networks: {current_pairwise_comparison}')

                # append the pairwise comparison to dict
                pairwise_comparison_order.append(current_pairwise_comparison)

                # compute the distance between centroid of network 1 and centroid of network 2 and append it to dict
                distance = yeo_cog_perm[n1] - yeo_cog_perm[n2]
                


                ### compute linear model for current pairwise between network dispersion (distance between centroids)
                print(f'lmer in progress...')

                # make a dataframe that will contain the data for linear mixed effects model (only change is the WN dispersion metric (1-7 networks)
                df = pd.DataFrame({'BN_dispersion': distance, 'sex': demographics_cleaned_final.Gender, 'age': demographics_cleaned_final.Age_in_Yrs, 'ICV': demographics_cleaned_final.FS_IntraCranial_Vol, 'Family_ID': demographics_cleaned_final.Family_ID, 'TwinStatus': demographics_cleaned_final.TwinStatus})

                # define model
                model = sm.MixedLM.from_formula(
                    formula = "BN_dispersion ~ 1 + sex + age + ICV", 
                    data = df, 
                    re_formula="1",
                    groups="Family_ID",
                    vc_formula={"TwinStatus": "0 + C(TwinStatus)"})  # random effect of TwinStatus nested in Family_ID

                # fit model
                results = model.fit()
                #results.summary()

                # save results
                BN_dispersion_perm_tval_sex_contrast_temp.append(results.tvalues[1])
                
    
    # append to permutation results
    BN_dispersion_perm_tval_sex_contrast.append(BN_dispersion_perm_tval_sex_contrast_temp)
    
    
    
    
### Clean results for export

print(f'\n---- Export results at /data/p_02667/sex_diff_gradients/results/HCP/WN and BN_dispersion_perm_tval_sex_contrast_null_distr ----')

# pack into arrays
WN_dispersion_perm_tval_sex_contrast = np.array(WN_dispersion_perm_tval_sex_contrast)
BN_dispersion_perm_tval_sex_contrast = np.array(BN_dispersion_perm_tval_sex_contrast)


# contains the null distribution of t values for the sex contrast on the Within Network dispersion per network (7)
WN_dispersion_perm_tval_sex_contrast_null_distr = {'visual': WN_dispersion_perm_tval_sex_contrast.T[0],
                                                   'sensory motor': WN_dispersion_perm_tval_sex_contrast.T[1], 
                                                   'dorsal attention': WN_dispersion_perm_tval_sex_contrast.T[2], 
                                                   'ventral attention': WN_dispersion_perm_tval_sex_contrast.T[3], 
                                                   'limbic': WN_dispersion_perm_tval_sex_contrast.T[4], 
                                                   'fronto parietal': WN_dispersion_perm_tval_sex_contrast.T[5], 
                                                   'DMN': WN_dispersion_perm_tval_sex_contrast.T[6]}


# contains the null distribution of t values for the sex contrast on the Between Network pairwise comparisons (21)
# column names is the pairwise comparison labels (turned into strings)
BN_dispersion_perm_tval_sex_contrast_null_distr = pd.DataFrame(BN_dispersion_perm_tval_sex_contrast, columns = [str(pair) for pair in pairwise_comparison_order])  


# export 
pd.DataFrame(WN_dispersion_perm_tval_sex_contrast_null_distr).to_csv(resdir_hcp+'WN_dispersion_perm_tval_sex_contrast_null_distr.csv', header = True, index = False)
BN_dispersion_perm_tval_sex_contrast_null_distr.to_csv(resdir_hcp+'BN_dispersion_perm_tval_sex_contrast_null_distr.csv', header = True, index = False)



    
    
