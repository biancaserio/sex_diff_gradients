# Sex differences in functional brain organization

### This is the repository for the publication:
Bianca Serio, Meike D. Hettwer, Lisa Wiersch, Giacomo Bignardi, Julia Sacher, Susanne Weis, Simon B. Eickhoff, & Sofie L. Valk, (2024) **Sex differences in functional cortical organization reflect differences in network topology rather than cortical morphometry**. _Nature Communications_, 15, 7714. https://doi.org/10.1038/s41467-024-51942-1.

Preprint version available [here](https://www.biorxiv.org/content/10.1101/2023.11.23.568437v1)

## Scripts

**1. Preparing the data**
- `calculate_fc_matrices_hcp_schaefer.ipynb` computes functional connectivity matrices in Schaefer 400 parcellation from the HCP BOLD timeseries (averaged across 4 resting state sessions) at the individual level
- `p1_main_function.ipynb` gets the functional connectivity matrices at the individual level
- `p1_main_microstructure.ipynb` gets the microstructural intensity profile matrices at the individual level
- `p1_geodesic_distance.ipynb` gets the geodesic distance of the functional connectivity profiles at the individual level
- `p1_connectivity_profiles_binary_matrices.py` computing binary matrices for the Chi Square test of independence contingency tables
- `p1_connectivity_profiles_strength_fc_sub_level.py` computes functional connectivity strength at the subject level for top 10% fucntional connections


**2. Main analyses**
- `p1_main.ipynb` computes and visualizes main analyses
- `p1_main.R` computes linear mixed effects models for main analyses
- `p1_connectivity_profiles_sex_diff.py`  runs Chi Square test of independence on contingency tables, testing for sex differences in the odds of connections belonging to the seed's top 10% functional connections at the individual level
- `p1_spin_permutation_lmer_within_between_network_dispersion.py` spin permutation test to construct empirical null distribution of beta values for within and between network dispersion analyses

**3. Functions**
- `p1_myfunctions.ipynb` contains functions used for main analyses


## Data
- Sample includes young adults (N=1000) from Human Connectome Project [(HCP)](https://db.humanconnectome.org/) S1200 release
- `Source Data.ods` contains the source data used to make all figures included in the publication


## Support
Please address any questions about the analyses or code to [Bianca Serio](mailto:serio@cbs.mpg.de)

---

### Research poster presented at:
- Annual Meeting of the Organization for Human Brain Mapping (OHBM), Montreal 2023
- Annual meeting of the Organization for the Study of Sex Differences (OSSD), Bergen 2024

![alt text](https://github.com/biancaserio/sex_diff_gradients/blob/master/Poster.png?raw=true)
