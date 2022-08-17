# HMM_epilepsy
MATLAB code to automatically detect interictal activity in timeseries data from epilepsy patients using Hidden Marov Modelling. 
MEG data were used from a 275-channel CTF system but any timeseries data could theoretically be used for model inference (including EEG or OPM-MEG).

HMM_TE_inference.m runs the model inference on multivariate MEG data (number of timepoints x number of channels) and saves the output for use with beamformer source reconstruction to localise voxels in the brain where the variance is greatest when the state is active.

HMM_state_source_localisation.m uses a regular LCMV beamformer together with the HMM output to create a 3D map of state variance (highlighting regions of the brain where virtual depth electorde (VDE) variance is high when the state is active) and a VDE timecourse from the peak voxel.



MORE MATLAB SCRIPTS TO FOLLOW INCLUDING COMPARISON WITH KURTOSIS MAPPING RESULTS
