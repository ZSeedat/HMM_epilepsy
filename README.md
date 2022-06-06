# HMM_epilepsy
MATLAB code to automatically detect interictal activity in timeseries data from epilepsy patients using Hidden Marov Modelling. 
MEG data were used from a 275-channel CTF system but any timeseries data could theoretically be used for model inference (including EEG or OPM-MEG).

HMM_TE_inference.m runs the model inference on multivariate MEG data (time x number of channels) and saves the output for use in the source localisations

MORE MATLAB SCRIPTS TO FOLLOW INC. SOURCE LOCALISATION (LCMV BEAMFORMER WITH STATE VARIANCE MAPPING) AND COMPARISON WITH KURTOSIS MAPPING RESULTS
