% This script utilises a beamformer function (not specified here, but you 
% could use e.g. FieldTrip) to return a 4mmx4mm map of state variance, for
% each state (as defined by the HMM, see script HMM_TE_inference.m).
% Data are filtered between 20 and 70Hz for weights estimation.
% Reconstructed virtual depth electrode (VDE) timeseries are from data 
% filtered 1-150Hz.
% Zelekha A. Seedat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear variables; clear global;
% Add paths to beamformer functions
addpath('/home/Zelekha/')
% Add fsl path and set environment
addpath '/opt/magres/fsl/Linux_64/fsl/etc/fslconf'
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
% Add filepath with the HMM_Master file in it
addpath(genpath('/net/huail/data_local/Zelekha/Project_2/HMM_MAR_Master'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select patient anatomical files
% Select patient
patient_num = '01';
% 2-minute datasets without large head movement artefacts or SQUID resets,
% there are 12 good runs in the case of this patient:
data_nums = {'02','03','04','05','06','07','08','09','10','11','12','16'};
data_base = '/home/Zelekha/Epilepsy/';

file = strcat(data_base,'patient_',patient_num,'/MRI/');
cd(file)
mrifile = dir(strcat(file,'*coreg.mri')); % Coregistered CTF .mri format
mrifilename = strcat(file,mrifile.name);
Anatfile = dir(strcat(file,'*coreg_brain.nii')); % Coregistered brain extracted using fsl BET
Anatfilename = strcat(file,Anatfile.name);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic parameters
mri_resolution = 0.859375;
f_samp = 1200;                   % MEG sample frequency
duration = 120;                  % MEG trial duration
% Set beamfromer parameters
mu = 0;                          % beamformer regularisation parameter
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load in the anatomical MRI and downsample it
% Load in nifti
brain_3d = cbiReadNifti(Anatfilename); % 256x256x256
NewDim = 64;
Newvox = 256/NewDim; % 4mm voxels
Down_anat = (zeros(NewDim,NewDim,NewDim));
for xdown = 1:NewDim
    for ydown = 1:NewDim
        for zdown = 1:NewDim
            xup = Newvox*xdown-(Newvox/2);
            yup = Newvox*ydown-(Newvox/2);
            zup = Newvox*zdown-(Newvox/2);
            localmat = brain_3d(xup-(Newvox/2 - 1):xup+(Newvox/2 - 1),yup-(Newvox/2 - 1):yup+(Newvox/2 - 1),zup-(Newvox/2 - 1):zup+(Newvox/2 - 1));
            Down_anat(xdown,ydown,zdown) = mean(mean(mean(localmat)));
        end
    end
    disp(sprintf('doing downsample of anatomical: slice %d',xdown));
end
downanat_filename =  strcat(file,'Downsampled_Anat.nii');
cbiWriteNifti(downanat_filename,Down_anat);
mask = Down_anat > 2000;
% Read in the transform from MRI space to CTF-MEG space
trans_struc = rdmrihead_new(mrifilename);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in the MEG data
% Read in each 2-minute run
for i = 1:length(data_nums)
    data_num = data_nums{i};
    filename_2070 = strcat(data_base,'patient_',patient_num,'/preproc_data_ds_',data_num,'_20_70.mat');
    filename_1150 = strcat(data_base,'patient_',patient_num,'/preproc_data_ds_',data_num,'_1_150.mat');
    save_to = strcat(data_base,'patient_',patient_num,'/HMM/ds_',data_num,'/');
    
    % Load preprocessed MEG data.
    % Data have been manually checked for motion artefact and filtered 
    % 20-70Hz for weights estimation:
    data_20_70 = load(filename_2070);
    % 1-150Hz for VDE reconstruction:
    data_1_150 = load(filename_1150);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HMM output
    % Model variables
    n_states = 5; % number of HMM states
    new_freq = 150; % the frequency of the downsampled data
    lags = 5; % window length over which the frequency content of the states was detected, here ~73ms
    T = ones(1,Ntrials).*(duration*f); % trials and durations needed for HMM inference
    
    % Set the HMM Options to be the same as in HMM_TE_inference.m
    % See OHBA's HMM-MAR GitHub wiki page for details of different options
    options = struct();
    options.K = n_states; % 5 states
    options.verbose = 1; 
    options.Fs = f; 
    options.order = 0; 
    options.embeddedlags = -lags:lags; % window span ~73ms
    options.zeromean = 1;
    options.covtype = 'full';
    options.useMEX = 1; 
    options.dropstates = 1; % can output fewer states than asked for where appropriate
    options.DirichletDiag = 10; % diagonal of prior of trans-prob matrix (default 10 anyway)
    options.useParallel = 0; 
    options.standardise = 1; % normalise
    options.pca = 50; % reduce dimensions
    options.downsample = new_freq; % downsample from 1200 to 150Hz
    options.onpower = 0; % Hilbert envelope of signal, not used
            
    % Downsample:
    f_bf = 600; % Beamforming frequency (600Hz rather than 150Hz gives it more data)
    data_20_70 = downsample(data_20_70',f_samp/f_bf)'; % 1200Hz to 600Hz
            
    % Load HMM output
    cd(save_to)
    load hmm.mat
    load Gamma.mat
    
    % Pad gamma to account for loss of timepoints during model inference
    Gamma = padGamma(Gamma,T,options);
    % Actual number of states outputted (in case any states were dropped)
    n_states = size(Gamma,2);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beamform and create 3D a map of state variance for each state
    for state = 1:n_states
        % Interpolate 150Hz state probability timecourse (Gamma) so we can map it onto 600Hz data
        Gamma_600 = interp(Gamma(:,state),size(data_20_70,2)/size(Gamma,1));
        % State defined as active when probability threshold exceeds 2/3
        state_mask = Gamma_600 > 2/3;
        
        % Compute data covariance, regularise and invert
        % Use all the data to create the covariance for weights calculation
        C = cov(data_20_70');
        noise = min(svd(C));
        Noise_Cr = noise.*eye(size(C));
        Cr = C + mu.*max(svd(C)).*eye(size(C));
        inv_Cr = inv(Cr);
                
        % Derive beamformer weights and VDE for each 4mm voxel
        SttVar = zeros(NewDim,NewDim,NewDim); % Create 3D map fo state variance
        Npoints = sum(mask(:)); % number of voxel locations
        pixelcount = 0;
        for xc = 1:NewDim
            for yc = 1:NewDim
            	for zc = 1:NewDim
                	if mask(xc,yc,zc) == 1
                        pixelcount = pixelcount+1;
                        % Coordinate transform
                        mri_coord = [256-(Newvox*xc-(Newvox/2))  256-(Newvox*yc-(Newvox/2))  256-(Newvox*zc-(Newvox/2))];
                        [ctfmeg_coord] = ctfmri2head(mri_coord,trans_struc.T,mri_resolution);
                        ctfmeg_coord = ctfmeg_coord./10;
                                
                        % Beamform to get weights estimation
                        [Weights, lead_fields] = Quick_beamformer(filename_2070,resource,Cr,inv_Cr,Noise_Cr,ctfmeg_coord(1),ctfmeg_coord(2),ctfmeg_coord(3));
                                
                        % Get VE timecourse
                        VE = (Weights'*data_1_150)./sqrt(Weights'*Noise_Cr*Weights);
                        % Get VE data when state is active
                        stt_data = VE.*state_mask';
                        stt_data(stt_data == 0) = [];
                        % Get non-state data:
                        non_stt_mask = ~state_mask;
                        non_stt_data = VE.*non_stt_mask';
                        non_stt_data(non_stt_data == 0) = [];
                        % Variance in state data for this voxel
                        SttVar(xc,yc,zc) = var(stt_data)/var(non_stt_data); % normalised by non-state variance
                        
                        % Print progress
                        if rem(pixelcount,100) == 0
                            clc; disp(sprintf('Done %f percent',100*(pixelcount/Npoints)));
                        end
                    end
                end
            end
        end
        
        % Save state variance map in .nii format
        cd(save_to)
        cbiWriteNifti(strcat('Var_map_state_',num2str(state),'.nii'),SttVar);
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Identify the voxel with peak state variance
        [v, l] = max(SttVar,[],'all','linear');
        [i1,i2,i3] = ind2sub(size(SttVar),l);
        sourcepos = [i1, i2, i3];
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate 1-150Hz VDE for peak location
        % x,y,z coordinates
        xpos = sourcepos(1); ypos = sourcepos(2); zpos = sourcepos(3);
        % Coordinate transform
        mri_coord = [256-(Newvox*xpos-(Newvox/2))  256-(Newvox*ypos-(Newvox/2))  256-(Newvox*zpos-(Newvox/2))];
        [ctfmeg_coord] = ctfmri2head(mri_coord,trans_struc.T,mri_resolution);
        ctfmeg_coord = ctfmeg_coord./10;

        % Beamform
        [Weights, lead_fields] = Quick_beamformer(filename_1150,resource,Cr,inv_Cr,Noise_Cr,ctfmeg_coord(1),ctfmeg_coord(2),ctfmeg_coord(3));
        % Get virtual depth electrode
        peak_VE = (Weights'*data_1_150)./sqrt(Weights'*Noise_Cr*Weights);

        % Save virtual depth electrodes and peak location
        cd(save_to);
        save(strcat('VDE_state_',num2str(state),'.mat'),'peak_VE', 'sourcepos');       
    end
end
