% Script to find epileptiform activity in SENSOR SPACE MEG data using a
% time-delay embedded Hidden Markov Model.
% Note that data were acquired on a 275-channel CTF-MEG system but in 
% theory any multivariate timeseries data would work. Data have been read
% into MATLAB and preprocessed before being loaded into this script.
% Data were acquired in 2-minute runs because paediatric patients struggle
% to remain still for long periods of time, which means that the position
% of the head relative to the MEG sensors will change. To avoid this
% corrupting the spatial encoding of state activity in the multivariate HMM 
% output, model inference is run on each 2-minute run independently.
% HMM toolbox downloaded from https://github.com/OHBA-analysis/HMM-MAR
% Zelekha A. Seedat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear variables; clear global;
% Add filepaths for HMM functions
addpath('/home/Zelekha/')
addpath(genpath('/home/Zelekha/HMM_MAR_Master/HMM-MAR-master'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Enter filenames
% Select patient:
patient_num = '01';
% 2-minute datasets without large head movement artefacts or SQUID resets,
% there are 12 good runs in the case of this patient:
data_nums = {'02','03','04','05','06','07','08','09','10','11','12','16'};

% Compute model inference for each run individually since there may be head
% motion between runs which could corrupt the spatial pattern of HMM
% epileptiform state activity across sensors. 
for i = 1:length(data_nums)
    data_num = data_nums{i};
    data_base = '/home/Zelekha/Data/Epilepsy/';
    filename = strcat(data_base,'patient_',patient_num,'/preproc_data_ds_',data_num,'_20_70.mat');
    save_to = strcat(data_base,'patient_',patient_num,'/HMM/ds_',data_num,'/');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Parameters
    f = 1200;                   % sample frequency
    duration = 120;             % trial duration
    Ntrials = 1;                % a single 2-minute long trial
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in the preprocessed MEG data
    % Data have been manually checked for motion artefact, filtered between
    % 20-70Hz and notch-filtered at 60Hz (mains frequency)
    data = load(filename);
    % Data will also be downsampled and have dimensionality reduced during
    % the model infeerence step - see options below.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Train HMM on all the channels (multivariate)
    % Model variables
    n_states = 5; % number of HMM states
    new_freq = 150; % the frequency of the downsampled data
    lags = 5; % window length over which the frequency content of the states is detected, here ~73ms
    T = ones(1,Ntrials).*(duration*f); % trials and durations needed for HMM inference
    
    % Set the HMM Options
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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run HMM function
    disp(['Initiating model inference for dataset ',data_nums{i}])
    [hmm, Gamma] = hmmmar(data,T,options);
    % Actual number of states outputted
    n_states = size(Gamma,2);
    disp(['Dataset ',data_nums{i},' output: ',num2str(n_states),' states'])
    
    % Save the output
    if exist(save_to)
        cd(save_to)
    else
        mkdir(save_to)
        cd(save_to)
    end
    
    save hmm.mat hmm
    save Gamma.mat Gamma
end
