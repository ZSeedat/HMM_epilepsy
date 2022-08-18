% Script to plot the HMM output so that a user can manually select the
% epileptiform state. 
% Both spatial maps and virtual depth electrodes are plotted as well as a
% state transition probability matrix which provides the likelihood of
% transitionsing from any one state to another (this could be useful where
% 2 epileptifrom states exist to see if one focus leads the other).
% Zelekha A. Seedat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear variables; clear global;
% Add filepath with the HMM_Master file in it
addpath(genpath('/home/Zelekha/HMM_MAR_Master'))
% Add path with freezeColors function on it
% This function enables overlay of functional data on structural MRI. You
% can download it from here: 
% https://uk.mathworks.com/matlabcentral/fileexchange/7943-freezecolors-unfreezecolors
addpath('/home/Zelekha')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select patient data
% Select patient:
patient_num = '01';
% 2-minute datasets without large head movement artefacts or SQUID resets,
% there are 12 good runs in the case of this patient:
data_nums = {'02','03','04','05','06','07','08','09','10','11','12','16'};

data_base = '/home/Zelekha/Data/Epilepsy/';
f_samp = 1200; % sample frequency
duration = 120; % trial duration
f_bf = 600; % 600Hz beaform frequency

% Get anatomical file:
anat_filepath = strcat(data_base,'patient_',patient_num,'/MRI/');

% Select the epileptiform state(s) for each dataset
for i = 1:length(data_nums) % Dataset number
    data_num = data_nums{i};
    save_to = strcat(data_base,'patient_',patient_num,'/HMM/ds_',data_num,'/');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load HMM output
    % Variables
    n_states = 5; % number of HMM states
    new_freq = 150; % the frequency of the downsampled data
    lags = 5; % window length over which the frequency content of the states was detected, here ~73ms
    Ntrials = 1; % single run, 120s long
    T = ones(1,Ntrials).*(duration*f_samp);
    
    % Set the HMM Options
    options = struct(); % Create options struct
    options.K = n_states;
    options.verbose = 1;
    options.Fs = f_samp;
    options.order = 0;
    options.embeddedlags = -lags:lags; % window span ~73ms
    options.zeromean = 1;
    options.covtype = 'full';
    options.useMEX = 1; % runs much faster
    options.dropstates = 1; % can output fewer states than asked for where appropriate
    options.DirichletDiag = 10; % diagonal of prior of trans-prob matrix (default 10 anyway)
    options.useParallel = 0;
    % Preprocessing parameters
    options.standardise = 1; % normalise
    options.pca = 50; % reduce dimensions
    options.downsample = new_freq; % downsample from 1200 to 300Hz
    options.onpower = 0; % Hilbert envelope of signal.
    
    % Load HMM output
    cd(save_to)
    load hmm.mat
    load Gamma.mat
    
    % Pad gamma to account for loss of timepoints during model inference
    Gamma = padGamma(Gamma,T,options);
    % Actual number of states outputted (in case any states were dropped)
    n_states = size(Gamma,2);
    
    % Get state transition probability matrix, i.e. the probability of
    % transitioning from any one state to any other
    trans_mat = getTransProbs(hmm);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot state varinace map
    % Load anatomical MRI to map state variance map onto
    cd (anat_filepath);
    anat_ds = niftiread('downsampled_brain.nii');
    
    Gamma_all = [];
    for state = 1:n_states
        cd(save_to);
        % Load in peak voxel
        load(strcat('VDE_state_',num2str(state),'.mat'));
        % Load state variance map
        stt_map = niftiread(['Var_map_state_',num2str(state),'.nii']);
        % Get peak voxel location
        [v, l] = max(stt_map,[],'all','linear');
        [i1,i2,i3] = ind2sub(size(stt_map),l);
        sourcepos = [i1, i2, i3];
        
        % Threshold the state variance map
        % Can alter the threshold value for visualisation
        thresh_stt_map = stt_map;
        thresh_stt_map(stt_map < 0.6*max(stt_map(:))) = nan;
        
        figure('Color','w','units','normalized','Position',[0.1 0.3 0.55 0.45]);
        % Axial slice with peak voxel visible
        subplot(2,4,1);
        I = squeeze(anat_ds(i1,:,:)); Ir = imrotate(I,90);
        h = pcolor(Ir); colormap('bone');
        caxis([min(anat_ds(:)) 0.4*max(anat_ds(:))]); freezeColors; hold on;
        set(h,'EdgeColor','none'); set(gca,'YDir','rev','YTickLabel',[],'XTickLabel',[])
        I = squeeze(thresh_stt_map(i1,:,:)); Ir = imrotate(I,90);
        h = pcolor(Ir); colormap('hot');
        caxis([min(thresh_stt_map(:)) max(thresh_stt_map(:))]); freezeColors;
        set(h,'EdgeColor','none'); set(gca,'YDir','rev','YTickLabel',[],'XTickLabel',[])
        
        % Coronal slice with peak voxel visible
        subplot(2,4,2);
        I = squeeze(anat_ds(:,i2,:)); Ir = imrotate(I,90);
        h = pcolor(Ir); colormap('bone');
        caxis([min(anat_ds(:)) 0.4*max(anat_ds(:))]); freezeColors; hold on;
        set(h,'EdgeColor','none'); set(gca,'YDir','rev','YTickLabel',[],'XTickLabel',[])
        I = squeeze(thresh_stt_map(:,i2,:)); Ir = imrotate(I,90);
        h = pcolor(Ir); colormap('hot');
        caxis([min(thresh_stt_map(:)) max(thresh_stt_map(:))]); freezeColors;
        set(h,'EdgeColor','none'); set(gca,'YDir','rev','YTickLabel',[],'XTickLabel',[])
        
        % Sagittal slice with peak voxel visible
        subplot(2,4,3);
        I = squeeze(anat_ds(:,:,i3)); Ir = imrotate(I,90);
        h = pcolor(Ir); colormap('bone');
        caxis([min(anat_ds(:)) 0.4*max(anat_ds(:))]); freezeColors; hold on;
        set(h, 'EdgeColor', 'none'); set(gca,'YDir','rev','YTickLabel',[],'XTickLabel',[])
        I = squeeze(thresh_stt_map(:,:,i3)); Ir = imrotate(I,90);
        h = pcolor(Ir); colormap('hot');
        caxis([min(thresh_stt_map(:)) max(thresh_stt_map(:))]); freezeColors;
        set(h, 'EdgeColor', 'none'); set(gca,'YDir','rev','YTickLabel',[],'XTickLabel',[])
        
        % Plot transition probability matrix
        subplot(2,4,4); imagesc(trans_mat); title('Transition Probabilities')
        xlabel('state at t+1'); ylabel('state at t'); set(gcf,'Color','w')
        h = colorbar; ylabel(h,'probability'); colormap('parula')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot peak voxel timecourse from state variance map       
        % Interpolate Gamma so we can map it onto 600Hz data
        Gamma_600 = interp(Gamma(:,state),f_bf/new_freq);
        state_mask = Gamma_600 > 2/3;
        % Get state data
        state_data = state_mask'.*peak_VE;
        state_data(state_data == 0) = nan; % for plotting purposes
        
        % Plot peak VDE
        x = (1:length(state_data))./f_bf;
        subplot(2,4,[5,6,7,8]);
        plot(x,peak_VE); hold on;
        plot(x,state_data); hold off;
        ylabel({'VDE';'A.U.'}); xlabel('Time, s'); title(['State ',num2str(state)]);
        set(gca,'FontSize',12); xlim([0 10])       
    end
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User defines epi state from graphs above
    epi_state = input('Please input the epilepsy state(s), use square brackets []. ');
    epi_states{i} = epi_state;
    close all;
end

% Save out epi states for further analysis
cd([save_to,'../'])
save('epi_states.mat','epi_states')



