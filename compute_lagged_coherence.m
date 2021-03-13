function compute_lagged_coherence(study_info, varargin)
% COMPUTE_LAGGED_COHERENCE - Computes lagged coherence in the C3 and C4 
% clusters across a range of frequencies and lags
%
% Syntax:  compute_lagged_coherence(study_info)
%
% Inputs:
%    study_info - structure containing the study information
%
% Example: 
%   compute_lagged_coherence(study_info)

% Parse optional arguments
defaults=struct();
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

% Frequencies to run on. We start at 5Hz because the trials aren't long
% enough to look at lower frequencies over long lags
foi=[5:.5:40];

% Lags to run on
lags=[2:.5:7];

% Electrode clusters to look in
clusters={'C3','C4'};
cluster_channels={{'E16', 'E20', 'E21', 'E22'},...
    {'E41', 'E49', 'E50', 'E51'}};

% Number of subjects
n_subjects=size(study_info.participant_info,1);

% Remove full fieldtrip path to avoid conflicts with EEGLab plugin
rmpath(genpath('C:\Users\jbonaiuto\Dropbox (Personal)\Toolboxes\fieldtrip-20190329'));

% Open EEGlab
[ALLEEG, EEG, CURRENTSET] = eeglab;

% Remove EEGlab fieldtrip plugin
rmpath(genpath('C:\Users\jbonaiuto\Documents\MATLAB\eeglab14_1_1b\plugins\Fieldtrip-lite20210304'));
% Add full fieldtrip
addpath('C:\Users\jbonaiuto\Dropbox (Personal)\Toolboxes\fieldtrip-20190329');
ft_defaults

lagged_coh=zeros(n_subjects,length(clusters),length(foi),length(lags));

for s=1:n_subjects
    
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s};
    
    % Path containing subject data
    subject_data_dir=fullfile(study_info.output_dir, subj_id, 'eeg');
    
    % Baseline and experimental epoch files
    base_fname=sprintf('%s_11_Epoch_Matched_CSD_baseline.set',subj_id);
    exp_fname=sprintf('%s_11_Epoch_Matched_CSD_experimental.set',subj_id);
    
    % Load data
    base_EEG=pop_loadset('filepath', subject_data_dir,...
        'filename', base_fname);        
    [ALLEEG, EEG] = eeg_store(ALLEEG, base_EEG, 1);        
    exp_EEG=pop_loadset('filepath', subject_data_dir,...
        'filename', exp_fname);  
    [ALLEEG, EEG] = eeg_store(ALLEEG, exp_EEG, 2);
    
    % Merge baseline and experimental EEG data
    merged_EEG=pop_mergeset(ALLEEG,1:2,0);

    % Data averaged within each cluster
    cluster_data=zeros(length(clusters), merged_EEG.pnts, merged_EEG.trials);

    for c_idx=1:length(clusters)

        % Channels in this cluster
        channels=cluster_channels{c_idx};

        % Find indices of cluster channels
        chan_idx=zeros(1,length(channels));
        for k=1:length(channels)
            chan_idx(k)=find(strcmp({merged_EEG.chanlocs.labels},channels{k}));
        end    

        % Average over channels in this cluster
        cluster_data(c_idx,:,:)=double(squeeze(mean(merged_EEG.data(chan_idx,:,:),1)));
    end
    
    % Get lagged coherence
    data=create_ft_data(clusters, cluster_data, merged_EEG.times,...
        merged_EEG.srate);
    
    for l_idx=1:length(lags)
        lag=lags(l_idx);

        % Configuration for frequency analysis
        cfg_F             = [];
        cfg_F.method      = 'wavelet';  % using Morlet wavelets
        cfg_F.width       = lag;        % number of wavelet cycles
                                        % according to lag - same as
                                        % Fourier transform of
                                        % Hanning-tapered signal
        cfg_F.output      = 'fourier';  % return the complex Fourier-
                                        % spectra because we need the phase
        cfg_F.keeptrials  = 'yes';      % return individual trials
        cfg_F.pad         = 'nextpow2'; % data padding
        
        % Configuration for lagged coherence
        cfg_LC            = [];                
        cfg_LC.lag        = lag;               % lag to look at
        cfg_LC.method     = 'laggedcoherence'; 
        cfg_LC.trialsets  = 'all';
        
        for f_idx = 1:length(foi);
            
            % Set freq range
            cfg_F.foi     = foi(f_idx);
            cfg_LC.foi    = foi(f_idx);
            
            % width of time windows for phase comparison (seconds)
            width         = cfg_LC.lag/cfg_F.foi;
            % width of time windows for phase comparison (data points)
            width_pts     = merged_EEG.srate*width;
            % Go from half window width after trial start to half window
            % width before trial end
            toi_start     = data.time{1}(1) + ceil(width_pts/2)/merged_EEG.srate;
            toi_stop      = data.time{1}(end) - ceil(width_pts/2)/merged_EEG.srate;
            cfg_F.toi     = toi_start:width:toi_stop;
            
            % Run frequency analysis
            freqout       = ft_freqanalysis(cfg_F,data);
            % Run lagged coherence
            lcoh = ft_connectivity_laggedcoherence(cfg_LC,freqout);
            lagged_coh(s,:,f_idx,l_idx)=lcoh.laggedcoh;            
        end        
    end    
end

save('../../../data/UMD/adult/deriv/lagged_coherence.mat','lags','foi','lagged_coh');


figure();
for c_idx=1:length(clusters)
    subplot(length(clusters),1,c_idx);
    contourf(lags,foi,squeeze(mean(lagged_coh(:,c_idx,:,:))),100,'linecolor','none');
    %set(gca,'clim',[-.08 .04]);
    colorbar();
    xlabel('Lag (cycles)');
    ylabel('Frequency (Hz)');
end
