% Open EEGlab
[ALLEEG, EEG, CURRENTSET] = eeglab;

% Init fieldtrip
ft_defaults;

study_info=init_umdadult_study_info();
subj_id=study_info.participant_info.participant_id{1};

% Path containing subject data
subject_data_dir=fullfile(study_info.deriv_dir, subj_id, 'eeg');

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

% Get C3 channels
c3_idx=find(strcmp(study_info.clusters,'C3'));
channels=study_info.cluster_channels{c3_idx};
merged_EEG = pop_select( merged_EEG,'channel',channels);

% Frequencies to run on. We start at 5Hz because the trials aren't long
% enough to look at lower frequencies over long lags
foi=[5:1:40];

% Lags to run on
lags=[2:.1:7];

% Create fieldtrip data structure
data=create_ft_data(merged_EEG);

% Configuration for frequency analysis
cfg_F             = [];
cfg_F.method      = 'mtmconvol';
cfg_F.taper       = 'hanning';
cfg_F.output      = 'fourier';
cfg_F.keeptrials  = 'yes';
cfg_F.pad         = 'nextpow2';

% Configuration for lagged coherence
cfg_LC            = [];            
cfg_LC.method     = 'laggedcoherence';
cfg_LC.trialsets  = 'all';
fs                = data.fsample;

lagged_coh=zeros(length(channels),length(foi),length(lags));

for f_idx = 1:length(foi);

    % Set freq range
    cfg_F.foi     = foi(f_idx);
    cfg_LC.foi    = foi(f_idx);

    for l_idx=1:length(lags)
        lag=lags(l_idx);
        cfg_LC.lag        = lag;
        cfg_F.width       = lag;

        % width of time windows for phase comparison (seconds)
        width         = cfg_F.width/cfg_F.foi;                        
        cfg_F.t_ftimwin = width;

        % half width of time window (seconds)
        halfwidth     = ceil(fs*width/2)/fs;

        % Go from half window width after trial start to half window
        % width before trial end
        toi_start     = data.time{1}(1) + halfwidth;
        toi_stop      = data.time{1}(end) - halfwidth;

        % Step size
        step          = ceil(fs*cfg_LC.lag/cfg_F.foi)/fs;
        cfg_F.toi     = toi_start:step:toi_stop;

        % Run frequency analysis
        freqout       = ft_freqanalysis(cfg_F,data);

        % Compute lagged coherence
        lcoh = ft_connectivity_laggedcoherence(cfg_LC,freqout);
        lagged_coh(:,f_idx,l_idx)=lcoh.laggedcoh;            
    end
end   

% Plot lagged coherence averaged within C3 electrode cluster
figure();
contourf(lags,foi,squeeze(nanmean(lagged_coh,1)),100, 'linecolor','none');    
colorbar();
ylabel('Frequency (Hz)');
xlabel('Lag (cycles)');