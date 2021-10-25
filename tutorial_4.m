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

% Sampling rate
srate=merged_EEG.srate;

% Get C3 channels
c3_idx=find(strcmp(study_info.clusters,'C3'));
channels=study_info.cluster_channels{c3_idx};
chan_idx=cellfun(@(x) find(strcmp({merged_EEG.chanlocs.labels},x)),...
    channels);

% Frequency range
foi=[19.25 22.75];

% Time stamps in each trial
all_times=merged_EEG.times;
n_times=length(all_times);

% Number of trials
trials=merged_EEG.trials;

% Compute amplitude in each channel of C3 cluster
cluster_amp=zeros(length(chan_idx),n_times,trials);

for i=1:length(chan_idx)
    chan_data=squeeze(merged_EEG.data(chan_idx(i),:,:));
                
    % Get amplitude using the filter-Hilbert method
    % Pad data to avoid edge artifacts
    % Pad with dc offset
    dc=squeeze(mean(chan_data));
    % Pad with 1s on either side
    p_chan_data=[repmat(dc, srate, 1); chan_data; repmat(dc, srate, 1)];

    % Bandpass filter in frequency range
    % 6th order two-pass Butterworth filter
    % reduce order in case of instability
    f_ch_data = ft_preproc_bandpassfilter(p_chan_data',...
        srate, foi, 6, 'but', 'twopass', 'reduce')';
    
    % Get amplitude from Hilbert transform
    ch_amp=abs(hilbert(f_ch_data));
    
    % Get rid of padding
    cluster_amp(i,:,:)=ch_amp(srate+1:srate+n_times,:);
end

% Average amplitude over cluster channels
cluster_amp=squeeze(mean(cluster_amp));
             
% Compute correlation between number of bursts and mean amplitude in each
% trial for a range of standard deviations above the median
sds=[.1:.1:3];
subj_corr=zeros(1,length(sds));
for sd_idx=1:length(sds)
    
    % Relative burst threshold
    thresh_sd=sds(sd_idx);

    % Compute absolute burst threshold
    threshold=median(cluster_amp(:))+(thresh_sd*std(cluster_amp(:)));

    % Number of bursts in each trial
    trial_bursts=zeros(1,trials);
    
    % Go through amplitude each trial 
    for t_idx=1:trials
        % 1 when amp is over threshold, 0 otherwise
        over_idx = cluster_amp(:,t_idx)>=threshold;
        % Change in over threshold gives threshold crossing
        thresh_crossing = find(diff(over_idx));
        trial_bursts(t_idx)=length(thresh_crossing);                
    end
    
    % Nonparametric correlation between nuber of bursts in each
    % trial and mean amplitude
    subj_corr(sd_idx)=corr(trial_bursts', mean(cluster_amp)',...
        'type', 'Spearman');
end

% Find SD that maximizes mean correlation
[max_corr,max_corr_idx]=max(subj_corr);
max_corr_sd=sds(max_corr_idx);
disp(sprintf('threshold=%.2f SDs, correlation=%.3f\n',max_corr_sd,...
    max_corr));

threshold=median(cluster_amp(:))+(max_corr_sd*std(cluster_amp(:)));

figure();
subplot(2,1,1);
hold all;
% Bins to compute amplitude density in
bins=linspace(0,100,101);
bins=bins(2:end);
% Plot amplitude density
density=ksdensity(cluster_amp(:),bins);
plot(bins,density);
plot([threshold threshold],ylim(),'r--');
xlabel('Amplitude (uV)');
ylabel('Density');
subplot(2,1,2);
hold all;
plot(sds,subj_corr);
plot([max_corr_sd max_corr_sd],ylim(),'r--');
xlabel('Threshold (SDs above median)');
ylabel('\rho');