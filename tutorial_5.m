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

% Relative threshold
rel_thresh=1.6;

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

% Average amplitude over all trials
cluster_amp=squeeze(mean(cluster_amp));
             
% Compute absolute burst threshold
threshold=median(cluster_amp(:))+(rel_thresh*std(cluster_amp(:)));

bursts=[];
% For each burst - trial it occured in
bursts.trial=[];
% For each burst - time within trial it occured
bursts.peak_time=[];
% For each burst - onset time within trial
bursts.onset_time=[];
% For each burst - offset time within trial
bursts.offset_time=[];
% Burst peak amplitude
bursts.peak_amp=[];

% Go through signal each trial 
for t_idx=1:trials

    % All times when amp is over threshold
    over_thresh_idx = cluster_amp(:,t_idx)>=threshold;
    % Change in threshold crossing
    over_thresh_diff = diff(over_thresh_idx);
    % All times when amp is over threshold and previous time is
    % under threshold
    all_start_idx=find(over_thresh_diff)+1;

    %% Process each burst
    for k=1:length(all_start_idx)
        start_idx=all_start_idx(k);

        % Find next time amplitude goes below threshold
        end_idx=find(cluster_amp(start_idx:end,t_idx)<threshold,1);
        
        % If it goes below threshold before the end of the trial
        if ~isempty(end_idx) && end_idx>2
            
            % Time when amplitude goes back below threshold
            end_idx=start_idx+end_idx-1;
        
            % Find peak amplitude
            burst_amp=cluster_amp(start_idx:end_idx,t_idx);
            [max_burst_amp,burst_peak_idx]=max(burst_amp);
            burst_peak_idx=burst_peak_idx+start_idx-1;

            % Save burst information
            bursts.trial(end+1)=t_idx;
            bursts.peak_time(end+1)=all_times(burst_peak_idx);
            bursts.onset_time(end+1)=all_times(start_idx);
            bursts.offset_time(end+1)=all_times(end_idx);
            bursts.peak_amp(end+1)=max_burst_amp;                        
        end
    end
end

figure();
subplot(2,1,1);
hold all;
for i=1:trials
    t_bursts=find(bursts.trial==i);
    plot(bursts.peak_time(t_bursts),i*ones(1,length(t_bursts)),'k.');
end
ylabel('Trial');
subplot(2,1,2);
plot(all_times,mean(cluster_amp,2));
xlabel('Time (ms)');