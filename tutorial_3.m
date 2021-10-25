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

% Get C3 channels
c3_idx=find(strcmp(study_info.clusters,'C3'));
channels=study_info.cluster_channels{c3_idx};
chan_idx=cellfun(@(x) find(strcmp({base_EEG.chanlocs.labels},x)),...
    channels);

% Get single trial
single_trial_data=squeeze(base_EEG.data(chan_idx,:,1));

% Frequency ranges identified
fois=[9 12; 19.25 22.75];
bands={'alpha','beta'};

% Time stamps in each trial
all_base_times=base_EEG.times;
n_times=length(all_base_times);

% Sampling rate
srate=base_EEG.srate;

% Filtered and amplitude data for each channel in alpha and beta bands
filtered_data=zeros(length(bands),length(chan_idx),n_times);
amp=zeros(length(bands),length(chan_idx),n_times);

for i=1:length(chan_idx)
    chan_data=squeeze(single_trial_data(i,:));
    
    % Get amplitude using the filter-Hilbert method
    % Pad data to avoid edge artifacts
    % Pad with dc offset
    dc=squeeze(mean(chan_data));
    % Pad with 1s on either side
    p_chan_data=[repmat(dc, 1, srate) chan_data repmat(dc, 1, srate)];

    for j=1:length(bands)
        foi=fois(j,:);
        
        % Bandpass filter in frequency range
        % 6th order two-pass Butterworth filter
        % reduce order in case of instability
        f_ch_data = ft_preproc_bandpassfilter(p_chan_data,...
            srate, foi, 6, 'but', 'twopass', 'reduce')';
    
        % Get amplitude from Hilbert transform
        ch_amp=abs(hilbert(f_ch_data));
    
        % Get rid of padding
        f_ch_data=f_ch_data(srate+1:srate+size(chan_data,2));
        ch_amp=ch_amp(srate+1:srate+size(chan_data,2));
    
        filtered_data(j,i,:)=f_ch_data;
        amp(j,i,:)=ch_amp;    
    end
end

% Average over cluster channels
single_trial_data=mean(single_trial_data);
filtered_data=squeeze(mean(filtered_data,2));
amp=squeeze(mean(amp,2));

figure();
hold all;
lgd_labels={'raw'};
plot(all_base_times,single_trial_data);
for i=1:length(bands)
    plot(all_base_times,filtered_data(i,:));
    lgd_labels{end+1}=bands{i};
    plot(all_base_times,amp(i,:));
    lgd_labels{end+1}=sprintf('%s amplitude',bands{i});
end
xlim(all_base_times([1 end]));
legend(lgd_labels);
xlabel('Time (ms)');
ylabel('Potential \muV');
