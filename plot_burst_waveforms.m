function plot_burst_waveforms(study_info, foi, woi, thresh_sd, varargin)
% PLOT_BURST_WAVEFORMS - Plot burst waveforms in the C3 and C4 clusters
%
% Syntax:  plot_burst_waveforms(study_info, foi, woi, thresh_sd)
%
% Inputs:
%    study_info - structure containing the study information
%    foi - frequency band of interest ([low high], Hz)
%    woi - time window of interest (to cut off filtering edge artifacts;
%        [low high], ms)
%    thresh_sd - number of standard deviations above the median amplitude
%        to set the threshold at 
%
% Optional keyword inputs:
%     burst_window - width of time window around burst to extract waveform
%         (ms; default = 100)
%
% Example: 
%   plot_burst_waveforms(study_info, [13 30], [-1000 1000], 1.5)

% Parse optional arguments
defaults=struct('burst_window',200,'align',true);
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

% Electrode clusters to look in
clusters={'C3','C4'};
cluster_channels={{'E16', 'E20', 'E21', 'E22'},...
    {'E41', 'E49', 'E50', 'E51'}};
cluster_colors={'b','r'};

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

% Mean burst waveforms for each subject
subj_waveforms=[];

for s=1:n_subjects
    
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s};
    
    % Path containing subject data
    subject_data_dir=fullfile(study_info.output_dir, subj_id, 'eeg');
    
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

    % Time stamps in each trial
    all_times=merged_EEG.times;
    
    for c_idx=1:length(clusters)

        % Channels in this cluster
        channels=cluster_channels{c_idx};

        % Find indices of cluster channels
        chan_idx=zeros(1,length(channels));
        for k=1:length(channels)
            chan_idx(k)=find(strcmp({merged_EEG.chanlocs.labels},channels{k}));
        end                       

        % Average over channels in cluster
        cluster_data=double(squeeze(mean(merged_EEG.data(chan_idx,:,:),1)));
        
        % Remove low frequency drift <2Hz
        % Twopass 4th order Butterworth highpass filter, reduce filter
        % order in case of instability
        cluster_data=ft_preproc_highpassfilter(cluster_data',...
            merged_EEG.srate, 2, 4, 'but', 'twopass', 'reduce')';
        
        % Extract bursts
        bursts=extract_bursts(cluster_data, all_times,...
            merged_EEG.srate, foi, thresh_sd,...
            'burst_window', params.burst_window);
        
        % Average burst waveforms within time window and baseline correct
        base_idx=(bursts.onset_time>=woi(1)) &...
            (bursts.offset_time<=woi(2));
        subj_mean_waveform=nanmean(bursts.waveform(base_idx,:),1);
        subj_mean_waveform=subj_mean_waveform-mean(subj_mean_waveform);

        % Resolve burst direction (oppositely oriented dipoles 
        % will have reversed polarities)
        mid_times=[-20 20];
        mid_idx=knnsearch(bursts.times',mid_times');
        mid_waveform=subj_mean_waveform(mid_idx(1):mid_idx(2));
        [extreme,extreme_idx]=max(abs(mid_waveform));
        if mid_waveform(extreme_idx)>0
            subj_mean_waveform=subj_mean_waveform.*-1;
        end

        % Save mean waveforms
        subj_waveforms(c_idx,s,:)=subj_mean_waveform;
    end
end

figure();
for c_idx=1:length(clusters)
    shadedErrorBar(bursts.times,...
        squeeze(mean(subj_waveforms(c_idx,:,:),2)),...
        squeeze(std(subj_waveforms(c_idx,:,:),[],2))./sqrt(n_subjects),...
        'LineProps', {'color',cluster_colors{c_idx}});    
end
xlim(bursts.times([1 end]));
xlabel('Time (ms)');
ylabel('Potential (uV)');
legend(clusters);

% for c_idx=1:length(clusters)
%     cluster_base_waveforms=squeeze(subj_base_waveforms(c_idx,:,:));
%     [wout,lags]=woody(cluster_base_waveforms',[],[],'woody','biased');
%     n_aligned=size(cluster_base_waveforms,2)-abs(min(lags))-max(lags);
% 
%     subj_aligned_waveforms=[];
%     for subj_idx=1:size(cluster_base_waveforms,1)
%         start_idx=1+(lags(subj_idx)-min(lags));
%         subj_aligned_waveforms(subj_idx,:)=cluster_base_waveforms(subj_idx,start_idx:start_idx+n_aligned-1);
%         aligned_times=base_bursts.times(start_idx:start_idx+n_aligned-1);
%     end
% end
% % figure();
% % subplot(2,1,1);
% % plot(aligned_times,subj_aligned_waveforms');
% % xlim(aligned_times([1 end]));
% % subplot(2,1,2);
% % plot(aligned_times,nanmean(subj_aligned_waveforms,1));
% % xlim(aligned_times([1 end]));
% 
% 
% 
% % %% GET WOODY ALIGNMENT FOR EACH BURST
% % 
% % % Extract the bursts %
% % cn=1;
% % burst_time_series=[];
% % burst_filtered_time_series=[];
% % burst_aligned_time_series=[];
% % burst_aligned_idx=[];
% % for trial_idx=1:size(burst_times,2)
% %     trial_burst_times=burst_times{trial_idx};
% %     for k=1:size(trial_burst_times,2)
% %         if ~isempty(trial_burst_times)
% %             %if trial_burst_times(k)>0
% %                 try
% %                     burst_time_idx=find(time==trial_burst_times(k));
% %                     if burst_time_idx-n_burst_extract_time_steps/2>=1 && burst_time_idx+n_burst_extract_time_steps/2<=size(cluster_data,1)
% %                         burst_time_series(cn,:)=cluster_data(burst_time_idx-round(n_burst_extract_time_steps/2):burst_time_idx+round(n_burst_extract_time_steps/2),trial_idx);
% %                         burst_filtered_time_series(cn,:)=filtered_cluster_data(burst_time_idx-round(n_burst_extract_time_steps/2):burst_time_idx+round(n_burst_extract_time_steps/2),trial_idx);
% %                         burst_idx(cn,1)=trial_idx;
% %                         burst_idx(cn,2)=burst_time_idx;
% %                         cn=cn+1;
% %                     end
% %                 end
% %             %end
% %         end
% %     end
% % end
% % 
% % %% frequency parameters
% % % min_freq =  study_params.beta_range(1);
% % % max_freq = study_params.beta_range(2);
% % % num_frex = 120;
% % % 
% % % % frequencies vector
% % % %frex = logspace(log10(min_freq), log10(max_freq), num_frex);
% % % frex = linspace(min_freq, max_freq, num_frex);
% % % 
% % % %% wavelet cycles - variable : min 4 max 10
% % % range_cycles = [ 4 10 ];
% % % % cycles vector
% % % cylvec = logspace(log10(range_cycles(1)),log10(range_cycles(end)), num_frex)./ (2*pi*frex);
% % % %% wavelet parameters
% % % wavtime = -1:1/Fs:1; % length of wavelet
% % % half_wave = (length(wavtime)-1)/2;
% % % nWave = length(wavtime);
% % % 
% % % for i=1:size(burst_time_series,1)
% % %     burst_tc=burst_time_series(i,:);
% % %     %% FFT parameters            
% % %     nData = length(burst_tc);
% % %     nConv = nWave + nData - 1;
% % %     dataX=fft(burst_tc,nConv);
% % %     burst_tf_data=zeros(length(frex),nData);
% % %     for fi=1:length(frex) % loop through all frequencies
% % % 
% % %         %% Create wavelate
% % %         wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*cylvec(fi)^2));
% % %         waveletX = fft(wavelet, nConv); % fft of wavelet
% % %         waveletX = waveletX ./ max(waveletX); % normalize fft of wavelet
% % % 
% % %         %% run convolution
% % %         data_conv = ifft(waveletX .* dataX);
% % %         data_conv = data_conv(half_wave+1:end-half_wave);                    
% % % 
% % %         %% compute power
% % %         amp = abs(data_conv);
% % %         burst_tf_data(fi,:)=amp;
% % %     end
% % %     
% % %     %[M,I]=max(burst_tf_data(:));
% % %     %[peak_freq,peak_time]=ind2sub(size(burst_tf_data),I);
% % %     %burst_peak_freqs(i)=frex(peak_freq);
% % %     [M,I]=max(burst_tf_data(:,round(size(burst_tf_data,2)/2)));
% % %     burst_peak_freqs(i)=frex(I);
% % %     
% % %     figure();
% % %     imagesc(1:size(burst_tf_data,2),frex,burst_tf_data);
% % %     set(gca,'ydir','normal');
% % %     colorbar();
% % %     disp('');
% % % end
% % if length(burst_time_series)>0
% %     if params.align
% %         [wout,lags]=woody(burst_time_series',[],[],'woody','biased');
% % 
% %         md=round(size(burst_time_series,2)/2);
% % 
% %         cn=1;
% %         for trial_idx=1:size(burst_time_series,1)
% %             %try
% %             if md+lags(trial_idx)-n_burst_time_steps/2>=1 && md+lags(trial_idx)+n_burst_time_steps/2<=size(burst_time_series,2)
% %                 burst_aligned_time_series(cn,:)=burst_time_series(trial_idx,md+lags(trial_idx)-round(n_burst_time_steps/2):md+lags(trial_idx)+round(n_burst_time_steps/2));
% % 
% %                 burst_aligned_idx(cn,1)=burst_idx(trial_idx,1);
% %                 burst_aligned_idx(cn,2)=burst_idx(trial_idx,2)+lags(trial_idx)+1;
% % 
% %                 cn=cn+1;
% %             end
% %         end
% %     else
% %         burst_aligned_time_series=burst_time_series;
% %         burst_aligned_idx=burst_idx;
% %     end
% end