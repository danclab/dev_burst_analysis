function plot_bursts_and_amplitude(study_info, foi, thresh_sd, varargin)
% PLOT_BURSTS_AND_AMPLITUDE - Plot distribution of burst timing and
% amplitude time course within frequency band in the C3 and C4 clusters
%
% Syntax:  plot_bursts_and_amplitude(study_info, foi, thresh_sd)
%
% Inputs:
%    study_info - structure containing the study information
%    foi - frequency band of interest ([low high], Hz)
%    thresh_sd - number of standard deviations above the median amplitude
%        to set the threshold at 
%
% Example: 
%   plot_bursts_and_amplitude(study_info, [13 30], 1.5)

% Parse optional arguments
defaults=struct('min_ntrials',5);
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

condition='execute';
cond_idx=find(strcmp(study_info.conditions,condition));
base_event=study_info.baseline_evts{cond_idx};
exp_event=study_info.exp_evts{cond_idx};

% time window of interest (to cut off filtering edge artifacts)
woi=[-1250 1250];

% Number of subjects
n_subjects=size(study_info.participant_info,1);

all_subj_base_bursts=struct('trial',{},'peak_time',{},'onset_time',{},...
    'offset_time',{},'peak_amp',{},'waveform',{},'times',{});
all_subj_exp_bursts=struct('trial',{},'peak_time',{},'onset_time',{},...
    'offset_time',{},'peak_amp',{},'waveform',{},'times',{});

% Baseline mean amplitude for each subject in each cluster
subj_amp_base=[];
% Experimental mean amplitude for each subject in each cluster
subj_amp_exp=[];
% Number of baseline trials per subject
ntrials_base=zeros(1,n_subjects);
% Number of experimental trials per subject
ntrials_exp=zeros(1,n_subjects);

for s=1:n_subjects
    
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s};
    
    % Path containing subject data
    subject_data_dir=fullfile(study_info.deriv_dir, subj_id, 'eeg');
    
    % Baseline and experimental epoch files
    base_fname=sprintf('%s_11_Epoch_Matched_CSD_baseline.set',subj_id);
    exp_fname=sprintf('%s_11_Epoch_Matched_CSD_experimental.set',subj_id);
    
    if exist(fullfile(subject_data_dir,base_fname),'file')==2 &&...
            exist(fullfile(subject_data_dir,exp_fname),'file')==2
        
        % Load data
        base_EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', base_fname);        
        exp_EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', exp_fname); 
            
        % Time stamps in each trial
        all_base_times=base_EEG.times;
        all_exp_times=exp_EEG.times;

        % Get condition trials
        base_trials=find(strcmp({base_EEG.event.type},base_event));
        exp_trials=find(strcmp({exp_EEG.event.type},exp_event));
                
        % Number of trials
        ntrials_base(s)=length(base_trials);
        ntrials_exp(s)=length(exp_trials);

        % If min number of trials in baseline and experimental epochs
        if ntrials_base(s)>=params.min_ntrials &&...
                ntrials_exp(s)>=params.min_ntrials
        
            % Process each cluster
            for c_idx=1:length(study_info.clusters)

                % Channels in this cluster
                channels=study_info.cluster_channels{c_idx};

                % Find indices of cluster channels
                chan_idx=zeros(1,length(channels));
                for k=1:length(channels)
                    chan_idx(k)=find(strcmp({exp_EEG.chanlocs.labels},...
                        channels{k}));
                end                       

                % Average over channels in cluster
                base_data=double(squeeze(mean(base_EEG.data(chan_idx,:,:),1)));
                exp_data=double(squeeze(mean(exp_EEG.data(chan_idx,:,:),1)));

                % Compute burst threshold
                threshold=compute_threshold(base_data, exp_data,...
                    all_base_times, all_exp_times, base_EEG.srate, foi,...
                    woi, thresh_sd);
                
                base_data=base_data(:,base_trials);
                exp_data=exp_data(:,exp_trials);
                
                % Extract bursts
                subj_base_bursts=extract_bursts(base_data, all_base_times,...
                    base_EEG.srate, foi, threshold);
                subj_exp_bursts=extract_bursts(exp_data, all_exp_times,...
                    exp_EEG.srate, foi, threshold);

                % Get amplitude
                [~, amp_base]=filter_hilbert(base_data, base_EEG.srate, foi);        
                [~, amp_exp]=filter_hilbert(exp_data, exp_EEG.srate, foi);        

                % Baseline-correct amplitude
                woi_idx=knnsearch(all_base_times',woi');
                mean_amp_base_cluster_data=mean(mean(amp_base(woi_idx(1):woi_idx(2),:)));
                subj_amp_base(s,c_idx,:)=nanmean(amp_base-mean_amp_base_cluster_data,2);
                subj_amp_exp(s,c_idx,:)=mean(amp_exp-mean_amp_base_cluster_data,2);

                % Save bursts
                all_subj_base_bursts(s,c_idx)=subj_base_bursts;
                all_subj_exp_bursts(s,c_idx)=subj_exp_bursts;
            end
        end
    end
end

% Plot each cluster
for c_idx=1:length(study_info.clusters)
    figure();
    
    % Raster plot of baseline bursts in all trials
    subplot(3,2,[1 3]);
    hold all
    trial_offset=0;
    for s=1:n_subjects
        if ~isempty(all_subj_base_bursts(s,c_idx))
            base_bursts=all_subj_base_bursts(s,c_idx);
            plot(base_bursts.peak_time,base_bursts.trial+trial_offset,'.k');
            trial_offset=trial_offset+ntrials_base(s);
            if s<n_subjects
                plot(woi,[trial_offset-.5 trial_offset-.5],'r--');
            end
        end
    end
    xlim(woi);
    ylim([0 trial_offset]);
    ylabel('Trial');

    % Raster plot of experimental bursts in all trials
    subplot(3,2,[2 4]);
    hold all
    trial_offset=0;
    for s=1:n_subjects
        if ~isempty(all_subj_exp_bursts(s,c_idx))
            exp_bursts=all_subj_exp_bursts(s,c_idx);
            plot(exp_bursts.peak_time,exp_bursts.trial+trial_offset,'.k');
            trial_offset=trial_offset+ntrials_exp(s);
            if s<n_subjects
                plot(woi,[trial_offset-.5 trial_offset-.5],'r--');
            end
        end
    end
    xlim(woi);
    ylim([0 trial_offset]);
    ylabel('Trial');

    % Bin and smoothing width
    binwidth=2;
    smoothw=30;
    
    % Compute bin centers
    base_bins=[all_base_times(1):binwidth:all_base_times(end)];
    base_bins=base_bins(1:end-1)+binwidth/2;
    exp_bins=[all_exp_times(1):binwidth:all_exp_times(end)];
    exp_bins=exp_bins(1:end-1)+binwidth/2;

    % Probability of a burst in each bin for each subject
    burst_base_prob=zeros(n_subjects,length(base_bins))*NaN;
    burst_exp_prob=zeros(n_subjects,length(exp_bins))*NaN;
    for s=1:n_subjects
        
        % Histogram of baseline burst times
        if ~isempty(all_subj_base_bursts(s,c_idx))
            base_bursts=all_subj_base_bursts(s,c_idx);
            [counts,~]=histc(base_bursts.peak_time,base_bins);
            if ~isempty(counts)
                % Compute mean bursts per bin and convolve
                burst_base_prob(s,:)=counts./ntrials_base(s);
                burst_base_prob(s,:)=filtfilt(ones(smoothw+1,1)/(smoothw+1),1,burst_base_prob(s,:));
            end
        end

        % Histogram of experimental burst times
        if ~isempty(all_subj_exp_bursts(s,c_idx))
            exp_bursts=all_subj_exp_bursts(s,c_idx);
            [counts,~]=histc(exp_bursts.peak_time,base_bins);
            if ~isempty(counts)
                % Compute mean bursts per bin and convolve
                burst_exp_prob(s,:)=counts./ntrials_exp(s);
                burst_exp_prob(s,:)=filtfilt(ones(smoothw+1,1)/(smoothw+1),1,burst_exp_prob(s,:));
            end
        end
    end

    % Cut off edges to remove bursts caused by filter artifacts
    base_b_idx=knnsearch(base_bins',woi');
    base_t_idx=knnsearch(all_base_times',woi');
    exp_b_idx=knnsearch(exp_bins',woi');
    exp_t_idx=knnsearch(all_exp_times',woi');

    % Mean and std error of baseline burst probability across subjects
    mean_base_prob=nanmean(burst_base_prob(:,base_b_idx(1):base_b_idx(2)));
    stderr_base_prob=nanstd(burst_base_prob(:,base_b_idx(1):base_b_idx(2)),[],1)./sqrt(n_subjects);
    % Mean and std error of experimental burst probability across subjects
    mean_exp_prob=nanmean(burst_exp_prob(:,exp_b_idx(1):exp_b_idx(2)));
    stderr_exp_prob=nanstd(burst_exp_prob(:,exp_b_idx(1):exp_b_idx(2)),[],1)./sqrt(n_subjects);
    
    % Mean and std error of baseline amplitude across subjects
    mean_base_amp=squeeze(nanmean(subj_amp_base(:,c_idx,base_t_idx(1):base_t_idx(2))));
    stderr_base_amp=squeeze(nanstd(subj_amp_base(:,c_idx,base_t_idx(1):base_t_idx(2))))./sqrt(n_subjects);
    % Mean and std error of experimental amplitude across subjects
    mean_exp_amp=squeeze(nanmean(subj_amp_exp(:,c_idx,exp_t_idx(1):exp_t_idx(2))));
    stderr_exp_amp=squeeze(nanstd(subj_amp_exp(:,c_idx,exp_t_idx(1):exp_t_idx(2))))./sqrt(n_subjects);

    % Burst probability y axis limits
    min_prob=min([mean_base_prob mean_exp_prob]-[stderr_base_prob stderr_exp_prob]);
    max_prob=max([mean_base_prob mean_exp_prob]+[stderr_base_prob stderr_exp_prob]);
    prob_ylim=[min_prob-.1*(max_prob-min_prob) max_prob+.1*(max_prob-min_prob)];

    % Ampltude y axis limits
    min_amp=min([mean_base_amp; mean_exp_amp]-[stderr_base_amp; stderr_exp_amp]);
    max_amp=max([mean_base_amp; mean_exp_amp]+[stderr_base_amp; stderr_exp_amp]);
    amp_ylim=[min_amp-.1*(max_amp-min_amp) max_amp+.1*(max_amp-min_amp)];

    % Plot baseline burst probability and amplitude
    subplot(3,2,5);
    ax=shadedErrorBaryy(base_bins(base_b_idx(1):base_b_idx(2)),...
        mean_base_prob, stderr_base_prob, 'b',...
        all_base_times(base_t_idx(1):base_t_idx(2)),...
        mean_base_amp,stderr_base_amp, 'r');
    xlim(ax(1),woi);
    ylim(ax(1),prob_ylim);
    xlim(ax(2),woi);
    ylim(ax(2),amp_ylim);
    xlabel('Time (ms)');
    ylabel(ax(1),'Burst probability');
    
    % Plot experimental burst probability and amplitude
    subplot(3,2,6);
    ax=shadedErrorBaryy(exp_bins(exp_b_idx(1):exp_b_idx(2)),...
        mean_exp_prob, stderr_exp_prob, 'b',...
        all_exp_times(exp_t_idx(1):exp_t_idx(2)),...
        mean_exp_amp, stderr_exp_amp, 'r');
    xlim(ax(1),woi);
    ylim(ax(1),prob_ylim);
    xlim(ax(2),woi);
    ylim(ax(2),amp_ylim);
    xlabel('Time (ms)');
    ylabel(ax(2),'Amplitude');
end