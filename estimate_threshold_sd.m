function max_corr_sds=estimate_threshold_sd(study_info, fois, varargin)
% ESTIMATE_THRESHOLD_SD - Estimate relative amplitude threshold (standard
% deviations above the median) for all participants in the C3 and C4
% electrode clusters. Computes the correlation between the number of bursts
% and mean amplitude across a range of relative thresholds
%
% Syntax:  estimate_threshold_sd(study_info, fois)
%
% Inputs:
%    study_info - structure containing the study information
%    fois - frequency band of interest for each cluster ([low_C3 high_C3; low_C4 high_C4], Hz)
%
% Example: 
%   estimate_threshold_sd(study_info, [13 30; 13 30])

% Parse optional arguments
defaults=struct('min_ntrials',5);
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

cluster_colors=cbrewer('qual','Dark2',7);
cluster_colors=cluster_colors(3:7,:);

% Range of standard deviations above the median to try
sds=[.1:.1:3];

% Number of subjects
n_subjects=size(study_info.participant_info,1);

% Correlations between burst count and amplitude for both clusters,
% all subjects, at each threshold level
subj_corr=zeros(length(study_info.clusters),n_subjects,length(sds)).*NaN;

% Open EEGlab
[ALLEEG, EEG, CURRENTSET] = eeglab;

subj_amp={};

for s=1:1:n_subjects
    
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
        [ALLEEG, EEG] = eeg_store(ALLEEG, base_EEG, 1);        
        exp_EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', exp_fname);        
        [ALLEEG, EEG] = eeg_store(ALLEEG, exp_EEG, 2);

        % If min number of trials in baseline and experimental epochs
        if base_EEG.trials>=params.min_ntrials &&...
                exp_EEG.trials>=params.min_ntrials
            
            % Merge baseline and experimental EEG data
            merged_EEG=pop_mergeset(ALLEEG,1:2,0);

            % Time stamps in each trial
            all_times=merged_EEG.times;
            n_times=length(all_times);
            
            % Number of trials
            trials=merged_EEG.trials;
            
            % Sampling rate
            srate=merged_EEG.srate;
            
            % Process each cluster
            for c_idx=1:length(study_info.clusters)
                
                % Get cluster channels
                channels=study_info.cluster_channels{c_idx};
                chan_idx=cellfun(@(x) find(strcmp({merged_EEG.chanlocs.labels},x)),...
                    channels);

                % Compute amplitude in each channel of C3 cluster
                cluster_amp=zeros(length(chan_idx),n_times,trials);

                for i=1:length(chan_idx)
                    chan_data=squeeze(merged_EEG.data(chan_idx(i),:,:));

                    [ch_filt, ch_amp]=filter_hilbert(chan_data, srate, fois(c_idx,:));
                    
                    % Get rid of padding
                    cluster_amp(i,:,:)=ch_amp;
                end

                % Average amplitude over cluster channels
                cluster_amp=squeeze(mean(cluster_amp));
                subj_amp{c_idx,s}=cluster_amp(:);
                    
                for sd_idx=1:length(sds)
                    thresh_sd=sds(sd_idx);

                    % Compute absolute burst threshold
                    threshold=median(cluster_amp(:))+(thresh_sd*std(cluster_amp(:)));

                    % Number of bursts in each trial
                    trial_bursts=zeros(1,trials);

                    % Go through amplitude each trial 
                    for t_idx=1:trials
                        
                        % All times when amp is over threshold
                        over_thresh_idx = cluster_amp(:,t_idx)>=threshold;
                        % Change in threshold crossing
                        over_thresh_diff = diff(over_thresh_idx);
                        % All times when amp is over threshold and previous time is
                        % under threshold
                        all_burst_start_idx=find(over_thresh_diff)+1;

                        %% Process each burst
                        for k=1:length(all_burst_start_idx)
                            burst_start_idx=all_burst_start_idx(k);

                            % Find next time amplitude goes below threshold
                            burst_end_idx=find(cluster_amp(burst_start_idx+1:end,t_idx)<threshold,1);

                            % If it goes below threshold before the end of the trial
                            if ~isempty(burst_end_idx)% && burst_end_idx>1
                               trial_bursts(t_idx)=trial_bursts(t_idx)+1;                
                            end
                        end
                    end

                    % Nonparametric correlation between nuber of bursts in each
                    % trial and mean amplitude
                    subj_corr(c_idx,s,sd_idx)=corr(trial_bursts', mean(cluster_amp)',...
                        'type', 'Spearman');
                    
                end
            end
        end
    end
end

% Compute sds that maximize mean correlation for each cluster
max_corr_idxs=zeros(1,length(study_info.clusters));
max_corr_sds=zeros(1,length(study_info.clusters));
for c_idx=1:length(study_info.clusters)
    % Average correlations over subjects for this cluster
    mean_corr=squeeze(nanmean(subj_corr(c_idx,:,:),2));

    % Find SD that maximizes mean correlation
    [max_corr,max_corr_idxs(c_idx)]=max(mean_corr);
    max_corr_sds(c_idx)=sds(max_corr_idxs(c_idx));
    fprintf('%s threshold=%.2f, correlation=%.3f\n',study_info.clusters{c_idx},...
        max_corr_sds(c_idx),max_corr);
end

% Plot correlations
figure();
subplot(1,2,1);
hold all
for c_idx=1:length(study_info.clusters)
    mean_corr=squeeze(nanmean(subj_corr(c_idx,:,:),2));
    stderr_corr=squeeze(nanstd(subj_corr(c_idx,:,:),[],2))/sqrt(n_subjects);
    shadedErrorBar(sds,mean_corr, stderr_corr,...
        'LineProps',{'color',cluster_colors(c_idx,:)});
end
yl=ylim();
yl(2)=1;
for c_idx=1:length(study_info.clusters)
    max_corr_sd=sds(max_corr_idxs(c_idx));
    plot([max_corr_sd max_corr_sd],yl,'--',...
        'Color',cluster_colors(c_idx,:));
end
ylim(yl);
legend(study_info.clusters)
xlabel('Threshold (SDs above median)');
ylabel('Correlation');

% Plot amplitude distribution
subplot(1,2,2);
hold all
% Bins to compute amplitude density in
bins=linspace(0,60,101);
bins=bins(2:end);
% Plot amplitude density
for s=1:n_subjects
     for c_idx=1:length(study_info.clusters)
        if ~isempty(subj_amp{c_idx,s})
            density=ksdensity(subj_amp{c_idx,s},bins);
            plot(bins,density,'color',cluster_colors(c_idx,:));            
        end
    end
end
for s=1:n_subjects
     for c_idx=1:length(study_info.clusters)
        if ~isempty(subj_amp{c_idx,s})
            threshold=median(subj_amp{c_idx,s})+(max_corr_sds(c_idx)*std(subj_amp{c_idx,s}));
            plot([threshold threshold],ylim(),'--','Color',cluster_colors(c_idx,:));
        end
     end
end
legend(study_info.clusters);
xlabel('Amplitude (uV)');
ylabel('Density');