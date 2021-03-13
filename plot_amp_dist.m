function plot_amp_dist(study_info, foi, woi, thresh_sd, varargin)
% PLOT_AMP_DIST - Plot distribution of amplitude within frequency band
% across all time points for each subject, and subject-specific thresholds
% in the C3 and C4 clusters
%
% Syntax:  plot_amp_dist(study_info, foi, woi, thresh_sd)
%
% Inputs:
%    study_info - structure containing the study information
%    foi - frequency band of interest ([low high], Hz)
%    woi - time window of interest (to cut off filtering edge artifacts;
%        [low high], ms)
%    thresh_sd - number of standard deviations above the median amplitude
%        to set the threshold at 
%
% Example: 
%   plot_amp_dist(study_info, [13 30], [-1000 1000], 1.5)

% Parse optional arguments
defaults=struct();
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

% Amplitude at each time point, for each subject in each cluster
subj_amp={};
% Max amplitude for each subject in each cluster
subj_max_amp=zeros(length(clusters),n_subjects);
% Subject-specific thresholds for each cluster
subj_thresholds=zeros(length(clusters),n_subjects);

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
    exp_EEG=pop_loadset('filepath', subject_data_dir,...
        'filename', exp_fname);        

    % Time stamps in each trial
    all_base_times=base_EEG.times;
    all_exp_times=exp_EEG.times;

    % Process each cluster
    for c_idx=1:length(clusters)
        
        % Channels in this cluster
        channels=cluster_channels{c_idx};

        % Find indices of cluster channels
        chan_idx=zeros(1,length(channels));
        for k=1:length(channels)
            chan_idx(k)=find(strcmp({exp_EEG.chanlocs.labels},channels{k}));
        end                       

        % Average over channels in cluster
        base_data=double(squeeze(mean(base_EEG.data(chan_idx,:,:),1)));
        exp_data=double(squeeze(mean(exp_EEG.data(chan_idx,:,:),1)));

        % Get amplitude
        [~, amp_base]=filter_hilbert(base_data, base_EEG.srate, foi);        
        [~, amp_exp]=filter_hilbert(exp_data, exp_EEG.srate, foi);   
        
        % Cut off edges to avoid edge effects
        base_woi_idx=knnsearch(all_base_times',woi');
        amp_base=amp_base(base_woi_idx(1):base_woi_idx(2),:);
        exp_woi_idx=knnsearch(all_exp_times',woi');
        amp_exp=amp_exp(exp_woi_idx(1):exp_woi_idx(2),:);
        
        % Compute burst threshold
        threshold=compute_threshold(base_data, exp_data, all_base_times,...
            all_exp_times, base_EEG.srate, foi, woi, thresh_sd);
        
        % Save amplitude in each time point, max amplitude, and threshold
        subj_amp{c_idx,s}=[reshape(amp_base,size(amp_base,1)*size(amp_base,2),1);...
            reshape(amp_exp,size(amp_exp,1)*size(amp_exp,2),1)];
        subj_max_amp(c_idx,s)=max(subj_amp{c_idx,s});
        subj_thresholds(c_idx,s)=threshold;
    end
end

figure();
hold all;
% Bins to compute amplitude density in
bin_max=max(subj_max_amp(:));
bins=linspace(0,bin_max,101);
bins=bins(2:end);
% Plot amplitude density
for s=1:n_subjects
     for c_idx=1:length(clusters)
        density=ksdensity(subj_amp{c_idx,s},bins);
        plot(bins,density,cluster_colors{c_idx});
    end
end
legend(clusters);
% Plot thresholds
for s=1:n_subjects
    for c_idx=1:length(clusters)
        plot([subj_thresholds(c_idx,s) subj_thresholds(c_idx,s)],ylim(),'--','color',cluster_colors{c_idx});
    end        
end
xlabel('Amplitude (uV)');
ylabel('Density');