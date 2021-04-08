function estimate_threshold_sd(study_info, foi, varargin)
% ESTIMATE_THRESHOLD_SD - Estimate relative amplitude threshold (standard
% deviations above the median) for all participants in the C3 and C4
% electrode clusters. Computes the correlation between the number of bursts
% and mean amplitude across a range of relative thresholds
%
% Syntax:  estimate_threshold_sd(study_info, foi,)
%
% Inputs:
%    study_info - structure containing the study information
%    foi - frequency band of interest ([low high], Hz)
%
% Example: 
%   estimate_threshold_sd(study_info, [13 30])

% Parse optional arguments
defaults=struct('min_ntrials',5);
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

cluster_colors={'b','r'};

% time window of interest (to cut off filtering edge artifacts)
woi=[-1250 1250];

% Range of standard deviations above the median to try
sds=[.1:.1:3];

% Number of subjects
n_subjects=size(study_info.participant_info,1);

% Correlations between burst count and amplitude for both clusters,
% all subjects, at each threshold level
subj_corr=zeros(length(study_info.clusters),n_subjects,length(sds)).*NaN;

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
        exp_EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', exp_fname);        

        % If min number of trials in baseline and experimental epochs
        if base_EEG.trials>=params.min_ntrials &&...
                exp_EEG.trials>=params.min_ntrials
            
            % Time stamps in each trial
            all_base_times=base_EEG.times;
            all_exp_times=exp_EEG.times;

            % Number of trials
            base_trials=base_EEG.trials;
            exp_trials=exp_EEG.trials;

            % Data averaged within each cluster
            base_cluster_data=get_cluster_data(study_info, base_EEG);                        
            exp_cluster_data=get_cluster_data(study_info, exp_EEG);                        
                            
            % Process each cluster
            for c_idx=1:length(study_info.clusters)
                
                % Data for this cluster
                base_data=squeeze(base_cluster_data(c_idx,:,:));
                exp_data=squeeze(exp_cluster_data(c_idx,:,:));
                
                % Get amplitude
                [~, amp_base]=filter_hilbert(base_data, base_EEG.srate,...
                    foi);        
                [~, amp_exp]=filter_hilbert(exp_data, exp_EEG.srate,...
                    foi);   

                % Cut off edges to avoid filter effects
                base_woi_idx=knnsearch(all_base_times',woi');
                amp_base=amp_base(base_woi_idx(1):base_woi_idx(2),:);
                exp_woi_idx=knnsearch(all_exp_times',woi');
                amp_exp=amp_exp(exp_woi_idx(1):exp_woi_idx(2),:);

                for sd_idx=1:length(sds)
                    thresh_sd=sds(sd_idx);

                    % Compute burst threshold
                    threshold=compute_threshold(base_data, exp_data,...
                        all_base_times, all_exp_times, base_EEG.srate,...
                        foi, woi, thresh_sd);

                    % Extract bursts
                    base_bursts=extract_bursts(base_data, all_base_times,...
                        base_EEG.srate, foi, threshold);
                    exp_bursts=extract_bursts(exp_data, all_exp_times,...
                        exp_EEG.srate, foi, threshold);

                    % Number of bursts in each baseline trial epoch
                    base_bursts_per_trial=zeros(1,base_trials);
                    for t_idx=1:base_trials

                        % Bursts in this trial within time window of interest (to
                        % avoid bursts caused by filter edge artifact)
                        trial_bursts=(base_bursts.trial==t_idx);
                        woi_bursts=(base_bursts.peak_time(trial_bursts)>=woi(1)) &...
                            (base_bursts.peak_time(trial_bursts)<=woi(2));

                        % Number of bursts in this trial
                        base_bursts_per_trial(t_idx)=length(find(woi_bursts));
                    end

                    % Number of bursts in each experimental trial epoch
                    exp_bursts_per_trial=zeros(1,exp_trials);
                    for t_idx=1:exp_trials

                        % Bursts in this trial within time window of interest (to
                        % avoid bursts caused by filter edge artifact)
                        trial_bursts=find(exp_bursts.trial==t_idx);
                        woi_bursts=(exp_bursts.peak_time(trial_bursts)>=woi(1)) &...
                            (exp_bursts.peak_time(trial_bursts)<=woi(2));

                        % Number of bursts in this trial
                        exp_bursts_per_trial(t_idx)=length(find(woi_bursts));
                    end

                    % Mean amplitude over whole trial (baseline and experimental
                    % epochs)
                    mean_amp_per_trial=mean([amp_base amp_exp]);
                    % Number of bursts over whole trial (baseline and experimental
                    % epochs)
                    bursts_per_trial=[base_bursts_per_trial exp_bursts_per_trial];

                    % Nonparametric correlation between nuber of bursts in each
                    % trial and mean amplitude
                    subj_corr(c_idx,s,sd_idx)=corr(bursts_per_trial',...
                        mean_amp_per_trial', 'type', 'Spearman');
                end        
            end
        end
    end
end

% Compute sds that maximize mean correlation for each cluster
max_corr_idxs=zeros(1,length(study_info.clusters));
for c_idx=1:length(study_info.clusters)
    
    % Average correlations over subjects for this cluster
    mean_corr=squeeze(nanmean(subj_corr(c_idx,:,:),2));
    
    % Find SD that maximizes mean correlation
    [max_corr,max_corr_idxs(c_idx)]=max(mean_corr);
    max_corr_sd=sds(max_corr_idxs(c_idx));
    fprintf('%s threshold=%.2f, correlation=%.3f\n',study_info.clusters{c_idx},...
        max_corr_sd,max_corr);
end

% Plot correlations
figure();
hold all
for c_idx=1:length(study_info.clusters)
    mean_corr=squeeze(nanmean(subj_corr(c_idx,:,:),2));
    stderr_corr=squeeze(nanstd(subj_corr(1,:,:),[],2))/sqrt(n_subjects);
    shadedErrorBar(sds,mean_corr, stderr_corr,...
        'LineProps',{'color',cluster_colors{c_idx}});
end
for c_idx=1:length(study_info.clusters)
    max_corr_sd=sds(max_corr_idxs(c_idx));
    plot([max_corr_sd max_corr_sd],[0 1],'--',...
        'Color',cluster_colors{c_idx});
end
legend(study_info.clusters)
xlabel('Threshold (SDs above median)');
ylabel('Correlation');