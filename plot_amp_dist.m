function plot_amp_dist(study_info, foi, varargin)
% PLOT_AMP_DIST - Plot distribution of amplitude within frequency band
% across all time points for each subject, and subject-specific thresholds
% in the C3 and C4 clusters
%
% Syntax:  plot_amp_dist(study_info, foi, woi, thresh_sd)
%
% Inputs:
%    study_info - structure containing the study information
%    foi - frequency band of interest ([low high], Hz)
%
% Example: 
%   plot_amp_dist(study_info, [13 30])

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

% Number of subjects
n_subjects=size(study_info.participant_info,1);

% Amplitude at each time point, for each subject in each cluster
subj_amp={};

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

        % If min number of trials in baseline and experimental epochs
        if base_EEG.trials>=params.min_ntrials &&...
                exp_EEG.trials>=params.min_ntrials
            
            % Time stamps in each trial
            all_base_times=base_EEG.times;
            all_exp_times=exp_EEG.times;

            % Data averaged within each cluster
            base_cluster_data=get_cluster_data(study_info, base_EEG);                        
            exp_cluster_data=get_cluster_data(study_info, exp_EEG);
            
            % Process each cluster
            for c_idx=1:length(study_info.clusters)
                
                % Data for this cluster
                base_data=squeeze(base_cluster_data(c_idx,:,:));
                exp_data=squeeze(exp_cluster_data(c_idx,:,:));

                % Get amplitude
                [~, amp_base]=filter_hilbert(base_data, base_EEG.srate, foi);        
                [~, amp_exp]=filter_hilbert(exp_data, exp_EEG.srate, foi);   

                % Cut off edges to avoid filter effects
                base_woi_idx=knnsearch(all_base_times',woi');
                amp_base=amp_base(base_woi_idx(1):base_woi_idx(2),:);
                exp_woi_idx=knnsearch(all_exp_times',woi');
                amp_exp=amp_exp(exp_woi_idx(1):exp_woi_idx(2),:);

                % Save amplitude in each time point, max amplitude, and threshold
                subj_amp{c_idx,s}=[reshape(amp_base,size(amp_base,1)*size(amp_base,2),1);...
                    reshape(amp_exp,size(amp_exp,1)*size(amp_exp,2),1)];
            end
        end
    end
end

figure();
hold all;
% Bins to compute amplitude density in
bins=linspace(0,100,101);
bins=bins(2:end);
% Plot amplitude density
for s=1:n_subjects
     for c_idx=1:length(study_info.clusters)
        if ~isempty(subj_amp{c_idx,s})
            density=ksdensity(subj_amp{c_idx,s},bins);
            plot(bins,density,cluster_colors{c_idx});
        end
    end
end
legend(study_info.clusters);
xlabel('Amplitude (uV)');
ylabel('Density');