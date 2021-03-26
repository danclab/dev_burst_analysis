function compute_psd(study_info, varargin)
% COMPUTE_PSD - Compute power spectrial density in the C3 and C4 clusters
%
% Syntax:  compute_psd(study_info, varargin)
%
% Inputs:
%    study_info - structure containing the study information
%
% Example: 
%   compute_psd(study_info, varargin)

% Parse optional arguments
defaults=struct('min_ntrials',5);
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

% Open EEGlab
[ALLEEG, EEG, CURRENTSET] = eeglab;

% Power spectra for each subject in each cluster
spectra=[];

% Number of subjects
n_subjects=size(study_info.participant_info,1);

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

            % Data averaged within each cluster
            cluster_data=zeros(length(study_info.clusters), merged_EEG.pnts,...
                merged_EEG.trials);

            for c_idx=1:length(study_info.clusters)

                % Channels in this cluster
                channels=study_info.cluster_channels{c_idx};

                % Find indices of cluster channels
                chan_idx=zeros(1,length(channels));
                for k=1:length(channels)
                    chan_idx(k)=find(strcmp({merged_EEG.chanlocs.labels},...
                        channels{k}));
                end    

                % Average over channels in this cluster
                cluster_data(c_idx,:,:)=double(squeeze(mean(merged_EEG.data(chan_idx,:,:),1)));
            end

            % Use a window size of twice the sampling rate with a 50%
            % overlap
            winsize=merged_EEG.srate*2;
            overlap=round(winsize/2);

            % Compute power spectral density for each cluster using Welch's
            % method
            [subj_spectra,frex,~,~,~] = spectopo(cluster_data,...
                merged_EEG.pnts, merged_EEG.srate, 'winsize', winsize,...
                'overlap', overlap, 'plot', 'off');

            spectra(end+1,:,:)=subj_spectra;    
        end
    end
end

save(fullfile(study_info.deriv_dir,'psd.mat'),'frex','spectra');
