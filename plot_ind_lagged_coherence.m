function plot_ind_lagged_coherence(varargin)
% PLOT_IND_LAGGED_COHERENCE - Plot individual lagged coherence
% in the C3 and C4 clusters
%
% Syntax:  plot_ind_lagged_coherence()
%
% Example: 
%   plot_ind_lagged_coherence()

% Parse optional arguments
defaults=struct();
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

study_infos={};
study_infos{1}=init_umd12m_study_info();
study_infos{2}=init_umdadult_study_info();

fig=figure();

for st=1:length(study_infos)
    study_info=study_infos{st};
    
    % Load lagged coherence
    load(fullfile(study_info.deriv_dir,'lagged_coherence.mat'));
    
    % Plot lag_idx at 2-2.5 and 3-3.5 cycles
    lag1_idx=find((lags>=2) & (lags<=2.5));
    lag2_idx=find((lags>=3) & (lags<=3.5));
    
    subj_id=study_info.participant_info.participant_id{1};
    subject_data_dir=fullfile(study_info.deriv_dir, subj_id, 'eeg');
    base_fname=sprintf('%s_11_Epoch_Matched_CSD_baseline.set',subj_id);
    base_EEG=pop_loadset('filepath', subject_data_dir,...
        'filename', base_fname);
    
    % Plot mean and std error of spectra in each cluster from 1 to 40Hz
    hold all;
    for c_idx=1:length(study_info.clusters)
        
        subplot(length(study_info.clusters),length(study_infos),(st-1)*length(study_info.clusters)+c_idx);
        hold all;
        
        % Get cluster channels
        channels=study_info.cluster_channels{c_idx};
        chan_idx=cellfun(@(x) find(strcmp({base_EEG.chanlocs.labels},x)),...
            channels);
        
        % Average over subjects
        subj_lagged_coh1=squeeze(nanmean(nanmean(lagged_coh(:,chan_idx,:,lag1_idx),4),2));
        subj_lagged_coh2=squeeze(nanmean(nanmean(lagged_coh(:,chan_idx,:,lag2_idx),4),2));
        
        % Plot mean PSD and lagged coherence at 2 cycles
        plot(foi,subj_lagged_coh1);
        plot(foi,subj_lagged_coh2,'--');
        xlabel('Frequency (Hz)');
        ylabel('Lagged coherence');
    end
end
