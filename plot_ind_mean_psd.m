function plot_ind_mean_psd(varargin)
% PLOT_MEAN_PSD - Plot individual power spectrial density
% in the C3 and C4 clusters
%
% Syntax:  plot_ind_mean_psd()
%
% Example: 
%   plot_ind_mean_psd()

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
    
    % Load power spectral densities
    load(fullfile(study_info.deriv_dir,'psd.mat'));

    subj_id=study_info.participant_info.participant_id{1};
    subject_data_dir=fullfile(study_info.deriv_dir, subj_id, 'eeg');    
    base_fname=sprintf('%s_11_Epoch_Matched_CSD_baseline.set',subj_id);
    base_EEG=pop_loadset('filepath', subject_data_dir,...
        'filename', base_fname);        

    % Plot spectra in each cluster from 1 to 40Hz
    hold all;
    for c_idx=1:length(study_info.clusters)

        subplot(length(study_infos),length(study_info.clusters),(st-1)*length(study_info.clusters)+c_idx);
        hold all;
        
        % Get cluster channels
        channels=study_info.cluster_channels{c_idx};
        chan_idx=cellfun(@(x) find(strcmp({base_EEG.chanlocs.labels},x)),...
            channels);

        % Average over subjects
        subj_spectra=squeeze(nanmean(periodic(:,chan_idx,:),2));
        
        % Plot PSD
        plot(frex, subj_spectra);
        xlabel('Frequency (Hz)');
        ylabel('log(Power uV^2)');
        
        xlim(frex([1 end]));
    end
end
