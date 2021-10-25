function plot_mean_psd_and_lagged_coherence(varargin)
% PLOT_MEAN_PSD_AND_LAGGED_COHERENCE - Plot mean power spectrial density 
% and lagged coherence in the C3 and C4 clusters
%
% Syntax:  plot_mean_psd_and_lagged_coherence(study_info)
%
% Inputs:
%    study_info - structure containing the study information
%
% Example: 
%   plot_mean_psd_and_lagged_coherence(study_info)

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
    
lc_lims=[0.05 0.3];
pow_lims=[-3 3.75; -1.5 5.25];
for st=1:length(study_infos)
    study_info=study_infos{st};
    
    % Load power spectral densities
    load(fullfile(study_info.deriv_dir,'psd.mat'));

    % Number of subjects
    n_subjects=size(spectra,1);

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

        subplot(length(study_info.clusters)*length(study_infos),3,(st-1)*length(study_info.clusters)*3+(c_idx-1)*3+1);

        % Get cluster channels
        channels=study_info.cluster_channels{c_idx};
        chan_idx=cellfun(@(x) find(strcmp({base_EEG.chanlocs.labels},x)),...
            channels);

        % Average over subjects
        mean_spectra=squeeze(nanmean(nanmean(periodic(:,chan_idx,:),2),1));
        stderr_spectra=squeeze(nanstd(nanmean(periodic(:,chan_idx,:),2),[],1))./sqrt(n_subjects);    
        mean_lagged_coh1=squeeze(nanmean(nanmean(nanmean(lagged_coh(:,chan_idx,:,lag1_idx),4),2),1));
        stderr_lagged_coh1=squeeze(nanstd(nanmean(nanmean(lagged_coh(:,chan_idx,:,lag1_idx),4),2),[],1))./sqrt(n_subjects);
        mean_lagged_coh2=squeeze(nanmean(nanmean(nanmean(lagged_coh(:,chan_idx,:,lag2_idx),4),2),1));
        stderr_lagged_coh2=squeeze(nanstd(nanmean(nanmean(lagged_coh(:,chan_idx,:,lag2_idx),4),2),[],1))./sqrt(n_subjects);

        % Plot mean PSD and lagged coherence at 2 cycles
        ax=shadedErrorBaryy(frex, mean_spectra, stderr_spectra, [19 135 255]./255,...
            foi,mean_lagged_coh1,stderr_lagged_coh1, [167 19 0]./255);
        xlabel('Frequency (Hz)');
        ylabel(ax(1),'log(Power uV^2)');
        ylabel(ax(2),'Lagged coherence');

        % Plot lagged coherence at 4 cycles
        axes(ax(2));
        hold all;
        shadedErrorBar(foi,mean_lagged_coh2,stderr_lagged_coh2,'LineProps',...
            {'color',[255 129 112]/255});

        % Plot beta band limits
        plot(study_info.beta_band(c_idx,[1 1]),lc_lims,'k--');
        plot(study_info.beta_band(c_idx,[2 2]),lc_lims,'k--');

        xlim(ax(1),foi([1 end]));
        xlim(ax(2),foi([1 end]));
        ylim(ax(1),pow_lims(st,:));
        ylim(ax(2),lc_lims);


        % Plot mean lagged coherence
        subplot(length(study_info.clusters)*length(study_infos),3,[(st-1)*length(study_info.clusters)*3+(c_idx-1)*3+2 (st-1)*length(study_info.clusters)*3+(c_idx-1)*3+3]);
        contourf(lags,foi,squeeze(nanmean(nanmean(lagged_coh(:,chan_idx,:,:),2),1)),100,...
            'linecolor','none');    
        set(gca,'clim',[0.05 0.25]);
        pos=get(gca,'position');
        colorbar();
        set(gca,'position',pos);
        ylabel('Frequency (Hz)');
        xlabel('Lag (cycles)');
    end
end
