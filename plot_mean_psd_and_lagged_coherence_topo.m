function plot_mean_psd_and_lagged_coherence_topo(varargin)
% PLOT_MEAN_PSD_AND_LAGGED_COHERENCE_TOPO - Plot mean power spectrial 
% density and lagged coherence topography
%
% Syntax:  plot_mean_psd_and_lagged_coherence_topo(study_info)
%
% Inputs:
%    study_info - structure containing the study information
%
% Example: 
%   plot_mean_psd_and_lagged_coherence_topo(study_info)

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
    
    subj_id=study_info.participant_info.participant_id{1};
    % Path containing subject data
    subject_data_dir=fullfile(study_info.deriv_dir, subj_id, 'eeg');
    % Baseline and experimental epoch files
    base_fname=sprintf('%s_11_Epoch_Matched_CSD_baseline.set',subj_id);
    EEG=pop_loadset('filepath', subject_data_dir, 'filename', base_fname);

    % Load power spectral densities
    load(fullfile(study_info.deriv_dir,'psd.mat'));

    % Load lagged coherence
    load(fullfile(study_info.deriv_dir,'lagged_coherence.mat'));

    % Plot lag_idx at 2-2.5 and 3-3.5 cycles
    lag1_idx=find((lags>=2) & (lags<=2.5));
    lag2_idx=find((lags>=3) & (lags<=3.5));

    foi_lims=[nanmean(study_info.alpha_band,1);nanmean(study_info.beta_band,1)];

    foi_lc1=[];
    foi_lc2=[];
    for f=1:size(foi_lims,1)
        flims=foi_lims(f,:);        
    end

    for f=1:size(foi_lims,1)
        flims=foi_lims(f,:);
        spec_freq_idx=find((frex>=flims(1)) & (frex<=flims(2)));
        ps=squeeze(nanmean(nanmean(periodic(:,:,spec_freq_idx),3),1)');
        subplot(size(foi_lims,1)*length(study_infos),3,(st-1)*size(foi_lims,1)*3+(f-1)*3+1);
        
        tmpEEG = eeg_emptyset;
        tmpEEG.chanlocs = EEG.chanlocs;
        tmpEEG.nbchan   = length(tmpEEG.chanlocs);
        tmpEEG.data     = zeros(tmpEEG.nbchan,1,1);
        tmpEEG.trials   = 1;
        tmpEEG.pnts     = 1;
        tmpEEG.xmin     = 0;
        tmpEEG.srate    = 1;
        tmpEEG.xmax     = 0;
        tmpEEG.data     = ps;
        tmpEEG = eeg_checkset(tmpEEG);
        ftEEG = eeglab2fieldtrip(tmpEEG, 'preprocessing', 'none');
        
        ftEEG.avg=ps;
        ftEEG.sig=ones(size(ftEEG.avg));
        
        cfg=[];
        cfg.rotate=90;
        cfg.elec = ftEEG.elec;
        lay = ft_prepare_layout(cfg);
        
        cfg              = [];
        cfg.colorbar     = 'yes';
        cfg.colorbartext = 'log(Power uV^2)';
        cfg.marker       = 'on';
        cfg.comment      = 'no';
        cfg.maskparameter = 'sig';
        cfg.style = 'straight';
        cfg.zlim=[min(ftEEG.avg(:)) max(ftEEG.avg(:))];
        cfg.layout=lay;
        
        ft_topoplotER(cfg, ftEEG);
                
        title(sprintf('Power, %.1f-%.1fHz',flims(1),flims(2)));
        
        lc_freq_idx=find((foi>=flims(1)) & (foi<=flims(2)));
        lc1=squeeze(nanmean(nanmean(nanmean(lagged_coh(:,:,lc_freq_idx,lag1_idx),3),4),1))';
        lc2=squeeze(nanmean(nanmean(nanmean(lagged_coh(:,:,lc_freq_idx,lag2_idx),3),4),1))';
        
        subplot(size(foi_lims,1)*length(study_infos),3,(st-1)*size(foi_lims,1)*3+(f-1)*3+2);
        tmpEEG = eeg_emptyset;
        tmpEEG.chanlocs = EEG.chanlocs;
        tmpEEG.nbchan   = length(tmpEEG.chanlocs);
        tmpEEG.data     = zeros(tmpEEG.nbchan,1,1);
        tmpEEG.trials   = 1;
        tmpEEG.pnts     = 1;
        tmpEEG.xmin     = 0;
        tmpEEG.srate    = 1;
        tmpEEG.xmax     = 0;
        tmpEEG.data     = lc1;
        tmpEEG = eeg_checkset(tmpEEG);
        ftEEG = eeglab2fieldtrip(tmpEEG, 'preprocessing', 'none');
        
        ftEEG.avg=lc1;
        ftEEG.sig=ones(size(ftEEG.avg));
        
        cfg=[];
        cfg.rotate=90;
        cfg.elec = ftEEG.elec;
        lay = ft_prepare_layout(cfg);
        
        cfg              = [];
        cfg.colorbar     = 'no';
        cfg.marker       = 'on';
        cfg.comment      = 'no';
        cfg.maskparameter = 'sig';
        cfg.style = 'straight';
        cfg.zlim=[0 max(lc1(:))];
        cfg.layout=lay;
        
        ft_topoplotER(cfg, ftEEG);
        
        title(sprintf('2 cycles, %.1f-%.1fHz',flims(1),flims(2)));
        
        subplot(size(foi_lims,1)*length(study_infos),3,(st-1)*size(foi_lims,1)*3+(f-1)*3+3);
        tmpEEG = eeg_emptyset;
        tmpEEG.chanlocs = EEG.chanlocs;
        tmpEEG.nbchan   = length(tmpEEG.chanlocs);
        tmpEEG.data     = zeros(tmpEEG.nbchan,1,1);
        tmpEEG.trials   = 1;
        tmpEEG.pnts     = 1;
        tmpEEG.xmin     = 0;
        tmpEEG.srate    = 1;
        tmpEEG.xmax     = 0;
        tmpEEG.data     = lc2;
        tmpEEG = eeg_checkset(tmpEEG);
        ftEEG = eeglab2fieldtrip(tmpEEG, 'preprocessing', 'none');
        
        ftEEG.avg=lc2;
        ftEEG.sig=ones(size(ftEEG.avg));
        
        cfg=[];
        cfg.rotate=90;
        cfg.elec = ftEEG.elec;
        lay = ft_prepare_layout(cfg);
        
        cfg              = [];
        cfg.colorbar     = 'yes';
        cfg.colorbartext = 'lagged coherence';
        cfg.marker       = 'on';
        cfg.comment      = 'no';
        cfg.maskparameter = 'sig';
        cfg.style = 'straight';
        cfg.zlim=[0 max(lc1(:))];
        cfg.layout=lay;
        
        ft_topoplotER(cfg, ftEEG);
        
        title(sprintf('3 cycles, %.1f-%.1fHz',flims(1),flims(2)));
        pos=get(gca,'Position');
        colorbar();
        set(gca,'Position',pos);
    end
end