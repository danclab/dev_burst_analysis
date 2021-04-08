function plot_mean_psd_and_lagged_coherence(study_info, varargin)
% PLOT_MEAN_PSD - Plot mean power spectrial density and lagged coherence
% in the C3 and C4 clusters
%
% Syntax:  plot_mean_psd_and_lagged_coherence(study_info)
%
% Inputs:
%    study_info - structure containing the study information
%
% Example: 
%   plot_mean_psd(study_info)

% Parse optional arguments
defaults=struct();
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

% Load power spectral densities
load(fullfile(study_info.deriv_dir,'psd.mat'));
% Look at frequencies from 5-40Hz because we don't have lagged coherence at
% lower frequencies (trials aren't long enough)
freq_idx=find((frex>=5) & (frex<=40));
frex=frex(freq_idx);
spectra=spectra(:,:,freq_idx);

% Number of subjects
n_subjects=size(spectra,1);

% Remove 1/f component
subj_psd_resids=zeros(length(study_info.clusters),n_subjects,length(frex));

for c_idx=1:length(study_info.clusters)
    for s_idx=1:n_subjects
        % Subject spectrum
        psd=squeeze(spectra(s_idx,c_idx,:));
        % 1/f
        oof=(1./frex);
        % Fit 1/f to spectrum
        lm_psd=fitlm(oof,psd);
        % Get residuals
        subj_psd_resids(c_idx,s_idx,:)=lm_psd.Residuals.Standardized;
    end
end
   
% Load lagged coherence
load(fullfile(study_info.deriv_dir,'lagged_coherence.mat'));

% Plot lag_idx at 2 and 4 cycles
lag1_idx=find(lags==2);
lag2_idx=find(lags==4);

% Max mean lagged coherence
grand_mean_lagged_coh=nanmean(lagged_coh);
grand_se_lagged_coh=nanstd(lagged_coh)/sqrt(n_subjects);
max_mean_lagged_coh=max(grand_mean_lagged_coh(:)+grand_se_lagged_coh(:));
grand_mean_spectra=nanmean(subj_psd_resids,2);
grandse_mean_spectra=nanstd(subj_psd_resids,[],2)./sqrt(n_subjects);
min_mean_spectra=min(grand_mean_spectra(:)+grandse_mean_spectra(:));
max_mean_spectra=max(grand_mean_spectra(:)+grandse_mean_spectra(:));

% Plot mean and std error of spectra in each cluster from 1 to 40Hz
figure();
hold all;
for c_idx=1:length(study_info.clusters)

    subplot(length(study_info.clusters),3,1+(c_idx-1)*3);
    % Average over subjects
    mean_spectra=squeeze(nanmean(subj_psd_resids(c_idx,:,:),2));
    stderr_spectra=squeeze(nanstd(subj_psd_resids(c_idx,:,:),[],2))./sqrt(n_subjects);    
    mean_lagged_coh1=squeeze(nanmean(lagged_coh(:,c_idx,:,lag1_idx)));
    stderr_lagged_coh1=squeeze(nanstd(lagged_coh(:,c_idx,:,lag1_idx)))./sqrt(n_subjects);
    mean_lagged_coh2=squeeze(nanmean(lagged_coh(:,c_idx,:,lag2_idx)));
    stderr_lagged_coh2=squeeze(nanstd(lagged_coh(:,c_idx,:,lag2_idx)))./sqrt(n_subjects);
    
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
    plot(study_info.beta_band([1 1]),[0 max_mean_lagged_coh],'k--');
    plot(study_info.beta_band([2 2]),[0 max_mean_lagged_coh],'k--');
        
    xlim(ax(1),frex([1 end]));
    xlim(ax(2),frex([1 end]));
    ylim(ax(1),[min_mean_spectra max_mean_spectra]);
    ylim(ax(2),[0 max_mean_lagged_coh]);
        
    
    % Plot mean lagged coherence
    subplot(length(study_info.clusters),3,[2+(c_idx-1)*3 3+(c_idx-1)*3]);
    contourf(lags,foi,squeeze(nanmean(lagged_coh(:,c_idx,:,:))),100,...
        'linecolor','none');    
    set(gca,'clim',[0 max_mean_lagged_coh]);
    pos=get(gca,'position');
    colorbar();
    set(gca,'position',pos);
    ylabel('Frequency (Hz)');
    xlabel('Lag (cycles)');


end
