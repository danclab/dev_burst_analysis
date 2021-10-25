function plot_mean_psd(study_info, varargin)
% PLOT_MEAN_PSD - Plot mean power spectrial density
% in the C3 and C4 clusters
%
% Syntax:  plot_mean_psd(study_info)
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

subj_id=study_info.participant_info.participant_id{1};
% Path containing subject data
subject_data_dir=fullfile(study_info.deriv_dir, subj_id, 'eeg');    
% Baseline and experimental epoch files
base_fname=sprintf('%s_11_Epoch_Matched_CSD_baseline.set',subj_id);
EEG=pop_loadset('filepath', subject_data_dir, 'filename', base_fname);  

% Load power spectral densities
load(fullfile(study_info.deriv_dir,'psd.mat'));

ch_space=[[EEG.chanlocs.X];[EEG.chanlocs.Y];[EEG.chanlocs.Z]]';
for i=1:3
    ch_space(:,i)=ch_space(:,i)-min(ch_space(:,i));
    ch_space(:,i)=ch_space(:,i)./max(ch_space(:,i));
end


fig=figure();
hold all
for ch=1:length(EEG.chanlocs)
    shadedErrorBar(frex,squeeze(mean(periodic(:,ch,:))),squeeze(std(periodic(:,ch,:)))/sqrt(size(periodic,1)),'LineProps',{'color',ch_space(ch,:)});
end
yl=ylim();
ylim([-2 yl(2)]);
xlabel('Frequency (Hz)');
ylabel('log(power)');
axes('Position',[.65 0.65 .2 .2]);
plot_sensors(zeros(1,length(EEG.chanlocs)),EEG.chanlocs,'electrodes','on',...
    'style','blank','emarker','o','electcolor',ch_space,'whitebk','on');
drawnow;
