function plot_mean_psd(study_info, varargin)
% PLOT_MEAN_PSD - Plot mean power spectrial density in the C3 and C4 
% clusters
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

% Open EEGlab
[ALLEEG, EEG, CURRENTSET] = eeglab;

% Electrode clusters to look in
clusters={'C3','C4'};
cluster_channels={{'E16', 'E20', 'E21', 'E22'},...
    {'E41', 'E49', 'E50', 'E51'}};
cluster_colors={'b','r'};

% Power spectra for each subject in each cluster
spectra=[];

% Number of subjects
n_subjects=size(study_info.participant_info,1);

for s=1:n_subjects
    
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s};
    
    % Path containing subject data
    subject_data_dir=fullfile(study_info.output_dir, subj_id,...
        'eeg');

    % Baseline and experimental epoch files
    base_fname=sprintf('%s_11_Epoch_Matched_CSD_baseline.set',subj_id);
    exp_fname=sprintf('%s_11_Epoch_Matched_CSD_experimental.set',subj_id);
    
    % Load data
    base_EEG=pop_loadset('filepath', subject_data_dir,...
        'filename', base_fname);        
    [ALLEEG, EEG] = eeg_store(ALLEEG, base_EEG, 1);        
    exp_EEG=pop_loadset('filepath', subject_data_dir,...
        'filename', exp_fname);  
    [ALLEEG, EEG] = eeg_store(ALLEEG, exp_EEG, 2);
    
    % Merge baseline and experimental EEG data
    merged_EEG=pop_mergeset(ALLEEG,1:2,0);

    % Data averaged within each cluster
    cluster_data=zeros(length(clusters), merged_EEG.pnts, merged_EEG.trials);

    for c_idx=1:length(clusters)

        % Channels in this cluster
        channels=cluster_channels{c_idx};

        % Find indices of cluster channels
        chan_idx=zeros(1,length(channels));
        for k=1:length(channels)
            chan_idx(k)=find(strcmp({merged_EEG.chanlocs.labels},channels{k}));
        end    

        % Average over channels in this cluster
        cluster_data(c_idx,:,:)=double(squeeze(mean(merged_EEG.data(chan_idx,:,:),1)));
    end
    
    % Use a window size of twice the sampling rate with a 50% overlap
    winsize=merged_EEG.srate*2;
    overlap=round(winsize/2);

    % Compute power spectral density for each cluster using Welch's method
    [subj_spectra,freqs,~,~,~] = spectopo(cluster_data,...
        merged_EEG.pnts, merged_EEG.srate, 'winsize', winsize,...
        'overlap', overlap, 'plot', 'off');
    
    spectra(end+1,:,:)=subj_spectra;    
end

save('../../../data/UMD/adult/deriv/psd.mat','freqs','spectra');

% Plot mean and std error of spectra in each cluster from 1 to 40Hz
freq_idx=(freqs>=1) & (freqs<=40);
figure();
hold all;
for c_idx=1:length(clusters)
    mean_spectra=squeeze(mean(spectra(:,c_idx,freq_idx)));
    stderr_spectra=squeeze(std(spectra(:,c_idx,freq_idx)))./sqrt(n_subjects);
    
    shadedErrorBar(freqs(freq_idx),mean_spectra,stderr_spectra,...
        'LineProps',{'color',cluster_colors{c_idx}});
end
legend(clusters);
xlabel('Frequency (Hz)');
ylabel('log(Power uV^2)');
