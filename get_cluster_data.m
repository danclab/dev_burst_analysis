function cluster_data=get_cluster_data(study_info, EEG)
% GET_CLUSTER_DATA - Averages over channels within each cluster
%
% Syntax:  cluster_data=get_cluster_data(study_info, EEG)
%
% Inputs:
%    study_info - structure containing the study information
%    EEG - EEGlab dataset
%
% Outputs:
%    cluster_data - data averaged over channels in each cluster (cluster x
%        time x trial)
%
% Example: 
%   cluster_data=get_cluster_data(study_info, EEG)

% Data averaged within each cluster
cluster_data=zeros(length(study_info.clusters), EEG.pnts, EEG.trials);

for c_idx=1:length(study_info.clusters)

    % Channels in this cluster
    channels=study_info.cluster_channels{c_idx};

    % Find indices of cluster channels
    chan_idx=zeros(1,length(channels));
    for k=1:length(channels)
        chan_idx(k)=find(strcmp({EEG.chanlocs.labels},channels{k}));
    end    

    % Average over channels in this cluster
    cluster_data(c_idx,:,:)=double(squeeze(mean(EEG.data(chan_idx,:,:),1)));
end