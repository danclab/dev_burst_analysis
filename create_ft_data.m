function data=create_ft_data(clusters, cluster_data, times, srate)
% CREATE_FT_DATA - Create data structure for fieldtrip from cluster data
%
% Syntax:  data=create_ft_data(clusters, cluster_data, times, srate)
%
% Inputs:
%    clusters - cell array of cluster names
%    cluster_data - data averaged over channels in each cluster (cluster x
%        time x trial)
%    times - timestamps for each trial (ms)
%    srate - sampling rate (Hz)
%
% Outputs:
%    filtered_data - bandpass-filtered data in the frequency range (time x
%        trials)
%    amp - amplitude in the frequency range (time x trials)
%
% Example: 
%   data=create_ft_data({'C3','C4'}, cluster_data, times, 250)

% Number of time points per trial
n_pts=size(cluster_data,2);
% Number of trials
n_trials=size(cluster_data,3);

% Initialize data structure
data=[];
% Channel labels given by cluster names
data.label=clusters;
% Sampling rate
data.fsample=srate;
data.trial={};
data.time={};
% Trial definition
data.cfg.trl=zeros(n_trials,3);
for t_idx=1:n_trials
    % Trial data
    data.trial{end+1}=cluster_data(:,:,t_idx);
    % Convert trial timestamps to seconds
    data.time{end+1}=times/1000.0;       
    % Create fake trial start time
    data.cfg.trl(t_idx,1)=(t_idx-1)*n_pts+1000;
    % Create fake trial end time
    data.cfg.trl(t_idx,2)=data.cfg.trl(t_idx,1)+n_pts-1;
    % Trial time zero
    data.cfg.trl(t_idx,3)=find(times==0);
end        
% Empty trialinfo
data.trialinfo=zeros(n_trials,1);