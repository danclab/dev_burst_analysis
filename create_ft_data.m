function ft_data=create_ft_data(EEG)
% CREATE_FT_DATA - Create data structure for fieldtrip from EEGlab data
%
% Syntax:  ft_data=create_ft_data(EEG)
%
% Inputs:
%    EEG - EEGLab-formatted data
%
% Outputs:
%    ft_data - FieldTrip-formatted data
%
% Example: 
%   ft_data=create_ft_data(EEG);

% Number of time points per trial
n_pts=size(EEG.data,2);
% Number of trials
n_trials=size(EEG.data,3);

% Initialize data structure
ft_data=[];
% Channel labels given by cluster names
ft_data.label={EEG.chanlocs.labels};
% Sampling rate
ft_data.fsample=EEG.srate;
ft_data.trial={};
ft_data.time={};
% Trial definition
ft_data.cfg.trl=zeros(n_trials,3);
for t_idx=1:n_trials
    % Trial data
    ft_data.trial{end+1}=EEG.data(:,:,t_idx);
    % Convert trial timestamps to seconds
    ft_data.time{end+1}=EEG.times/1000.0;       
    % Create fake trial start time
    ft_data.cfg.trl(t_idx,1)=(t_idx-1)*n_pts+1000;
    % Create fake trial end time
    ft_data.cfg.trl(t_idx,2)=ft_data.cfg.trl(t_idx,1)+n_pts-1;
    % Trial time zero
    ft_data.cfg.trl(t_idx,3)=find(EEG.times==0);
end        
% Empty trialinfo
ft_data.trialinfo=zeros(n_trials,1);