function bursts=extract_bursts(data, all_times, srate, foi, threshold,...
    varargin)
% EXTRACT_BURSTS - Extracts burst from the data as amplitude within a
% frequency band exceeding a threshold
%
% Syntax:  bursts = extract_bursts(data, all_times, srate, foi, threshold);
%
% Inputs:
%    data - trial data (time x trials)
%    all_times - timestamps (ms)
%    srate - sampling rate (Hz)
%    foi - frequency band of interest ([low high], Hz)
%    threshold - amplitude threshold
%
% Outputs:
%    bursts - structure containing data for each burst:
%        trial - trial burst occurred in
%        peak_time - time of burst peak (ms)
%        onset_time - onset of burst (ms)
%        offset_time - offset of burst (ms)
%
% Example: 
%   bursts = extract_bursts(exp_data, all_exp_times, 250, [13 30], threshold);

% Parse optional arguments
defaults=struct();
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

% Get amplitude
[~, amp_data]=filter_hilbert(data, srate, foi);

bursts=[];
% For each burst - trial it occured in
bursts.trial=[];
% For each burst - time within trial it occured
bursts.peak_time=[];
% For each burst - onset time within trial
bursts.onset_time=[];
% For each burst - offset time within trial
bursts.offset_time=[];
% Burst peak amplitude
bursts.peak_amp=[];

% Go through signal each trial 
for t_idx=1:size(data,2)

    % All times when amp is over threshold
    over_thresh_idx = amp_data(:,t_idx)>=threshold;
    % Change in threshold crossing
    over_thresh_diff = diff(over_thresh_idx);
    % All times when amp is over threshold and previous time is
    % under threshold
    all_burst_start_idx=find(over_thresh_diff)+1;

    %% Process each burst
    for k=1:length(all_burst_start_idx)
        burst_start_idx=all_burst_start_idx(k);

        % Find next time amplitude goes below threshold
        burst_end_idx=find(amp_data(burst_start_idx:end,t_idx)<threshold,1);
        
        % If it goes below threshold before the end of the trial
        if ~isempty(burst_end_idx)
            
            % Time when amplitude goes back below threshold
            burst_end_idx=burst_start_idx+burst_end_idx-1;
        
            % Find peak amplitude
            burst_amp=amp_data(burst_start_idx:burst_end_idx,t_idx);
            [max_burst_amp,burst_peak_idx]=max(burst_amp);
            burst_peak_idx=burst_peak_idx+burst_start_idx-1;

            % Save burst information
            bursts.trial(end+1)=t_idx;
            bursts.peak_time(end+1)=all_times(burst_peak_idx);
            bursts.onset_time(end+1)=all_times(burst_start_idx);
            bursts.offset_time(end+1)=all_times(burst_end_idx);
            bursts.peak_amp(end+1)=max_burst_amp;                        
        end
    end
end
