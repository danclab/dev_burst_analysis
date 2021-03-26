function threshold=compute_threshold(base_data, exp_data, all_base_times,...
    all_exp_times, srate, foi, woi, thresh_sd, varargin)
% COMPUTE_THRESHOLD - Compute absolute amplitude threshold for a single
% participant within a single electrode or electrode cluster. Combines
% baseline and experimental trials and computes median + thresh_sd*stddev...
% of amplitude
%
% Syntax:  threshold = compute_threshold(base_data, exp_data,...
%                          all_base_times, all_exp_times, srate, foi,...
%                          woi, thresh_sd)
%
% Inputs:
%    base_data - baseline trial data (time x trials)
%    exp_data - experimental trial data (time x trials)
%    all_base_times - timestamps in baseline trials (ms)
%    all_exp_times - timestamps in experimental trials (ms)
%    srate - sampling rate (Hz)
%    foi - frequency band of interest ([low high], Hz)
%    woi - time window of interest (to cut off filtering edge artifacts;
%        [low high], ms)
%    thresh_sd - number of standard deviations above the median amplitude
%        to set the threshold at 
%
% Outputs:
%    threshold - absolute threshold for defining burst - median +
%        thresh_sd*std dev of amplitude in all time points / trials
%
% Example: 
%   threshold = compute_threshold(base_data, exp_data, all_base_times,...
%                   all_exp_times, 250, [13 30], [-1000 1000], 1.5)

% Parse optional arguments
defaults=struct('plot', false);
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

% Get amplitude using the filter-Hilbert method
[~, amp_base_data]=filter_hilbert(base_data, srate, foi);
[~, amp_exp_data]=filter_hilbert(exp_data, srate, foi);

% Cut off edges to remove filter artifacts
base_woi_idx=knnsearch(all_base_times',woi');
amp_base_data=amp_base_data(base_woi_idx(1):base_woi_idx(2),:);

% Cut off edges to remove filter artifacts
exp_woi_idx=knnsearch(all_exp_times',woi');
amp_exp_data=amp_exp_data(exp_woi_idx(1):exp_woi_idx(2),:);

% Compute threshold
amp_all_data=[amp_base_data amp_exp_data];
threshold=median(amp_all_data(:))+(thresh_sd*std(amp_all_data(:)));
