function [filtered_data, amp]=filter_hilbert(data, srate, foi)
% FILTER_HILBERT - Compute amplitude within a frequency band using the
% filter-Hilbert method
%
% Syntax:  [filtered_data, amp]=filter_hilbert(data, srate, foi)
%
% Inputs:
%    data - trial data (time x trials)
%    srate - sampling rate (Hz)
%    foi - frequency band of interest ([low high], Hz)
%
% Outputs:
%    filtered_data - bandpass-filtered data in the frequency range (time x
%        trials)
%    amp - amplitude in the frequency range (time x trials)
%
% Example: 
%   [filtered_data, amp]= filter_hilbert(exp_data, 250, [13 30])

% Bandpass filter in frequency range
% 6th order two-pass Butterworth filter
% reduce order in case of instability
filtered_data = ft_preproc_bandpassfilter(data', srate, foi, 6, 'but',...
    'twopass', 'reduce')';

% Get amplitude from Hilbert transform
amp=abs(hilbert(filtered_data));
