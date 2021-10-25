%function bursts=extract_bursts(data, all_times, srate, foi, threshold,...
%    varargin)
function export_bursts(study_info, fois, thresh_sds, varargin)
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
defaults=struct('min_ntrials',5);
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

% Number of subjects
n_subjects=size(study_info.participant_info,1);

fid1=fopen(fullfile(study_info.deriv_dir,'bursts.csv'),'w');
fprintf(fid1, 'Subject,Condition,Trial,Epoch,Cluster,Time,Onset,Offset,Amp\n');
fid2=fopen(fullfile(study_info.deriv_dir,'nbursts.csv'),'w');
fprintf(fid2, 'Subject,Condition,Trial,Epoch,Cluster,Count\n');

for s=1:n_subjects
    
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s};
    
    % Path containing subject data
    subject_data_dir=fullfile(study_info.deriv_dir, subj_id, 'eeg');
    
    % Baseline and experimental epoch files
    base_fname=sprintf('%s_11_Epoch_Matched_CSD_baseline.set',subj_id);
    exp_fname=sprintf('%s_11_Epoch_Matched_CSD_experimental.set',subj_id);
    
    if exist(fullfile(subject_data_dir,base_fname),'file')==2 &&...
            exist(fullfile(subject_data_dir,exp_fname),'file')==2
        
        % Load data
        base_EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', base_fname);        
        exp_EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', exp_fname); 
            
        % Time stamps in each trial
        all_base_times=base_EEG.times;
        all_exp_times=exp_EEG.times;

        for cond_idx=1:length(study_info.conditions)
            condition=study_info.conditions{cond_idx};

            base_event=study_info.baseline_evts{cond_idx};
            exp_event=study_info.exp_evts{cond_idx};
            
            % Get condition trials
            base_trials=find(strcmp({base_EEG.event.type},base_event));
            exp_trials=find(strcmp({exp_EEG.event.type},exp_event));

            % If min number of trials in baseline and experimental epochs
            if length(base_trials)>=params.min_ntrials &&...
                    length(exp_trials)>=params.min_ntrials

                % Process each cluster
                for c_idx=1:length(study_info.clusters)     
                    cluster=study_info.clusters{c_idx};

                    % Get cluster channels
                    channels=study_info.cluster_channels{c_idx};
                    chan_idx=cellfun(@(x) find(strcmp({base_EEG.chanlocs.labels},x)),...
                        channels);

                    % Compute amplitude in each channel of C3 cluster
                    cluster_base_amp=zeros(length(chan_idx),length(all_base_times),length(base_trials));
                    cluster_exp_amp=zeros(length(chan_idx),length(all_exp_times),length(exp_trials));

                    for i=1:length(chan_idx)
                        % Data for this cluster and condition
                        base_data=squeeze(base_EEG.data(chan_idx(i),:,base_trials));
                        exp_data=squeeze(exp_EEG.data(chan_idx(i),:,exp_trials));

                        % Get amplitude
                        [~, ch_base_amp]=filter_hilbert(base_data, base_EEG.srate, fois(c_idx,:));
                        [~, ch_exp_amp]=filter_hilbert(exp_data, exp_EEG.srate, fois(c_idx,:));

                        cluster_base_amp(i,:,:)=ch_base_amp;
                        cluster_exp_amp(i,:,:)=ch_exp_amp;
                    end

                    % Average amplitude over cluster channels
                    cluster_base_amp=squeeze(mean(cluster_base_amp));
                    cluster_exp_amp=squeeze(mean(cluster_exp_amp));
                    
                    all_amp=[cluster_base_amp cluster_exp_amp];

                    % Compute absolute burst threshold
                    threshold=median(all_amp(:))+(thresh_sds(c_idx)*std(all_amp(:)));

                    % Extract bursts
                    extract_epoch_bursts(subj_id, condition, cluster, 'baseline',...
                        base_trials, all_base_times, cluster_base_amp, threshold, fid1, fid2);
                    extract_epoch_bursts(subj_id, condition, cluster, 'exp',...
                        exp_trials, all_exp_times, cluster_exp_amp, threshold, fid1, fid2);                    
                end
            end
        end
    end
end             
fclose(fid1);
fclose(fid2);

end


function extract_epoch_bursts(subj_id, condition, cluster, epoch, trials,...
    times, cluster_amp, threshold, fid1, fid2)
    for t_idx=1:size(cluster_amp,2)

        trial=trials(t_idx);
        
        % All times when amp is over threshold
        over_thresh_idx = cluster_amp(:,t_idx)>=threshold;
        % Change in threshold crossing
        over_thresh_diff = diff(over_thresh_idx);
        % All times when amp is over threshold and previous time is
        % under threshold
        all_burst_start_idx=find(over_thresh_diff)+1;

        %% Process each burst
        n_bursts=0;
        for k=1:length(all_burst_start_idx)
            burst_start_idx=all_burst_start_idx(k);
            burst_start_time=times(burst_start_idx);
            
            % Find next time amplitude goes below threshold
            burst_end_idx=find(cluster_amp(burst_start_idx+1:end,t_idx)<threshold,1);

            if isempty(burst_end_idx) || burst_end_idx>1
                
                % If it goes below threshold before the end of the trial
                if ~isempty(burst_end_idx)% && burst_end_idx>1
                    % Time when amplitude goes back below threshold
                    burst_end_idx=burst_start_idx+burst_end_idx;                
                    burst_end_time=times(burst_end_idx);
                else
                    burst_end_idx=size(cluster_amp,1);
                    burst_end_time=NaN;
                end
                % Find peak amplitude
                burst_amp=cluster_amp(burst_start_idx:burst_end_idx,t_idx);
                [max_burst_amp,burst_peak_idx]=max(burst_amp);
                % Peak time
                burst_peak_idx=burst_peak_idx+burst_start_idx-1;
                burst_peak_time=times(burst_peak_idx);

                % Save burst information
                fprintf(fid1, '%s,%s,%d,%s,%s,%.3f,%.3f,%.3f,%.3f\n',...
                    subj_id,...
                    condition,...
                    trial,...
                    epoch,...
                    cluster,...
                    burst_peak_time,...
                    burst_start_time,...
                    burst_end_time,...
                    max_burst_amp);
                n_bursts=n_bursts+1;
            end
        end
        fprintf(fid2,'%s,%s,%d,%s,%s,%d\n',...
            subj_id,...
            condition,...
            trials(t_idx),...
            epoch,...
            cluster,...
            n_bursts);
    end
end