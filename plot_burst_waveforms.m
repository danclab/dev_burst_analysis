function plot_burst_waveforms(study_info, foi, thresh_sd, varargin)
% PLOT_BURST_WAVEFORMS - Plot burst waveforms in the C3 and C4 clusters
%
% Syntax:  plot_burst_waveforms(study_info, foi, woi, thresh_sd)
%
% Inputs:
%    study_info - structure containing the study information
%    foi - frequency band of interest ([low high], Hz)
%    thresh_sd - number of standard deviations above the median amplitude
%        to set the threshold at 
%
% Optional keyword inputs:
%     burst_window - width of time window around burst to extract waveform
%         (ms; default = 100)
%
% Example: 
%   plot_burst_waveforms(study_info, [13 30], 1.5)

% Parse optional arguments
defaults=struct('align',true,'min_ntrials',5);
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

% time window of interest (to cut off filtering edge artifacts)
woi=[-1250 1250];

% Burst window = 3 cycles of FOI center frequency
burst_window=3.5*1/mean(foi(:))*1000;

% Number of subjects
n_subjects=size(study_info.participant_info,1);

% Open EEGlab
[ALLEEG, EEG, CURRENTSET] = eeglab;
%ft_defaults;

% Mean burst waveforms for each subject
subj_waveforms=zeros(n_subjects,round(burst_window/2)+1).*NaN;

for s=1:n_subjects
    
    % Get subject ID from study info
    subj_id=study_info.participant_info.participant_id{s};
    
    % Path containing subject data
    subject_data_dir=fullfile(study_info.deriv_dir, subj_id, 'eeg');
    
    base_fname=sprintf('%s_11_Epoch_Matched_CSD_baseline.set',subj_id);
    exp_fname=sprintf('%s_11_Epoch_Matched_CSD_experimental.set',subj_id);
    
    if exist(fullfile(subject_data_dir,base_fname),'file')==2 &&...
            exist(fullfile(subject_data_dir,exp_fname),'file')==2
        
        % Load data
        base_EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', base_fname);        
        [ALLEEG, EEG] = eeg_store(ALLEEG, base_EEG, 1);        
        exp_EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', exp_fname);  
        [ALLEEG, EEG] = eeg_store(ALLEEG, exp_EEG, 2);

        % If min number of trials in baseline and experimental epochs
        if base_EEG.trials>=params.min_ntrials &&...
                exp_EEG.trials>=params.min_ntrials
        
            % Merge baseline and experimental EEG data
            merged_EEG=pop_mergeset(ALLEEG,1:2,0);      

            % Time stamps in each trial
            all_times=merged_EEG.times;

            % All subject bursts
            all_subj_waveforms=[];
            
            for c_idx=1:length(study_info.clusters)

                % Channels in this cluster
                channels=study_info.cluster_channels{c_idx};

                % Find indices of cluster channels
                chan_idx=zeros(1,length(channels));
                for k=1:length(channels)
                    chan_idx(k)=find(strcmp({merged_EEG.chanlocs.labels},...
                        channels{k}));
                end                       

                % Average over channels in cluster
                cluster_data=double(squeeze(mean(merged_EEG.data(chan_idx,:,:),1)));

                % Extract bursts
                bursts=extract_bursts(cluster_data, all_times,...
                    merged_EEG.srate, foi, thresh_sd,...
                    'burst_window', burst_window);

                % Average burst waveforms within time window
                base_idx=find((bursts.onset_time>=woi(1)) &...
                    (bursts.offset_time<=woi(2)));
                all_subj_waveforms(end+1:end+length(base_idx),:)=bursts.waveform(base_idx,:);
            end
            
            subj_mean_waveform=nanmean(all_subj_waveforms,1);
            subj_mean_waveform=subj_mean_waveform-mean(subj_mean_waveform);

            % Resolve burst direction (oppositely oriented dipoles 
            % will have reversed polarities)
            mid_times=[-20 20];
            mid_idx=knnsearch(bursts.times',mid_times');
            mid_waveform=subj_mean_waveform(mid_idx(1):mid_idx(2));
            [extreme,extreme_idx]=max(abs(mid_waveform));
            if mid_waveform(extreme_idx)>0
                subj_mean_waveform=subj_mean_waveform.*-1;
            end

            % Save mean waveforms
            subj_waveforms(s,:)=subj_mean_waveform;
        end
    end
end

nn_subj_waveforms=subj_waveforms(~all(isnan(subj_waveforms),2),:);
[wout,lags]=woody(nn_subj_waveforms',[],[],'woody','biased');
subj_aligned_waveforms=zeros(size(nn_subj_waveforms)).*NaN;
for i=1:size(nn_subj_waveforms,1)
    if(lags(i)>0)
        subj_aligned_waveforms(i,1:end-(lags(i)-1))=nn_subj_waveforms(i,lags(i):end);
    elseif(lags(i)<0)
        subj_aligned_waveforms(i,(lags(i)*-1)+1:end)=nn_subj_waveforms(i,lags(i)*(-1)+1:end);
    else
        subj_aligned_waveforms(i,:)=nn_subj_waveforms(i,:);
    end
end

figure();
 subplot(1,2,1);
shadedErrorBar(bursts.times,...
    nanmean(subj_waveforms),...
    nanstd(subj_waveforms)./sqrt(n_subjects),...
    'LineProps', {'color','k'});    
xlim(bursts.times([1 end]));
xlabel('Time (ms)');
ylabel('Potential (uV)');

subplot(1,2,2);
shadedErrorBar(bursts.times,...
    nanmean(subj_aligned_waveforms),...
    nanstd(subj_aligned_waveforms)./sqrt(n_subjects),...
    'LineProps', {'color','k'});    
xlim(bursts.times([1 end]));
xlabel('Time (ms)');
ylabel('Potential (uV)');
