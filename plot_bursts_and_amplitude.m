function plot_bursts_and_amplitude(varargin)
% PLOT_BURSTS_AND_AMPLITUDE - Plot distribution of burst timing and
% amplitude time course within frequency band in the C3 and C4 clusters
% during the execution condition
%
% Syntax:  plot_bursts_and_amplitude()
%
% Example:
%   plot_bursts_and_amplitude()

% Parse optional arguments
defaults=struct('min_ntrials',5);
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

% Load study infos
study_infos={};
study_infos{1}=init_umd12m_study_info();
study_infos{2}=init_umdadult_study_info();

fig=figure();

for st=1:length(study_infos)
    study_info=study_infos{st};
    
    % Plot execute condition
    condition='execute';
    cond_idx=find(strcmp(study_info.conditions,condition));
    base_event=study_info.baseline_evts{cond_idx};
    exp_event=study_info.exp_evts{cond_idx};

    % time window of interest (to cut off filtering edge artifacts)
    woi=[-1250 1250];

    % Number of subjects
    n_subjects=size(study_info.participant_info,1);

    % Baseline mean amplitude for each subject in each cluster
    subj_amp_base=[];
    % Experimental mean amplitude for each subject in each cluster
    subj_amp_exp=[];
    % Baseline mean burst rate for each subject in each cluster
    subj_rate_base=[];
    % Experimental mean burst rate for each subject in each cluster
    subj_rate_exp=[];
    % Number of baseline trials per subject
    ntrials_base=zeros(1,n_subjects);
    % Number of experimental trials per subject
    ntrials_exp=zeros(1,n_subjects);

    % Load burst data
    bursts=readtable(fullfile(study_info.deriv_dir,'bursts.csv'));
    % Bin and smoothing width
    binwidth=10;
    smoothw=25;

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

            % Compute bin centers
            base_bins=[all_base_times(1):binwidth:all_base_times(end)];
            base_bins=base_bins(1:end-1)+binwidth/2;
            exp_bins=[all_exp_times(1):binwidth:all_exp_times(end)];
            exp_bins=exp_bins(1:end-1)+binwidth/2;

            % Get condition trials
            base_trials=find(strcmp({base_EEG.event.type},base_event));
            exp_trials=find(strcmp({exp_EEG.event.type},exp_event));

            % Number of trials
            ntrials_base(s)=length(base_trials);
            ntrials_exp(s)=length(exp_trials);

            % If min number of trials in baseline and experimental epochs
            if ntrials_base(s)>=params.min_ntrials &&...
                    ntrials_exp(s)>=params.min_ntrials

                % Process each cluster
                for c_idx=1:length(study_info.clusters)
                    cluster=study_info.clusters(c_idx);

                    channels=study_info.cluster_channels{c_idx};
                    chan_idx=cellfun(@(x) find(strcmp({base_EEG.chanlocs.labels},x)),...
                        channels);

                    % Compute amplitude in each channel of luster
                    cluster_base_amp=zeros(length(chan_idx),length(all_base_times),length(base_trials));
                    cluster_exp_amp=zeros(length(chan_idx),length(all_exp_times),length(exp_trials));

                    for i=1:length(chan_idx)
                        chan_base_data=squeeze(base_EEG.data(chan_idx(i),:,base_trials));
                        [~, ch_base_amp]=filter_hilbert(chan_base_data, base_EEG.srate, study_info.beta_band(c_idx,:));
                        cluster_base_amp(i,:,:)=ch_base_amp;

                        chan_exp_data=squeeze(exp_EEG.data(chan_idx(i),:,exp_trials));
                        [~, ch_exp_amp]=filter_hilbert(chan_exp_data, exp_EEG.srate, study_info.beta_band(c_idx,:));
                        cluster_exp_amp(i,:,:)=ch_exp_amp;
                    end

                    % Average amplitude over cluster channels
                    cluster_base_amp=squeeze(mean(cluster_base_amp));
                    cluster_exp_amp=squeeze(mean(cluster_exp_amp));

                    % Baseline-correct amplitude
                    amp_woi_idx=knnsearch(all_base_times',woi');
                    cluster_base_amp=cluster_base_amp(amp_woi_idx(1):amp_woi_idx(2),:);
                    cluster_exp_amp=cluster_exp_amp(amp_woi_idx(1):amp_woi_idx(2),:);
                    mean_amp_base_cluster_data=mean(mean(cluster_base_amp));
                    subj_amp_base(s,c_idx,:)=nanmean(cluster_base_amp-mean_amp_base_cluster_data,2);
                    subj_amp_exp(s,c_idx,:)=nanmean(cluster_exp_amp-mean_amp_base_cluster_data,2);

                    % Compute burst rate baseline
                    rate_base=zeros(length(base_bins),length(base_trials));
                    for t=1:length(base_trials)
                        rows=strcmp(bursts.Subject,subj_id) & strcmp(bursts.Epoch,'baseline') & (bursts.Trial==base_trials(t)) & strcmp(bursts.Cluster,cluster);
                        [n,bin]=histc(bursts.Time(rows),base_bins);
                        rate_base(:,t)=n;
                    end
                    rate_base=rate_base*(1000/binwidth);
                    
                    % Smooth baseline burst rate
                    smoothed_rate_base=zeros(size(rate_base));
                    w=gausswin(smoothw);
                    for t=1:length(base_trials)
                        smoothed_rate_base(:,t)=filtfilt(w,1,squeeze(rate_base(:,t)));
                    end
                    rate_woi_idx=knnsearch(base_bins',woi');
                    smoothed_rate_base=smoothed_rate_base(rate_woi_idx(1):rate_woi_idx(2),:);
                    
                    % Compute experimental burst rate
                    rate_exp=zeros(length(exp_bins),length(exp_trials));
                    for t=1:length(exp_trials)
                        rows=strcmp(bursts.Subject,subj_id) & strcmp(bursts.Epoch,'exp') & (bursts.Trial==exp_trials(t)) & strcmp(bursts.Cluster,cluster);
                        [n,bin]=histc(bursts.Time(rows),exp_bins);
                        rate_exp(:,t)=n;
                    end
                    rate_exp=rate_exp*(1000/binwidth);
                    
                    % Smooth experimental burst rate
                    smoothed_rate_exp=zeros(size(rate_exp));
                    w=gausswin(smoothw);
                    for t=1:length(exp_trials)
                        smoothed_rate_exp(:,t)=filtfilt(w,1,squeeze(rate_exp(:,t)));
                    end
                    smoothed_rate_exp=smoothed_rate_exp(rate_woi_idx(1):rate_woi_idx(2),:);
                    
                    % Baseline correct
                    mean_rate_base_cluster_data=mean(mean(smoothed_rate_base));
                    subj_rate_base(s,c_idx,:)=nanmean(smoothed_rate_base-mean_rate_base_cluster_data,2);
                    subj_rate_exp(s,c_idx,:)=nanmean(smoothed_rate_exp-mean_rate_base_cluster_data,2);                    
                end
            end
        end
    end

    % Plot each cluster
    for c_idx=1:length(study_info.clusters)
        cluster=study_info.clusters{c_idx};

        % Raster plot of baseline bursts in all trials
        subplot(3*length(study_infos),...
            length(study_info.clusters)*2,...
            [(st-1)*length(study_infos)*3*length(study_info.clusters)+(c_idx-1)*length(study_info.clusters)+1 
            (st-1)*length(study_infos)*3*length(study_info.clusters)+(c_idx-1)*length(study_info.clusters)+5]);
        hold all
        trial_offset=0;
        for s=1:n_subjects
            subj_id=study_info.participant_info.participant_id{s};
            rows=strcmp(bursts.Subject,subj_id) & strcmp(bursts.Epoch,'baseline') & strcmp(bursts.Cluster,cluster);
            if any(rows)
                plot(bursts.Time(rows),bursts.Trial(rows)+trial_offset,'.k');
                trial_offset=trial_offset+ntrials_base(s);
                if s<n_subjects
                    plot(woi,[trial_offset-.5 trial_offset-.5],'r--');
                end
            end
        end
        xlim(woi);
        ylim([0 trial_offset]);
        ylabel('Trial');

        % Raster plot of experimental bursts in all trials
        subplot(3*length(study_infos),...
            length(study_info.clusters)*2,...
            [(st-1)*length(study_infos)*3*length(study_info.clusters)+(c_idx-1)*length(study_info.clusters)+2 
            (st-1)*length(study_infos)*3*length(study_info.clusters)+(c_idx-1)*length(study_info.clusters)+6]);
        hold all
        trial_offset=0;
        for s=1:n_subjects
            subj_id=study_info.participant_info.participant_id{s};
            rows=strcmp(bursts.Subject,subj_id) & strcmp(bursts.Epoch,'exp') & strcmp(bursts.Cluster,cluster);
            if any(rows)
                plot(bursts.Time(rows),bursts.Trial(rows)+trial_offset,'.k');
                trial_offset=trial_offset+ntrials_exp(s);
                if s<n_subjects
                    plot(woi,[trial_offset-.5 trial_offset-.5],'r--');
                end
            end
        end
        xlim(woi);
        ylim([0 trial_offset]);
        ylabel('Trial');

        % Mean and std error of baseline burst rate across subjects
        mean_base_rate=squeeze(nanmean(subj_rate_base(:,c_idx,:)));
        stderr_base_rate=squeeze(nanstd(subj_rate_base(:,c_idx,:),[],1))./sqrt(n_subjects);
        % Mean and std error of experimental burst rate across subjects
        mean_exp_rate=squeeze(nanmean(subj_rate_exp(:,c_idx,:)));
        stderr_exp_rate=squeeze(nanstd(subj_rate_exp(:,c_idx,:),[],1))./sqrt(n_subjects);

        % Mean and std error of baseline amplitude across subjects
        mean_base_amp=squeeze(nanmean(subj_amp_base(:,c_idx,:)));
        stderr_base_amp=squeeze(nanstd(subj_amp_base(:,c_idx,:)))./sqrt(n_subjects);
        % Mean and std error of experimental amplitude across subjects
        mean_exp_amp=squeeze(nanmean(subj_amp_exp(:,c_idx,:)));
        stderr_exp_amp=squeeze(nanstd(subj_amp_exp(:,c_idx,:)))./sqrt(n_subjects);

        % Plot baseline burst probability and amplitude
        subplot(3*length(study_infos),...
            length(study_info.clusters)*2,...
            (st-1)*length(study_infos)*3*length(study_info.clusters)+(c_idx-1)*length(study_info.clusters)+9);
        ax=shadedErrorBaryy(base_bins(rate_woi_idx(1):rate_woi_idx(2)),...
            mean_base_rate, stderr_base_rate, 'b',...
            all_base_times(amp_woi_idx(1):amp_woi_idx(2)),...
            mean_base_amp,stderr_base_amp, 'r');
        xlim(ax(1),woi);
        set(ax(1), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
        ylim(ax(1),[-120 60]);
        xlim(ax(2),woi);
        set(ax(2), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
        ylim(ax(2),[-2.5 1.25]);
        xlabel('Time (ms)');
        ylabel(ax(1),'Burst rate');

        % Plot experimental burst probability and amplitude
        subplot(3*length(study_infos),...
            length(study_info.clusters)*2,...
            (st-1)*length(study_infos)*3*length(study_info.clusters)+(c_idx-1)*length(study_info.clusters)+10);
        ax=shadedErrorBaryy(exp_bins(rate_woi_idx(1):rate_woi_idx(end)),...
            mean_exp_rate, stderr_exp_rate, 'b',...
            all_exp_times(amp_woi_idx(1):amp_woi_idx(2)),...
            mean_exp_amp, stderr_exp_amp, 'r');
        xlim(ax(1),woi);
        set(ax(1), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
        ylim(ax(1),[-120 60]);
        xlim(ax(2),woi);
        set(ax(2), 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
        ylim(ax(2),[-2.5 1.25]);
        xlabel('Time (ms)');
        ylabel(ax(2),'Amplitude');
    end
end