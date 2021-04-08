function plot_single_trial_data(study_info)

% Get subject ID from study info
subj_id=study_info.participant_info.participant_id{1};
    
% Path containing subject data
subject_data_dir=fullfile(study_info.deriv_dir, subj_id, 'eeg');
    
% Baseline and experimental epoch files
base_fname=sprintf('%s_11_Epoch_Matched_CSD_baseline.set',subj_id);
exp_fname=sprintf('%s_11_Epoch_Matched_CSD_experimental.set',subj_id);
    
% Load data
base_EEG=pop_loadset('filepath', subject_data_dir,...
    'filename', base_fname);        
exp_EEG=pop_loadset('filepath', subject_data_dir,...
    'filename', exp_fname);        

% Time stamps in each trial
all_base_times=base_EEG.times;
all_exp_times=exp_EEG.times;

% Sampling rate
srate=base_EEG.srate;

% Data averaged within each cluster
base_cluster_data=get_cluster_data(study_info, base_EEG);                        
exp_cluster_data=get_cluster_data(study_info, exp_EEG);
            
% C3 cluster data
base_data=squeeze(base_cluster_data(1,:,:));
exp_data=squeeze(exp_cluster_data(1,:,:));

% Get amplitude using the filter-Hilbert method
[filt_base_beta_data, amp_base_beta_data]=filter_hilbert(base_data, srate,...
    study_info.beta_band);
[filt_exp_beta_data, amp_exp_beta_data]=filter_hilbert(exp_data, srate,...
    study_info.beta_band);
[filt_base_mu_data, amp_base_mu_data]=filter_hilbert(base_data, srate,...
    study_info.mu_band);
[filt_exp_mu_data, amp_exp_mu_data]=filter_hilbert(exp_data, srate,...
    study_info.mu_band);

figure();
subplot(1,2,1);
hold all;
plot(all_base_times,base_data(:,1));
plot(all_base_times,filt_base_beta_data(:,1));
plot(all_base_times,amp_base_beta_data(:,1));
plot(all_base_times,filt_base_mu_data(:,1));
plot(all_base_times,amp_base_mu_data(:,1));
xlim(all_base_times([1 end]));
legend('raw','beta','beta amplitude','mu','mu amplitude');
xlabel('Time (ms)');
ylabel('Potential \muV');

subplot(1,2,2);
hold all;
plot(all_exp_times,exp_data(:,1))
plot(all_exp_times,filt_exp_beta_data(:,1));    
plot(all_exp_times,amp_exp_beta_data(:,1));
plot(all_exp_times,filt_exp_mu_data(:,1));    
plot(all_exp_times,amp_exp_mu_data(:,1));
xlim(all_exp_times([1 end]));
xlabel('Time (ms)');
ylabel('Potential \muV');
