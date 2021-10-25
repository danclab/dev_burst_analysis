% Open EEGlab
[ALLEEG, EEG, CURRENTSET] = eeglab;

study_info=init_umdadult_study_info();
subj_id=study_info.participant_info.participant_id{1};

% Path containing subject data
subject_data_dir=fullfile(study_info.deriv_dir, subj_id, 'eeg');

% Baseline and experimental epoch files
base_fname=sprintf('%s_11_Epoch_Matched_CSD_baseline.set',subj_id);
exp_fname=sprintf('%s_11_Epoch_Matched_CSD_experimental.set',subj_id);
    
% Load data
base_EEG=pop_loadset('filepath', subject_data_dir,...
    'filename', base_fname);        
[ALLEEG, EEG] = eeg_store(ALLEEG, base_EEG, 1);        
exp_EEG=pop_loadset('filepath', subject_data_dir,...
    'filename', exp_fname);  
[ALLEEG, EEG] = eeg_store(ALLEEG, exp_EEG, 2);

% Merge baseline and experimental EEG data
merged_EEG=pop_mergeset(ALLEEG,1:2,0);

% Get C3 channels
c3_idx=find(strcmp(study_info.clusters,'C3'));
channels=study_info.cluster_channels{c3_idx};
chan_idx=cellfun(@(x) find(strcmp({base_EEG.chanlocs.labels},x)),...
    channels);

% Use a window size (in datapoints) equal to the sampling rate (so 1s) with 
% a 50% overlap
winsize=merged_EEG.srate;
overlap=round(winsize/2);

% Compute power spectral density for each channel using Welch's
% method
[spectra,frex,~,~,~] = spectopo(merged_EEG.data, merged_EEG.pnts,...
    merged_EEG.srate, 'winsize', winsize, 'overlap', overlap,...
    'plot', 'off', 'freqfac',2);

% Get frequencies from 0.5 to 40
freq_idx=find((frex>=0.5) & (frex<=40));
frex=frex(freq_idx);
spectra=spectra(:,freq_idx);

% Subject spectrum
psd=mean(spectra(chan_idx,:));

% 1/f
oof=log10(1./frex);
% Fit 1/f to spectrum
lm_psd=fitlm(oof,psd,'RobustOpts','on');
% Get residuals
resids=lm_psd.Residuals.Raw;
    
% Find peaks and sort in descending order
[pks,locs]=findpeaks(resids,frex);    
[sorted_pks,sort_idx]=sort(pks,1,'descend');
sorted_locs=locs(sort_idx);

% Peak freqs and FWHMs
pk_freqs=zeros(1,length(sorted_pks));
fwhms=zeros(1,length(sorted_pks));

% For each peak - fit Gaussian
for p_idx=1:length(sorted_pks)

    % Peak frequency and residual power
    pk_freq=sorted_locs(p_idx);
    pk_idx=find(frex==pk_freq);
    pk_pow=resids(pk_idx);
    
    % Find FWHM
    % Index of frequency with half peak power on the left
    l_idx=find(resids(1:pk_idx)<pk_pow*.5);
    l_freq=frex(l_idx(end));
    % Index of frequency with half peak power on right
    r_idx=find(resids(pk_idx+1:end)<pk_pow*.5);
    r_freq=frex(pk_idx-1+r_idx(1));
    % Sum the two to compute fwhm
    fwhm=(r_freq-pk_freq)+(pk_freq-l_freq);

    pk_freqs(p_idx)=pk_freq;
    fwhms(p_idx)=fwhm;
end

% Plot power spectra
fig=figure();
hold all
plot(frex,psd);
plot(frex,lm_psd.Coefficients.Estimate(1)+...
    lm_psd.Coefficients.Estimate(2).*log10(1./frex));

for p_idx=1:length(sorted_pks)
    % Band limits
    l_freq=pk_freqs(p_idx)-fwhms(p_idx)*.5;
    r_freq=pk_freqs(p_idx)+fwhms(p_idx)*.5;
    disp(sprintf('peak freq=%.2fHz, fwhm=%.2fHz, range=%.2f-%.2fHz',...
        pk_freqs(p_idx), fwhms(p_idx), l_freq, r_freq));
    plot(pk_freqs(p_idx),psd(frex==pk_freqs(p_idx)),'or');
    plot([l_freq l_freq],ylim(),'--','color','k');
    plot([r_freq r_freq],ylim(),'--','color','k');
end    
legend('PSD','Aperiodic component');    
xlabel('Frequency (Hz)');
ylabel('log(power)');
    
    