function fois=determine_fois(study_info)

subj_id=study_info.participant_info.participant_id{1};
% Path containing subject data
subject_data_dir=fullfile(study_info.deriv_dir, subj_id, 'eeg');
% Baseline and experimental epoch files
base_fname=sprintf('%s_11_Epoch_Matched_CSD_baseline.set',subj_id);
EEG=pop_loadset('filepath', subject_data_dir, 'filename', base_fname);

% Load power spectral densities
load(fullfile(study_info.deriv_dir,'psd.mat'));

% Number of subjects
n_subjects=size(spectra,1);

% Look at residuals from 5-40Hz because we don't have lagged coherence at
% lower frequencies (trials aren't long enough)
freq_idx=find(frex>=5);
frex=frex(freq_idx);
periodic=periodic(:,:,freq_idx);

fois={};
for c=1:length(study_info.clusters)
    
    % Get cluster channels
    channels=study_info.cluster_channels{c_idx};
    chan_idx=cellfun(@(x) find(strcmp({EEG.chanlocs.labels},x)),...
        channels);
    
    mean_resids=squeeze(nanmean(nanmean(periodic(:,chan_idx,:),2),1));
    
    % Find peaks and sort in descending order
    [pks,locs]=findpeaks(mean_resids,frex);
    [sorted_pks,sort_idx]=sort(pks,1,'descend');
    sorted_locs=locs(sort_idx);
    
    mean_resids=mean_resids-min(mean_resids);
    
    c_fois=[];
    % For each peak - fit Gaussian
    for p_idx=1:length(sorted_pks)
        
        % Peak frequency and residual power
        pk_freq=sorted_locs(p_idx);
        pk_idx=find(frex==pk_freq);
        pk_pow=mean_resids(pk_idx);
        
        % Find fwhm
        l_idx=find(mean_resids(1:pk_idx)<pk_pow*.5);
        r_idx=find(mean_resids(pk_idx+1:end)<pk_pow*.5);
        if length(l_idx) && length(r_idx)
            l_freq=frex(l_idx(end));
            r_freq=frex(pk_idx-1+r_idx(1));
            fwhm=(r_freq-pk_freq)+(pk_freq-l_freq);
        elseif length(l_idx)
            l_freq=frex(l_idx(end));
            fwhm=2*(pk_freq-l_freq);
        elseif length(r_idx)
            r_freq=frex(pk_idx-1+r_idx(1));
            fwhm=2*(r_freq-pk_freq);
        end
        
        % Band limits
        l_freq=pk_freq-fwhm*.5;
        r_freq=pk_freq+fwhm*.5;
        if l_freq>=frex(1) && r_freq<=frex(end)
            c_fois(end+1,:)=[l_freq r_freq];
            disp(sprintf('%s - %s: peak freq=%.2fHz, fwhm=%.2fHz, range=%.2f-%.2fHz', study_info.age, study_info.clusters{c}, pk_freq, fwhm, l_freq, r_freq));
        end
        
    end
    fois{c}=c_fois;
end
