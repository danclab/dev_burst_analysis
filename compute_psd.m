function compute_psd(study_info, varargin)
% COMPUTE_PSD - Compute power spectral density in all channels and remove
% aperiodic component
%
% Syntax:  compute_psd(study_info, varargin)
%
% Inputs:
%    study_info - structure containing the study information
%
% Example: 
%   compute_psd(study_info, varargin)

% Parse optional arguments
defaults=struct('min_ntrials',5);
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

% Open EEGlab
[ALLEEG, EEG, CURRENTSET] = eeglab;

% Power spectra for each subject in each channel
spectra=[];
% Periodic spectra for each subject in each channel
periodic=[];
% Subjects included
subjects={};

% Number of subjects
n_subjects=size(study_info.participant_info,1);

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
        [ALLEEG, EEG] = eeg_store(ALLEEG, base_EEG, 1);        
        exp_EEG=pop_loadset('filepath', subject_data_dir,...
            'filename', exp_fname);  
        [ALLEEG, EEG] = eeg_store(ALLEEG, exp_EEG, 2);

        % If min number of trials in baseline and experimental epochs
        if base_EEG.trials>=params.min_ntrials &&...
                exp_EEG.trials>=params.min_ntrials
            
            % Merge baseline and experimental EEG data
            merged_EEG=pop_mergeset(ALLEEG,1:2,0);

            % Use a window size (in datapoints) equal to the sampling rate
            % (so 1s) with a 50% overlap
            winsize=merged_EEG.srate;
            overlap=round(winsize/2);

            % Compute power spectral density for each channel using Welch's
            % method
            [subj_spectra,frex,~,~,~] = spectopo(merged_EEG.data,...
                merged_EEG.pnts, merged_EEG.srate, 'winsize', winsize,...
                'overlap', overlap, 'plot', 'off', 'freqfac',2);

            % Get frequencies from 0.5 to 40
            freq_idx=find((frex>=0.5) & (frex<=40));
            frex=frex(freq_idx);
            subj_spectra=subj_spectra(:,freq_idx);

            ch_resids=zeros(size(subj_spectra));
            for i=1:size(subj_spectra,1)
                % Channel spectrum
                chan_psd=subj_spectra(i,:);
                % 1/f
                oof=log10(1./frex);
                % Fit 1/f to spectrum
                lm_psd=fitlm(oof,chan_psd,'RobustOpts','on');
                % Get residuals
                ch_resids(i,:)=lm_psd.Residuals.Raw;                
            end
            
            spectra(end+1,:,:)=subj_spectra;    
            periodic(end+1,:,:)=ch_resids;
            subjects{end+1}=subj_id;
        end
    end
end

save(fullfile(study_info.deriv_dir,'psd.mat'),'frex','spectra',...
    'periodic','subjects');
