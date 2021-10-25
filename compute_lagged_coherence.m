function compute_lagged_coherence(study_info, varargin)
% COMPUTE_LAGGED_COHERENCE - Computes lagged coherence in all channels
% across a range of frequencies and lags
%
% Syntax:  compute_lagged_coherence(study_info)
%
% Inputs:
%    study_info - structure containing the study information
%
% Example: 
%   compute_lagged_coherence(study_info)

% Parse optional arguments
defaults=struct('min_ntrials',5);
params=struct(varargin{:});
for f=fieldnames(defaults)',
    if ~isfield(params, f{1})
        params.(f{1})=defaults.(f{1});
    end
end

% Frequencies to run on. We start at 5Hz because the trials aren't long
% enough to look at lower frequencies over long lags
foi=[5:1:40];

% Lags to run on
lags=[2:.1:7];

% Number of subjects
n_subjects=size(study_info.participant_info,1);

% Open EEGlab
[ALLEEG, EEG, CURRENTSET] = eeglab;

% Lagged coherence for each subject in each channel (64) over all
% frequencies and lags
lagged_coh=zeros(n_subjects,64,length(foi),length(lags)).*NaN;

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

            % Create fieldtrip data structure
            data=create_ft_data(merged_EEG);
            
            % Configuration for frequency analysis
            cfg_F             = [];
            cfg_F.method      = 'mtmconvol';
            cfg_F.taper       = 'hanning';
            cfg_F.output      = 'fourier';
            cfg_F.keeptrials  = 'yes';
            cfg_F.pad         = 'nextpow2';
            
            % Configuration for lagged coherence
            cfg_LC            = [];            
            cfg_LC.method     = 'laggedcoherence';
            cfg_LC.trialsets  = 'all';
            fs                = data.fsample;
            
            for f_idx = 1:length(foi);
                
                % Set freq range
                cfg_F.foi     = foi(f_idx);
                cfg_LC.foi    = foi(f_idx);
                
                for l_idx=1:length(lags)
                    lag=lags(l_idx);
                    cfg_LC.lag        = lag;
                    cfg_F.width       = lag;
                   
                    % width of time windows for phase comparison (seconds)
                    width         = cfg_F.width/cfg_F.foi;                        
                    cfg_F.t_ftimwin = width;
                    
                    % half width of time window (seconds)
                    halfwidth     = ceil(fs*width/2)/fs;

                    % Go from half window width after trial start to half
                    % window width before trial end
                    toi_start     = data.time{1}(1) + halfwidth;
                    toi_stop      = data.time{1}(end) - halfwidth;
                    
                    % Step size
                    step          = ceil(fs*cfg_LC.lag/cfg_F.foi)/fs;
                    cfg_F.toi     = toi_start:step:toi_stop;
                    
                    % Run frequency analysis
                    freqout       = ft_freqanalysis(cfg_F,data);
                    
                    % Compute lagged coherence
                    lcoh = ft_connectivity_laggedcoherence(cfg_LC,freqout);
                    lagged_coh(s,:,f_idx,l_idx)=lcoh.laggedcoh;            
                end
            end    
        end
    end
    save(fullfile(study_info.deriv_dir,'lagged_coherence.mat'),'lags',...
        'foi','lagged_coh');
end

