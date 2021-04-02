Developmental Transient Burst EEG Analysis
=======================

Analysis code for analyzing frequency-specific transient burst events in developmental EEG data.

> H Rayson, R Debnath, S Alavizadeh, N Fox, PF Ferrari, JJ Bonaiuto<br>
> **Detection and Analysis of Cortical Beta Bursts in Developmental EEG Data**<br>

## Requirements:

* EEGlab v14.1.1: https://sccn.ucsd.edu/eeglab/index.php
* FieldTrip v20190329: https://www.fieldtriptoolbox.org/
* shadedErrorBar: https://fr.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar?
* shadedErrorBaryy: https://uk.mathworks.com/matlabcentral/fileexchange/42578-shaded-error-bar-yy

## Usage

    % Initialize study information
    infant_study_info=init_umd12m_study_info();
    adult_study_info=init_umdadult_study_info();
    
    % Compute lagged coherence
    compute_lagged_coherence(infant_study_info);
    compute_lagged_coherence(adult_study_info);
    
    % Compute power spectral density
    compute_psd(infant_study_info);
    compute_psd(adult_study_info);
    
    % Plot mean PSD and lagged coherence
    plot_mean_psd_and_lagged_coherence(infant_study_info);
    plot_mean_psd_and_lagged_coherence(adult_study_info);
    
    % Plot single trial data for adult showing filter Hilbert method
    plot_single_trial_data(adult_study_info);
    
    % Plot amplitude distribution within identified beta band
    plot_amp_dist(infant_study_info, infant_study_info.beta_band);
    plot_amp_dist(adult_study_info, adult_study_info.beta_band);
    
    % Estimate the optimal relative threshold for beta
    estimate_threshold_sd(infant_study_info, infant_study_info.beta_band);
    estimate_threshold_sd(adult_study_info, adult_study_info.beta_band);
    
    % Plot beta bursts and beta amplitude using the optimal relative threshold
    plot_bursts_and_amplitude(infant_study_info, infant_study_info.beta_band, infant_study_info.beta_thresh_sd);
    plot_bursts_and_amplitude(adult_study_info, adult_study_info.beta_band, adult_study_info.beta_thresh_sd);
