study_infos={};
study_infos{1}=init_umd12m_study_info();
study_infos{2}=init_umdadult_study_info();

for st=1:length(study_infos)
    study_info=study_infos{st};
    
    compute_psd(study_info);
    plot_mean_psd(study_info);
    compute_lagged_coherence(study_info);
    
    c_fois=determine_fois(study_info);
    alphas=[];
    betas=[];
    for c=1:length(c_fois)
        fois=c_fois{c};
        sorted_fois=sortrows(fois,1);
        alphas(end+1,:)=sorted_fois(1,:);
        betas(end+1,:)=sorted_fois(2,:);
    end
    study_info.alpha_band=alphas;
    study_info.beta_band=betas;
    
    study_info.beta_thresh_sd=estimate_threshold_sd(study_info, study_info.beta_band);
    
    export_bursts(study_info, study_info.beta_band, study_info.beta_thresh_sd);  
    
    plot_ind_mean_psd(study_info);
    plot_ind_lagged_coherence(study_info);
end

plot_bursts_and_amplitude();
plot_mean_psd_and_lagged_coherence();
plot_mean_psd_and_lagged_coherence_topo();
