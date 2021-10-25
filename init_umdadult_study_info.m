%% Individualized_Info
function study_info=init_umdadult_study_info()

study_info=[];
study_info.age='adult';

%% Specify Paths
% Locate data on hard drive
study_info.data_dir = '../../data/adult';
study_info.deriv_dir = '../../data/UMD/adult/deriv';

%% Subject list
study_info.participant_info=readtable([fullfile(study_info.data_dir, 'participants.tsv')],...
    'FileType','text','Delimiter','\t','TreatAsEmpty',{''});

% Electrode clusters to look in
study_info.clusters={'C3','C4'};
study_info.cluster_channels={
    {'E16', 'E20', 'E21', 'E22'},...
    {'E41', 'E49', 'E50', 'E51'}};

% Conditions and events
study_info.conditions={'observe','execute'};
study_info.baseline_evts={'OBBASE_P','EXBASE_P'};
study_info.exp_evts={'FTGO_P','FTGE_P'};

study_info.beta_band=[19.5 25.5; 18.75 24.25];
study_info.beta_thresh_sd=[1.1 1.2];
study_info.alpha_band=[9.5 12.5; 9.5 12.5];
