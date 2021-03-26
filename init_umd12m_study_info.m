%% Individualized_Info
function study_info=init_umd12m_study_info()

study_info=[];

%% Specify Paths
% Locate data on hard drive
study_info.data_dir = '../../data/UMD/12m/bids';
study_info.deriv_dir = '../../data/UMD/12m/deriv';

%% Subject list
study_info.participant_info=readtable([fullfile(study_info.data_dir, 'participants.tsv')],...
    'FileType','text','Delimiter','\t','TreatAsEmpty',{''});

% Electrode clusters to look in
study_info.clusters={'C3','C4'};
study_info.cluster_channels={{'E16', 'E20', 'E21', 'E22'},...
    {'E41', 'E49', 'E50', 'E51'}};

% Conditions and events
study_info.conditions={'observe','execute'};
study_info.baseline_evts={'OBBASE_P','EXBASE_P'};
study_info.exp_evts={'FTGO_P','FTGE_P'};

study_info.beta_band=[12 18];
study_info.beta_thresh_sd=1.25;
study_info.mu_band=[6 10];

% Open EEGlab
%[ALLEEG, EEG, CURRENTSET] = eeglab;

% Init fieldtrip
%ft_defaults;