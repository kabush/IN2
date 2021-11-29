%% ========================================
%% Load in path data
load('proj.mat');

if(proj.flag.clean_build)

    %% Create project directory
    disp(['Base directories']);

    %% Create project directory
    disp(['Creating fresh sub-directories']);

    %% Create all top-level directories
%     eval(['! rm -rf ',proj.path.data,proj.path.analysis.name]);
%     eval(['! mkdir ',proj.path.data,proj.path.analysis.name]);
    
%     eval(['! rm -rf ',proj.path.data,proj.path.mvpa.name]);
%     eval(['! mkdir ',proj.path.data,proj.path.mvpa.name]);

%     eval(['! rm -rf ',proj.path.data,proj.path.haufe.name]);
%     eval(['! mkdir ',proj.path.data,proj.path.haufe.name]);
    
    disp(['Clearing tmp']);
    eval(['! rm ',proj.path.code,'tmp/*']);

end