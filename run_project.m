%% ------------------------------------------------------------
%% Clean up matlab environment
matlab_reset;

tic

%% ------------------------------------------------------------
%% Link all source code
addpath(genpath('./source/'));

%% ------------------------------------------------------------
%%  Initialize the projects directories and parameters.
init_project;

%% ----------------------------------------
%% Clear and reconstruct the project data folder
clean_project;

%% ------------------------------------------------------------
%% Apply Affect Models to Resting State Data
% mvpa_fmri_rest_gm_cls;

%% ------------------------------------------------------------ 
%% Calculate (RST) Stimuli Beta-Series
% calc_fmri_rest_beta;

%% ------------------------------------------------------------ 
%% Compute Dynamics
% dynamics_fmri_rest_gm;

%% ------------------------------------------------------------ 
%% MVPA Regression of 2nd/1st Deriv (RST)
% mvpa_fmri_rest_2drv_rgr;  
% mvpa_fmri_rest_1drv_rgr;  

%% ------------------------------------------------------------ 
%% Analysis 2nd and 1st Deriv Predictions (RST)
% analyze_fmri_rest_mvpa_2drv;
% analyze_fmri_rest_mvpa_1drv;

%% ------------------------------------------------------------ 
%% Haufe Encoding of 1st & 2nd Deriv Hyperplanes
% haufe_fmri_rest_2drv_rgr_perm_v;
% haufe_fmri_rest_2drv_rgr_perm_a;
% haufe_fmri_rest_1drv_rgr_perm_v;
% haufe_fmri_rest_1drv_rgr_perm_a;

%% ------------------------------------------------------------ 
%% Affect Simulations
% simulate_fmri_rest_gm_v;
% simulate_fmri_rest_gm_a;

toc
