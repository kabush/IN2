%%========================================
%%========================================
%%
%% Keith Bush, PhD (2021)
%% Univ. of Arkansas for Medical Sciences
%% Brain Imaging Research Center (BIRC)
%%
%%========================================
%%========================================

%% Load in path data
load('proj.mat');

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger(['Intra-subject LOOCV MVPA of 1st Deriv  '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mvpa.fmri_rest_1drv_rgr]);
    eval(['! rm -rf ',proj.path.mvpa.fmri_rest_1drv_rgr]);
    disp(['Creating ',proj.path.mvpa.fmri_rest_1drv_rgr]);
    eval(['! mkdir ',proj.path.mvpa.fmri_rest_1drv_rgr]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

rho_v_all = [];
rho_a_all = [];

%% ----------------------------------------
%% iterate over study subjects
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    try

        %% Load 1st deriv targets
        load([proj.path.dyn.rest,subj_study,'_',name,'_prds.mat']);
        
        %% Reshape targets & indices
        
        dh_indx_v = reshape(prds.v_indx.dh',1,prod(size(prds.v_indx.dh)));
        dh_trgs_v = reshape(prds.v_dcmp.dh',1,prod(size(prds.v_dcmp.dh)));
        dh_indx_a = reshape(prds.a_indx.dh',1,prod(size(prds.a_indx.dh)));
        dh_trgs_a = reshape(prds.a_dcmp.dh',1,prod(size(prds.a_dcmp.dh)));
        
        %% Clean-up indices that were zeroed during beta-series
        %% 1 eroded during deriv calcs (erode 4 more on front and 9
        %% more on tail)
        dh_indx_v = dh_indx_v(5:(end-9));
        dh_trgs_v = dh_trgs_v(5:(end-9));
        dh_indx_a = dh_indx_a(5:(end-9));
        dh_trgs_a = dh_trgs_a(5:(end-9));
        
        %% Only use data without NaN
        v_nan_mask = isnan(dh_trgs_v);
        dh_indx_v(v_nan_mask) = [];
        dh_trgs_v(v_nan_mask) = [];
        
        a_nan_mask = isnan(dh_trgs_a);
        dh_indx_a(a_nan_mask) = [];
        dh_trgs_a(a_nan_mask) = [];
        
        %% Save out NaN mask for later use in Simulation step
        
        save([proj.path.mvpa.fmri_rest_1drv_rgr,subj_study,'_',name,'_v_nan_mask.mat'],'v_nan_mask');
        save([proj.path.mvpa.fmri_rest_1drv_rgr,subj_study,'_',name,'_a_nan_mask.mat'],'a_nan_mask');

        %% Load gray matter mask 
        gm_nii = load_nii([proj.path.mri.gm_mask,'group_gm_mask.nii']);
        mask = double(gm_nii.img);
        brain_size=size(mask);
        mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
        in_brain=find(mask==1);  
        
        %% Load beta-series
        base_nii = load_nii([proj.path.betas.fmri_rest_beta,subj_study,'_',name,'_lss.nii']);
        brain_size = size(base_nii.img);
        
        %% Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4));
        
        %% Concatenate the MASKED base image
        all_img = base_img(in_brain,:)';

        %% Subselect data for which 1st derivs are calculated
        dh_img_v = all_img(dh_indx_v,:);
        dh_img_a = all_img(dh_indx_a,:);
        
        %% Peform quality check of generated features
%         qlty = check_gm_img_qlty(dh_img);
        qlty.ok = 1;

        if(qlty.ok)

            % organize predictions
            prds = struct();
            
            tic
            % Model valence
            [out,trg,mdl,stats] = regress_intra_loocv(dh_img_v,dh_trgs_v,proj.param.mvpa.kernel);            
            prds.v.out = out;
            prds.v.trg = trg;
            prds.v.mdl = mdl;
            prds.v.stats = stats;

            disp([' v rho=',num2str(stats.rho)]);
            rho_v_all = [rho_v_all,stats.rho];

            % Model arousal
            [out,trg,mdl,stats] = regress_intra_loocv(dh_img_a,dh_trgs_a,proj.param.mvpa.kernel);            
            prds.a.out = out;
            prds.a.trg = trg;
            prds.a.mdl = mdl;
            prds.a.stats = stats;
            
            disp([' a rho=',num2str(stats.rho)]);
            rho_a_all = [rho_a_all,stats.rho];

            toc
            
            % save predictions
            save([proj.path.mvpa.fmri_rest_1drv_rgr,subj_study,'_',name,'_prds.mat'],'prds');

        end
        
    catch
        logger(['  -MVPA Error: possible missing beta series'],proj.path.logfile);
    end

end