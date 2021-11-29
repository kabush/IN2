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

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.betas.fmri_rest_beta]);
    eval(['! rm -rf ',proj.path.betas.fmri_rest_beta]);
    disp(['Creating ',proj.path.betas.fmri_rest_beta]);
    eval(['! mkdir ',proj.path.betas.fmri_rest_beta]);
end

%% ----------------------------------------
%% load all subjs
subjs = load_subjs(proj);

%% Initialize log section
logger(['************************************************'],proj.path.logfile);
logger(['Calculating fMRI beta-series (REST) of ',num2str(numel(subjs)),' subjects'],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Identify subjs with good models 
cnt = 0;
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    try
        
        % Test for data completeness
        rest_path = [proj.path.mri.mri_clean,subj_study,'_',name,'/rest/run1/'];
        censor_name = [subj_study,'.',name,'.rest.run1.censor.1D'];
        censor = load([rest_path,censor_name]);
        
        sum_censor = sum(censor);
        
        if (sum_censor>proj.param.rest.n_pseudo+29)
        
            cv_svm_path = [proj.path.mvpa.fmri_ex_gm_mdl];
            cv_v_svm_name = [subj_study,'_',name,'_v_model.mat'];
            cv_a_svm_name = [subj_study,'_',name,'_a_model.mat'];

            v_model = load([cv_svm_path,cv_v_svm_name]);
            a_model = load([cv_svm_path,cv_a_svm_name]);

            cnt = cnt + 1;
            good_subjs{cnt} = subjs{i};
            
        else
            
            logger(['*** Too censored to run REST ***: ',subj_study,'_',name],proj.path.logfile);
            
        end
        
    catch
        disp(['   ',subj_study,':',name,', Model error']);
    end

end
subjs = good_subjs; %%Subjects are now correctly sized

%% ----------------------------------------
%% iterate over study subjects
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    %% debug
    logger([subj_study,':',name],proj.path.logfile);

    %% Set-up data paths
    tmp_path = [proj.path.code,'tmp/'];
    rest_path = [proj.path.mri.mri_clean,subj_study,'_',name,'/rest/run1/'];

    %% Subject processed resting state data files
    rest_name = [subj_study,'.',name,'.rest.run1.scaled.resid'];
    censor_name = [subj_study,'.',name,'.rest.run1.censor.1D'];
    motion_name = [subj_study,'.',name,'.rest.run1.motion.1D'];
    motion_sqr_name = [subj_study,'.',name,'.rest.run1.motion.square.1D'];
    motion_pre_t_name = [subj_study,'.',name,'.rest.run1.motion_pre_t.1D'];
    motion_pre_t_sqr_name = [subj_study,'.',name,'.rest.run1.motion_pre_t_square.1D'];

    %% Pre clean-up tmp 
    eval(['! rm ',tmp_path,'*']);

    %% Copy rest data to tmp (& simplify name)
    eval(['! cp ',rest_path,rest_name,'+tlrc.BRIK ',tmp_path,'rest+tlrc.BRIK']);
    eval(['! cp ',rest_path,rest_name,'+tlrc.HEAD ',tmp_path,'rest+tlrc.HEAD']);
    eval(['! cp ',rest_path,censor_name,' ',tmp_path,'censor.1D']);
    eval(['! cp ',rest_path,motion_name,' ',tmp_path,'motion.1D']);
    eval(['! cp ',rest_path,motion_sqr_name,' ',tmp_path,'motion_square.1D']);
    eval(['! cp ',rest_path,motion_pre_t_name,' ',tmp_path,'motion_pre_t.1D']);
    eval(['! cp ',rest_path,motion_pre_t_sqr_name,' ',tmp_path,'motion_pre_t_square.1D']);

    %% Locally rename project params
    N_trs = proj.param.mri.n_trs_rest;
    N_sample = proj.param.rest.n_pseudo; %stimulus sample
    N_iter = proj.param.rest.n_resample;  %resampling of stimuli
    N_trans = proj.param.rest.n_trs_trans; %start volumes trimmed
    N_tail = proj.param.rest.n_trs_tail; %end volumes trimmed

    %% Load mask
    gm_nii = load_nii([proj.path.mri.gm_mask,'group_gm_mask.nii']);
    Nvox = prod(size(gm_nii.img));

    %% Initialize storage
    mu_img = zeros(Nvox,proj.param.mri.n_trs_rest);
    mu_cnt = zeros(Nvox,proj.param.mri.n_trs_rest);

    %% Estimate beta and affect multiple times
    for j=1:N_iter
        
        %% Sample a stimulus timing set (50% of total TRs)
        stim_ids = randsample((N_trans+1):(N_trs-N_tail),N_sample);
        stim_times = stim_ids * proj.param.mri.TR;
        
        %% Write stimulutus timings to file
        fid = fopen([tmp_path,'stim_times.1D'],'w');
        fprintf(fid,'%5.2f\n',stim_times);
        fclose(fid);
      
        %% Execute beta-series regression
        eval(['! ./source/mvpa/mvpa_3dlss']);       

        %% Load beta-series
        path = [tmp_path,'rest_lss.nii'];
        base_nii = load_nii(path);        
        brain_size = size(base_nii.img);
        
        %% Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);        
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4));

        mu_img(:,stim_ids) = mu_img(:,stim_ids)+base_img;
        mu_cnt(:,stim_ids) = mu_cnt(:,stim_ids)+1;
        
    end %(of iterations, index: j)

    %% Compute average beta series at rest
    mdl_seq = (N_trans+1):(N_trs-N_tail);
    mu_rest = 0*mu_img;
    mu_rest(:,mdl_seq) = mu_img(:,mdl_seq)./mu_cnt(:,mdl_seq); % Point-wise avg of values

    %% Build nifti format of beta series from mask
    gm_nii = load_nii([proj.path.mri.gm_mask,'group_gm_mask.nii']);
    mu_rest_nii = build_beta_nii_from_gm_mask(mu_rest,gm_nii,1:Nvox);
    
    %% save nifti (match naming format for ex & in variants)
    save_nii(mu_rest_nii,[proj.path.betas.fmri_rest_beta,subj_study,'_',name,'_lss.nii']);

end % (of subjs, index: i)