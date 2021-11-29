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
logger(['Predicting Resting State Affect (V & A)         '],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.mvpa.fmri_rest_gm_cls]);
    eval(['! rm -rf ',proj.path.mvpa.fmri_rest_gm_cls]);
    disp(['Creating ',proj.path.mvpa.fmri_rest_gm_cls]);
    eval(['! mkdir ',proj.path.mvpa.fmri_rest_gm_cls]);
end

%% ----------------------------------------
%% load all subjs
subjs = load_subjs(proj);

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

    %% Iteration based data from regression
    id_rest = zeros(N_trs,N_iter);

    %% Will become avg values (used to keep running sum also)
    cv_mu_v_rest = zeros(N_trs,1);
    cv_mu_a_rest = zeros(N_trs,1);
    slf_mu_v_rest = zeros(N_trs,1);
    slf_mu_a_rest = zeros(N_trs,1);

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
        
        %% ============================================================
        %% For every loocv subject's model estimate the valence
        %% ============================================================
        
        %% store predicted affect (CV models)
        itr_cv_v = zeros(numel(stim_ids),numel(subjs)-1);
        itr_cv_a = zeros(numel(stim_ids),numel(subjs)-1);

        %% store predicted affect (self models)
        itr_slf_v = zeros(numel(stim_ids),1);
        itr_slf_a = zeros(numel(stim_ids),1);

        %% loop over all CV subjects
        cnt = 0;
        for k=1:numel(subjs)
            
            %% extract subject info
            cv_subj_study = subjs{k}.study;
            cv_subj_id = subjs{k}.id;
            cv_name = subjs{k}.name;
            
            cv_svm_path = [proj.path.mvpa.fmri_ex_gm_mdl];
            cv_v_svm_name = [cv_subj_study,'_',cv_name,'_v_model.mat'];
            cv_a_svm_name = [cv_subj_study,'_',cv_name,'_a_model.mat'];

            cv_mask_path = [proj.path.mri.gm_mask];
            cv_mask_name =  [cv_subj_study,'.',cv_name,'.gm.nii'];

            %% Copy cv subj SVM to tmp (& simplify name)
            eval(['! cp ',cv_svm_path,cv_v_svm_name,' ',tmp_path,'cv_v_model.mat']);
            eval(['! cp ',cv_svm_path,cv_a_svm_name,' ',tmp_path,'cv_a_model.mat']);
                
            %% Copy subj GM mask to tmp (& simplify name)
            eval(['! cp ',cv_mask_path,cv_mask_name,' ',tmp_path,'gm.nii']);
            
            %% Load gray matter mask 
            gm_nii = load_nii([tmp_path,'gm.nii']);
            mask = double(gm_nii.img);
            brain_size=size(mask);
            mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
            in_brain=find(mask==1);  
            
            %% Concatenate the MASKED base image
            all_img = base_img(in_brain,:)';
            
            %% Predict valence|arousal quantities REST
            load([tmp_path,'cv_v_model.mat']);
            load([tmp_path,'cv_a_model.mat']);
            
            [cls hd] = predict(v_model,all_img);
            cv_v_tilde = hd(:,2);
            
            [cls hd] = predict(a_model,all_img);
            cv_a_tilde = hd(:,2);
            
            %% Inter-subject predictions
            if((strcmp(cv_subj_study,subj_study) ~=0 & strcmp(cv_name,name) ...
                == 0) | strcmp(cv_subj_study,subj_study) == 0)
                
                disp(['   ',cv_subj_study,':',subj_study,':',cv_name,':',name,':',num2str(j)]);
                
                cnt = cnt+1;
                
                %% Store predictions
                itr_cv_v(:,cnt) = cv_v_tilde;
                itr_cv_a(:,cnt) = cv_a_tilde;
                
                % %% ----------------------------------------
                % %% *** DEBUG: View Incremental Predictions ***
                % mu = mean(itr_cv_v(:,1:cnt),2);
                % sig=0;
                % if(cnt>3)
                %     sig = std(itr_cv_v(:,1:cnt)');
                % end
                % 
                % plot(mu,'r-','LineWidth',2);
                % hold on;
                % plot(mu+1.96*(sig/sqrt(cnt))','b-');
                % plot(mu-1.96*(sig/sqrt(cnt))','b-');
                % hold off;
                % drawnow;
                % %% ----------------------------------------

                
                %% Self predictions
            else
                
                disp(['   ',cv_subj_study,':',subj_study,':',cv_name,':',name,':',num2str(j),':self']);
                
                %% Store predictions
                itr_slf_v = cv_v_tilde;
                itr_slf_a = cv_a_tilde;
                
            end
            
        end % (of CV subjects, index: k)
        
        %% Combined CV predicted affect for these pseudo stims
        mu_itr_cv_v = mean(itr_cv_v,2);
        mu_itr_cv_a = mean(itr_cv_a,2);
        
        %% store the ids involved in these pseudo stims
        id_rest(stim_ids,j) = ones(numel(stim_ids),1);
        
        %% build total affect scores (via running sums)
        cv_mu_v_rest(stim_ids) = cv_mu_v_rest(stim_ids)+mu_itr_cv_v;
        cv_mu_a_rest(stim_ids) = cv_mu_a_rest(stim_ids)+mu_itr_cv_a;
        slf_mu_v_rest(stim_ids) = slf_mu_v_rest(stim_ids)+itr_slf_v;
        slf_mu_a_rest(stim_ids) = slf_mu_a_rest(stim_ids)+itr_slf_a;
        
        %% delete computed beta series
        eval(['! rm ',tmp_path,'rest_lss.nii']);
        
    end %(of iterations, index: j)
    
    %% Use id_count to compute averages scores|states
    id_count = sum(id_rest,2);
    good_ids = find(id_count>0);
    for j=1:numel(good_ids)
        id = good_ids(j);
        cnt = id_count(id);
        cv_mu_v_rest(id) = cv_mu_v_rest(id)/cnt;
        cv_mu_a_rest(id) = cv_mu_a_rest(id)/cnt;
        slf_mu_v_rest(id) = slf_mu_v_rest(id)/cnt;
        slf_mu_a_rest(id) = slf_mu_a_rest(id)/cnt;
    end
    
    %% Save out data
    save([proj.path.mvpa.fmri_rest_gm_cls,subj_study,'_',name,'_cv_mu_v_rest.mat'],'cv_mu_v_rest');
    save([proj.path.mvpa.fmri_rest_gm_cls,subj_study,'_',name,'_cv_mu_a_rest.mat'],'cv_mu_a_rest');
    save([proj.path.mvpa.fmri_rest_gm_cls,subj_study,'_',name,'_slf_mu_v_rest.mat'],'slf_mu_v_rest');
    save([proj.path.mvpa.fmri_rest_gm_cls,subj_study,'_',name,'_slf_mu_a_rest.mat'],'slf_mu_a_rest');
    
    %% Post clean-up tmp 
    eval(['! rm ',tmp_path,'*']);
    
end % (of subjs, index: i)