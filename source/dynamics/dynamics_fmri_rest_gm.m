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
logger(['*************************************************'],proj.path.logfile);
logger(['Computing REST affect dynamics'],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.dyn.rest]);
    eval(['! rm -rf ',proj.path.dyn.rest]);
    disp(['Creating ',proj.path.dyn.rest]);
    eval(['! mkdir ',proj.path.dyn.rest]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% Transform beta-series into affect series {v,a}
for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log processing of subject
    logger([subj_study,'_',name],proj.path.logfile);

    data_exist = 0;
    try
        
        %% Load gray matter mask 
        gm_nii = load_nii([proj.path.mri.gm_mask,subj_study,'.',name,'.gm.nii']);
        mask = double(gm_nii.img);
        brain_size=size(mask);
        mask = reshape(mask,brain_size(1)*brain_size(2)*brain_size(3),1);
        in_brain=find(mask==1);  
        
        %% Load beta-series
        path = [proj.path.betas.fmri_rest_beta,subj_study,'_',name,'_lss.nii'];
        base_nii = load_nii(path);
        brain_size = size(base_nii.img);

        %% Data is present
        data_exist = 1;

        %% Vectorize the base image
        base_img = vec_img_2d_nii(base_nii);
        base_img = reshape(base_img,brain_size(1)*brain_size(2)*brain_size(3),brain_size(4));
        
        %% Concatenate the MASKED base image
        subj_img = base_img(in_brain,:)';
        
        %% Perform quality
%         qlty = check_gm_img_qlty(subj_img);
        qlty.ok = 1;
        
    catch
        logger(['   -mask or beta-series does not exist'],proj.path.logfile);
    end
    
    if(qlty.ok & data_exist)
    
        %%Extract only usable betas for modeling
        rest_id = 1:proj.param.mri.n_trs_rest;
        % rest_id = proj.param.rest.n_trs_trans:(proj.param.mri.n_trs_rest-proj.param.rest.n_trs_tail);
        subj_img = subj_img(rest_id,:);

        %% Initialize the prediction structure of this subject
        prds = struct();
        prds.v_hd = zeros(numel(rest_id),1);
        prds.a_hd = zeros(numel(rest_id),1);
        
        mdl_exist = 0;
        try
            %% Load SVM models
            load([proj.path.mvpa.fmri_ex_gm_mdl,subj_study,'_',name,'_v_model.mat']);
            load([proj.path.mvpa.fmri_ex_gm_mdl,subj_study,'_',name,'_a_model.mat']);
            mdl_exist = 1;
        catch
            logger(['   -model does not exist.'],proj.path.logfile);
        end
            
        if(mdl_exist==1)

            %% ----------------------------------------
            %% predict REST task using EX-based models
            for j=1:numel(rest_id)
                
                %% valence
                [tst_predict,hd] = predict(v_model,subj_img(j,:));
                prds.v_hd(j) = hd(2);
                
                %% arousal
                [tst_predict,hd] = predict(a_model,subj_img(j,:));
                prds.a_hd(j) = hd(2);
                
            end
            
            %% ----------------------------------------
            %% Count NaN Volumes
            prds.v_pre_nan_cnt = sum(isnan(prds.v_hd));
            prds.a_pre_nan_cnt = sum(isnan(prds.a_hd));
            
            logger([num2str(prds.v_pre_nan_cnt), ' pre-outlier NaN volumes'],proj.path.logfile);
            
            %% ----------------------------------------
            %% Outlier detection based on (Cox, 2002), replaces with NaN
            % valence
            v_median = nanmedian(prds.v_hd);
            v_mad = nanmedian(abs(prds.v_hd-v_median));
            v_categorical = (1-norminv(0.05/length(prds.v_hd)))*((pi/2)^(1/2));
            v_low_bound = v_median-v_categorical*v_mad;
            v_upper_bound = v_median+v_categorical*v_mad;

            v_outlier_mask = zeros(length(prds.v_hd),1);
            for z=1:length(prds.v_hd)
                if (prds.v_hd(z)<v_low_bound) || (prds.v_hd(z)>v_upper_bound)
                   v_outlier_mask(z) = 1; 
                end
            end
            v_outlier_ids = find(v_outlier_mask);
%             v_outlierness = -log10((1-normcdf(abs(prds.v_hd-v_median)./(v_mad*((pi/2)^(1/2))))));
            prds.v_hd(v_outlier_ids) = NaN;

            % arousal
            a_median = nanmedian(prds.a_hd);
            a_mad = nanmedian(abs(prds.a_hd-a_median));
            a_categorical = (1-norminv(0.05/length(prds.a_hd)))*((pi/2)^(1/2));
            a_low_bound = a_median-a_categorical*a_mad;
            a_upper_bound = a_median+a_categorical*a_mad;

            a_outlier_mask = zeros(length(prds.a_hd),1);
            for z=1:length(prds.a_hd)
                if (prds.a_hd(z)<a_low_bound) || (prds.a_hd(z)>a_upper_bound)
                   a_outlier_mask(z) = 1; 
                end
            end
            a_outlier_ids = find(a_outlier_mask);
%             a_outlierness = -log10((1-normcdf(abs(prds.a_hd-a_median)./(a_mad*((pi/2)^(1/2))))));
            prds.a_hd(a_outlier_ids) = NaN;
            
            %% Count num of removed outlier volumes
            prds.v_out_nan_cnt = numel(v_outlier_ids);
            prds.a_out_nan_cnt = numel(a_outlier_ids);
            
            logger([num2str(prds.v_out_nan_cnt), ' valence outlier NaN volumes'],proj.path.logfile);
            logger([num2str(prds.a_out_nan_cnt), ' arousal outlier NaN volumes'],proj.path.logfile);
            
            %% ----------------------------------------
            %% decompose predicted trajectories (& derivs)
            
            %% valence
            [prds.v_dcmp,prds.v_indx] = decompose_rest(proj,prds.v_hd');
            
            %% arousal
            [prds.a_dcmp,prds.a_indx] = decompose_rest(proj,prds.a_hd');
            
            logger('   -success',proj.path.logfile);
            
            % debug
            figure(99)
            plot(prds.v_dcmp.h(1,3:end-2));
            hold on;
            plot(prds.v_dcmp.dh(1,2:end-1));
            plot(prds.v_dcmp.d2h(1,1:end));
            hold off;              
            drawnow
            
            %% Save out prediction structure
            save([proj.path.dyn.rest,subj_study,'_',name,'_prds.mat'],'prds');

        end
            
    else
        logger('   -failed quality check',proj.path.logfile);
    end
    
end