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
logger([' Global Perm. Testing of MVPA Hyperplanes 1drv a'],proj.path.logfile);
logger(['************************************************'],proj.path.logfile);

%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.haufe.fmri_rest_1drv_rgr_a]);
    eval(['! rm -rf ',proj.path.haufe.fmri_rest_1drv_rgr_a]);
    disp(['Creating ',proj.path.haufe.fmri_rest_1drv_rgr_a]);
    eval(['! mkdir ',proj.path.haufe.fmri_rest_1drv_rgr_a]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% iterate over permuations
Nperm = proj.param.haufe.npermute;
Nloop = Nperm + 1; %% first loop is true model structure
Nchunk = proj.param.haufe.chunk;

%% storage for group Haufe 
grp_haufe_a = zeros(172800,Nloop);

%% permutation significance levels
alpha05 = 0.05;
alpha01 = 0.01;
alpha001 = 0.001;

for i = 1:Nloop

    tic

    %%storage for group haufe 
    all_haufe_a_wts = zeros(172800,numel(subjs));
    all_haufe_a_mask = zeros(172800,numel(subjs));

    qlty_vec = zeros(1,numel(subjs));

    %% ----------------------------------------
    %% iterate over study subjects
    for j = 1:numel(subjs)
        
        %% extract subject info
        subj_study = subjs{j}.study;
        name = subjs{j}.name;
        id = subjs{j}.id;
        
        %% debug
        logger([subj_study,':',name,':i=',num2str(i)],proj.path.logfile);
%         disp([subj_study,':',name,':i=',num2str(i)]);
        
        try

            %% Load 1st deriv targets
            load([proj.path.dyn.rest,subj_study,'_',name,'_prds.mat']);
            
            %% Reshape targets & indices
            dh_indx = reshape(prds.a_indx.dh',1,prod(size(prds.a_indx.dh))); 
            dh_trgs_a = reshape(prds.a_dcmp.dh',1,prod(size(prds.a_dcmp.dh)));

            %% Only use data without NaN
            a_nan_mask = isnan(dh_trgs_a);
            dh_indx(a_nan_mask) = [];
            dh_trgs_a(a_nan_mask) = [];
            
            %% Clean up indices that were zeroed during beta-series
            %% 2 eroded during deriv calcs (erode 3 more on front and 8
            %% more on tail)
            dh_indx = dh_indx(5:(end-9));
            dh_trgs_a = dh_trgs_a(5:(end-9));
            
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
            dh_img = all_img(dh_indx,:);
            
            %% Peform quality check of generated features
%             qlty = check_gm_img_qlty(dh_img);
            qlty.ok = 1;

            if(qlty.ok)
                
                %% ----------------------------------------
                %% AROUSAL 1st Deriv
                label_ids = 1:numel(dh_trgs_a);

                %% Only first iteration is structure (remaining
                %% loops are permutations, therefore randomize labels
                if(i>1)
                    label_ids = randsample(label_ids,numel(label_ids));
                end

                %% Fit classifier
                mdl = fitrsvm(dh_img,dh_trgs_a(label_ids),'KernelFunction',proj.param.mvpa.kernel);
                
                %% Construct Arousal Haufe tranform
                wts = mdl.Beta;
                haufe_a_wts = zscore(fast_haufe(dh_img,wts,Nchunk));
                all_haufe_a_wts(in_brain,j) = haufe_a_wts;
                all_haufe_a_mask(in_brain,j) = 1;
                
                qlty_vec(j)=1;
                
            end
            
        catch
            logger(['  -Haufe Error'],proj.path.logfile);
        end

    end

    %% ----------------------------------------
    %% Extract Quality fits        
    qlty_ids = find(qlty_vec==1);
    qlty_n = numel(qlty_ids);
    
    qlty_haufe_a_wts = all_haufe_a_wts(:,qlty_ids);
    qlty_haufe_a_mask = all_haufe_a_mask(:,qlty_ids);
    
    %% ----------------------------------------
    %% Group Mean Haufe transforms (1 saved per permutation)
    
    %% Group arousal Haufe
    ahf_sum = sum(qlty_haufe_a_mask,2);
    row_ids_a = find(ahf_sum>(qlty_n/2));

    % Distribution version
    grp_haufe_a_dist = zeros(size(grp_haufe_a,1),qlty_n);
    grp_haufe_a_dist(row_ids_a,:) = qlty_haufe_a_wts(row_ids_a,:);

    % Mean version
    grp_haufe_a(row_ids_a,i) = mean(qlty_haufe_a_wts(row_ids_a,:),2);

    % T-score version
    if(i==1)
        grp_haufe_a_tstat = 0*grp_haufe_a(:,1);
        for k=1:numel(row_ids_a)
            [h p ci stat] = ttest(qlty_haufe_a_wts(row_ids_a(k),:));
            grp_haufe_a_tstat(row_ids_a(k),1)=stat.tstat;
        end
    end

    %% ----------------------------------------
    %% Do permutation test given samples available

    if(i>1)

        save([proj.path.haufe.fmri_rest_1drv_rgr_a,'grp_haufe_a_n=',num2str(i-1),'_of_N=',num2str(Nperm),'.mat'],'grp_haufe_a');

        if(i>2)
            eval(['! rm ',proj.path.haufe.fmri_rest_1drv_rgr_a,'grp_haufe_a_n=',num2str(i-2),'_of_N=',num2str(Nperm),'.mat']);
        end
    
        %% ----------------------------------------
        %% Arousal
        sig_ids_05_a = [];
        sig_ids_01_a = [];
        sig_ids_001_a = [];
        
        for k=1:numel(row_ids_a)
            
            % Count extrem random samples
            Next = 0;
            if(grp_haufe_a(row_ids_a(k),1)>0)
                Next = numel(find(grp_haufe_a(row_ids_a(k),2:i)>grp_haufe_a(row_ids_a(k),1)));
            else
                Next = numel(find(grp_haufe_a(row_ids_a(k),2:i)<grp_haufe_a(row_ids_a(k),1)));
            end
            
            % Do 2-sided tests
            if(Next<round((alpha05/2)*i))
                sig_ids_05_a = [sig_ids_05_a,row_ids_a(k)];
            end
            
            if(Next<round((alpha01/2)*i))
                sig_ids_01_a = [sig_ids_01_a,row_ids_a(k)];
            end
            
            if(Next<round((alpha001/2)*i))
                sig_ids_001_a = [sig_ids_001_a,row_ids_a(k)];
            end
            
        end

        % ----------------------------------------
        % Save out: mean encoding of group gray-matter voxels
        if(numel(row_ids_a)>0)

            mu_a_haufe_nii = build_nii_from_gm_mask(grp_haufe_a_tstat(row_ids_a,1),gm_nii,row_ids_a);
            save_nii(mu_a_haufe_nii,[proj.path.haufe.fmri_rest_1drv_rgr_a,'mu_haufe_a_N=',num2str(Nperm),'.nii']);

        end

        % ----------------------------------------
        % Save out: mean encoding of permstrap sign. (p<0.05) group
        % gray-matter voxels
        if(numel(sig_ids_05_a)>0)
            mu_perm_a_haufe_nii = build_nii_from_gm_mask(grp_haufe_a_tstat(sig_ids_05_a,1),gm_nii,sig_ids_05_a);
            save_nii(mu_perm_a_haufe_nii,[proj.path.haufe.fmri_rest_1drv_rgr_a,'mu_perm_haufe_a_N=',num2str(Nperm),'_05.nii']);
        end

        % ----------------------------------------
        % Save out: mean encoding of permstrap sign. (p<0.01) group
        % gray-matter voxels
        if(numel(sig_ids_01_a)>0)
            mu_perm_a_haufe_nii = build_nii_from_gm_mask(grp_haufe_a(sig_ids_01_a,1),gm_nii,sig_ids_01_a);
            save_nii(mu_perm_a_haufe_nii,[proj.path.haufe.fmri_rest_1drv_rgr_a,'mu_perm_haufe_a_N=',num2str(Nperm),'_01.nii']);
        end

        % ----------------------------------------
        % Save out: mean encoding of permstrap sign. (p<0.001) group gray-matter voxels
        if(numel(sig_ids_001_a)>0)
             mu_perm_a_haufe_nii = build_nii_from_gm_mask(grp_haufe_a(sig_ids_001_a,1),gm_nii,sig_ids_001_a);
             save_nii(mu_perm_a_haufe_nii,[proj.path.haufe.fmri_rest_1drv_rgr_a,'mu_perm_haufe_a_N=',num2str(Nperm),'_001.nii']);
        end

    else
        
        % ----------------------------------------
        % Save out: wts of encodign for power analysis
        disp('saving distribution***************')
        save([proj.path.haufe.fmri_rest_1drv_rgr_a,'grp_haufe_distr.mat'],'grp_haufe_a_dist');
    end

    toc

end