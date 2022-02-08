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
logger([' Simulate REST affect dynamics VALENCE'],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure for fMRI betas
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.dyn.sim_rest_v]);
    eval(['! rm -rf ',proj.path.dyn.sim_rest_v]);
    disp(['Creating ',proj.path.dyn.sim_rest_v]);
    eval(['! mkdir ',proj.path.dyn.sim_rest_v]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% ----------------------------------------
%% persistent data storage for analysis
running_all_err_free = [];
b_vals = [];

%% ----------------------------------------
%% loop simulation
loops = 30; % key param

for loop = 1:loops
    
    % log iteration
    logger(['Simulation Iteration ',num2str(loop)],proj.path.logfile);
    
    %% data storage
    all_err_rplc = [];
    all_err_free = [];
    bad_count = 0;

    %% ----------------------------------------
    %% Transform beta-series into affect series {v}
    for i = 1:numel(subjs)
        
        %% extract subject info
        subj_study = subjs{i}.study;
        name = subjs{i}.name;
        id = subjs{i}.id;

        try
            
            nsteps = [];
            rmse = [];

            %% Load affective state predictions
            load([proj.path.dyn.rest,subj_study,'_',name,'_prds.mat'],'prds');
            h_indx = 6:215;
            h_trg = prds.v_hd(h_indx);

            clear prds;

            %% Load 2drv (VAL) predictions
            load([proj.path.mvpa.fmri_rest_2drv_rgr,subj_study,'_',name,'_prds.mat'],'prds');
            d2h_out = prds.v.out;
            d2h_trg = prds.v.trg;
            d2h_indx = 6:215; %% equates to 1-210;
            
            clear prds;

            %% Load 2drv NaN mask
            load([proj.path.mvpa.fmri_rest_2drv_rgr,subj_study,'_',name,'_v_nan_mask.mat'],'v_nan_mask');
            v_nan_mask_2drv = v_nan_mask;

            %% Add back NaN to equalize vector lengths
            for ii=1:length(v_nan_mask_2drv)
                if (v_nan_mask_2drv(ii) && ii==1)
                    d2h_out = [NaN; d2h_out(ii:end)];
                    d2h_trg = [NaN; d2h_trg(ii:end)];
                elseif (v_nan_mask_2drv(ii) && ii==210)
                    d2h_out = [d2h_out; NaN];
                    d2h_trg = [d2h_trg; NaN];
                elseif v_nan_mask_2drv(ii)
                    d2h_out = [d2h_out(1:(ii-1)); NaN; d2h_out(ii:end)];
                    d2h_trg = [d2h_trg(1:(ii-1)); NaN; d2h_trg(ii:end)];
                end
            end

            %% Load 1drv (VAL) predictions
            load([proj.path.mvpa.fmri_rest_1drv_rgr,subj_study,'_',name,'_prds.mat'],'prds');
            dh_out = prds.v.out;
            dh_trg = prds.v.trg;
            dh_indx = 6:215;  %% equates to 1-210;

            %% Load 1drv NaN mask
            load([proj.path.mvpa.fmri_rest_1drv_rgr,subj_study,'_',name,'_v_nan_mask.mat'],'v_nan_mask');
            v_nan_mask_1drv = v_nan_mask;

            %% Add back NaN to equalize vector lengths
            for ii=1:length(v_nan_mask_1drv)
                if (v_nan_mask_1drv(ii) && ii==1)
                    dh_out = [NaN; dh_out(ii:end)];
                    dh_trg = [NaN; dh_trg(ii:end)];
                elseif (v_nan_mask_1drv(ii) && ii==210)
                    dh_out = [dh_out; NaN];
                    dh_trg = [dh_trg; NaN];
                elseif v_nan_mask_1drv(ii)
                    dh_out = [dh_out(1:(ii-1)); NaN; dh_out(ii:end)];
                    dh_trg = [dh_trg(1:(ii-1)); NaN; dh_trg(ii:end)];
                end
            end

            %% Collect reduced vectors
            indx = 7:214; %% prune ends

            rdx_d2h_trg = d2h_trg;
            rdx_d2h_out = d2h_out; %% should be same b/c limiting factor

            rdx_dh_trg =[]; 
            rdx_dh_out =[]; 
            rdx_h_trg =[]; 

            for j= 1:numel(indx)

                % gather reduced 1st deriv trgs/prds
                rdx_dh_trg = [rdx_dh_trg,dh_trg(find(dh_indx==indx(j)))];
                rdx_dh_out = [rdx_dh_out,dh_out(find(dh_indx==indx(j)))];

                % gather reduced state trgs
                rdx_h_trg = [rdx_h_trg,h_trg(find(h_indx==indx(j)))];

            end

            %% Logic section to check for potential for valid simulation (no NaN)
            good_indx = [];
            max_Nstep = [];

            for j = 1:(numel(indx)-1)
                %% First check current index for NaNs
                if (isnan(rdx_h_trg(j))~=1 && isnan(rdx_dh_out(j))~=1 && isnan(rdx_d2h_out(j))~=1)
                    for Nstep = 1:6
                        %% Make sure checking Nsteps would not exceed matrix dimensions
                        if (Nstep==2 && j>(numel(indx)-2))
                            Nstep=Nstep-1;
                            break
                        elseif (Nstep==3 && j>(numel(indx)-3))
                            Nstep=Nstep-1;
                            break
                        elseif (Nstep==4 && j>(numel(indx)-4))
                            Nstep=Nstep-1;
                            break
                        elseif (Nstep==5 && j>(numel(indx)-5))
                            Nstep=Nstep-1;
                            break
                        elseif (Nstep==6 && j>(numel(indx)-6))
                        Nstep=Nstep-1;
                        break
                        end
                        %% Make sure there are no NaNs when moving Nsteps forward
                        if (isnan(rdx_h_trg(j+Nstep))~=1 && isnan(rdx_d2h_out(j+Nstep))~=1)
                            continue
                        else
                            Nstep=Nstep-1;
                            break
                        end
                    end
                    %% Check we can go at least 1 Nstep forward
                    if Nstep >=1
                        good_indx = [good_indx, j];
                        max_Nstep = [max_Nstep, Nstep];
                    end
                end
            end

            %% Setup data vectors
            err_rplc = NaN(numel(good_indx),1);
            sham_err_rplc = NaN(numel(good_indx),1);
            real_err_rplc = NaN(numel(good_indx),1);

            err_free = NaN(numel(good_indx),max(max_Nstep));
            sham_err_free = NaN(numel(good_indx),max(max_Nstep));
            real_err_free = NaN(numel(good_indx),max(max_Nstep));

            %% Iterate through good indices to conduct simulation
            for ii = 1:numel(good_indx)

                %% Integrate with replacement (TRUE)
                h = rdx_h_trg(good_indx(ii));
                dh = rdx_dh_out(good_indx(ii));

                h_new = h+dh;            % predict state
                h_trg = rdx_h_trg(good_indx(ii)+1);  % get true state
                real_err_rplc(ii,1) = h_trg-h_new;  % compute error

                %% Integrate with replacement (SHAM)
                sham_indx = randsample(1:numel(good_indx),1);

                h = rdx_h_trg(good_indx(ii));
                dh = rdx_dh_out(good_indx(sham_indx));

                h_new = h+dh;            % predict state
                h_trg = rdx_h_trg(good_indx(ii)+1);  % get true state
                sham_err_rplc(ii,1) = h_trg-h_new;  % compute error

                %% Integration length (KEY VARIABLE)
                Nstep = max_Nstep(ii);

                %% Integrate without replacement (TRUE)
                h = rdx_h_trg(good_indx(ii));      %known state
                dh = rdx_dh_out(good_indx(ii));    %estimate dh
                d2h = rdx_d2h_out(good_indx(ii));  %estimate d2h

                dh_sim = [];

                for k=1:Nstep
                    dh_sim = [dh_sim,dh];
                    h = h+dh;
                    dh = dh+d2h;    
                    d2h = rdx_d2h_out(good_indx(ii)+k);  %estimate d2h

                    %Conduct RMSE analysis
                    if loop == 1
                        mse_current = sqrt(immse(dh_sim, rdx_dh_out(good_indx(ii):(good_indx(ii)+k-1))));
                        rmse = [rmse; mse_current];
                        nsteps = [nsteps; k];
                    end

                    h_new = h;
                    h_trg = rdx_h_trg(good_indx(ii)+k);  % get true state
                    real_err_free(ii,k) = h_trg-h_new;  % compute error

                end

                %% Integrate without replacement (SHAM)
                h = rdx_h_trg(good_indx(ii));      %known state

                sham_indx = randsample(1:numel(good_indx),1);
                dh = rdx_dh_out(good_indx(sham_indx));

                sham_indx = randsample(1:numel(good_indx),1);
                d2h = rdx_d2h_out(good_indx(sham_indx)); 

                for k=1:Nstep

                    h = h+dh;
                    dh = dh+d2h;

                    sham_indx = randsample(1:numel(good_indx),1);
                    d2h = rdx_d2h_out(good_indx(sham_indx));  %estimate d2h

                    h_trg = rdx_h_trg(good_indx(ii)+k);  % get true state
                    sham_err_free(ii,k) = h_trg-h;  % compute error

                end

            end

            %% Data averaging for analysis

            err_rplc = (median(real_err_rplc.^2,'omitnan')-median(sham_err_rplc.^2,'omitnan'));
            all_err_rplc = [all_err_rplc,err_rplc];

            err_free = (median(real_err_free(1:numel(good_indx),:).^2,'omitnan')-median(sham_err_free(1:numel(good_indx),:).^2,'omitnan'));
            all_err_free = [all_err_free;err_free];
            
            %% MSE analysis of simulation via iteratively reweighted least-squares regression
            
            if loop==1
                [b,stat] = robustfit(nsteps,rmse);
                b_vals = [b_vals, b(2)];
            end

        catch
            % log subject
            logger([subj_study,'_',name],proj.path.logfile);
            logger(['   -error in fmri rest simulation'],proj.path.logfile);
            bad_count = bad_count + 1;
        end
        
    end
    
    % Group-level stats for MSE analysis
    if loop==1
        mean_b_val=mean(b_vals);
        [h_val,p_val,ci_val,stats_val]=ttest(b_vals);
    end
        
    running_all_err_free = cat(3,running_all_err_free,all_err_free);
    
end

%% Calculate means over simulation
mean_all_err_free = mean(running_all_err_free,3);

%% Plot Figure
close all

figure(1)
set(gcf,'color','w');

plot(mean_all_err_free(1,:),'color',proj.param.plot.very_light_grey,'LineWidth',1);
hold on;
for i=2:size(mean_all_err_free,1)
    plot(mean_all_err_free(i,:),'color',proj.param.plot.very_light_grey,'LineWidth',1);
end

lo_bound = [];
hi_bound = [];
for i=1:size(mean_all_err_free,2)
    [h p ci stat] = ttest(mean_all_err_free(:,i));
    lo_bound = [lo_bound,ci(1)];
    hi_bound = [hi_bound,ci(2)];
end

plot(lo_bound,'r-','LineWidth',2);
plot(hi_bound,'r-','LineWidth',2);
plot([1,6],[0,0],'k--');
plot(mean(mean_all_err_free),'r-','LineWidth',3);

hold off;

ylim([-6,1]);
set(gca,'Layer','top');

% logger(['comparing deriv vs sham'],proj.path.logfile);
% [p,h,stats] = ranksum(all_err_free(:,1),sham_err_free(:,1));
% Nsmall = min(numel(all_err_free(:,1)),numel(sham_err_free(:,1)));
% r = stats.zval/sqrt(Nsmall)
% logger(['***EFFECT*** r(z/sqrt(nsmall))=',num2str(r)],proj.path.logfile);
logger(['Bad Subject Count =',num2str(bad_count)],proj.path.logfile);

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'valence_simulation.png' -r300  
eval(['! mv ',proj.path.code,'valence_simulation.png ', ...
      proj.path.fig]);