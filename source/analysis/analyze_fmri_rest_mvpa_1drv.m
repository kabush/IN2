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
logger([' Analyzing 1st Deriv Predictions                 '],proj.path.logfile);
logger(['*************************************************'],proj.path.logfile);

%% ----------------------------------------
%% Set-up Directory Structure
if(proj.flag.clean_build)
    disp(['Removing ',proj.path.analysis.analyze_fmri_rest_mvpa_1drv]);
    eval(['! rm -rf ',proj.path.analysis.analyze_fmri_rest_mvpa_1drv]);
    disp(['Creating ',proj.path.analysis.analyze_fmri_rest_mvpa_1drv]);
    eval(['! mkdir ',proj.path.analysis.analyze_fmri_rest_mvpa_1drv]);
end

%% ----------------------------------------
%% load subjs
subjs = load_subjs(proj);

%% Valence
logger(['----------------------------------------'], ...
       proj.path.logfile);
logger(['Valence'],proj.path.logfile);

%% Set figure
figure(3)
set(gcf,'color','w');

%% ----------------------------------------
%% scatter the underlying rest
predictors = [];
measures = [];
subjects = [];

sig_cnt = 1;
non_cnt = 1;

sig_subjs = {};
non_subjs = {};

b_all = []; %% for power
sex = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log analysis of subject
    logger([subj_study,'_',name],proj.path.logfile);

    try

        %% Load IN trajectory structures
        load([proj.path.mvpa.fmri_rest_1drv_rgr,subj_study,'_',name,'_prds.mat']);

        if(isfield(prds,'v'))
            
            out = prds.v.out;
            trg = prds.v.trg;

            %% build data for group GLMM
            predictors = [predictors;double(out)];
            measures = [measures;double(trg)];
            subjects = [subjects;repmat(i,numel(out),1)];
            
            %% scatter plot specific points        
            scatter(out,trg,10,'MarkerFaceColor', ...
                    proj.param.plot.white,'MarkerEdgeColor', ...
                    proj.param.plot.light_grey);
            hold on;
            
            %% build individual subject structures
            subj = struct();
            subj.study = subj_study;
            subj.name = name;
            
            subj.out = out;
            subj.trg = trg;

            demo = readtable(['/raw/bush/demo/',subj_study,'.csv']);
            id = find(strcmp(demo.ID,name)~=0);
            sex = [sex,demo.Type(id)];
            
            [b stat] = robustfit(out,trg);
            b_all = [b_all,b(2)]; %% for power

            subj.b1 = b(2); % slope
            subj.b0 = b(1); % intercept
            subj.p1 = stat.p(2); %slope
            subj.p0 = stat.p(1); %intercept
 
            %% sort subjects by significance
             if(subj.p1<0.05)
                 sig_subjs{sig_cnt} = subj;
                 sig_cnt = sig_cnt + 1;
             else
                 non_subjs{non_cnt} = subj;
                 non_cnt = non_cnt + 1;
             end
             
        else
            disp(['  -Could not find v_dcmp for: ',subj_study,'_', ...
                  name],proj.path.logfile);
        end
        
     catch
         % do nothing
        logger(['  -Could not find/load prds for: ',subj_study,'_', ...
                name],proj.path.logfile);
    end

end

%% ----------------------------------------
%% identify max/min x-range|y-rang

xmin = -3;
xmax = 3;
ymin = -1.5;
ymax = 1.5;
vseq = linspace(xmin,xmax);

%% ----------------------------------------
%% overlay the optimal outcome
% plot(vseq,vseq,'k:','LineWidth',2);

%% ----------------------------------------
%% overlay the individual 1st deriv plot
for i =1:numel(non_subjs)
    plot(non_subjs{i}.out,non_subjs{i}.out*non_subjs{i}.b1+ ...
         non_subjs{i}.b0,'Color',proj.param.plot.light_grey, ...
         'LineWidth',1);
end

for i =1:numel(sig_subjs)
    plot(sig_subjs{i}.out,sig_subjs{i}.out*sig_subjs{i}.b1+ ...
         sig_subjs{i}.b0,'Color',proj.param.plot.dark_grey, ...
         'LineWidth',2);
end

%% ----------------------------------------
%% Group GLMM fit
tbl = table(measures,predictors,subjects,'VariableNames', ...
            {'measures','predictors','subjects'});
mdl = fitlme(tbl,['measures ~ 1 + predictors + (1|subjects) + ' ...
                  '(predictors-1|subjects)']);
[~,~,FE] = fixedEffects(mdl);

%% Save out model
save([proj.path.analysis.analyze_fmri_rest_mvpa_1drv,'valence_1drv_mdl.mat'],'mdl');

%% ----------------------------------------
%% overlay the group 1 deriv plot
y_hat = FE.Estimate(1) + FE.Estimate(2)*vseq;
plot(vseq,y_hat,'r-','LineWidth',3);

%% ----------------------------------------
%% compute effect size
SS_res=sum((mdl.residuals).^2);
SS_tot=sum((measures-mean(measures)).^2);
Rsqr = 1-(SS_res/SS_tot);
Fsqr = Rsqr/(1-Rsqr);
logger(['Rsqr=',num2str(Rsqr)],proj.path.logfile);
logger(['Fsqr=',num2str(Fsqr)],proj.path.logfile);

%% ----------------------------------------
%% format figure
xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;
set(gca,'Layer','top');

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'valence_mvpa_rest_1drv.png' -r300  
eval(['! mv ',proj.path.code,'valence_mvpa_rest_1drv.png ', ...
      proj.path.fig]);

%% ----------------------------------------
%% Output summary

% fixed effects
ge = num2str(FE.Estimate(2));
gep = num2str(FE.pValue(2));

% single subject results
n_sig = numel(sig_subjs);
n_tot = numel(sig_subjs)+numel(non_subjs);
ns = num2str(n_sig);
nt = num2str(n_tot);
ssp = num2str(100*(n_sig/n_tot));

n_tot
n_sig

logger(['----------------------------------------'], ...
       proj.path.logfile);
logger(['Valence Statistical Summary'],proj.path.logfile);
logger(['  -Group effect: ',ge,', p=',gep],proj.path.logfile);
logger(['  -Percent sign. subjs. (p<0.05): ',ssp,'% (',ns,'/',nt, ...
        ')'],proj.path.logfile);

b_m = b_all(find(sex==1));
b_f = b_all(find(sex==2));
p = ranksum(b_m,b_f);
disp(['ranksum of b_all sex diffs, p=',num2str(p)]);

%% Arousal
logger(['----------------------------------------'], ...
       proj.path.logfile);
logger(['Arousal'],proj.path.logfile);

%% Set figure
figure(4)
set(gcf,'color','w');

%% ----------------------------------------
%% scatter the underlying rest
predictors = [];
measures = [];
subjects = [];

sig_cnt = 1;
non_cnt = 1;

sig_subjs = {};
non_subjs = {};

b_all = []; %% for power
sex = [];

for i = 1:numel(subjs)

    %% extract subject info
    subj_study = subjs{i}.study;
    name = subjs{i}.name;
    id = subjs{i}.id;

    % log analysis of subject
    logger([subj_study,'_',name],proj.path.logfile);

    try

        %% Load IN trajectory structures
        load([proj.path.mvpa.fmri_rest_1drv_rgr,subj_study,'_',name,'_prds.mat']);

        if(isfield(prds,'a'))
            
            out = prds.a.out;
            trg = prds.a.trg;

            %% build data for group GLMM
            predictors = [predictors;double(out)];
            measures = [measures;double(trg)];
            subjects = [subjects;repmat(i,numel(out),1)];
            
            %% scatter plot specific points        
            scatter(out,trg,10,'MarkerFaceColor', ...
                    proj.param.plot.white,'MarkerEdgeColor', ...
                    proj.param.plot.light_grey);
            hold on;
            
            %% build individual subject structures
            subj = struct();
            subj.study = subj_study;
            subj.name = name;
            
            subj.out = out;
            subj.trg = trg;

            demo = readtable(['/raw/bush/demo/',subj_study,'.csv']);
            id = find(strcmp(demo.ID,name)~=0);
            sex = [sex,demo.Type(id)];
            
            [b stat] = robustfit(out,trg);
            b_all = [b_all,b(2)]; %% for power

            subj.b1 = b(2); % slope
            subj.b0 = b(1); % intercept
            subj.p1 = stat.p(2); %slope
            subj.p0 = stat.p(1); %intercept
 
            %% sort subjects by significance
             if(subj.p1<0.05)
                 sig_subjs{sig_cnt} = subj;
                 sig_cnt = sig_cnt + 1;
             else
                 non_subjs{non_cnt} = subj;
                 non_cnt = non_cnt + 1;
             end
             
        else
            disp(['  -Could not find a_dcmp for: ',subj_study,'_', ...
                  name],proj.path.logfile);
        end
        
     catch
         % do nothing
        logger(['  -Could not find/load prds for: ',subj_study,'_', ...
                name],proj.path.logfile);
    end

end

%% ----------------------------------------
%% identify max/min x-range|y-rang

xmin = -3;
xmax = 3;
ymin = -1.5;
ymax = 1.5;
vseq = linspace(xmin,xmax);

%% ----------------------------------------
%% overlay the optimal outcome
% plot(vseq,vseq,'k:','LineWidth',2);

%% ----------------------------------------
%% overlay the individual 1st deriv plot
for i =1:numel(non_subjs)
    plot(non_subjs{i}.out,non_subjs{i}.out*non_subjs{i}.b1+ ...
         non_subjs{i}.b0,'Color',proj.param.plot.light_grey, ...
         'LineWidth',1);
end

for i =1:numel(sig_subjs)
    plot(sig_subjs{i}.out,sig_subjs{i}.out*sig_subjs{i}.b1+ ...
         sig_subjs{i}.b0,'Color',proj.param.plot.dark_grey, ...
         'LineWidth',2);
end

%% ----------------------------------------
%% Group GLMM fit
tbl = table(measures,predictors,subjects,'VariableNames', ...
            {'measures','predictors','subjects'});
mdl = fitlme(tbl,['measures ~ 1 + predictors + (1|subjects) + ' ...
                  '(predictors-1|subjects)']);
[~,~,FE] = fixedEffects(mdl);

%% Save out model
save([proj.path.analysis.analyze_fmri_rest_mvpa_1drv,'arousal_1drv_mdl.mat'],'mdl');

%% ----------------------------------------
%% overlay the group 1 deriv plot
y_hat = FE.Estimate(1) + FE.Estimate(2)*vseq;
plot(vseq,y_hat,'r-','LineWidth',3);

%% ----------------------------------------
%% compute effect size
SS_res=sum((mdl.residuals).^2);
SS_tot=sum((measures-mean(measures)).^2);
Rsqr = 1-(SS_res/SS_tot);
Fsqr = Rsqr/(1-Rsqr);
logger(['Rsqr=',num2str(Rsqr)],proj.path.logfile);
logger(['Fsqr=',num2str(Fsqr)],proj.path.logfile);

%% ----------------------------------------
%% format figure
xlim([xmin,xmax]);
ylim([ymin,ymax]);

hold off;
fig = gcf;
ax = fig.CurrentAxes;
ax.FontSize = proj.param.plot.axisLabelFontSize;
set(gca,'Layer','top');

%% ----------------------------------------
%% explot hi-resolution figure
export_fig 'arousal_mvpa_rest_1drv.png' -r300  
eval(['! mv ',proj.path.code,'arousal_mvpa_rest_1drv.png ', ...
      proj.path.fig]);

%% ----------------------------------------
%% Output summary

% fixed effects
ge = num2str(FE.Estimate(2));
gep = num2str(FE.pValue(2));

% single subject results
n_sig = numel(sig_subjs);
n_tot = numel(sig_subjs)+numel(non_subjs);
ns = num2str(n_sig);
nt = num2str(n_tot);
ssp = num2str(100*(n_sig/n_tot));

n_tot
n_sig

logger(['----------------------------------------'], ...
       proj.path.logfile);
logger(['Arousal Statistical Summary'],proj.path.logfile);
logger(['  -Group effect: ',ge,', p=',gep],proj.path.logfile);
logger(['  -Percent sign. subjs. (p<0.05): ',ssp,'% (',ns,'/',nt, ...
        ')'],proj.path.logfile);

b_m = b_all(find(sex==1));
b_f = b_all(find(sex==2));
p = ranksum(b_m,b_f);
disp(['ranksum of b_all sex diffs, p=',num2str(p)]);