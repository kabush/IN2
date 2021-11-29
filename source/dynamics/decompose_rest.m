function [dcmp,indx] = decompose_rest(proj,hd)

%% ----------------------------------------
%% Calculate dynamcs from RST

%% Length params
Ntot = numel(hd);

%% Setup decomposition struct
dcmp = struct();    
indx = struct();

%%----------------------------------------
%% Construct RST trajectories

%% Construct plan base
dcmp.h = hd;
indx.h = 1:numel(hd);

%%Construct plant 1st derivative
dcmp.dh = zeros(1,numel(2:(Ntot-1)));
indx.dh = zeros(1,numel(2:(Ntot-1)));
for i = 2:(Ntot-1)
    dcmp.dh(i-1) = (dcmp.h(i+1)-dcmp.h(i-1))/2;
    indx.dh(i-1) = i;
end

%%Construct plant 2nd derivative
dcmp.d2h = zeros(1,numel(3:(Ntot-2)));
indx.d2h = zeros(1,numel(3:(Ntot-2)));
for i = 3:(Ntot-2)
    dcmp.d2h(i-2) = (dcmp.dh(i)-dcmp.dh(i-2))/2;
    indx.d2h(i-2) = i;
end

%%Construct plant 3rd derivative
dcmp.d3h = zeros(1,numel(4:(Ntot-3)));
indx.d3h = zeros(1,numel(4:(Ntot-3)));
for i = 4:(Ntot-3)
    dcmp.d3h(i-3) = (dcmp.d2h(i-1)-dcmp.d2h(i-3))/2;
    indx.d3h(i-3) = i;
end