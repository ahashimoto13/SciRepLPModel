%%% December 2020
%%%%%%% Based on the distances from the actual find point (badreps_Dh_ics),
%%%%%%% this script finds the best fitting behaviors based on the E-dist.
%%%%%%% It takes out one ic, calculates the avg beh for the rest of the ics
%%%%%%% and saves this behavior. This beh will be used to run the simulation
%%%%%%% for the ic that was left out.

clearvars; close all; clc
% tic
load('beh_dist_6.mat')
LPtype = 'hiker';
load('mapdim20_updatethresh.mat')
ind_ic_neq_fdt = map.goodsim88;
ics = length(ind_ic_neq_fdt);
reps = 500;

%% leave one out analysis
tic
%%%% chosen ic
icf = ind_ic_neq_fdt;
%%%%
newics = ics -1;          %% minus one for the one out

for iicf = 1:length(icf)
    iicf
    icout = icf(iicf);
    ics_chop = ind_ic_neq_fdt(ind_ic_neq_fdt~=icout);
    
    bestbeh = cell(newics,1);
    allbp = zeros(newics,6);
    allE = zeros(newics,1);
    E = zeros(length(probs),newics);
    dippf = zeros(newics,1);
    
    for iic = 1:newics
        ic = ics_chop(iic);
        load(['../sims500t/sim_',LPtype,'_ic',num2str(ic),'_t100.mat'],'pointdata');
        final_point = map.findxy(ic,:);
        ipp = map.icsxy(ic,:);
        
        % calculate distance from find to ipp for weighting
        dippf(iic) = pdist([final_point;ipp]);
        
        % energy distance calculation
        A = zeros(length(probs),1);
        B = zeros(length(probs),1);
        nreps = zeros(length(probs),1);
        for iprob = 1:length(probs)
            X = zeros(reps,2);
            Y = zeros(reps,2);
            for irep = 1:reps
                x = pointdata{iprob,irep}(2,1);
                y = pointdata{iprob,irep}(2,2);
                X(irep,:) = [x, y];
                Y(irep,:) = final_point;
            end
            inan = isnan(X);        %% takes out the Nan where agent went off map
            X(inan(:,1),:) = [];
            Y(inan(:,1),:) = [];
            nreps(iprob) = length(X);
            Ap = pdist2(X,Y);
            Aps = sum(Ap(:,1));       %% same as manually: sum(sqrt(sum((X-Y).^2,2)))
            A(iprob) = (1/nreps(iprob))*Aps;
            Bp = 2*sum(pdist(X));      %% same as sum(squareform(pdist(X)),'all'), or sqrt((X(1,1)-X(2,1))^2+(X(1,2)-X(2,2))^2)
            B(iprob) = (1/nreps(iprob)^2)*Bp;
            E(iprob,iic) = 2*A(iprob) - B(iprob);
            
        end
        %%% sort by lowest value of E statistic
        E(nreps<98,iic) = nan;                  %%% threshold for number of good reps
        [Ev,Iprob] = sort(E(:,iic));
        
        %%% prob dist and E values for top fits
        bestbeh{iic} = [Iprob,Ev];
        allbp(iic,:) = probs(Iprob(1),:);
        allE(iic) = Ev(1);
        
    end
    %%%%%%% weight 1/E^L (no threshold) %%%%%%%%%
    LE = 0.5;
    allwt_EL = allE.^(-LE);
    sumbeh_EL = sum(allbp.*allwt_EL,1);
    best_EL = sumbeh_EL/norm(sumbeh_EL,1)
    
    %%%% best behavioral profile for all ICs - weighted by 1/(E/d)^L
    allwt_EdL = (allE./dippf).^(-LE);
    sumbeh_EdL = sum(allbp.*allwt_EdL,1);
    best_EdL = sumbeh_EdL/norm(sumbeh_EdL,1)
    
    save(['fits500/bestbeh_icf',num2str(icout),'.mat'],'dippf','bestbeh','allbp','allE','E','LE','best_EL','best_EdL','allwt_EL','allwt_EdL','ics_chop')

end
toc

main_LOO2

