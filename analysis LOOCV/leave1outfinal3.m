%%% November
%%%%%%% Based on the distances from the actual find point (badreps_Dh_ics),
%%%%%%% this script finds the best fitting behaviors based on the t-values.
%%%%%%% It takes out one ic, calculates the avg beh for the rest of the ics
%%%%%%% and saves this behavior. This beh will be used to run the simulation
%%%%%%% for the ic that was left out.

clearvars; close all; clc
set(0,'DefaultFigureWindowStyle','docked');

tic
probs = 1;
LPtype = 'hiker';
load('mapdim20_updatethresh')
ind_ic_neq_fdt = map.goodsim88;
ics = length(ind_ic_neq_fdt);
reps = 500;

%% leave one out analysis
tic
%%%% chosen ic
icf = ind_ic_neq_fdt;
L1outbp = zeros(ics,6);
finalE = zeros(ics,1);
E = zeros(ics,1);
for iic = 1:ics
    ic = ind_ic_neq_fdt(iic);
    load(['simsfinal/sim_hiker_ic',num2str(ic),'_t100.mat'],'pointdata');
    load(['fits500/bestbeh_icf',num2str(ic),'.mat'],'best_EdL');
    
    final_point = map.findxy(ic,:);
    
    % energy distance calculation
    X = zeros(reps,2);
    Y = zeros(reps,2);
    for irep = 1:reps
        x = pointdata{1,irep}(2,1);
        y = pointdata{1,irep}(2,2);
        X(irep,:) = [x, y];
        Y(irep,:) = final_point;
    end
    inan = isnan(X);        %% takes out the Nan where agent went off map
    X(inan(:,1),:) = [];
    Y(inan(:,1),:) = [];
    nreps = length(X);
    Ap = pdist2(X,Y);
    Aps = sum(Ap(:,1));       %% same as manually: sum(sqrt(sum((X-Y).^2,2)))
    A = (1/nreps)*Aps;
    Bp = 2*sum(pdist(X));      %% same as sum(squareform(pdist(X)),'all'), or sqrt((X(1,1)-X(2,1))^2+(X(1,2)-X(2,2))^2)
    B = (1/nreps^2)*Bp;
    E(iic) = 2*A - B;
    
    
    %%% sort by lowest value of E statistic
    if nreps<90
        flagreps = 1
    else
    end
    
    %%% prob dist and E values for left out IC (test set)
    L1outbp(iic,:) = best_EdL;
    finalE(iic) = E(iic);
    
end
% %%%%%%% weight 1/E^L (no threshold) %%%%%%%%%
% LE = 0.5;
% allwt_EL = finalE.^(-LE);
% sumbeh_EL = sum(L1outbp.*allwt_EL,1);
% best_EL = sumbeh_EL/norm(sumbeh_EL,1)

save(['LOOCV results/leave1outEdistresultsfinal.mat'],'finalE','L1outbp')

toc
%% box plots
% clearvars;
close all
load(['LOOCV results/allweightedbesthiker_edist500t.mat'],'E','bestbeh','allE')
load('LOOCV results/leave1outEdistresultsfinal.mat')
bestprob = zeros(length(bestbeh),1);
for ii = 1:length(bestbeh)
    bestprob(ii) = bestbeh{ii}(1,1);
end
figure, boxchart(E,'MarkerStyle','+')
hold on
plot(finalE,'rp')
plot(allE,'gp')
legend('Edist values for all probs','Edist from weighted average','Edist for best prob from original sort')
xlabel('initial conditions')
ylabel('Energy distance for all probs')
title('Comparison of energy distance values for leave 1 out analysis')
set(gcf,'PaperPosition',[0,0,11,8],'paperorientation','landscape');
print('-dpdf',['plots/leave1out_boxplotsfinal.pdf'])

%% Percentile graphs
close all
load(['LOOCV results/allweightedbesthiker_edist500t.mat'],'E','bestbeh','allE')
load('LOOCV results/leave1outEdistresultsfinal.mat')
Eadd = [E; finalE'];
[Eadds, Ia] = sort(Eadd);
Ew = zeros(ics,1);
for kk = 1:ics
    Ew(kk) = find(Ia(:,kk)==463);
end
Emed = median(E,'omitnan');
perc = (1-(Ew-1)/462)*100;


for ij = 1:ics
    ic = map.goodsim88(ij)
    figure(ij+1)
    plot(Eadds(:,ij),'.-'), hold on
    plot(Ew(ij),finalE(ij),'rp')
    xline(Ew(ij),'-.',[num2str(perc(ij),'%.1f'),' percentile'])
    yline(Emed(ij),'r-.','median - 50th percentile')
    xlabel('462 + wtavg probs sorted from best to worst')
    ylabel('E values')
    title(['IC ',num2str(ic),' [',map.key{ic},']'])
    set(gca,'fontsize',14)
    set(gcf,'PaperPosition',[0,0,10.5,8],'paperorientation','landscape');
    print('-dpdf',['plots/p_ic',num2str(ij),'.pdf'])
end

figure
for ij = 1:ics
    ic = map.goodsim88(ij)
    subplot(5,13,ij)
    plot(Eadds(:,ij),'.-'), hold on
    plot(Ew(ij),finalE(ij),'rp')
    xline(Ew(ij),'-.',[num2str(perc(ij),'%.1f'),' percentile'])
    yline(Emed(ij),'r-.')
    %     xlabel('462 + wtavg probs sorted from best to worst')
    %     ylabel('E values')
    %     title(['IC ',num2str(ic),' [',map.key{ic},']'])
    %     set(gca,'fontsize',14)
    axis square
    hold off
end
set(gcf,'PaperPosition',[0,0,10.5,8],'paperorientation','landscape');
print('-dpdf',['plots/allperc500.pdf'])
% perc = prctile(Eadd(:,1),1:100);
% value = finalE;
% [c index] = min(abs(perc'-value(1)));
% x = index+1;
%
% function x = comp_percentile(datas,value)
% perc = prctile(datas,1:100);
% [c index] = min(abs(perc'-value));
% x = index+1;
% end

%% bar chart of all the percentiles
nty = sum(perc>=90)/ics*100;
ntyf = sum(perc>=95)/ics*100;
sf = sum(perc>=75)/ics*100;
fty = sum(perc>=50)/ics*100;
figure,
bar(perc,'EdgeColor','none'), hold on
yline(90,'r-.',['90th percentile (',num2str(nty),'%)'])
yline(95,'k-.',['95th percentile (',num2str(ntyf),'%)'])
yline(75,'k-.',['75th percentile (',num2str(sf),'%)'])
yline(50,'k-.',['50th percentile (',num2str(fty),'%)'])


xlabel('initial conditions')
title('Percentiles of E-stat for Weighted Average Hiker Profile from LOO Analysis')
set(gcf,'PaperPosition',[0,0,11,8],'paperorientation','landscape');
% print('-dpdf',['plots/percentilebarfinal.pdf'])

%% hist chart of all the percentiles
figure,
histogram(perc,10), hold on
% yline(90,'k-.','90th percentile (34%)')
xlabel('percentiles'), ylabel('number of ICs in percentile range')
title('Histogram - Percentiles of E-stat for Weighted Average Hiker Profile from LOO Analysis')
set(gcf,'PaperPosition',[0,0,11,8],'paperorientation','landscape');
print('-dpdf',['plots/percentilehistfinal.pdf'])

