% find the energy statistic for the distance between the closest
% sim point and the find point from Bob's data for all 462 beh and reps
clearvars; close all; clc
set(0,'DefaultFigureWindowStyle','docked');

tic
% take out ics where ipp=find
load('beh_dist_6.mat')
load('mapdim20_updatethresh.mat')
ind_ic_neq_fdt = map.goodsim88;
reps = 500;
LPtype = 'hiker';
ics = length(ind_ic_neq_fdt);
simname = {'0','1','2','3','4','500t'};
isim = 6;
%% calculate the best fits based on energy distance - not leaving any out
% for isim = 1:5
bestbeh = cell(ics,1);
allbp = zeros(ics,6);
allE = zeros(ics,1);
E = zeros(length(probs),ics);
dippf = zeros(ics,1);
for iic = 1:ics
    iic
    ic = ind_ic_neq_fdt(iic);
    load(['../sims',simname{isim},'/sim_',LPtype,'_ic',num2str(ic),'_t100.mat'],'pointdata');
    final_point = map.findxy(ic,:);
    ipp = map.icsxy(ic,:);
    
    % calculate distance from find to ipp for weighting
    dippf(iic) = pdist([final_point;ipp]);
    
    % energy distance calculation
    A = zeros(length(probs),1);
    B = zeros(length(probs),1);
    nreps = zeros(length(probs),1);
    for iprob = 1:length(probs)
        [iic iprob]
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
    
    %     [Ev,Iprob] = sort(E(:,iic));
    %     Esort = E;
    %     Esort(nreps<48) = nan;                  %%% threshold for number of good reps
    %     [Evthresh,Iprobthresh] = sort(Esort(:,iic));
    %     aux = [Iprob(1:15), nreps(Iprob(1:15))];
    %     aux1 = [aux1, aux];
    %     % prob dist and E values for top fits
    %     bestbeht{iic} = [Iprobthresh,Evthresh];
    %     allbpt(iic,:) = probs(Iprobthresh(1),:);
    %     allEt(iic) = Evthresh(1);
    
    %%% prob dist and E values for top fits
    bestbeh{iic} = [Iprob,Ev];
    allbp(iic,:) = probs(Iprob(1),:);
    allE(iic) = Ev(1);
    
end
save(['best fit data/allbesthiker_edist',simname{isim},'.mat'],'bestbeh','allbp','allE','E','dippf')
% end
% save(['sims/allbesthiker_edist500t.mat'],'bestbeh','allbp','allE','E','dippf')

%%
% load('sims/allbesthiker_edist500t.mat')
% for isim = 1:5
load(['best fit data/allbesthiker_edist',simname{isim},'.mat'])

%%%% best behavioral profile for all ICs - weighted by 1/E^L
LE = .5;
allwt_EL = allE.^(-LE);
sumbeh_EL = sum(allbp.*allwt_EL,1);
best_EL = sumbeh_EL/norm(sumbeh_EL,1)
%%%% best behavioral profile for all ICs - weighted by 1/(E/d)^L
ipps = map.icsxy(ind_ic_neq_fdt,:);
finds = map.findxy(ind_ic_neq_fdt,:);   % dippf = sqrt(sum((ipps-finds).^2,2))
ndippf = sqrt(sum((ipps-finds).^2,2));
allwt_EdL = (allE./ndippf).^(-LE);
sumbeh_EdL = sum(allbp.*allwt_EdL,1);
best_EdL = sumbeh_EdL/norm(sumbeh_EdL,1)
save(['best fit data/allweightedbesthiker_edist',simname{isim},'.mat'],'ndippf','bestbeh','allbp','allE','E','LE','best_EL','best_EdL','allwt_EL','allwt_EdL')
% end
%% find the closest reps

clpts = cell(ics,1);
for iic = 1:ics
    ic = ind_ic_neq_fdt(iic);
    load(['../sims',simname{isim},'/sim_',LPtype,'_ic',num2str(ic),'_t100.mat'],'pointdata');
    findxy = map.findxy(ic,:);
    for iprob = 1:462
        reppts = zeros(reps,2);
        for irep = 1:reps
            reppts(irep,:) = pointdata{iprob,irep}(2,1:2);
        end
        [k,dist] = dsearchn(reppts,findxy);
        clpts{iic}(iprob,:) = reppts(k,:);
    end
end
save(['best fit data/allclosestpoints',simname{isim},'.mat'],'clpts')
%% find the time of the closest reps

clpts = cell(ics,1);
tclpts = cell(ics,1);
cprep = cell(ics,1);
reps = 500;
for iic = 1:ics
    ic = ind_ic_neq_fdt(iic);
    load(['../sims',simname{isim},'/sim_',LPtype,'_ic',num2str(ic),'_t100.mat'],'pointdata');
    findxy = map.findxy(ic,:);
    for iprob = 1:462
        reppts = zeros(reps,2);
        treppts = zeros(reps,1);
        for irep = 1:reps
            reppts(irep,:) = pointdata{iprob,irep}(2,1:2);
            treppts(irep,:) = pointdata{iprob,irep}(4,1);
        end
        [k,dist] = dsearchn(reppts,findxy);
        clpts{iic}(iprob,:) = reppts(k,:);
        tclpts{iic}(iprob,:) = treppts;
        cprep{iic}(iprob,:) = k;
    end
end
save('best fit data/closestpoints_time500t.mat','clpts','tclpts','cprep')


%% plots of traj points for best fit behaviors
load('mapdim20_updatethresh')
load(['best fit data/allweightedbesthiker_edist',simname{isim},'.mat'])
load(['best fit data/allclosestpoints',simname{isim},'.mat'])
% load('edist/besthiker_edist.mat')
ind_ic_neq_fdt = map.goodsim88;
col = [0.1 0.2 0.2; 0.4660 0.6740 0.1880; 0 0.7 0.7];
reps = 500;
for iic = 1:ics
    figure(iic)
    ic = ind_ic_neq_fdt(iic)
    bestp = bestbeh{iic}(1,1);
    load(['../sims',simname{isim},'/sim_',LPtype,'_ic',num2str(ic),'_t100.mat']);
    otherp = [2, 5, 10];
    for iiprob = 1:3
        iprob = bestbeh{iic}(otherp(iiprob),1);
        for ii = 1:reps
            cpx = pointdata{iprob,ii}(2,1);
            cpy = pointdata{iprob,ii}(2,2);
            pl(iiprob,ii) = scatter(cpx,cpy,50,'o','markerfacecolor',col(iiprob,:),'markeredgecolor',col(iiprob,:));
            alpha(pl(iiprob,ii),0.4)
            hold on
        end
        pc(iiprob) = plot(clpts{iic}(iprob,1),clpts{iic}(iprob,2),'p','markerfacecolor',col(iiprob,:),'markersize',12,'MarkerEdgeColor','k',...
            'Displayname','Closest Rep');
        hold on
        
    end
    pl(1).DisplayName = ['2nd ',num2str(bestbeh{iic}(2,1)),' ',mat2str(probs(bestbeh{iic}(2,1),:),2)];
    pl(2).DisplayName = ['5th ',num2str(bestbeh{iic}(5,1)),' ',mat2str(probs(bestbeh{iic}(5,1),:),2)];
    pl(3).DisplayName = ['10th ',num2str(bestbeh{iic}(10,1)),' ',mat2str(probs(bestbeh{iic}(10,1),:),2)];
    for ii = 1:reps
        xb = pointdata{bestp,ii}(2,1);
        yb = pointdata{bestp,ii}(2,2);
        pb = plot(xb,yb,'o','markerfacecolor',[0.8500 0.3250 0.0980],'markeredgecolor','k','markersize',7,...
            'Displayname',['1st ',num2str(bestp),' ',mat2str(probs(bestp,:),2)]);
        hold on
    end
    pc(4) = plot(clpts{iic}(bestp,1),clpts{iic}(bestp,2),'p','markerfacecolor',[0.8500 0.3250 0.0980],'markersize',12,'MarkerEdgeColor','k',...
        'Displayname','1st best beh closest rep');
    pc(1).DisplayName = '2nd best beh closest rep';
    pc(2).DisplayName = '5th best beh closest rep';
    pc(3).DisplayName = '10th best beh closest rep';
    legend([pb,pl(1),pl(2),pl(3),pc],'location','northeastoutside')
    
    findxy = map.findxy(ic,:);
    plot(findxy(1),findxy(2),'p','markerfacecolor','y','markeredgecolor','k',...
        'linewidth',1.5,'markersize',15,'displayname','find')
    icxy = map.icsxy(ic,:);
    plot(icxy(1),icxy(2),'gs','markerfacecolor','g','markeredgecolor','k','markersize',10,'displayname','ipp')
    axis square, xlim([0 map.dim(1)]), ylim([0 map.dim(1)]), box on
    if strcmp(map.mobdomain{ic},'Dry')
        sgtitle(['closest pt - ic',num2str(ic),' [',map.key{ic},'], beh ',num2str(bestp),' [RW RT DT SP VE BT], E = ',num2str(allE(iic)),', w_{E/d}=',num2str(allwt_EdL(iic)),' (Dry)'])
    else
        sgtitle(['closest pt - ic',num2str(ic),' [',map.key{ic},'], beh ',num2str(bestp),' [RW RT DT SP VE BT], E = ',num2str(allE(iic)),', w_{E/d}=',num2str(allwt_EdL(iic)),' (Temperate)'])
    end
    hold off
    set(gca,'fontname','times','fontsize',14)
    set(gcf,'PaperPosition',[0,0,10.5,8],'paperorientation','landscape');
    print('-dpdf',['plots/traj_ic',num2str(ic),'.pdf'])
    close all
end

%% plots of final traj points for best fit behaviors - sanity check
close all
ind_ic_neq_fdt = map.goodsim88;
col = [0.1 0.2 0.2; 0.4660 0.6740 0.1880; 0 0.7 0.7];
% reps = 100;
for iic = 1:ics
    figure(iic)
    ic = ind_ic_neq_fdt(iic)
    bestp = bestbeh{iic}(1,1);
    load(['../sims',simname{isim},'/sim_',LPtype,'_ic',num2str(ic),'_t100.mat']);
    otherp = [2, 5, 10];
    for iiprob = 1:3
        iprob = bestbeh{iic}(otherp(iiprob),1);
        for ii = 1:reps
            cpx = pointdata{iprob,ii}(3,1);
            cpy = pointdata{iprob,ii}(3,2);
            pl(iiprob,ii) = scatter(cpx,cpy,50,'o','markerfacecolor',col(iiprob,:),'markeredgecolor',col(iiprob,:));
            alpha(pl(iiprob,ii),0.4)
            hold on
        end
        pc(iiprob) = plot(clpts{iic}(iprob,1),clpts{iic}(iprob,2),'p','markerfacecolor',col(iiprob,:),'markersize',12,'MarkerEdgeColor','k',...
            'Displayname','Closest Rep');
        hold on
        
    end
    pl(1).DisplayName = ['2nd ',num2str(bestbeh{iic}(2,1)),' ',mat2str(probs(bestbeh{iic}(2,1),:),2)];
    pl(2).DisplayName = ['5th ',num2str(bestbeh{iic}(5,1)),' ',mat2str(probs(bestbeh{iic}(5,1),:),2)];
    pl(3).DisplayName = ['10th ',num2str(bestbeh{iic}(10,1)),' ',mat2str(probs(bestbeh{iic}(10,1),:),2)];
    for ii = 1:reps
        xb = pointdata{bestp,ii}(3,1);
        yb = pointdata{bestp,ii}(3,2);
        pb = plot(xb,yb,'o','markerfacecolor',[0.8500 0.3250 0.0980],'markeredgecolor','k','markersize',7,...
            'Displayname',['1st ',num2str(bestp),' ',mat2str(probs(bestp,:),2)]);
        hold on
    end
    pc(4) = plot(clpts{iic}(bestp,1),clpts{iic}(bestp,2),'p','markerfacecolor',[0.8500 0.3250 0.0980],'markersize',12,'MarkerEdgeColor','k',...
        'Displayname','1st best beh closest rep');
    pc(1).DisplayName = '2nd best beh closest rep';
    pc(2).DisplayName = '5th best beh closest rep';
    pc(3).DisplayName = '10th best beh closest rep';
    legend([pb,pl(1),pl(2),pl(3),pc],'location','northeastoutside')
    
    findxy = map.findxy(ic,:);
    plot(findxy(1),findxy(2),'p','markerfacecolor','y','markeredgecolor','k',...
        'linewidth',1.5,'markersize',15,'displayname','find')
    icxy = map.icsxy(ic,:);
    plot(icxy(1),icxy(2),'gs','markerfacecolor','g','markeredgecolor','k','markersize',10,'displayname','ipp')
    axis square, xlim([0 map.dim(1)]), ylim([0 map.dim(1)]), box on
    if strcmp(map.mobdomain{ic},'Dry')
        sgtitle(['final pt - ic',num2str(ic),' [',map.key{ic},'], beh ',num2str(bestp),' [RW RT DT SP VE BT], E = ',num2str(allE(iic)),', w_{E/d}=',num2str(allwt_EdL(iic)),' (Dry)'])
    else
        sgtitle(['final pt - ic',num2str(ic),' [',map.key{ic},'], beh ',num2str(bestp),' [RW RT DT SP VE BT], E = ',num2str(allE(iic)),', w_{E/d}=',num2str(allwt_EdL(iic)),' (Temperate)'])
    end
    hold off
    set(gca,'fontname','times','fontsize',14)
    set(gcf,'PaperPosition',[0,0,10.5,8],'paperorientation','landscape');
    print('-dpdf',['plots/fintraj_ic',num2str(ic),'.pdf'])
    close all
end
%% plots of traj points for best fit behaviors - sorted by weight (least to most)
% load('mapdim20_updatethresh')
% load('allweightedbesthiker_edist.mat')
% load('allclosestpoints.mat')
% [Ewsort,IndEw] = sort(allwt_EdL);
% ind_ic_neq_fdt = map.goodsim88;
% col = [0.1 0.2 0.2; 0.4660 0.6740 0.1880; 0 0.7 0.7];
% % reps = 100;
% for iic = 1:ics
%     figure(iic)
%     ic = ind_ic_neq_fdt(IndEw(iic))
%     bestp = bestbeh{IndEw(iic)}(1,1);
%     load(['../newsims/sim_',LPtype,'_ic',num2str(ic),'_t100.mat']);
%     otherp = [2, 5, 10];
%     for iiprob = 1:3
%         iprob = bestbeh{IndEw(iic)}(otherp(iiprob),1);
%         for ii = 1:reps
%             cpx = pointdata{iprob,ii}(2,1);
%             cpy = pointdata{iprob,ii}(2,2);
%             pl(iiprob,ii) = scatter(cpx,cpy,50,'o','markerfacecolor',col(iiprob,:),'markeredgecolor',col(iiprob,:));
%             alpha(pl(iiprob,ii),0.4)
%             hold on
%         end
%         pc(iiprob) = plot(clpts{IndEw(iic)}(iprob,1),clpts{IndEw(iic)}(iprob,2),'p','markerfacecolor',col(iiprob,:),'markersize',12,'MarkerEdgeColor','k',...
%             'Displayname','Closest Rep');
%         hold on
%
%     end
%     pl(1).DisplayName = ['2nd ',num2str(bestbeh{IndEw(iic)}(2,1)),' ',mat2str(probs(bestbeh{IndEw(iic)}(2,1),:),2)];
%     pl(2).DisplayName = ['5th ',num2str(bestbeh{IndEw(iic)}(5,1)),' ',mat2str(probs(bestbeh{IndEw(iic)}(5,1),:),2)];
%     pl(3).DisplayName = ['10th ',num2str(bestbeh{IndEw(iic)}(10,1)),' ',mat2str(probs(bestbeh{IndEw(iic)}(10,1),:),2)];
%     for ii = 1:reps
%         xb = pointdata{bestp,ii}(2,1);
%         yb = pointdata{bestp,ii}(2,2);
%         pb = plot(xb,yb,'o','markerfacecolor',[0.8500 0.3250 0.0980],'markeredgecolor','k','markersize',7,...
%             'Displayname',['1st ',num2str(bestp),' ',mat2str(probs(bestp,:),2)]);
%         hold on
%     end
%     pc(4) = plot(clpts{IndEw(iic)}(bestp,1),clpts{IndEw(iic)}(bestp,2),'p','markerfacecolor',[0.8500 0.3250 0.0980],'markersize',12,'MarkerEdgeColor','k',...
%         'Displayname','1st best beh closest rep');
%     pc(1).DisplayName = '2nd best beh closest rep';
%     pc(2).DisplayName = '5th best beh closest rep';
%     pc(3).DisplayName = '10th best beh closest rep';
%     legend([pb,pl(1),pl(2),pl(3),pc],'location','northeastoutside')
%
%     findxy = map.findxy(ic,:);
%     plot(findxy(1),findxy(2),'p','markerfacecolor','y','markeredgecolor','k',...
%         'linewidth',1.5,'markersize',15,'displayname','find')
%     icxy = map.icsxy(ic,:);
%     plot(icxy(1),icxy(2),'gs','markerfacecolor','g','markeredgecolor','k','markersize',10,'displayname','ipp')
%     axis square, xlim([0 map.dim(1)]), ylim([0 map.dim(1)]), box on
%     if strcmp(map.mobdomain{ic},'Dry')
%         sgtitle(['closest pt - ic',num2str(ic),' [',map.key{ic},'], beh ',num2str(bestp),' [RW RT DT SP VE BT], E = ',num2str(allE(IndEw(iic))),', w_{E/d}=',num2str(allwt_EdL(IndEw(iic))),' (Dry)'])
%     else
%         sgtitle(['closest pt - ic',num2str(ic),' [',map.key{ic},'], beh ',num2str(bestp),' [RW RT DT SP VE BT], E = ',num2str(allE(IndEw(iic))),', w_{E/d}=',num2str(allwt_EdL(IndEw(iic))),' (Temperate)'])
%     end
%     hold off
%     set(gca,'fontname','times','fontsize',14)
%     set(gcf,'PaperPosition',[0,0,10.5,8],'paperorientation','landscape');
%     print('-dpdf',['plots/all200/sorted/wtraj_ic',num2str(iic),'.pdf'])
%     close all
% end

%% plots of traj points for best fit behaviors - sorted by weight (most to least)
close all
load('mapdim20_updatethresh')
load(['best fit data/allweightedbesthiker_edist',simname{isim},'.mat'])
load(['best fit data/allclosestpoints',simname{isim},'.mat'])
[Ewsort,IndEw] = sort(allwt_EdL);
IndEw = flip(IndEw);
ind_ic_neq_fdt = map.goodsim88;
col = [0.1 0.2 0.2; 0.4660 0.6740 0.1880; 0 0.7 0.7];
for iic = 1:ics
    figure(iic)
    ic = ind_ic_neq_fdt(IndEw(iic))
    bestp = bestbeh{IndEw(iic)}(1,1);
    load(['../sims',simname{isim},'/sim_',LPtype,'_ic',num2str(ic),'_t100.mat']);
    otherp = [2, 5, 10];
    for iiprob = 1:3
        iprob = bestbeh{IndEw(iic)}(otherp(iiprob),1);
        for ii = 1:reps
            cpx = pointdata{iprob,ii}(2,1);
            cpy = pointdata{iprob,ii}(2,2);
            pl(iiprob,ii) = scatter(cpx,cpy,50,'o','markerfacecolor',col(iiprob,:),'markeredgecolor',col(iiprob,:));
            alpha(pl(iiprob,ii),0.4)
            hold on
        end
        pc(iiprob) = plot(clpts{IndEw(iic)}(iprob,1),clpts{IndEw(iic)}(iprob,2),'p','markerfacecolor',col(iiprob,:),'markersize',12,'MarkerEdgeColor','k',...
            'Displayname','Closest Rep');
        hold on
        
    end
    pl(1).DisplayName = ['2nd ',num2str(bestbeh{IndEw(iic)}(2,1)),' ',mat2str(probs(bestbeh{IndEw(iic)}(2,1),:),2)];
    pl(2).DisplayName = ['5th ',num2str(bestbeh{IndEw(iic)}(5,1)),' ',mat2str(probs(bestbeh{IndEw(iic)}(5,1),:),2)];
    pl(3).DisplayName = ['10th ',num2str(bestbeh{IndEw(iic)}(10,1)),' ',mat2str(probs(bestbeh{IndEw(iic)}(10,1),:),2)];
    for ii = 1:reps
        xb = pointdata{bestp,ii}(2,1);
        yb = pointdata{bestp,ii}(2,2);
        pb = plot(xb,yb,'o','markerfacecolor',[0.8500 0.3250 0.0980],'markeredgecolor','k','markersize',7,...
            'Displayname',['1st ',num2str(bestp),' ',mat2str(probs(bestp,:),2)]);
        hold on
    end
    pc(4) = plot(clpts{IndEw(iic)}(bestp,1),clpts{IndEw(iic)}(bestp,2),'p','markerfacecolor',[0.8500 0.3250 0.0980],'markersize',12,'MarkerEdgeColor','k',...
        'Displayname','1st best beh closest rep');
    pc(1).DisplayName = '2nd best beh closest rep';
    pc(2).DisplayName = '5th best beh closest rep';
    pc(3).DisplayName = '10th best beh closest rep';
    legend([pb,pl(1),pl(2),pl(3),pc],'location','northeastoutside')
    
    findxy = map.findxy(ic,:);
    plot(findxy(1),findxy(2),'p','markerfacecolor','y','markeredgecolor','k',...
        'linewidth',1.5,'markersize',15,'displayname','find')
    icxy = map.icsxy(ic,:);
    plot(icxy(1),icxy(2),'gs','markerfacecolor','g','markeredgecolor','k','markersize',10,'displayname','ipp')
    axis square, xlim([0 map.dim(1)]), ylim([0 map.dim(1)]), box on
    if strcmp(map.mobdomain{ic},'Dry')
        sgtitle(['closest pt - ic',num2str(ic),' [',map.key{ic},'], beh ',num2str(bestp),' [RW RT DT SP VE BT], E = ',num2str(allE(IndEw(iic))),', w_{E/d}=',num2str(allwt_EdL(IndEw(iic))),' (Dry)'])
    else
        sgtitle(['closest pt - ic',num2str(ic),' [',map.key{ic},'], beh ',num2str(bestp),' [RW RT DT SP VE BT], E = ',num2str(allE(IndEw(iic))),', w_{E/d}=',num2str(allwt_EdL(IndEw(iic))),' (Temperate)'])
    end
    hold off
    set(gca,'fontname','times','fontsize',14)
    set(gcf,'PaperPosition',[0,0,10.5,8],'paperorientation','landscape');
    print('-dpdf',['plots/all',simname{isim},'/sorted/btraj_ic',num2str(iic),'.pdf'])
    close all
end

plotbehaviors

toc
