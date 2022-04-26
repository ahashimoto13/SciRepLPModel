%%%%% December 19, 2020
%%%%% This runs the best fit behavior from ebestfitbeh_ics2 on the IC that
%%%%% was left out of that analysis. The results from this simulation will
%%%%% be used to find the t statistics for these best fits
clearvars; close all; clc
set(0,'DefaultFigureWindowStyle','normal');
tic

% load map data - ics, find, body coords...
load('mapdim20_updatethresh.mat')

% limits of map
ll = 1;
LLx = map.dim(2);
LLy = map.dim(1);
LL = [LLx, LLy, ll];

% initial points
ind_ic_neq_fdt = map.goodsim88;
ics = map.icsxy(ind_ic_neq_fdt,:);
finds = map.findxy(ind_ic_neq_fdt,:);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load behavior distribution
% icfall = [1, 6, 16, 27, 39, 45, 52, 60, 66, 72];
% icf = icfall(10);

icf = ind_ic_neq_fdt;
probs = zeros(1,6);

parameters; % load reps, N, ics_num, and lpt

% instead of mobility, set time step in hours
simT = 100;
T = ts*simT;

%% Run the simulation

for iic = 1:length(icf)
    iic
    icout = icf(iic);
    load(['fits500/bestbeh_icf',num2str(icout),'.mat'],'best_EdL');
    probs = best_EdL;
    
    pointdata = cell(size(probs,1),reps); % cell array where each entry is a T(variable) x 3 matrix (broken down for each IC now)     %%%%%%%%%%%%%%%%%%
    initial_point = ics(iic,:); % initial starting point
    find_point = finds(iic,:);
    
    map_data = load(['../KoesterHikerMaps20/',map.loadfilenames{icout}]);
    map_data = parallel.pool.Constant(map_data);
    
    for iprob = 1:size(probs,1)
        [iic iprob]
        
        % define the behavior profile and length of simulation
        p_behavior = probs(iprob,:);
        
        parfor irep = 1:reps
            
            [closestpts,endpts,time_closest] = run_replicate(initial_point, find_point, map_data.Value, T, p_behavior, alpha, LL);
            pointdata{iprob,irep} = [initial_point;closestpts;endpts;time_closest];
            
        end
        save(['simsfinal/sim_hiker_ic',num2str(icout),'_t',num2str(simT),'.mat'],'pointdata','probs','-v7.3')
        
    end
    save(['simsfinal/sim_hiker_ic',num2str(icout),'_t',num2str(simT),'.mat'],'pointdata','probs','-v7.3')
    
end



toc

leave1outfinal3

%%
% close all; clearvars
% % set(0,'DefaultFigureWindowStyle','docked');
%
% load('mapdim10')
% load('leave1outresults.mat')
%
% for iicf = 1%:length(val.icleftout)
%     icf = val.icleftout(iicf);
%     figure(iicf)
%     load(['sim_hiker_icf',num2str(icf),'.mat'])
%     for iprob = 1
%         for ii = 1:100
%             x = all_data_ic{iprob,ii}(:,1);
%             y = all_data_ic{iprob,ii}(:,2);
%             plot(x(end),y(end),'bo','markerfacecolor','b','markersize',10),
%             hold on
%             plot(x,y)
%             %             for jj = 1:length(x)
%             %                 plot(x(jj),y(jj),'m.','markerfacecolor','m','markersize',5), hold on
%             %             end
%         end
%
%     end
%     findxy = map.findxy(icf,:);
%     plot(findxy(1),findxy(2),'kp','markerfacecolor','k','markersize',10)
%     icxy = map.icsxy(icf,:);
%     plot(icxy(1),icxy(2),'go','markerfacecolor','g','markersize',10)
%     axis square
%     title(['hiker - ic left out ',num2str(icf)])
%     hold off
%     set(gca,'fontname','times','fontsize',14)
% %     set(gcf,'PaperPosition',[0,0,8,8],'paperorientation','portrait');
% %     print('-dpdf',['plots/traj_icf',num2str(icf),'.pdf'])
% end