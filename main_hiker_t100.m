%%%%% February 18, 2021 - Runs LP model without mobility for hiker Koester ICs
clearvars; close all; clc
set(0,'DefaultFigureWindowStyle','normal');
tic

% load map data - ics, find points
load('mapdim20_updatethresh')

% limits of map
ll = 1;
LLx = map.dim(2);
LLy = map.dim(1);
LL = [LLx, LLy, ll];

% initial points
ind_ic_neq_fdt = map.goodsim88;
ics = map.icsxy(ind_ic_neq_fdt,:);
finds = map.findxy(ind_ic_neq_fdt,:);

%% load behavior distribution and parameters
load('beh_dist_6.mat');

parameters; % load reps, N, ics_num, and lpt

% instead of mobility, set time step in hours
simT = 100;
T = ts*simT;

%% Run the simulation

for iic = 1:size(ics,1)
    ic = ind_ic_neq_fdt(iic);
    pointdata = cell(size(probs,1),reps); % cell array where each entry is a T(variable) x 3 matrix (broken down for each IC now)     %%%%%%%%%%%%%%%%%%
    initial_point = ics(iic,:); % initial starting point
    find_point = finds(iic,:);
%     mypool = parpool ('local', 4);
    map_data = load(['KoesterHikerMaps20/',map.loadfilenames{ic}]);
    map_data = parallel.pool.Constant(map_data);
    
    for iprob = 1:size(probs,1)
        [iic iprob]
        
        % define the behavior profile and length of simulation
        p_behavior = probs(iprob,:);
        
        parfor irep = 1:reps
            %             [x, y, behavior] = run_replicate(initial_point, find_point, map_data, T, p_behavior, alpha, LL);
            %             all_data_ic{iprob,irep} = [x, y, behavior,];
            
            [closestpts,endpts,time_closest] = run_replicate(initial_point, find_point, map_data.Value, T, p_behavior, alpha, LL);
            pointdata{iprob,irep} = [initial_point;closestpts;endpts;time_closest];
            
        end
        save(['sims500t/sim_hiker_ic',num2str(ic),'_t',num2str(simT),'.mat'],'pointdata','-v7.3')
        
    end
    save(['sims500t/sim_hiker_ic',num2str(ic),'_t',num2str(simT),'.mat'],'pointdata','-v7.3')
%     clear pointdata map_data
%     delete(gcp);
end


toc

