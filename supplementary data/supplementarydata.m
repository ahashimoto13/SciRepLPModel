%%%% Supplementary Data for Scientific Reports Paper May 2021
% Create csv files of model simulation data:
% 1) SimulationResults: 65 ICs with closest point lat/lon and time for best fit behavior
% 2) BehavioralProfiles: Best fit behavior distribution for all 65 ICS
% 3) IPPandFindLocation: 65 ICs with IPP and find location
clearvars; close all; clc

% load map indices
load('../mapdim20_updatethresh.mat')
ind_ic_neq_fdt = map.goodsim88;
ics = length(ind_ic_neq_fdt);
reps = 500;

% load best fitting behaviors
load('../analysis energy dist/best fit data/allweightedbesthiker_edist500t.mat','bestbeh','allbp')


%% 1) 65 ICs with closest point lat/lon and time for best fit behavior
k = 0;  % row index for nested loop
closestpoints = zeros(ics*reps,5);
for iic = 1:ics
    ic = ind_ic_neq_fdt(iic);
    load(['../sims500t/sim_hiker_ic',num2str(ic),'_t100.mat'],'pointdata');
    bestp = bestbeh{iic}(1,1);              % best behavior by edist on closest pts
    
    % convert to lon (x) and lat (y)
    ycrds = linspace(map.latlim(ic,1),map.latlim(ic,2),map.dim(1));   %% switch 1 and 2 to flip y-axis
    xcrds = linspace(map.lonlim(ic,1),map.lonlim(ic,2),map.dim(1));
    for irep = 1:reps
        k = k+1;
        cpx = pointdata{bestp,irep}(2,1);
        cpy = pointdata{bestp,irep}(2,2);
        cplon = xcrds(cpx);
        cplat = ycrds(cpy);
        cpt = pointdata{bestp,irep}(4,1)/850;
        closestpoints(k,:) = [iic,irep,cplat,cplon,cpt];
    end
end
T1 = array2table(closestpoints);
T1.Properties.VariableNames(1:5) = {'incident_index','rep','closestpt_lat','closestpt_lon','closestpt_time_hr'};

%% 2) Best fit behavior distribution for all 65 ICS
behaviors = [(1:65)',allbp];
T2 = array2table(behaviors);
T2.Properties.VariableNames(1:7) = {'incident_index','RW','RT','DT','SP','VE','BT'};

%% 3) 65 ICs with IPP and find location
mapinfo = [(1:65)',map.ics(ind_ic_neq_fdt,:),map.find(ind_ic_neq_fdt,:)];
T3 = array2table(mapinfo);
T3.Properties.VariableNames(1:5) = {'incident_index','IPP_lat','IPP_lon','find_lat','find_lon'};

%% Save as csv files
writetable(T1,'SimulationResults.csv')
writetable(T2,'BehavioralProfiles.csv')
writetable(T3,'IPPandFindLocation.csv')