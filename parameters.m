% behavior combos
N = 1;
% replicates
reps = 500;

% choose the LPTs to run {'hiker';'child10to12';'child1to6';'dementia';'hunter';'angler';'child13to15';'despond';'snowboard';'worker'};
lpt = [1];  % hiker, hunter, angler, snowboard

% number of lost person types
lpts = numel(lpt);

% smoothing parameter
alpha = 0.55;

% walking speed - ts time steps per hour
ts = 850;   % speed(6.67m/6s) * conversion(3600s/6.67m) = ts/hr

% name list
savename = {'hiker'};
% full list below
% savename = {'hiker','hunter','angler','snowboard','worker','child10to12','child1to6','dementia','despond','child13to15'};
