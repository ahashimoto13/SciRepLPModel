# Guide to running scripts for the Scientific Reports Lost Person Model 
 [2022ScientificReportsPaper.pdf](https://github.com/ahashimoto13/SciRepLPModel/files/8558714/2022ScientificReportsPaper.pdf)



## Running the simulations
### Creating all the maps for the ICs (Python)
Note: Use [AGSTools](https://git.caslab.ece.vt.edu/hlarkin3/ags_grabber) for updated map scripts
1.	Open up Pycharm and navigate to the `importmap` project folder
2.	Open `feature_set.py` and enter in all the desired ICs
3.	Set the `extent` of the maps
4.	Run `feature_set.py` and the output will be `BW_LFandInac_Zelev_[IC lat, IC lon].mat`

### Run all possible behavior probs for 65 ICs (MATLAB)
1.	Run `main_hiker_t100.m` (uses `parameters.m`, `randsmpl.m`, `run_replicate.m`, `beh_dist_6.mat`, `mapdim20_updatethresh.mat`)
    -	Open `parameters.m` and adjust number of `probs` and `replicates`
    -	Adjust the simulation time variable `simT` if necessary
2.	Simulations save to a folder (`\sims500t`)

 
### Running the analysis on the simulations
***PART 1: Energy Distance*** 

**Goal**: finds the energy statistic for the distance between the closest sim point and the find point from ISRID data for all 462 behaviors and reps
1.	Navigate to folder `\analysis energy dist`
2.	Run `bestbeh_edist_clpts.m` (uses `beh_dist_6.mat`, `mapdim20_updatethresh.mat`)
    -	Calculates the best fits for all the simulations based on the energy dist
    -	Finds the best beh profile for all ICs using the weight, time of the closest reps
    -	Plots the trajectory points of the best fit behaviors
    -	Plots the final trajectory points of the best fit behaviors
    -	Plots these points sorted by weight
    -	Saves the best fit data in `\analysis energy dist\best fit data`
    -	**Output**: `allbesthiker_edist500t.mat`, `allclosestpoints500t.mat`, `allweightedbesthiker_edist500t.mat`
3.	Run `plotbehaviors.m`
    -	Plots bar charts of behavior probability dists with weights, energy distances, and distances from the IPP to find
    -	Plots the sensitivity analysis for L 
    -	Sorts the ICs by weights and behaviors
    -	Plots the weights versus distances
    -	Used to create some of the figures (those have been extricated to new scripts)

***PART 2: Leave-one-out cross validation***  

**Goal**: Based on the distances from the actual find point, this script finds the best fitting behaviors based on the E-dist. It takes out one IC, calculates the average behavior using the weights for the rest of the ICs and saves this behavior. This behavior will be used to run the simulation for the IC that was left out.
1.	Navigate to folder `\analysis LOOCV`
2.	Run `bestbeh_LOO1.m` (uses `beh_dist_6.mat`, `mapdim20_updatethresh.mat`)
    -	Uses the simulation data to run a LOOCV and saves the best fits for each IC
    -	Saves the best fits (`best_EdL`) in folder `\analysis LOOCV\fits500\bestbeh_icf#`
3.	Run `main_LOO2.m` (uses `parameters.m`, `randsmpl.m`, `run_replicate.m`, `beh_dist_6.mat`, `mapdim20_updatethresh.mat`)
    -	Runs the best fit behavior from `bestbeh_LOO1.m` on the IC that was left out of that analysis. The results from this simulation will used to find the percentiles
    -	Can change the number of simulation hours (`100`) and reps (`500`)
    -	Saves the output of the runs to folder `\analysis LOOCV\simsfinal`
4.	Run `leave1outfinal3.m`
    -	Using the output of the simulations run with the probabilities found from leaving one out, this script finds the best behavior by energy distance for each IC
    -	By comparing the left-out energy distances to the original energy distances, we can find the percentiles 
    - Plots box plots, bar charts, and a histogram
