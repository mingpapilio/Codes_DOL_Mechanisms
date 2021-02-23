This directory contains the Matlab scripts and data that were used to generate Figure4, Figure5 and FigureS5.

Figure4.m generates the 4 panels of Figure4. filamentgrow2.m must be in the same directory. 

simrun.m takes a given parameter combination given by input.csv and runs independent evolutionary simulations, storing the output in the data folder which must be in the same directory. We have done this to create .csv files for all of the parameter combinations considered, for 5 replicates each. Each .csv file is called "track_X_Y_Z.csv", where X is the value of phi, Y is the value of eta and Z is the index of the independent replicate. Filamentgrow.m must be in the same directory. 

Figure5AandB.m takes the output of simrun.m that are stored in the data folder and plots heat maps of the resulting evolved trait values (Figures 5A, 5B, S5A and S5B). This script outputs qmid.csv, smid.sv, dmid.csv and vmid.csv which are the average evolved trait values at each parameter combination across all independent replicated.

Figure5C takes the output of Figures5AandB.m and produces Figure5C. Filamentgrow3.m must be in the same directory. 

