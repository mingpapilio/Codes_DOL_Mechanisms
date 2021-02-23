# Codes_DOL_Mechanisms
This repository is for open access to source codes and data for the manuscript **The evolution of mechanisms to divide labour in microorganisms**

In ***analytical_and_numerical_scripts*** are the files for generating Fig.3 and Fig.S1-3. The enclosed scripts do not label axes or results. For this information please refer to the figures as they occur in the main text and supplementary information. Please include WK.m in the same folder in order to run: FiguresS1BandC.m, FigureS2C.m and FigureS2D.m.

In ***cyanobacteria_simulations*** are the files for generating Fig.4-5 and Fig.S5.

In ***simulation_codes*** are the source codes of Fig.S4. Because the random number generators files (dSFMT) are required for running, please **do not** change the relative location of these files. **PLEASE READ** the notes at the beginning of each file, containing the key parameters and notes.

In ***simulation_data*** are the raw data of Fig.S4. All levelplot data are saved as *level_sumary.csv*, you may use the files in *R_files* to replot the figures. If you wish to generate your own *level_summary.csv*, please assign the folders as in *simulation_data/Fig4_levelplot* and run *cord_ss.c* for each group size. The *merge_summaries.R* would help you create the *level_summary.csv* file. 

Please email *ming.liu@zoo.ox.ac.uk* if there is any problem, thanks! (Ming Liu)
