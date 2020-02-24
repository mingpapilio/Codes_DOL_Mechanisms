# Codes_DOL_Mechanisms
This repository is for open access to source codes and data for the manuscript **Dividing labour in social organisms: coordinated or random specialisation?**.

In *simulation_codes* folder, because the random number generators files (dSFMT) are required for running, please **DO NOT** change the relative location of these files. **PLEASE READ** the notes at the beginning of each file, containing the key parameters and notes.

In *simulation_data* folder, all levelplot data are saved as *level_sumary.csv*, you may use the files in *R_files* to replot the figures. 

If you wish to generate your own *level_summary.csv*, please assign the folders as in *simulation_data/Fig4_levelplot* and run *cord_ss.c* for each group size. The *merge_summaries.R* would help you create the *level_summary.csv* file. 

Please email *ming.liu@zoo.ox.ac.uk* if there is any problems, thanks! (Ming Liu)
