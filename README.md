# prob_dG
This repository has a Python script to calculate the total probability and deltaG between basins in the PMF

The script will read the grid file from PLUMED and calculate the total probability for each basin by integrating the Boltzmann distribution.

Then it will calculate the total probability over all energy values, normalize the total probabilities for each basin, and calculate the delta_G for each basin with respect to a basin (in the script, with respect to basin B).

Finally, it will calculate the mean, standard deviation, and error of the bias-corrected values, printing all to a txt file. 
