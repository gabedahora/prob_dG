import numpy as np
import scipy.integrate as spi

kb = 0.0019872041  # Boltzmann constant in kcal/(mol*K)

# Define the Boltzmann distribution
def boltzmann(E, T):
    A = 1.0  # Assume A=1 for simplicity; in practice, you'd need to calculate this
    return A * np.exp(-E / (kb * T))

# Load the data
data = np.genfromtxt('grid*.dat') # enter the grid file name

# Split the data into x, y, z
x = data[:, 0]
y = data[:, 1]
z = data[:, 2] * (-0.239006)  # Convert energy from kJ/mol to kcal/mol

# Define the ranges for each basin
basinA = ((2.6 <= x) & (x <= 3.2) & (-16 <= y) & (y <= -14))
basinB = ((2.6 <= x) & (x <= 3.3) & (-53 <= y) & (y <= -48))
basinC = ((2.6 <= x) & (x <= 3.2) & (-35 <= y) & (y <= -21))

# Get the energy values for each basin
E_basinA = z[basinA]
E_basinB = z[basinB]
E_basinC = z[basinC]

# Define the temperature
T = 300  # Adjust this value as needed

# Shift the energy scale so that the lowest energy is zero
#E_min = min(E_basinA.min(), E_basinB.min(), E_basinC.min())
#E_basinA -= E_min
#E_basinB -= E_min
#E_basinC -= E_min

# Calculate the total probability for each basin by integrating the Boltzmann distribution
total_prob_basinA = np.sum(boltzmann(E_basinA, T))
total_prob_basinB = np.sum(boltzmann(E_basinB, T))
total_prob_basinC = np.sum(boltzmann(E_basinC, T))

# Calculate the total probability over all energy values
total_prob_all = np.sum(boltzmann(z, T))

# Normalize the total probabilities
total_prob_basinA /= total_prob_all
total_prob_basinB /= total_prob_all
total_prob_basinC /= total_prob_all

# Calculate delta_G for each basin with respect to basin B
delta_G_A_B = -kb * T * np.log(total_prob_basinA / total_prob_basinB)
delta_G_C_B = -kb * T * np.log(total_prob_basinC / total_prob_basinB)

# Initialize lists to store the resampled and bias-corrected delta G values
delta_G_values_A_B = []
delta_G_values_C_B = []

delta_G_values_A_B_corrected = []
delta_G_values_C_B_corrected = []

# Define the number of iterations for the resampling
n_iterations = 1000

# Perform the resampling
for _ in range(n_iterations):
    # Resample the energy values for each basin
    E_basinA_resampled = np.random.choice(E_basinA, size=len(E_basinA), replace=True)
    E_basinB_resampled = np.random.choice(E_basinB, size=len(E_basinB), replace=True)
    E_basinC_resampled = np.random.choice(E_basinC, size=len(E_basinC), replace=True)

    # Calculate the total probability for each basin by integrating the Boltzmann distribution
    total_prob_basinA_resampled = np.sum(boltzmann(E_basinA_resampled, T))
    total_prob_basinB_resampled = np.sum(boltzmann(E_basinB_resampled, T))
    total_prob_basinC_resampled = np.sum(boltzmann(E_basinC_resampled, T))

    # Calculate delta_G for each basin with respect to basin B
    delta_G_A_B_resampled = -kb * T * np.log((total_prob_basinA_resampled + 1e-10) / (total_prob_basinB_resampled + 1e-10))
    delta_G_C_B_resampled = -kb * T * np.log((total_prob_basinC_resampled + 1e-10) / (total_prob_basinB_resampled + 1e-10))

    # Store the resampled delta G values
    delta_G_values_A_B.append(delta_G_A_B_resampled)
    delta_G_values_C_B.append(delta_G_C_B_resampled)

# Calculate the bias for each basin with respect to basin B
bias_A_B = np.mean(delta_G_values_A_B) - delta_G_A_B
bias_C_B = np.mean(delta_G_values_C_B) - delta_G_C_B

# Correct the resampled delta G values for bias
delta_G_values_A_B_corrected = delta_G_values_A_B - bias_A_B
delta_G_values_C_B_corrected = delta_G_values_C_B - bias_C_B

# Calculate the mean, standard deviation, and error of the bias-corrected values
delta_G_mean_corrected_A_B = np.mean(delta_G_values_A_B_corrected)
delta_G_std_corrected_A_B = np.std(delta_G_values_A_B_corrected)
delta_G_error_corrected_A_B = delta_G_std_corrected_A_B / np.sqrt(n_iterations)

delta_G_mean_corrected_C_B = np.mean(delta_G_values_C_B_corrected)
delta_G_std_corrected_C_B = np.std(delta_G_values_C_B_corrected)
delta_G_error_corrected_C_B = delta_G_std_corrected_C_B / np.sqrt(n_iterations)

# Write the results to a text file
with open('NOSHIFT_TOTAL_NORM_results.txt', 'w') as f:
    f.write(f"Total probability (population) of basin A: {total_prob_basinA}\n")
    f.write(f"Total probability (population) of basin B: {total_prob_basinB}\n")
    f.write(f"Total probability (population) of basin C: {total_prob_basinC}\n")
    f.write(f"Delta G from basin A to B: {delta_G_A_B} kcal/mol\n")
    f.write(f"Delta G from basin C to B: {delta_G_C_B} kcal/mol\n")
    f.write(f"Mean Delta G from bias-corrected bootstrapping (A to B): {delta_G_mean_corrected_A_B} kcal/mol\n")
    f.write(f"Standard Deviation of Delta G from bias-corrected bootstrapping (A to B): {delta_G_std_corrected_A_B} kcal/mol\n")
    f.write(f"Error in Delta G from bias-corrected bootstrapping (A to B): {delta_G_error_corrected_A_B} kcal/mol\n")
    f.write(f"Mean Delta G from bias-corrected bootstrapping (C to B): {delta_G_mean_corrected_C_B} kcal/mol\n")
    f.write(f"Standard Deviation of Delta G from bias-corrected bootstrapping (C to B): {delta_G_std_corrected_C_B} kcal/mol\n")
    f.write(f"Error in Delta G from bias-corrected bootstrapping (C to B): {delta_G_error_corrected_C_B} kcal/mol\n")
