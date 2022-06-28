import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import fnmatch

# Take user input
# Ask if the data should be normalise by the peak total emission and absorption
normalise = input('Would you like to normalise the data? (y/n) ')
# Ask for the file extension of the output plots
file_ext = input("Which file extension? (Include '.') ")

# Define the search directory where the spectra files are housed
searchdir = os.listdir('emission_absorption/fixed_ion_list/list/')
# Define a search pattern which finds all relevant files within the search directory
pattern = 'plotspecemission_*.txt'

# Initialise a files array which holds all file names in the search directory, matching the search pattern
files = []
# Loop through each file in the search directory and add all that match the search pattern to the files array
for entry in searchdir:
    if fnmatch.fnmatch(entry, pattern):
        files.append(entry)

# Loop through each file
for file in files:
    # Print the file name currently being worked on
    print(file)
    # Second index selects species
    # Every species has an emission and absorption spectrum
    # Define the header, which contains the column labels
    header = np.loadtxt('emission_absorption/fixed_ion_list/list/' + file, max_rows=1, delimiter=',', dtype=str)
    # Each spectra is then the numerical values of each column
    spectra = np.loadtxt('emission_absorption/fixed_ion_list/list/' + file, skiprows=1, delimiter=',')

    # Initialise empty dictionaries, which will contain the emission and absorption spectra of each species
    emission = {}
    absorption = {}
    # Initialise an empty array that will hold the keys for the dictionaries
    e_keys = []
    a_keys = []
    # These counters keep track of local indices for the dictionaries. The global loop variable is for the whole spectrum array, but the emission and absorption indices are different
    e_counter = 0
    a_counter = 0
    # Create an array of zeroes equal to the length of the spectra. This will be used to determine the total spectrum across all species
    total_e = np.zeros(1000)
    total_a = np.zeros(1000)
    # Now loop through each item (species) in the header
    for item in range(1,len(header)):
        # If the first letter of the header item is 'e', then the column refers to an emission spectra
        if header[item][0] == 'e':
            # We append the header name to keys, removing everything expect the species name i.e. 'Fe I'
            e_keys.append(header[item][17:])
            # Next append the spectral data to the relevant dictionary
            # The 1e13 is a scaling factor only
            emission[e_counter] = spectra[:,item]*1e13
            # Reassign the dictionary key with the string associated with the species
            emission[e_keys[-1]] = emission.pop(e_counter)
            # Increase the local emission index variable
            e_counter += 1
            # Add to the total emission spectra
            total_e += emission[e_keys[-1]]
        # Run the same process as above for the absorption spectrum
        elif header[item][0] == 'a':
            a_keys.append(header[item][19:])
            absorption[a_counter] = -spectra[:,item]*1.e13
            absorption[a_keys[-1]] = absorption.pop(a_counter)
            a_counter += 1
            total_a += absorption[a_keys[-1]]
        else:
            None

    # If the answer at the beggining of the script is yes the normalisation, run the following block of code. This will normalise the spectra based on the maximum flux for the total emission and absorption.
    if normalise == 'y':
        # Loop through both the emission and absorption spectrum for each species
        for (key_e, key_a) in zip(e_keys, a_keys):
            # Check that we don't divide by zero
            if np.max(total_e) != 0 and np.min(total_a) != 0:
                # Normalise the data
                emission[key_e] = emission[key_e]/np.max(total_e)
                absorption[key_a] = -(absorption[key_a]/np.min(total_a))

    # Normalise the total emission/absorption spectrum
    if np.max(total_e) != 0 and np.min(total_a) != 0:
        total_e = total_e/np.max(total_e)
        total_a = -(total_a/np.min(total_a))

    # Begin the plotting
    fig, ax = plt.subplots()
    # Plot the total emission and absorption spectra
    ax.plot(spectra[:,0], total_e, label = None, color='k', linewidth=1)
    ax.plot(spectra[:,0], total_a, label = None, color='k', linewidth=1)

    colour = list(plt.get_cmap('tab20')(np.linspace(0, 1.0, 20)))

    # Create a stackplot of the emission and absorption contribution from each species
    # ax.stackplot(spectra[:,0], emission.values(), labels = e_keys, colors = ['tab:blue', 'tab:red', 'tab:green', 'purple', 'tab:orange', 'tab:pink', 'tab:gray', 'gold', 'tab:cyan', 'darkblue', 'darkgreen', 'maroon', 'mediumvioletred', 'saddlebrown', 'darkslategrey', 'bisque', 'yellow'])
    # ax.stackplot(spectra[:,0], absorption.values(), colors = ['tab:blue', 'tab:red', 'tab:green', 'purple', 'tab:orange', 'tab:pink', 'tab:gray', 'gold', 'tab:cyan', 'darkblue', 'darkgreen', 'maroon', 'mediumvioletred', 'saddlebrown', 'darkslategrey', 'bisque', 'yellow'])
    ax.stackplot(spectra[:,0], emission.values(), labels = e_keys, colors = colour)
    ax.stackplot(spectra[:,0], absorption.values(), colors = colour)
    # Make the plot look nice
    plt.xlim(2500, 11000)
    plt.legend(ncol=3, loc=4)
    plt.xlabel(r'Wavelength [$\AA$]')
    plt.ylabel(r'F$_{\lambda}$ at 1 Mpc [$10^{-13}$ erg/s/cm$^2$/$\AA$]')
    plt.title('Emission/Absorption Spectrum ' + file[17:-4])
    if normalise == 'y':
        plt.ylabel('Normalised Flux')
        plt.title('Emission/Absorption Spectrum (Normalised) ' + file[17:-4])

    # print(file[17:-4])
    # input()

    if normalise == 'y':
        if file_ext == '.pdf':
            save_file = 'emission_absorption_' + file[17:-4] + '_normalised' + file_ext
        else:
            save_file = 'images/emission_absorption_' + file[17:-4] + '_normalised' + file_ext
    else:
        save_file = 'emission_absorption_' + file[17:-4] + file_ext

    plt.savefig('emission_absorption/fixed_ion_list/list/' + save_file)

    # Clear and close the figure in preparation for the next timestep
    fig.clf()
    plt.close()

print('Done!')
