import numpy as np
import matplotlib.pyplot as plt
import os
import fnmatch

# Define the search directory where the spectra files are housed
searchdir = os.listdir('emission_absorption/fixed_ion_list/list1/')
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
    header = np.loadtxt('emission_absorption/fixed_ion_list/list1/' + file, max_rows=1, delimiter=',', dtype=str)
    # Each spectra is then the numerical values of each column
    spectra = np.loadtxt('emission_absorption/fixed_ion_list/list1/' + file, skiprows=1, delimiter=',')

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
    ax.stackplot(spectra[:,0], emission.values(), labels = e_keys, colors = colour)
    ax.stackplot(spectra[:,0], absorption.values(), colors = colour)
    # Make the plot look nice
    plt.xlim(2500, 11000)
    plt.legend(ncol=3, loc=4)
    plt.xlabel(r'Wavelength [$\AA$]')
    plt.ylabel(r'F$_{\lambda}$ at 1 Mpc [$10^{-13}$ erg/s/cm$^2$/$\AA$]')
    plt.title('Emission/Absorption Spectrum ' + file[17:-4])

    plt.savefig('emission_absorption_' + file[17:-4] + '_normalised.pdf')
    plt.savefig('emission_absorption_' + file[17:-4] + '_normalised.png')


    # Clear and close the figure in preparation for the next timestep
    fig.clf()
    plt.close()

print('Done!')
