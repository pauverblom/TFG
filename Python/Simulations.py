import numpy as np
import matplotlib.pyplot as plt
import os

# Define the values for corona, angle, and alpha
corona_values = ['1', '10', '50']
angle_values = ['0.00', '1.57', '3.14']
alpha_values = ['0.0', '0.1', '1.0']

# Base directory for the data files
base_dir = '/Users/pauverdeguer/TFG/MATLAB/'


# Iterate over the alpha values to create separate figures
fig, axs = plt.subplots(len(alpha_values), len(corona_values), figsize=(10, 10))
# Iterate over the angle and corona values to create subplots
for i, alpha in enumerate(alpha_values):
    for j, corona in enumerate(corona_values):
        # Construct the filename
    
        
        icecubefilepath = os.path.join(base_dir, 'ngc1068_spectrum_95.txt')
        fermifilepath = os.path.join(base_dir, 'FermiData.txt')
        # Load the data from the CSV file
        icecubedata = np.loadtxt(icecubefilepath, delimiter='\t', skiprows=1)
        fermidatapoints = np.loadtxt(fermifilepath, delimiter=' ', max_rows=5)
        fermidataupperlims = np.loadtxt(fermifilepath, delimiter=' ', skiprows=5)

        fermidatapointsx = fermidatapoints[:, 0]
        fermidatapointsy = fermidatapoints[:, 2] * 6.242e11 # Convert to eV from erg
        fermidatapointsyerr = fermidatapoints[:, 3] * 6.242e11 # Convert to eV from erg
        fermidataupperlimsx = fermidataupperlims[:, 0]
        fermidataupperlimsy = fermidataupperlims[:, 2] * 6.242e11 # Convert to eV from erg

        icecubex = icecubedata[:, 0] * 1e9 # Convert to eV
        icecubey = icecubedata[:, 1] * 1e12 # Convert to eV
        icecubelowerlimit = icecubedata[:, 2] * 1e12
        icecubeupperlimit = icecubedata[:, 3] * 1e12
        ax = axs[i, j]

        ax.errorbar(fermidatapointsx, fermidatapointsy, yerr=fermidatapointsyerr, marker='.', linestyle = ' ', color='k', label='Fermi')
        ax.errorbar(fermidataupperlimsx, fermidataupperlimsy, yerr=0.1*fermidataupperlimsy, uplims=True, linestyle=' ', color='k')
        ax.loglog(icecubex, icecubey, label='IceCube', color='b')
        ax.fill_between(icecubex, icecubelowerlimit, icecubeupperlimit, color='b', alpha=0.2)
        ax.loglog(icecubex, icecubelowerlimit, color='b', linestyle='-', alpha=0.5)
        ax.loglog(icecubex, icecubeupperlimit, color='b', linestyle='-', alpha=0.5)

        for angle in angle_values:
            datafilename = f'data/corona_{corona}_angle_{angle}_alpha_{alpha}.csv'
            datafilepath = os.path.join(base_dir, datafilename)
            data = np.loadtxt(datafilepath, delimiter=',', skiprows=1)
            data = data[data[:, 1] != 0]
            x = data[:, 0]
            y = data[:, 1]
            if angle == '0.00':
                ax.loglog(x, y, marker='o', linestyle=' ', markerfacecolor='none', color='darkorange', markersize=4, label=rf'$\theta=${angle}')  # Brighter orange
            elif angle == '1.57':
                ax.loglog(x, y, marker='s', linestyle=' ', markerfacecolor='none', color='crimson', markersize=4, label=rf'$\theta=${angle}')  # Vibrant red-orange
            else:
                ax.loglog(x, y, marker='^', linestyle=' ', markerfacecolor='none', color='royalblue', markersize=4, label=rf'$\theta=${angle}')  # Steel blue

        # Filter out rows where the y value (second column) is 0
        
            

        # Assuming the data has two columns: x and y
        

        # Create a log-log plot in the appropriate subplot
        
        
        ax.set_title(rf'$R_c/R_s=${corona}, $\alpha=${alpha}')
        ax.set_xlabel(r'$E$ [eV]')
        ax.set_ylabel(r'$E^2 dN/dE$ [eV cm$^{-2}$ s$^{-1}$]')
        
        # Adjust layout and set the figure size to square dimensions
# Adjust layout
fig.set_size_inches(8, 8)

plt.tight_layout()        
plt.subplots_adjust(bottom=0.11)  # Add space for the legend

# Create a single legend below all subplots
handles, labels = axs[0, 0].get_legend_handles_labels()
fig.legend(handles, ['IceCube Data', r'$\theta_{\max}=0$', r'$\theta_{\max}=\pi/2$', r'$\theta_{\max}=\pi$', 'Fermi Data'], loc='lower center', fontsize=10, ncol=len(labels), bbox_to_anchor=(0.5, 0))

plt.savefig('/Users/pauverdeguer/TFG/LaTeX/Figures/simulations_plot.pdf')
