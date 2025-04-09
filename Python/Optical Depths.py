import numpy as np
import matplotlib.pyplot as plt

# Load the file and store its contents as columns
file_path = '/Users/pauverdeguer/TFG/MATLAB/optical_depths_as_function_of_Rc.txt'
data = np.loadtxt(file_path, delimiter=',')

# Access columns
energies = data[:, 0]
Rc_1 = data[:, 1]
Rc_10 = data[:, 2]
#Rc_50 = data[:, 3]
Rc_100 = data[:, 4]

fig, ax = plt.subplots(figsize=(10, 6))

# Set axis limits
ax.set_xlim(1e2, 1e14)
ax.set_ylim(5e-14, 5e8)

# Plot the data with new colors
ax.loglog(energies, Rc_1, label='Rc = Rs', color='darkviolet', linewidth=3)
ax.loglog(energies, Rc_10, label='Rc = 10 Rs', color='mediumpurple', linewidth=3)
#ax.loglog(energies, Rc_50, label='Rc = 50 Rs', color='navy', linewidth=3)
ax.loglog(energies, Rc_100, label='Rc = 100 Rs', color='darkslategray', linewidth=3)


# Add horizontal line at y = 1
ax.axhline(y=1, color='crimson', linestyle='-', linewidth=2, label=r'$\tau_{\gamma\gamma}=1$')


# Add labels and title
ax.set_xlabel(r'$E_{\gamma}$ [eV]', fontsize=18)
ax.set_ylabel(r'$\tau_{\gamma\gamma}$', fontsize=18)

# Change the size of the numbers in the axes
ax.tick_params(axis='both', which='major', labelsize=18)

# Add grid
ax.grid(True, which='both', linestyle='-', linewidth=0.5)

ax.legend(fontsize=18, frameon=True)

# Adjust layout
plt.tight_layout()

# Show the plot
plt.savefig('/Users/pauverdeguer/TFG/LaTeX/Figures/OpticalDepths.pdf')

plt.show()
