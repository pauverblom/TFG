import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker

# Define the file path
file_path = '/Users/pauverdeguer/TFG/MATLAB/storedvariables.csv'

# Initialize lists to store the data
energy = []
interaction_length = []
x_position = []
y_position = []
z_position = []

# Read the CSV file
with open(file_path, 'r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header row
    for row in reader:
        energy.append(float(row[0]))
        interaction_length.append(float(row[1]))
        x_position.append(float(row[2]))
        y_position.append(float(row[3]))
        z_position.append(float(row[4]))



# Convert lists to numpy arrays
energy = np.array(energy)
interaction_length = np.array(interaction_length)
x_position = np.array(x_position)
y_position = np.array(y_position)
z_position = np.array(z_position)

# Create figure with constrained layout
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(2, 2, height_ratios=[3, 2])

# Define colormap
cmap = cm.get_cmap('spring')

# Normalize by step index
norm = mcolors.Normalize(vmin=0, vmax=len(x_position) - 1)

# Generate colors based on step indices
colors = [tuple(1 - np.array(cmap(norm(i))[:3])) for i in range(len(x_position))]

# 3D Plot - Larger and centered at the top
ax3 = fig.add_subplot(gs[0, :], projection='3d')

for i in range(len(x_position)-2):
    ax3.plot(
        [x_position[i], x_position[i + 1]],
        [y_position[i], y_position[i + 1]],
        [z_position[i], z_position[i + 1]],
        color=colors[i],
    )

# Draw a sphere of radius R = 20 * 1.4e-6
R = 20 * 1.4e-6
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x_sphere = R * np.outer(np.cos(u), np.sin(v))
y_sphere = R * np.outer(np.sin(u), np.sin(v))
z_sphere = R * np.outer(np.ones(np.size(u)), np.cos(v))

ax3.plot_surface(x_sphere, y_sphere, z_sphere, color='gray', edgecolor='gray', alpha=0.2, linewidth=0.1)

# Add text for corona
ax3.text(0, 0, -R * 1.3, 'Edge of Corona', color='gray', fontsize=10, ha='center',
         bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.3'))

# Adjust orientation
ax3.view_init(elev=23, azim=46, roll=0)

# Zoom in on the corona region
zoom_factor = 0.9  # Increase to zoom more
ax3.set_xlim(-R * zoom_factor, R * zoom_factor)
ax3.set_ylim(-R * zoom_factor, R * zoom_factor)
ax3.set_zlim(-R * zoom_factor, R * zoom_factor)

# Add an arrow for the last position (outside graph)
arrow_length = R * 1.2  # Length of the arrow
direction = np.array([x_position[-2] - x_position[-2], 
                      y_position[-1] - y_position[-2], 
                      z_position[-1] - z_position[-2]])
direction = direction / np.linalg.norm(direction)  # Normalize the direction vector
ax3.quiver(x_position[-2], y_position[-2], z_position[-2],  
           direction[0] * arrow_length, 
           direction[1] * arrow_length,
           direction[2] * arrow_length, 
           color='blue', alpha=0.8, arrow_length_ratio=0.2, linewidth=1.5, edgecolor='blue')

        # Disable the scientific notation offset text for all axes
        # and force scientific notation for each tick label


# Ensure equal scaling
ax3.set_box_aspect([1, 1, 1])
ax3.set_title('3D Position of the Particle')

# Define the formatter function to multiply by 1000
def format_scaled(value, tick_number):
    # Multiply the value by 1000 and format it
    # Using scientific notation with 1 decimal place for consistency
    return f'{value * 1000000:.0f}'

# Apply the custom formatter to each axis
ax3.xaxis.set_major_formatter(mticker.FuncFormatter(format_scaled))
ax3.yaxis.set_major_formatter(mticker.FuncFormatter(format_scaled))
ax3.zaxis.set_major_formatter(mticker.FuncFormatter(format_scaled))

# Update axis labels to reflect the new scale (assuming original unit was mpc, now pc)
ax3.set_xlabel(r'x position [$\mu$pc]')
ax3.set_ylabel(r'y position [$\mu$pc]')
ax3.set_zlabel(r'z position [$\mu$pc]')


# Energy Plot - Bottom left
ax1 = fig.add_subplot(gs[1, 0])
for i in range(len(energy) - 1):
    ax1.plot([i, i + 1], [energy[i], energy[i + 1]], color=colors[i])
ax1.set_yscale('log')
ax1.set_title('Energy vs Step')
ax1.set_xlabel('Step #')
ax1.set_ylabel('Energy [eV]')
ax1.set_aspect(1.0 / ax1.get_data_ratio(), adjustable='box')
ax1.grid(True, which='major', linestyle='--', linewidth=0.5)

# Interaction Length Plot - Bottom right
ax2 = fig.add_subplot(gs[1, 1])
for i in range(len(interaction_length) - 1):
    ax2.plot([i, i + 1], [interaction_length[i], interaction_length[i + 1]], color=colors[i])
ax2.axhline(y=20 * 1.4e-6, color='gray', linestyle='dotted')
ax2.text(len(interaction_length) * 0.5, 21 * 1.4e-6, 'Corona Radius', color='gray', fontsize=10, ha='center', va='bottom')
ax2.set_yscale('log')
ax2.set_title('Interaction Length vs Step')
ax2.set_xlabel('Step #')
ax2.set_ylabel('Interaction Length [pc]')
ax2.set_aspect(1.0 / ax2.get_data_ratio(), adjustable='box')
ax2.grid(True, which='major', linestyle='--', linewidth=0.5)

# Adjust layout and show
plt.tight_layout()
plt.savefig('LaTeX/Figures/example_random_walk.pdf', dpi=300)
plt.show()