import numpy as np
import pandas as pd
import joypy
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter

# Set style
sns.set_theme(style="whitegrid")

# Parameters
alphas = [0.01, 0.1, 1, 5, 10]
num_samples = 10000000
colors = plt.cm.plasma(np.linspace(0, 1, len(alphas)))
# Generate data
data = pd.DataFrame()
for alpha in alphas:
    u = np.random.rand(num_samples)
    distance = u ** (1 / alpha)
    temp_df = pd.DataFrame({'distance': distance, 'alpha': f"α = {alpha}"})
    data = pd.concat([data, temp_df])

# Ensure the data is sorted by 'alpha' before plotting
data['alpha'] = pd.Categorical(data['alpha'], categories=[f"α = {alpha}" for alpha in alphas], ordered=True)
data = data.sort_values(by='alpha')

# Create ridge plot with percentage y-axis and reduced overlap
fig, axes = joypy.joyplot(
    data,
    by='alpha',
    column='distance',
    bins=60,
    hist=True,
    linecolor='black',
    linewidth=0.5,
    fade=False,
    figsize=(10, 6),
    x_range=[0, 1],
    grid=True,  # Disable the default grid
    xlabelsize=12,
    ylabelsize=12,
    alpha=0.8,
)

# Custom function to format y-ticks as percentages
def percentage_formatter(x, pos):
    return f"{x:.0f}%" if x >= 1 else f"{x:.1f}%"

# Convert y-axis to percentage scale and normalize histogram data
for ax, color in zip(axes, colors):
    ax.set_yscale('log')

    ax.set_ylim(0.011, 100)  # Avoid log(0) and set the lower limit to 0.1%

    # Apply the custom formatter for the y-axis
    ax.yaxis.set_major_formatter(FuncFormatter(percentage_formatter))
    
    # Set y-ticks explicitly to include 0.1, 1, 10, and 100
    ax.set_yticks([0.1, 1, 10, 100])
    
    # Make y-ticks bigger
    ax.tick_params(axis='y', labelsize=14)
    #ax.y
    # Get the current y-values (counts) and x-values (bin edges)
    for patch in ax.patches:
        patch.set_facecolor(color)
    
    # Normalize the y-values to add up to 100%
    heights = [patch.get_height() for patch in ax.patches]
    total_height = sum(heights)
    if total_height != 0:  # Prevent division by zero
        normalized_heights = [height / total_height * 100 for height in heights]
    else:
        normalized_heights = [0 for _ in heights]  # Handle zero height case
    
    # Apply the normalized heights back to the patches
    for patch, normalized_height in zip(ax.patches, normalized_heights):
        patch.set_height(normalized_height)

    # Add a less intrusive grid
    ax.set_xticks([0, 0.5, 1])  # Explicit x-ticks
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.7)

# Customize plot
plt.xlabel(r'$r\,/\,R_c$', fontsize=18)
plt.xticks([0, 0.5, 1], ['0',r'$R_c/2$', r'$R_c$'], fontsize=14)
handles = [plt.Rectangle((0, 0), 0.1, 0.1, color=c) for c in colors]
labels = [f"α = {alpha}" for alpha in alphas]
plt.legend(handles, labels, loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=12)
# Add a manual text label for the y-axis
plt.tight_layout()  # Adjust padding to prevent clipping
plt.subplots_adjust(left=0.1, right=0.85)  # Adjust left and right margins
plt.text(-0.11, 0.5, 'Probability Density (%)', fontsize=18, rotation='vertical', va='center', ha='center')
# Move legend to the side



# Adjust layout to prevent clipping of legend
plt.savefig('/Users/pauverdeguer/TFG/LaTeX/Figures/InitialDistributions.pdf')

# Show plot
plt.show()

