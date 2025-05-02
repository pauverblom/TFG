import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches

# Load the CSV
data = pd.read_csv('/Users/pauverdeguer/TFG/Python/parameter_combinations_sorted.csv')

# Keep only the top N lowest MSE combinations
top_data = data[data['mse']<=0.3]

# Optionally round parameters to reduce noise
alpha_step = 1/3
corona_step = 2

top_data['alpha_binned'] = np.round(top_data['alpha'] / alpha_step) * alpha_step
top_data['corona_binned'] = np.round((top_data['corona'] - 1) / corona_step) * corona_step + 1

# Count how often each (alpha, corona) pair appears
freq = top_data.groupby(['alpha_binned', 'corona_binned']).size().reset_index(name='count')

# Pivot to create a 2D matrix for heatmap
heatmap_data = freq.pivot(index='corona_binned', columns='alpha_binned', values='count').fillna(0)

# Plot the heatmap
plt.figure(figsize=(8, 8))
ax = sns.heatmap(heatmap_data, cmap="RdYlGn", linewidths=0 ,linecolor='lightgray',cbar_kws={'orientation': 'horizontal', 'pad':0.12, 'aspect': 35})
ax.set_xticklabels(np.round(heatmap_data.columns.values,1), fontsize=12, rotation=60)
#ax.set_yticklabels(np.round(heatmap_data.rows.values, 1), fontsize=12)


for i in range(heatmap_data.shape[0]):
    for j in range(heatmap_data.shape[1]):
        value = int(heatmap_data.iloc[i, j])
        if value >= 24:
            rect = patches.Rectangle((j, i), 1, 1, facecolor='none', edgecolor='black', lw=2)
            ax.text(j + 0.5, i + 0.5, value, color='black', ha='center', va='center', fontsize=9)
            ax.add_patch(rect)

ax.invert_yaxis()
# Adjust subplot parameters to minimize space below the colorbar
# The value might need tuning depending on figure size and other elements
# Increase the bottom margin slightly to accommodate the colorbar and x-label
plt.xlabel(r"$\alpha$", fontsize = 16)
plt.ylabel(r"$R_c/R_s$", fontsize = 16)
plt.legend(title='MSE < 0.30', loc='best', title_fontsize=18)
plt.tight_layout()
plt.subplots_adjust(bottom=-0.1) 

plt.show()