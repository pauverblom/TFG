import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
from sklearn.decomposition import PCA
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
import csv
import seaborn as sns
import matplotlib.patches as patches
from sklearn.base import clone


base_dir = '/Users/pauverdeguer/TFG/MATLAB/data'

# --- Parameters ---
# Generate corona_values and alpha_values from filenames in base_dir
corona_values = set()
alpha_values = set()



for filename in os.listdir(base_dir):
    if filename.startswith("corona_") and filename.endswith(".csv"):
        try:
            parts = filename.split('_')
            corona = float(parts[1])
            alpha = float(parts[-1].replace('.csv', ''))
            corona_values.add(corona)
            alpha_values.add(alpha)
        except ValueError:
            continue

corona_values = sorted(corona_values)
alpha_values = sorted(alpha_values)
theta_values = [np.pi]
num_points = 100


# --- Plot heatmap of existing combinations ---

# Create a set to store existing combinations
existing_combinations = set()



# Scan the directory for files and extract combinations
for filename in os.listdir(base_dir):
    if filename.startswith("corona_") and filename.endswith(".csv"):
        try:
            parts = filename.split('_')
            corona = float(parts[1])
            alpha = float(parts[-1].replace('.csv', ''))
            existing_combinations.add((corona, alpha))
        except ValueError:
            continue

# Unpack existing combinations into separate lists
corona_vals, alpha_vals = zip(*existing_combinations) if existing_combinations else ([], [])

plt.figure(figsize=(10, 8))
plt.scatter(alpha_vals, corona_vals, color='blue', s=60, edgecolors='k')  # s: marker size

plt.xlabel(r"$\alpha$", fontsize=14)
plt.ylabel(r"$R_c/R_s$", fontsize=14)
plt.title("Existing Combinations of Alpha and Corona", fontsize=16)

#plt.xticks(sorted(alpha_values), rotation=45, fontsize=10)
#plt.yticks(sorted(corona_values), fontsize=10)

plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
#plt.show()
        


external_data_dir = '/Users/pauverdeguer/TFG/MATLAB'

# --- Determine common log10(x) grid ---

# --- Parameters ---
num_points = 200
tail_length = 10
tail_dx = 0.02  # step size in log(x) for extrapolated tail
tail_samples_for_slope = 3

logx_min, logx_max = [], []
X, Y = [], []

# --- First pass: scan all files to find logx range ---
for corona in corona_values:
    for alpha in alpha_values:
        for theta in theta_values:
            filename = f'corona_{corona:.2f}_angle_{theta:.2f}_alpha_{alpha:.2f}.csv'
            filepath = os.path.join(base_dir, filename)
            try:
                data = np.loadtxt(filepath, delimiter=',', skiprows=1)
                x, y = data[:, 0], data[:, 1]
                valid = (x > 0) & (y > 0)
                logx = np.log10(x[valid])
                logx_min.append(logx.min())
                logx_max.append(logx.max())
            except:
                continue

common_logx_grid = np.linspace(min(logx_min), max(logx_max), num_points)

# --- Second pass: load, interpolate, and augment ---
for corona in corona_values:
    for alpha in alpha_values:
        for theta in theta_values:
            filename = f'corona_{corona:.2f}_angle_{theta:.2f}_alpha_{alpha:.2f}.csv'
            filepath = os.path.join(base_dir, filename)
            try:
                data = np.loadtxt(filepath, delimiter=',', skiprows=1)
                x, y = data[:, 0], data[:, 1]
                valid = (x > 0) & (y > 0)
                if valid.sum() < tail_samples_for_slope:
                    continue

                logx = np.log10(x[valid])
                logy = np.log10(y[valid])

                logy_interp = np.interp(common_logx_grid, logx, logy, right=-4, left=-4)

                X.append([corona, alpha, theta])
                Y.append(logy_interp)

            except Exception as e:
                print(f"Could not load file: {filepath} ({e})")
                continue

X, Y = np.array(X), np.array(Y)
extended_logx_grid = common_logx_grid  # same for all samples
print("Loaded", X.shape[0], "samples. Each with", Y.shape[1], "points.")


# --- Train/test split ---
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

""" THIS IS FOR GRID SEARCH, UNCOMMENT TO USE
# --- Define pipeline with MLPRegressor ---
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('mlp', MLPRegressor())
])

# --- Define hyperparameter grid ---
param_grid = {
    'mlp__hidden_layer_sizes': [(100,100)],
    'mlp__alpha': [1e-3],
    'mlp__learning_rate_init': [1e-4],
    'mlp__solver': ['adam', 'lbfgs'],
    'mlp__early_stopping': [True],
    'mlp__max_iter': [5000],
}
# --- Grid search ---
search = GridSearchCV(
    pipeline,
    param_grid,
    cv=3,
    scoring='neg_mean_squared_error',
    verbose=10,
    n_jobs=-1,
)
print("Starting GridSearchCV...")
search.fit(X_train, Y_train)
print("Grid search done.")

print("Best Parameters Found:")
print(search.best_params_)
# --- Use best model ---
best_model = search.best_estimator_

# --- Evaluate on test set ---
Y_pred_test = best_model.predict(X_test)
test_mse = mean_squared_error(Y_test, Y_pred_test)
print("Best MLP Test MSE:", test_mse)
"""
# --- Define pipeline with MLPRegressor ---
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('mlp', MLPRegressor())
])



#print("Starting Training...")
#model.fit(X_train, Y_train)

# --- Evaluate on test set ---
#Y_pred_test = model.predict(X_test)
#test_mse = mean_squared_error(Y_test, Y_pred_test)
#print("Test MSE:", test_mse)

# --- Load IceCube and Fermi data ---
try:
    icecube_path = os.path.join(external_data_dir, 'ngc1068_spectrum_95.txt')
    fermi_path = os.path.join(external_data_dir, 'FermiData.txt')

    # IceCube
    icecube = np.loadtxt(icecube_path, delimiter='\t', skiprows=1)
    icecubex = icecube[:, 0] * 1e9
    icecubey = icecube[:, 1] * 1e12
    ic_lower = icecube[:, 2] * 1e12
    ic_upper = icecube[:, 3] * 1e12

    # Fermi
    fermidata_points = np.loadtxt(fermi_path, delimiter=' ', max_rows=5)
    fermidata_ulims = np.loadtxt(fermi_path, delimiter=' ', skiprows=5)
    fpx = fermidata_points[:, 0]
    fpy = fermidata_points[:, 2] * 6.242e11
    fperr = fermidata_points[:, 3] * 6.242e11
    fulx = fermidata_ulims[:, 0]
    fuly = fermidata_ulims[:, 2] * 6.242e11
    fulerr = 0.2 * fuly
except Exception as e:
    print(f"Error loading Fermi or IceCube data: {e}")



# --- Fit to Fermi data ---

def get_mse_on_fermi_with_ul(x_pred, y_pred, fermipx, fermipy):
    # Interpolate y_pred to match the x-coordinates of Fermi points
    y_pred_interp = np.interp(fermipx, x_pred, y_pred)

    # Calculate the mean squared error
    mse = np.mean((y_pred_interp - fermipy) ** 2)

    return mse


# --- Grid search over corona and alpha values ---
#corona_grid = np.sort(np.concatenate((np.linspace(1, 1000, 1750), corona_values)))
corona_grid = np.sort(np.concatenate((np.linspace(1, 1000, 500), corona_values)))
alpha_grid = np.sort(np.concatenate((np.linspace(0, 10, 200), alpha_values)))
# 1) prepare storage
n_runs = 1
mse_runs = np.zeros((n_runs, len(alpha_grid), len(corona_grid)))

# 2) loop over seeds
for seed in range(n_runs):
    # clone your “best” pipeline and give it a new random_state
    model = clone(pipeline)
    model.set_params(mlp__hidden_layer_sizes=(100,100),
                    mlp__activation='tanh',
                    mlp__alpha=1e-3,
                    mlp__learning_rate_init=0.001,
                    mlp__solver='lbfgs',
                    mlp__early_stopping=True,
                    mlp__max_iter=5000,
                    mlp__random_state=seed)
    model.fit(X_train, Y_train)
    print(f"Training completed for seed {seed + 1}/{n_runs}.")
    # 3) for each (i,j) in your α vs Rc grid, predict & compute mse
    for i, alpha in enumerate(alpha_grid):
        for j, corona in enumerate(corona_grid):
            # build your test‐flux spectrum X_new → Y_log_pred → y_pred
            X_new = np.array([[corona, alpha, np.pi]])
            Y_log = model.predict(X_new)[0]
            y_pred = 10**Y_log

            # compute MSE vs the Fermi points
            mse_runs[seed, i, j] = get_mse_on_fermi_with_ul(
                10**extended_logx_grid,
                y_pred,
                np.concatenate((fpx, fulx)),
                np.concatenate((fpy, fuly))
            )
        print(f"Seed {seed + 1}/{n_runs}: Completed row {i + 1}/{len(alpha_grid)}.")
# 4) average (or median) across runs
mean_mse_grid = mse_runs.mean(axis=0)
log_mse_grid  = np.log10(mean_mse_grid + 1e-12)

# --- Identify the best combinations (lowest log(MSE)) ---
# You could also tweak this part to find multiple minima if necessary (e.g., top N)
best_combinations = []
for j in range(log_mse_grid.shape[1]):
    min_mse_index = np.argmin(log_mse_grid[:, j])
    best_combinations.append((min_mse_index, j))

# Find the overall best combination
overall_best_index = np.unravel_index(np.argmin(log_mse_grid), log_mse_grid.shape)
best_alpha_value = alpha_grid[overall_best_index[0]]
best_corona_value = corona_grid[overall_best_index[1]]

# --- Create the contour plot ---
plt.figure(figsize=(8, 8))
alpha_mesh, corona_mesh = np.meshgrid(corona_grid, alpha_grid)

contour = plt.contourf(
    corona_mesh, alpha_mesh, log_mse_grid,
    levels=250, cmap='RdYlGn_r'
)
cb = plt.colorbar(contour, location='top', pad=0.02)
cb.ax.tick_params(labelsize=14)
cb.set_label('log(MSE)', fontsize=16, labelpad=10)

# Optional: Draw contours on top of log-scale data
CS = plt.contour(
    corona_mesh, alpha_mesh, log_mse_grid,
    levels=9, colors='black', linestyles='--'
)
plt.clabel(CS, inline=True, fmt='%.1f' , inline_spacing=10, fontsize=12, use_clabeltext=True, manual=True)

# --- Overlay lines for the best combinations ---
# Extract the corresponding corona and alpha values for the best combinations
best_corona_values = [corona_grid[j] for i, j in best_combinations]
best_alpha_values = [alpha_grid[i] for i, j in best_combinations]

# Plot the best combinations as a line
#plt.plot(best_alpha_values, best_corona_values, '-', color='magenta', linewidth=2, label='Best Combinations')

# Highlight the overall best combination
#plt.scatter(best_alpha_value, best_corona_value, color='black', s=100, label='Best Combination Overall', zorder=5)

plt.xlabel(r"$\alpha$", fontsize=16)
plt.ylabel(r"$R_c/R_s$", fontsize=16)
rect = patches.Rectangle((j, i), 1, 1, facecolor='none', edgecolor='red', lw=2)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
#plt.legend()
plt.tight_layout()
#plt.savefig('LaTeX/Figures/10x60_contour_plot.pdf')
plt.show()


# --- Sort the best parameters by their MSE ---
#sorted_best_combinations = sorted(
#    best_combinations,
#    key=lambda idx: mse_grid[idx[0], idx[1]]
#)

#--- Dump all the parameters sorted by their mse into a file ---
output_file = '/Users/pauverdeguer/TFG/LaTeX/parameter_combinations_sorted.csv'
# Create a list of tuples (alpha, corona, mse)
results = []
for i, alpha in enumerate(alpha_grid):
    for j, corona in enumerate(corona_grid):
        mse = mean_mse_grid[i, j]
        results.append((alpha, corona, mse))

# Sort the results by mse (ascending)
results.sort(key=lambda x: x[2])

# Save to a CSV file
#with open(output_file, mode='w', newline='') as file:
#    writer = csv.writer(file)
#    writer.writerow(['alpha', 'corona', 'mse'])
#    writer.writerows(results)

print(f"Parameter combinations sorted by MSE have been saved to {output_file}")


# --- Selected best combinations as provided ---
# Format: (Rank, alpha, corona, MSE)
# Select specific ranks (e.g., 1, 1000, 5000, ...)
specific_ranks = [1, 1000, 3500, 5000, 10000, 20000, 50000, 300000, 350000]

# Get the flat indices of sorted MSE values
sorted_indices = np.argsort(mean_mse_grid, axis=None)

# Convert flat indices to (i, j) coordinates
ij_pairs = list(zip(*np.unravel_index(sorted_indices, mean_mse_grid.shape)))

# Select only the requested ranks
selected = [(rank, alpha_grid[i], corona_grid[j], mean_mse_grid[i, j])
            for rank, (i, j) in enumerate(ij_pairs, start=1)
            if rank in specific_ranks]

# --- Selected best combinations as provided ---
# Format: (Rank, alpha, corona, MSE)         
selected = [
    (1, 0.60,13, 0.2400858689715596),
    (1000, 0.75,13.57,0.42411418277845125),
    (3500, 1.96,1.0,1.887580770184489),
    (5000, 0.10,57.55,3.763363992000076),
    (10000, 5.38,9.57,54.07225123750154),
    (20000, 8.79,16.99,380.23983585055873),
    (50000, 9.10,78.11,2493.1851860660013),
    (300000, 8.94,383.69,6506.684374790608),
    (430000, 0.0,624.73,13320.48396355941)
]
# --- Compute predictions for each selected combo
# Here we will store the predictions in a list (one per selected combo)
predictions = []  # Each entry will be the y_pred corresponding to the extended x grid

for rank, alpha, corona, mse in selected:
    # Build the model input; note that here we fix theta = pi
    X_new = np.array([[corona, alpha, np.pi]])
    Y_log_pred = model.predict(X_new)[0]
    y_pred = 10 ** Y_log_pred  # convert log-flux back to flux if needed
    predictions.append(y_pred)
    
# --- Create a grid of subplots for each selected combination
num_plots = len(selected)
ncols = int(np.ceil(np.sqrt(num_plots)))
nrows = int(np.ceil(num_plots / ncols))
fig, axs = plt.subplots(nrows, ncols, figsize=(15, 15))
axs = axs.flatten()  # Easier indexing

# --- Loop over the selected combinations and plot each in its own subplot
for idx, (sel_tuple, y_pred) in enumerate(zip(selected, predictions)):
    rank, alpha, corona, mse = sel_tuple
    ax = axs[idx]
    ax.legend(title=f"#{rank}\n log(MSE) = {np.log(mse):.1f}", fontsize=6, loc='lower left')

    # Plot Fermi data
    
    # Plot IceCube data
    icecubefilepath = os.path.join("MATLAB", 'ngc1068_spectrum_95.txt')
    icecubedata = np.loadtxt(icecubefilepath, delimiter='\t', skiprows=1)
    icecubex = icecubedata[:, 0] * 1e9  # Convert to eV
    icecubey = icecubedata[:, 1] * 1e12 # e.g., conversion factor (adjust if needed)
    icecubelowerlimit = icecubedata[:, 2] * 1e12
    icecubeupperlimit = icecubedata[:, 3] * 1e12
    ax.loglog(icecubex, icecubey, label='IceCube', color='purple')
    ax.fill_between(icecubex, icecubelowerlimit, icecubeupperlimit, color='purple', alpha=0.2)
    ax.loglog(icecubex, icecubelowerlimit, color='purple', linestyle='-', alpha=0.5)
    ax.loglog(icecubex, icecubeupperlimit, color='purple', linestyle='-', alpha=0.5)
    
    # Plot original simulation data
    datafilename = f'corona_{corona:.2f}_angle_3.14_alpha_{alpha:.2f}.csv'
    datafilepath = os.path.join("MATLAB/data", datafilename)
    
    # Plot prediction (the y_pred curve is defined over the extended x-grid)
    ax.loglog(10 ** extended_logx_grid, y_pred, label='Prediction', color='r', linewidth=2)
    
    # Titles and labels
    ax.set_title(rf'$R_c/R_s$ = {corona:.1f}, $\alpha$ = {alpha:.1f}', fontsize=10)
    ax.set_xlabel(r'$E$ [eV]')
    ax.set_ylabel(r'$E^2 dN/dE$ [eV cm$^{-2}$ s$^{-1}$]')
    ax.errorbar(fpx, fpy, yerr=fperr, marker='.', linestyle=' ', color='k', label='Fermi Data')
    ax.errorbar(fulx, fuly, yerr=0.1 * fuly, uplims=True, linestyle=' ', color='k')
    try:
        data = np.loadtxt(datafilepath, delimiter=',', skiprows=1)
        # Remove any zero flux points to avoid log-scale issues:
        data = data[data[:, 1] != 0]
        x_orig = data[:, 0]
        y_orig = data[:, 1]
        ax.loglog(x_orig, y_orig, 'o', markerfacecolor='none', color='royalblue',
                  markersize=4, label='Original Simulation')
    except Exception as e:
        print(f"Could not load file: {datafilepath} ({e})")
    
    ax.set_xlim(1e4, 1e14)
    ax.grid(True, which='major', linestyle='--', linewidth=0.5)
    ax.set_yticks([1e-3, 1e-1, 1e1, 1e3])
    ax.tick_params(which='minor', bottom=False, left=False)
    
fig.set_size_inches(8, 8)

plt.tight_layout()        
plt.subplots_adjust(bottom=0.11)  # Add space for the legend

# Create a single legend below all subplots
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, ['IceCube Data', 'Model Prediction', 'Simulated Flux', 'Fermi Data'], loc='lower center', fontsize=10, ncol=len(labels), bbox_to_anchor=(0.5, 0))

#plt.savefig('/Users/pauverdeguer/TFG/LaTeX/Figures/selected_simulations_plot.pdf')
plt.show()

