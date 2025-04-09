import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.multioutput import MultiOutputRegressor
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
from scipy.signal import savgol_filter
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.multioutput import MultiOutputRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV




# --- Parameters ---
corona_values = [1, 2, 3, 5, 8, 10, 12, 15, 17, 18, 20, 22, 25, 30, 33, 35, 40, 45, 50, 60, 70, 80, 90, 100, 500, 1000]
alpha_values = [0, 0.01, 0.1, 0.5, 0.8, 1, 1.5, 3, 4, 6, 8, 10]
theta_values = [np.pi]
num_points = 100



selected_corona = 9
selected_alpha = 1.5
selected_theta = np.pi

base_dir = '/Users/pauverdeguer/TFG/MATLAB/data'
external_data_dir = '/Users/pauverdeguer/TFG/MATLAB'

# --- Determine common log10(x) grid ---
logx_min, logx_max = [], []

for corona in corona_values:
    for alpha in alpha_values:
        for theta in theta_values:
            filename = f'corona_{corona}_angle_{theta:.2f}_alpha_{alpha:.2f}.csv'
            filepath = os.path.join(base_dir, filename)
            try:
                data = np.loadtxt(filepath, delimiter=',', skiprows=1)
                x, y = data[:, 0], data[:, 1]
                valid = (x > 0) & (y > 0)
                logx = np.log10(x[valid])
                logx_min.append(logx.min())
                logx_max.append(logx.max())
            except:
                print(f"Could not load file: {filepath}")
                continue

common_logx_grid = np.linspace(min(logx_min), max(logx_max), num_points)

# --- Build dataset ---
X, Y = [], []

for corona in corona_values:
    for alpha in alpha_values:
        for theta in theta_values:
            filename = f'corona_{corona}_angle_{theta:.2f}_alpha_{alpha:.2f}.csv'
            filepath = os.path.join(base_dir, filename)
            try:
                data = np.loadtxt(filepath, delimiter=',', skiprows=1)
                x, y = data[:, 0], data[:, 1]
                valid = (x > 0) & (y > 0)
                if valid.sum() < 2:
                    continue
                logx = np.log10(x[valid])
                logy = np.log10(y[valid])
                logy_interp = np.interp(common_logx_grid, logx, logy)
                X.append([corona, alpha, theta])
                Y.append(logy_interp)
            except:
                print(f"Could not load file: {filepath}")
                continue

X, Y = np.array(X), np.array(Y)
print("Loaded", X.shape[0], "samples")


# --- Train/test split ---
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

# --- Neural network pipeline ---

# --- Emphasize the sample with corona == 17 and alpha == 0 ---
important_corona = 17
important_alpha = 0
important_theta = np.pi  # or use theta_values[0]

# Find indices of the important sample(s)
important_indices = np.where(
    (X_train[:, 0] == important_corona) &
    (X_train[:, 1] == important_alpha) &
    (X_train[:, 2] == important_theta)
)[0]

# Repeat the important sample N times to give it more weight
N = 40  # Adjust N to give more or less emphasis
if important_indices.size > 0:
    X_important = np.repeat(X_train[important_indices], N, axis=0)
    Y_important = np.repeat(Y_train[important_indices], N, axis=0)

    # Add to training set
    X_train = np.vstack([X_train, X_important])
    Y_train = np.vstack([Y_train, Y_important])

    # Shuffle the data
    from sklearn.utils import shuffle
    X_train, Y_train = shuffle(X_train, Y_train, random_state=42)

    print(f"Boosted importance of sample: corona={important_corona}, alpha={important_alpha} ({N}x)")
else:
    print("Important sample not found in training set.")


# --- Define pipeline with MLPRegressor ---
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('mlp', MLPRegressor(max_iter=5000, early_stopping=True))
])

# --- Define hyperparameter grid ---
param_grid = {
    'mlp__hidden_layer_sizes': [(400)],
    'mlp__activation': ['relu'],
    'mlp__alpha': [1e-3],
    'mlp__learning_rate_init': [0.001],
    'mlp__solver': ['lbfgs'],
    'mlp__batch_size': ['auto']
}

# --- Grid search ---
search = GridSearchCV(
    pipeline,
    param_grid,
    cv=2,
    scoring='neg_mean_squared_error',
    verbose=2,
    n_jobs=-1
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

# --- Predict for new parameters ---
X_new = np.array([[selected_corona, selected_alpha, selected_theta]])
Y_new_log = best_model.predict(X_new)[0]
y_new = 10**Y_new_log
Y_new_smooth = savgol_filter(Y_new_log, window_length=15, polyorder=5)
y_smooth = 10**Y_new_smooth


# --- Plot original simulation for comparison ---
# filename = f'corona_{selected_corona}_angle_{selected_theta:.2f}_alpha_{selected_alpha:.2f}.csv'
# filepath = os.path.join(base_dir, filename)
# try:
#     data = np.loadtxt(filepath, delimiter=',', skiprows=1)
#     x, y = data[:, 0], data[:, 1]
#     valid = (x > 0) & (y > 0)
#     plt.plot(x[valid], y[valid], 'o', label="Original Simulated", color='blue')
# except Exception as e:
#     print(f"Could not load original data: {e}")

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
    # Interpolate model prediction
    logx_pred = np.log10(x_pred)
    logy_pred = np.log10(y_pred)
    
    # Fermi detections (points)
    logfermipx = np.log10(fermipx)
    logfermipy = np.log10(fermipy)

    interp_points = np.interp(fermipx, x_pred, y_pred)
    
    mse = np.mean((fermipy - interp_points) ** 2)
    return mse


# --- Grid search over corona and alpha values ---
corona_grid = np.linspace(1, 100, 99)
alpha_grid = np.linspace(0, 2, 250)
mse_grid = np.zeros((len(alpha_grid), len(corona_grid)))

for i, alpha in enumerate(alpha_grid):
    for j, corona in enumerate(corona_grid):
        X_new = np.array([[corona, alpha, np.pi]])
        Y_log_pred = best_model.predict(X_new)[0]
        y_pred = 10 ** Y_log_pred
        mse = get_mse_on_fermi_with_ul(10**common_logx_grid, y_pred, np.concatenate((fpx, fulx)), np.concatenate((fpy, fuly)))
        if (alpha == 1 and corona == 10) or (alpha == 0 and corona == 17):
            print(f"Alpha: {alpha}, Corona: {corona}, MSE: {mse}")
        mse_grid[i, j] = mse

    # --- Create contour plot ---
plt.figure(figsize=(8, 9))
alpha_mesh, corona_mesh = np.meshgrid(corona_grid, alpha_grid)

contour = plt.contourf(
    corona_mesh, alpha_mesh, mse_grid,
    levels=250, cmap='RdYlGn_r'
)
cb = plt.colorbar(contour,location='top', pad=0.02)
cb.set_label('Mean Squared Error', fontsize=16, labelpad = 10)

# Draw contour lines for clarity
CS = plt.contour(
    corona_mesh, alpha_mesh, mse_grid,
    levels=5, colors='black', linestyles='--'
)
plt.clabel(CS, inline=True, fmt='%.1f', inline_spacing=10, fontsize=12, use_clabeltext=True)

plt.xlabel(r"$\alpha$", fontsize = 16)
plt.ylabel(r"$R_c/R_s$", fontsize = 16)

plt.grid(True)
plt.tight_layout()
plt.show()


# --- Plot predictions and original data for different values ---
plt.figure(figsize=(10, 6))

for corona, alpha in [(5,0),(17,0),(5, 0.5), (10, 1.0), (20, 1.5)]:
    X_new = np.array([[corona, alpha, np.pi]])
    Y_log_pred = best_model.predict(X_new)[0]
    y_pred = 10 ** Y_log_pred

    # Plot model predictions
    plt.plot(10**common_logx_grid, y_pred, label=f"Prediction (corona={corona}, alpha={alpha})")

    # Plot original data if available
    filename = f'corona_{corona}_angle_{np.pi:.2f}_alpha_{alpha:.2f}.csv'
    filepath = os.path.join(base_dir, filename)
    try:
        data = np.loadtxt(filepath, delimiter=',', skiprows=1)
        x, y = data[:, 0], data[:, 1]
        valid = (x > 0) & (y > 0)
        plt.plot(x[valid], y[valid], 'o', label=f"Original Data (corona={corona}, alpha={alpha})")
    except Exception as e:
        print(f"Could not load original data for corona={corona}, alpha={alpha}: {e}")
    
    plt.errorbar(fpx, fpy, yerr=fperr, marker='o', linestyle=' ', color='k', label='Fermi Data', linewidth=2, markersize=8)
    plt.errorbar(fulx, fuly, marker='o', yerr=fulerr, uplims=True, linestyle=' ', color='k', linewidth=2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Model Predictions vs Original Data")
    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.show()