import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from scipy.stats import t

# Base directory for the data files
base_dir = '/Users/pauverdeguer/TFG/MATLAB/'

corona = 25
angle = '3.14'
alpha = 0.4

plt.figure(figsize=(10, 7))

# Load IceCube and Fermi data
icecubefilepath = os.path.join(base_dir, 'ngc1068_spectrum_95.txt')
fermifilepath = os.path.join(base_dir, 'FermiData.txt')

icecubedata = np.loadtxt(icecubefilepath, delimiter='\t', skiprows=1)
fermidatapoints = np.loadtxt(fermifilepath, delimiter=' ', max_rows=5)
fermidataupperlims = np.loadtxt(fermifilepath, delimiter=' ', skiprows=5)

fermidatapointsx = fermidatapoints[:, 0]
fermidatapointsy = fermidatapoints[:, 2] * 6.242e11
fermidatapointsyerr = fermidatapoints[:, 3] * 6.242e11
fermidataupperlimsx = fermidataupperlims[:, 0]
fermidataupperlimsy = fermidataupperlims[:, 2] * 6.242e11

icecubex = icecubedata[:, 0] * 1e9
icecubey = icecubedata[:, 1] * 1e12
icecubelowerlimit = icecubedata[:, 2] * 1e12
icecubeupperlimit = icecubedata[:, 3] * 1e12

plt.errorbar(fermidatapointsx, fermidatapointsy, yerr=fermidatapointsyerr, marker='o', linestyle=' ', color='k', label='Fermi Data', linewidth=2, markersize=8)
plt.errorbar(fermidataupperlimsx, fermidataupperlimsy, marker='o', yerr=0.2*fermidataupperlimsy, uplims=True, linestyle=' ', color='k', linewidth=2)
plt.loglog(icecubex, icecubey, label='IceCube', color='purple')
plt.fill_between(icecubex, icecubelowerlimit, icecubeupperlimit, color='purple', alpha=0.2)
plt.loglog(icecubex, icecubelowerlimit, color='purple', linestyle='-', alpha=0.5)
plt.loglog(icecubex, icecubeupperlimit, color='purple', linestyle='-', alpha=0.5)

# Load your simulation data
datafilename = f'data/corona_{corona:.2f}_angle_{float(angle):.2f}_alpha_{alpha:.2f}.csv'
datafilepath = os.path.join(base_dir, datafilename)
data = np.loadtxt(datafilepath, delimiter=',', skiprows=1)
data = data[data[:, 1] != 0]

x = data[:, 0]
y = data[:, 1]



# Fit a polynomial in log-log space
logx = np.log10(x)
logy = np.log10(y)
deg = 5  # Degree of the polynomial (adjustable)

# Fit the polynomial
coeffs, cov = np.polyfit(logx, logy, deg=deg, cov=True)
p = np.poly1d(coeffs)

# Predict y and calculate confidence intervals
logy_fit = p(logx)
y_fit = 10**logy_fit

# Calculate confidence intervals
alpha_conf = 0.05  # 95% confidence
n = len(logy)
p_err = np.sqrt(np.diag(cov))
t_val = t.ppf(1.0 - alpha_conf / 2.0, n - deg - 1)

# Standard error of the fit
fit_stderr = np.sqrt(np.sum((logy - logy_fit)**2) / (n - deg - 1))
logy_upper = logy_fit + t_val * fit_stderr
logy_lower = logy_fit - t_val * fit_stderr

y_upper = 10**logy_upper
y_lower = 10**logy_lower

# Sort for plotting
sorted_idx = np.argsort(x)
plt.loglog(x[sorted_idx], y_fit[sorted_idx], color='darkblue', label='Polynomial Fit')
plt.fill_between(x[sorted_idx], y_lower[sorted_idx], y_upper[sorted_idx], color='darkblue', alpha=0.3, label='95% Confidence')

# Final plot adjustments
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title(rf'$R_c/R_s=${corona}, $\alpha=${alpha}, $\theta_\max=\pi$', fontsize=20, pad=24)
plt.xlabel(r'$E$ [eV]', fontsize=16)
plt.ylabel(r'$E^2 dN/dE$ [eV cm$^{-2}$ s$^{-1}$]', fontsize=16)
plt.tight_layout()
plt.legend(loc='best', fontsize=12)
plt.grid(True)

plt.loglog(x, y, marker='o', linestyle=' ', color='royalblue', markersize=8, label='Simulated Flux')
plt.savefig('/Users/pauverdeguer/TFG/LaTeX/Figures/example_simulation_plot.pdf')
plt.show()
