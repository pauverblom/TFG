import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, LogLocator

k = 8.67e-5
A = 2.72556e11
T = 239.083


data = np.array([
    [3.67282e-4, 1.00000e-4], [4.75883e-4, 3.01565e-4], [7.58578e-4, 1.99053e-3],
    [8.86135e-4, 3.27691e-3], [1.83021e-3, 5.45900e-2], [2.63027e-3, 2.26814e-1],
    [3.98107e-3, 1.16684e0], [5.72137e-3, 5.59014e0], [9.12011e-3, 3.20003e1],
    [1.53109e-2, 2.18880e2], [1.88365e-2, 4.30547e2], [2.93427e-1, 4.00952e2],
    [3.25462e-1, 1.89824e2], [3.60994e-1, 6.75923e1], [3.80189e-1, 4.90593e1],
    [4.00405e-1, 2.58447e1], [4.00405e-1, 1.06114e1], [4.21697e-1, 4.67846e0],
    [4.44120e-1, 1.92090e0], [4.92606e-1, 3.73392e-1], [5.18800e-1, 1.07381e-1],
    [5.75440e-1, 2.40682e-2], [6.06038e-1, 5.59014e-3], [6.38263e-1, 1.99053e-3],
    [6.72202e-1, 5.14446e-4], [7.07946e-1, 1.70592e-4], [7.07946e-1, 1.00000e-4],
    [7.41310e-3, 1.31388e1], [1.18168e-2, 8.07638e1], [1.69824e-2, 3.12497e2],
    [1.41661e-2, 1.53309e2], [1.03813e-2, 5.08377e1], [8.22243e-3, 2.01430e1],
    [5.85464e-4, 6.36977e-4], [4.65050e-3, 2.29522e0], [3.07256e-3, 4.79085e-1]
])

data2 = np.array([
    [2137962.089502241, 0.0001036248506253716], [1830206.1063110605, 0.0006600664043541983],
    [1487647.3740795061, 0.010611412465291561], [1341219.935405375, 0.04734321337998158],
    [1090184.4923851313, 0.17677574671401622], [933254.3007969924, 0.6600664043541976],
    [841395.1416451999, 1.4447536969037158], [758577.5750291852, 2.1374546404928227],
    [616595.001861486, 5.792774169758715], [475882.7909888171, 17.468999855573255],
    [429042.18923888745, 27.75229920291265], [314412.6420258812, 56.56884604299137],
    [242661.00950824242, 96.50194851573197], [137246.0961007571, 218.88005318453503],
    [77624.71166286927, 347.7259588841705], [48696.75251658651, 430.54705065022387],
    [30549.21113215522, 479.08477682064785]
])
# Set up plot with log-log scale

fig, ax = plt.subplots(figsize=(12, 5))

# Generate energy range (similar to Mathematica's ϵ range)
epsilon = np.logspace(-4, 8, 1000)  # 10^-4 to 10^8

# Define the l1, l2, l3 functions (replace with your actual functions)
def l1(epsilon):
    return np.piecewise(
        epsilon,
        [(epsilon > 1e-4) & (epsilon < 0.043181)],
        [lambda x: 10**9.7 * x**4, np.nan]
    )

def l2(epsilon):
    return np.piecewise(
        epsilon,
        [(epsilon > 0.043181) & (epsilon < 256268.75096)],
        [lambda x: 10**3.9 * x**(-0.25), np.nan]
    )

def l3(epsilon):
    return np.piecewise(
        epsilon,
        [(epsilon > 256268.75096) & (epsilon < 1e8)],
        [lambda x: 10**35 * x**(-6), np.nan]
    )

# Define Planck's Law fit function (replace with your actual fit function)
def planck(epsilon, a, c):
    return a * epsilon**5 / (np.exp(epsilon / (k * c)) - 1)

# 1. Plot l1, l2, l3 functions (replace with your actual functions)
ax.loglog(epsilon, l1(epsilon), color='red', linewidth=3)
ax.loglog(epsilon, l2(epsilon), color='red', linewidth=3)
ax.loglog(epsilon, l3(epsilon), color='red', linewidth=3)

# 2. Plot Planck's Law Fit (replace with your actual fit function)
ax.loglog(epsilon, planck(epsilon, A, T), color='gray', linestyle='--', linewidth=2)

# 3. Plot data points (replace with your actual data)
ax.plot(data[:, 0], data[:, 1], 'o', color='black', markersize=5)
ax.plot(data2[:, 0], data2[:, 1], '^', color='black', markersize=5)

# Set axis limits
ax.set_xlim(9e-5, 1e8)
ax.set_ylim(5e-5, 5e4)

# Configure grid lines
ax.grid(True, which='major', linestyle='--', alpha=0.7)
ax.set_axisbelow(True)

# Custom tick formatting


# Ensure log scale is correctly applied
ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=7))
ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=15))

ax.xaxis.set_minor_locator(LogLocator(base=10.0, numticks=10))
ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs = np.arange(2,10),numticks=150))


# Styling
ax.set_xlabel('$E$ [eV]', fontsize=18)
ax.set_ylabel('$E^2 \Phi$ [eV s$^{-1}$ cm$^{-2}$]', fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=16, width=2, length=6)

# Custom legend
legend_elements = [
    plt.Line2D([0], [0], color='red', lw=3, label=r'$E^2 \phi$'),
    plt.Line2D([0], [0], color='gray', linestyle='--', lw=1.5, label="Planck's Law Fit"),
    plt.Line2D([0], [0], marker='o', color='black', lw=0, markersize=5, label='Photon background \n (starburst)'),
    plt.Line2D([0], [0], marker='^', color='black', lw=0, markersize=5, label='Photon background \n (corona)')
]



ax.legend(bbox_to_anchor = (1.02, 0.75), handles=legend_elements, loc='upper left', fontsize=16, frameon=True, borderaxespad=0)

# Set axis color and thickness
for spine in ax.spines.values():
    spine.set_color('black')
    spine.set_linewidth(1)

plt.tight_layout()
plt.savefig('/Users/pauverdeguer/TFG/LaTeX/Figures/Background_SED_AGN.pdf')
