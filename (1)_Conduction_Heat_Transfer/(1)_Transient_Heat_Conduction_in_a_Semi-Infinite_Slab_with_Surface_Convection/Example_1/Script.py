import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

# --- Range of variables ---
tau = np.linspace(0.01, 3, 30)   # tau range
eta = np.linspace(0, 5, 20)      # eta range

# Generate x (eta) and t (tau) matrices for 3D plot
x, t = np.meshgrid(eta, tau)

# Define theta function
def theta(x, t):
    return erfc(0.5 * x / t) - np.exp(x + t**2) * erfc(0.5 * x / t + t)

# --- Figure 1: 3D mesh surface plot ---
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
ax1.plot_surface(x, t, theta(x, t), cmap='viridis')

ax1.set_xlabel(r'$\eta$')
ax1.set_ylabel(r'$\tau$')
ax1.set_zlabel(r'$\theta$')
ax1.set_title("Transient Heat Conduction in a Semi-Infinite Solid")

# --- Figure 2: 2D plots for different eta values ---
eta_vals = np.arange(0, 6, 1)           # eta = 0,1,2,3,4,5
tau_vals = np.linspace(0.01, 4, 40)

plt.figure()
for k in eta_vals:
    thet = theta(k, tau_vals)
    plt.plot(tau_vals, thet, 'k-')
    plt.text(0.92*4, 1.02*thet[-1], r'$\eta = {}$'.format(k))

plt.xlabel(r'$\tau$')
plt.ylabel(r'$\theta$')
plt.title("Temperature Distribution vs Time at Different $\eta$")
plt.grid(True)

plt.show()