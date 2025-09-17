import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.special import j0, j1   # Bessel functions J0 and J1

# -------------------------------------------------------------------
# Helper function: FindZeros
# -------------------------------------------------------------------
def find_zeros(fun, Nroot, x, w):
    """
    Locate up to Nroot zeros of a function fun(x, w)
    within the interval defined by vector x.
    """
    f = fun(x, w)
    indx = np.where(f[:-1] * f[1:] < 0)[0]   # sign changes

    L = len(indx)
    if L < Nroot:
        print(f"Warning: Requested {Nroot} roots, but only {L} sign changes detected.")
        Nroot = L

    Rt = np.zeros(Nroot)
    for k in range(Nroot):
        bracket = (x[indx[k]], x[indx[k] + 1])
        Rt[k] = brentq(lambda z: fun(z, w), *bracket)

    return Rt

# -------------------------------------------------------------------
# Problem setup
# -------------------------------------------------------------------
Bi = 0.5           # Biot number (dimensionless)
Nroot = 15         # Number of roots

# Define transcendental equation
def cylinder_roots(x, Bi):
    return x * j1(x) - Bi * j0(x)

# Find roots
r = find_zeros(cylinder_roots, Nroot, np.linspace(0, 50, 200), Bi)

# -------------------------------------------------------------------
# Time domain
tau = np.linspace(0, 1.5, 20)     # tau [1 x Nt]
t, rt = np.meshgrid(tau, r)       # t, rt [Nroot x Nt]

# Transient exponential decay
Fn = np.exp(-t * rt**2)

# Coefficients for each term
cn = 2 * j1(r) / (r * (j0(r)**2 + j1(r)**2))   # [Nroot]

# Broadcast coefficients
ccn = np.outer(cn, np.ones(len(tau)))          # [Nroot x Nt]

# Multiply with exponential terms
pro = ccn * Fn                                 # [Nroot x Nt]

# Radial position
rstar = np.linspace(0, 1, 20)                  # [1 x Nr]
R, rx = np.meshgrid(rstar, r)                  # [Nroot x Nr]

# Bessel function evaluated at scaled radial positions
Jo = j0(rx * R)                                # [Nroot x Nr]

# Temperature distribution (matrix multiply)
the = Jo.T @ pro                               # [Nr x Nt]

# -------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------
rr, tt = np.meshgrid(rstar, tau)   # [Nt x Nr]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(rr, tt, the.T, cmap='viridis')

ax.set_xlabel(r'$\xi$')   # Dimensionless radial position
ax.set_ylabel(r'$\tau$')  # Dimensionless time
ax.set_zlabel(r'$\theta$')# Dimensionless temperature
ax.view_init(elev=-34, azim=49.5)

plt.show()