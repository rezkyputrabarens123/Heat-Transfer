import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, solve_bvp

# ---------------------------
# Parameters
# ---------------------------
Bi = 0.1      # Biot number
Tr = 0.55     # Reference temperature ratio (right Dirichlet)
Sigma = 1.0   # Source term (same as MATLAB 'Sigma')

# Spatial and temporal discretization to match MATLAB
xi = np.linspace(0.0, 1.0, 41)    # 41 spatial nodes between 0 and 1
tau_plot = np.linspace(0.0, 1.0, 101)  # times to evaluate (same as tau in MATLAB)

Nx = xi.size
dx = xi[1] - xi[0]

# Initial condition function (pdeIC in MATLAB)
def pdeIC(x):
    return 1.0 - 0.45 * x

# ----------------------------------------------------------------------
# Method-of-lines for PDE: u_t = u_xx + Sigma
# Boundary conditions:
#   at x=0 -> Robin: u_x(0) = Bi * u(0)
#   at x=1 -> Dirichlet: u(1) = Tr
# ----------------------------------------------------------------------

def compute_u_xx(u):
    """
    Compute second spatial derivative u_xx with boundary conditions:
      left:  u_x(0) = Bi * u0   (Robin)
      right: u(N-1) = Tr       (Dirichlet, node pinned)
    u: array length Nx
    returns: u_xx array length Nx
    """
    u_xx = np.zeros_like(u)

    # interior points (central difference) for i = 1 .. Nx-2
    if Nx > 2:
        u_xx[1:-1] = (u[2:] - 2*u[1:-1] + u[0:-2]) / dx**2

    # Left boundary: Robin u_x(0) = Bi * u0
    # Use ghost node u_{-1} such that (u1 - u_{-1})/(2dx) = Bi * u0
    u_ghost_left = u[1] - 2.0 * dx * Bi * u[0]
    u_xx[0] = (u[1] - 2.0 * u[0] + u_ghost_left) / dx**2

    # Right boundary: pinned Dirichlet node (we'll set du_dt[-1] = 0 later).
    # We can set u_xx[-1] to zero (unused because du_dt[-1] will be forced to 0).
    u_xx[-1] = 0.0

    return u_xx

def pde_rhs(t, u):
    """
    Right-hand side for the semi-discretized PDE (MOL).
    We avoid modifying the integrator's input array in-place by working on a copy.
    """
    u_work = u.copy()

    # Enforce Dirichlet boundary at the right node for derivative computations
    u_work[-1] = Tr

    u_xx = compute_u_xx(u_work)

    du_dt = u_xx + Sigma  # PDE: u_t = u_xx + Sigma

    # For right boundary node, pin it: ensure time derivative is zero (Dirichlet)
    du_dt[-1] = 0.0

    return du_dt

# Initial state vector (ensure float)
u0 = pdeIC(xi).astype(float)

# Integrate in time and request solution at tau_plot points
sol = solve_ivp(fun=pde_rhs, t_span=(tau_plot[0], tau_plot[-1]), y0=u0,
                t_eval=tau_plot, method='BDF', atol=1e-8, rtol=1e-6)

# sol.y shape: (Nx, Nt). We want theta as (Nt, Nx)
theta = sol.y.T  # shape (Nt, Nx)
tau = sol.t      # times corresponding to rows of theta

# -----------------------
# Figure 1: transient response at selected positions
# -----------------------
z = np.arange(0.0, 1.0 + 1e-12, 0.25)  # [0, 0.25, 0.5, 0.75, 1.0]

plt.figure(1)
for k, z_k in enumerate(z, start=1):
    # find index in xi closest to z_k
    kk = np.argmin(np.abs(xi - z_k))
    plt.plot(tau, theta[:, kk], 'k-')

    # text labels similar to MATLAB's logic
    if k == 1:
        plt.text(0.5, 1.02 * theta[-1, kk], r'$\xi = 0.0 \; \text{and} \; 0.25$')
    elif k > 2:
        plt.text(0.5, theta[-1, kk] + 0.02, rf'$\xi = {xi[kk]:.2f}$')

plt.axis([0, 1, 0.5, 1.0])
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\theta$')
plt.title('Transient Response at Selected Positions')
plt.grid(False)

# -----------------------
# Figure 2: spatial profiles at different times and steady state
# -----------------------

# find time index where theta(:, 1) (xi=0 column) is minimum (MATLAB used min over times for xi=0)
thmin = np.min(theta[:, 0])
imin = np.argmin(theta[:, 0])  # index into tau where theta at xi=0 is minimum

plt.figure(2)
# plot tau=0 initial
plt.plot(xi, theta[0, :], 'k-', linewidth=2)

# plot at tau = tau[1] (MATLAB used the second time entry) if available
if theta.shape[0] >= 2:
    plt.plot(xi, theta[1, :], 'k')

# plot at time of minimum at xi=0
plt.plot(xi, theta[imin, :], 'k:')

# -----------------------
# Steady-state solution using solve_bvp (equivalent to bvp4c)
# ODE system:
#   y1' = y2
#   y2' = -Sigma
# BCs: y2(0) - Bi*y1(0) = 0
#      y1(1) = Tr
# -----------------------
def barode(x, y):
    # y has shape (2, len(x))
    dydx = np.zeros_like(y)
    dydx[0, :] = y[1, :]                 # y1' = y2
    dydx[1, :] = -Sigma * np.ones_like(x)   # y2' = -Sigma
    return dydx

def barbc(ya, yb):
    return np.array([ya[1] - Bi * ya[0],
                     yb[0] - Tr])

# initial mesh and guess (MATLAB used bvpinit with guess [1 1])
x_guess = np.linspace(0.0, 1.0, 20)
y_guess = np.zeros((2, x_guess.size))
y_guess[0, :] = 1.0  # theta initial guess
y_guess[1, :] = 1.0  # theta' initial guess

bvp_sol = solve_bvp(barode, barbc, x_guess, y_guess, tol=1e-6, max_nodes=1000)

# Evaluate steady-state solution on a finer grid
x_plot = np.linspace(0.0, 1.0, 100)
y_plot = bvp_sol.sol(x_plot)  # shape (2, len(x_plot))

# plot steady-state theta (y_plot[0,:])
plt.plot(x_plot, y_plot[0, :], 'k--')

# finalize figure 2 formatting
plt.xlabel(r'$\xi$')
plt.ylabel(r'$\theta$')

# Build legend entries matching plotted curves (4 curves plotted)
legend_entries = [
    r'$\tau = 0$ (Initial condition)'
]
if theta.shape[0] >= 2:
    legend_entries.append(r'$\tau = {:.3f}$'.format(tau[1]))
else:
    legend_entries.append(r'$\tau =$ N/A')
legend_entries.append(r'$\tau = {:.3f}$ (Minimum at ξ = 0)'.format(tau[imin]))
legend_entries.append(r'Steady state')

plt.legend(legend_entries, loc='lower left')

plt.title('Spatial Profiles at Different Times')
plt.grid(False)
plt.show()