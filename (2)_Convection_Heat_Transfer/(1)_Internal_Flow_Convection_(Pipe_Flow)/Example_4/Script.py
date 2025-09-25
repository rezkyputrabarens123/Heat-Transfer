import numpy as np
from scipy.integrate import solve_ivp, trapezoid
import matplotlib.pyplot as plt

# -------------------------------
# Parameters
# -------------------------------
Tw = 40       # Wall temperature [°C] (for constant Tw condition)
qw = 10       # Wall heat flux [W/m^2] (for constant q condition)
Re = 40       # Reynolds number (laminar flow, Re < 2300)
Pr = 5        # Prandtl number
R = 0.01      # Pipe radius [m]
L = 0.5       # Pipe length [m]
k = 0.6       # Thermal conductivity of fluid [W/m-K]

Rt = 401      # Number of radial grid points
zt = 50       # Number of axial grid points
dxi = 1 / (Rt - 1)   # Radial step size (normalized radius)

# Dimensionless coordinates
# xi: Dimensionless radial coordinate (0=center, 1=wall)
# zeta: Dimensionless axial coordinate (0=inlet, 1=outlet)
xi = np.linspace(0, 1, Rt)
zeta = np.linspace(0, 1, zt)

# Inlet fluid temperature
T_inlet = 20.0

# ----------------------------------------------------
# Case 1: Constant Wall Temperature (Isothermal)
# ----------------------------------------------------

def pde_system_T(zeta, T_inner, Tw, Re, Pr, R, L, dxi, xi_inner):
    """
    Defines the system of ODEs for the constant wall temperature case.
    The system is solved for the inner nodes (from center to node Rt-2).
    """
    # Reconstruct the full temperature profile including the wall
    T = np.append(T_inner, Tw)
    
    # Initialize the spatial derivative term (laplacian)
    laplacian = np.zeros_like(T_inner)
    
    # Centerline (i=0) using L'Hopital's rule for the cylindrical Laplacian
    laplacian[0] = 4 * (T[1] - T[0]) / (dxi**2)
    
    # Interior nodes (i=1 to Rt-2)
    for i in range(1, len(T_inner)):
        d2T_dxi2 = (T[i+1] - 2*T[i] + T[i-1]) / (dxi**2)
        dT_dxi = (T[i+1] - T[i-1]) / (2 * dxi)
        laplacian[i] = d2T_dxi2 + (1 / xi_inner[i]) * dT_dxi
        
    # Coefficient 'c' from the PDE
    c = Re * Pr * R / L * (1 - xi_inner**2)
    
    # dT/dzeta = laplacian / c
    dTdzeta = laplacian / c
    
    return dTdzeta

# Initial condition (uniform temperature at the inlet for inner nodes)
T0_inner = np.full(Rt - 1, T_inlet)
xi_inner = xi[:-1]

# Solve the ODE system
sol_T_obj = solve_ivp(
    pde_system_T,
    [0, 1],
    T0_inner,
    t_eval=zeta,
    args=(Tw, Re, Pr, R, L, dxi, xi_inner),
    method='BDF' # stiff solver
)

# Reconstruct the full solution matrix including the wall temperature
solT = np.hstack((sol_T_obj.y.T, np.full((zt, 1), Tw)))


# ----------------------------------------------------
# Case 2: Constant Wall Heat Flux
# ----------------------------------------------------

def pde_rhs_F(zeta, T, qw, R, k, dxi, xi):
    """
    Defines the Right-Hand Side (RHS) of the DAE system for the constant flux case.
    This function computes the spatial derivative term (laplacian) for all nodes.
    """
    laplacian = np.zeros_like(T)
    
    # Centerline (i=0)
    laplacian[0] = 4 * (T[1] - T[0]) / (dxi**2)
    
    # Interior nodes (i=1 to Rt-2)
    for i in range(1, Rt - 1):
        d2T_dxi2 = (T[i+1] - 2*T[i] + T[i-1]) / (dxi**2)
        dT_dxi = (T[i+1] - T[i-1]) / (2 * dxi)
        laplacian[i] = d2T_dxi2 + (1 / xi[i]) * dT_dxi

    # Wall (i=Rt-1) using a ghost point to enforce the flux boundary condition
    # Ghost point T_ghost = T[Rt] is derived from dT/dxi = -qw*R/k
    T_ghost = T[-2] - 2 * dxi * qw * R / k
    d2T_dxi2_wall = (T_ghost - 2*T[-1] + T[-2]) / (dxi**2)
    dT_dxi_wall = (T_ghost - T[-2]) / (2*dxi)
    laplacian[-1] = d2T_dxi2_wall + (1 / xi[-1]) * dT_dxi_wall
    
    return laplacian

# Initial condition (uniform temperature at the inlet)
T0_full = np.full(Rt, T_inlet)

# The coefficient 'c' becomes the diagonal of the mass matrix 'M'.
# M is singular because c=0 at the wall (xi=1), creating a DAE system.
c_vec = Re * Pr * R / L * (1 - xi**2)
M = np.diag(c_vec)

# Solve the DAE system
sol_F_obj = solve_ivp(
    pde_rhs_F,
    [0, 1],
    T0_full,
    t_eval=zeta,
    args=(qw, R, k, dxi, xi),
    method='BDF', # Stiff solver required for DAEs
    mass_matrix=M
)

solF = sol_F_obj.y.T

# -------------------------------
# Post-Processing and Analysis
# -------------------------------
NuT = np.zeros(zt)
NuF = np.zeros(zt)

for i in range(zt):
    # --- Constant Wall Temperature ---
    # Bulk mean temperature
    integrand_T = xi * (1 - xi**2) * solT[i, :]
    TmT = 4 * trapezoid(integrand_T, x=xi)
    # Nusselt number
    # Nu = -2 * (dT/dxi)_wall / (Tw - Tm)
    dTdxi_wall_T = (solT[i, -1] - solT[i, -2]) / dxi
    if abs(TmT - Tw) > 1e-6:
        NuT[i] = -2 * dTdxi_wall_T / (TmT - Tw)
    else:
        NuT[i] = NuT[i-1] if i > 0 else 0

    # --- Constant Wall Heat Flux ---
    # Bulk mean temperature
    integrand_F = xi * (1 - xi**2) * solF[i, :]
    TmF = 4 * trapezoid(integrand_F, x=xi)
    # Wall temperature at current z location
    Tw_local = solF[i, -1]
    # Nusselt number
    # Nu = (qw * D) / (k * (Tw_local - Tm))
    if abs(Tw_local - TmF) > 1e-6:
        NuF[i] = (qw * 2 * R) / (k * (Tw_local - TmF))
    else:
        NuF[i] = NuF[i-1] if i > 0 else 0

# Final dimensionless temperature profile at outlet (zeta = 1)
TmT_outlet = TmT
TmF_outlet = TmF
Tw_local_outlet = solF[-1, -1]

ThT = (solT[-1, :] - Tw) / (TmT_outlet - Tw)
ThF = (solF[-1, :] - Tw_local_outlet) / (TmF_outlet - Tw_local_outlet)


# -------------------------------
# Plotting
# -------------------------------
plt.style.use('seaborn-v0_8-whitegrid')

# Chart 1: Temperature Distribution (Outlet Profile)
plt.figure(1, figsize=(8, 6))
plt.plot(xi, ThT, 'k-', label='Constant wall temperature')
plt.plot(xi, ThF, 'k--', label='Constant wall heat flux')
plt.xlabel('$\\xi$ (r/R)', fontsize=12)
plt.ylabel('$\\theta$ (Dimensionless Temperature)', fontsize=12)
plt.title('Outlet Temperature Distribution ($\\zeta=1$)', fontsize=14)
plt.legend()
plt.tight_layout()

# Chart 2: Nusselt Number along Pipe
plt.figure(2, figsize=(8, 6))
plt.plot(zeta, NuT, 'k-', label='Constant wall temperature')
plt.plot(zeta, NuF, 'k--', label='Constant wall heat flux')
plt.xlabel('$\\zeta$ (z/L)', fontsize=12)
plt.ylabel('Nu (Nusselt number)', fontsize=12)
plt.title('Local Nusselt Number Along Pipe', fontsize=14)
plt.ylim([0, 6])
plt.legend()
plt.tight_layout()

plt.show()