import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

def natural_convection_along_a_heated_plate():
    """
    Solves and plots the boundary layer equations for natural convection 
    along a vertical heated plate using the Blasius-type formulation 
    with a BVP solver.
    """
    # --- Problem Parameters ---
    pr_numbers = [0.07, 0.7, 7]  # Prandtl numbers to consider
    eta_max = [11, 8, 8]        # Maximum eta (similarity variable) for each Pr
    x_lim = [10, 5, 5]          # x-axis limits for plotting
    y_lim_a = [2, 0.8, 0.5]     # y-axis limits for stream function/velocity/shear
    
    # Loop through each Prandtl number
    for i, pr in enumerate(pr_numbers):
        
        # --- Solve the Boundary Value Problem ---
        
        # 1. Define the initial mesh (eta values)
        # An initial mesh of 5 points is sufficient.
        eta_init = np.linspace(0, eta_max[i], 5)
        
        # 2. Define the initial guess for the solution variables [y1, y2, y3, y4, y5]
        # A guess of all zeros works well. The shape must be (num_equations, num_mesh_points)
        y_init = np.zeros((5, eta_init.size))
        
        # 3. Solve the BVP using SciPy's solve_bvp
        # We use lambda functions to pass the Prandtl number `pr` to our ODE and BC functions.
        sol = solve_bvp(
            lambda eta, y: odes(eta, y, pr),
            lambda ya, yb: bc(ya, yb, pr),
            eta_init,
            y_init
        )
        
        # --- Prepare Data for Plotting ---
        
        # Create a fine grid of eta points for a smooth plot
        eta_plot = np.linspace(0, eta_max[i], 300)
        
        # Evaluate the solution on the fine grid
        y_plot = sol.sol(eta_plot)
        
        # --- Create Plots ---
        fig = plt.figure(i + 1)
        fig.suptitle(f'Natural Convection for Prandtl Number (Pr) = {pr}', fontsize=14)
        
        # Plot 1: Stream function, velocity, and shear
        ax1 = plt.subplot(2, 1, 1)
        ax1.plot(eta_plot, y_plot[0], '-.', color='k', label=r'Stream function, $f = y_1$')
        ax1.plot(eta_plot, y_plot[1], '-', color='k', label=r'Velocity, $df/d\eta = y_2$')
        ax1.plot(eta_plot, y_plot[2], '--', color='k', label=r'Shear, $d^2f/d\eta^2 = y_3$')
        
        ax1.set_xlim(0, x_lim[i])
        ax1.set_ylim(-0.2, y_lim_a[i])
        ax1.set_xlabel(r'$\eta$')
        ax1.set_ylabel(r'$y_1, y_2, y_3$')
        ax1.legend()
        ax1.grid(True, linestyle=':', alpha=0.6)
        
        # Plot 2: Temperature and heat flux
        ax2 = plt.subplot(2, 1, 2)
        ax2.plot(eta_plot, y_plot[3], '-', color='k', label=r'Temperature, $T^* = y_4$')
        ax2.plot(eta_plot, y_plot[4], '--', color='k', label=r'Heat flux, $dT^*/d\eta = y_5$')
        
        ax2.set_xlim(0, x_lim[i])
        ax2.set_ylim(-1.2, 1.0)
        ax2.set_xlabel(r'$\eta$')
        ax2.set_ylabel(r'$y_4, y_5$')
        ax2.legend()
        ax2.grid(True, linestyle=':', alpha=0.6)
        
        plt.tight_layout(rect=[0, 0, 1, 0.95]) # Adjust layout to prevent title overlap

    plt.show()

def odes(eta, y, pr):
    """
    Defines the system of first-order Ordinary Differential Equations (ODEs).
    The system is derived from the boundary layer equations.
    y[0] = f
    y[1] = f'
    y[2] = f''
    y[3] = T*
    y[4] = T*'
    """
    return np.vstack([
        y[1],                              # dy1/deta = y2 (f')
        y[2],                              # dy2/deta = y3 (f'')
        -3 * y[0] * y[2] + 2 * y[1]**2 - y[3], # dy3/deta = -3*f*f'' + 2*(f')^2 - T*
        y[4],                              # dy4/deta = y5 (T*')
        -3 * pr * y[0] * y[4]             # dy5/deta = -3*Pr*f*T*'
    ])

def bc(ya, yb, pr):
    """
    Defines the boundary conditions for the BVP.
    'ya' represents the solution at the start of the interval (eta = 0).
    'yb' represents the solution at the end of the interval (eta = inf).
    The function should return the residuals, which the solver will attempt to make zero.
    """
    return np.array([
        ya[0],      # y1(0) = f(0) = 0
        ya[1],      # y2(0) = f'(0) = 0
        ya[3] - 1,  # y4(0) = T*(0) = 1
        yb[1],      # y2(inf) = f'(inf) = 0
        yb[3]       # y4(inf) = T*(inf) = 0
    ])

# --- Run the main function ---
if __name__ == "__main__":
    natural_convection_along_a_heated_plate()