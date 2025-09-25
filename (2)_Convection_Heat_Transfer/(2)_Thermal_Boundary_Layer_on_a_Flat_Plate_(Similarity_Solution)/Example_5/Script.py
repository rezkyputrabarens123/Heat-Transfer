import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

def Blasius_Formulation():
    """
    Python equivalent of the MATLAB Blasius Formulation code.
    Solves the coupled Blasius momentum and thermal boundary layer equations
    for different Prandtl numbers.
    """
    Pr = [0.07, 0.7, 7.0]  # Prandtl numbers
    eta_max = [15, 8, 8]   # Maximum eta for each Pr
    x_lim = [15, 5, 5]     # x-axis limit for plotting

    for k in range(3):
        # Initial guess for BVP solution
        eta_init = np.linspace(0, eta_max[k], 8)
        y_init = np.zeros((5, eta_init.size))

        # Solve boundary value problem
        # The `fun` and `bc` functions are defined below.
        sol = solve_bvp(lambda eta, y: blasius_ode(eta, y, Pr[k]),
                        lambda ya, yb: blasius_bc(ya, yb),
                        eta_init, y_init)

        # Evaluate solution on a fine eta grid
        eta = np.linspace(0, eta_max[k], 100)
        y = sol.sol(eta)

        fig, axs = plt.subplots(2, 1, num=k + 1)
        fig.suptitle(f'Prandtl Number (Pr) = {Pr[k]}')

        # --- Plot momentum boundary layer results ---
        ax1 = axs[0]
        ax1.plot(eta, y[0, :], '-.k', label='Stream function, $f = y_1$')
        ax1.plot(eta, y[1, :], '-k', label='Velocity, $df/d\eta = y_2$')
        ax1.plot(eta, y[2, :], '--k', label='Shear, $d^2f/d\eta^2 = y_3$')
        ax1.set_xlabel('$\eta$')
        ax1.set_ylabel('$y_1, y_2, y_3$')
        ax1.set_title('Momentum Boundary Layer')
        ax1.legend()
        ax1.set_xlim([0, x_lim[k]])
        ax1.set_ylim([0, 2])
        ax1.grid(True)

        # --- Plot thermal boundary layer results ---
        ax2 = axs[1]
        ax2.plot(eta, y[3, :], '-k', label='Temperature, $T^* = y_4$')
        ax2.plot(eta, y[4, :], '--k', label='Heat flux, $dT^*/d\eta = y_5$')
        ax2.set_xlabel('$\eta$')
        ax2.set_ylabel('$y_4, y_5$')
        ax2.set_title('Thermal Boundary Layer')
        ax2.legend()
        ax2.set_xlim([0, x_lim[k]])
        ax2.set_ylim([0, 2])
        ax2.grid(True)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.show()

# --- ODE system function for solve_bvp ---
def blasius_ode(eta, y, Pr):
    """
    Defines the system of ODEs.
    y[0]: f, y[1]: f', y[2]: f''
    y[3]: T*, y[4]: T*'
    """
    f, fp, fpp, T_star, T_star_p = y
    
    dydeta = np.vstack([
        fp,         # dy1/d(eta) = y2
        fpp,        # dy2/d(eta) = y3
        -0.5 * f * fpp,   # dy3/d(eta) = -0.5 * y1 * y3 (Momentum equation)
        T_star_p,   # dy4/d(eta) = y5
        -Pr * 0.5 * f * T_star_p # dy5/d(eta) = -Pr * 0.5 * y1 * y5 (Energy equation)
    ])
    return dydeta

# --- Boundary conditions function for solve_bvp ---
def blasius_bc(ya, yb):
    """
    Defines the boundary conditions.
    ya: solution at eta = 0
    yb: solution at eta = eta_max
    """
    return np.array([
        ya[0],      # f(0) = 0
        ya[1],      # f'(0) = 0
        ya[3],      # T*(0) = 0
        yb[1] - 1,  # f'(inf) = 1
        yb[3] - 1   # T*(inf) = 1
    ])

if __name__ == '__main__':
    Blasius_Formulation()