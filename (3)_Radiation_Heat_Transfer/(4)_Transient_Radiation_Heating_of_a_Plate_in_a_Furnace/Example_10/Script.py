import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

def transient_radiation_heating():
    """
    This script models the transient heating of a plate inside a furnace due to radiation.
    It determines the heat transfer rate (Q) required to bring the plate to a 
    specific target temperature (Te) at a specific time (th).
    """
    # --- Define physical constants ---
    # These constants are specific to the problem setup (e.g., thermal mass, 
    # surface area, Stefan-Boltzmann constant).
    P1 = 1.67e-5  # Constant for wall temperature change rate, related to heat input Q.
    P2 = 8.8e-14   # Constant for wall temperature change rate, related to radiation from the plate.
    P3 = 6.3e-13   # Constant for plate temperature change rate, related to radiation from the wall.

    # --- Define problem parameters and initial conditions ---
    Qguess = 1e5   # Initial guess for the total heat transfer rate (W) to the furnace wall.
    Te = 1100      # Target plate temperature (K).
    th = 600       # Target time (s) to reach the target temperature.
    tend = 660     # Total duration of the simulation (s).
    Two = 300      # Initial wall temperature (K).
    Tpo = 300      # Initial plate temperature (K).

    # --- Part 1: Find the required heat transfer rate Q ---
    # Use 'root_scalar' to find the value of Q that makes the plate temperature at time 'th' 
    # equal to 'Te'. This is equivalent to finding the root of the function QGen.
    
    # Bundle additional arguments for the QGen function
    args_for_QGen = (Te, th, Two, Tpo, P1, P2, P3)
    
    # Find the root (the correct Q value)
    solution = root_scalar(QGen, args=args_for_QGen, x0=Qguess, method='secant')
    if not solution.converged:
        print("Warning: Root-finding for Q did not converge.")
    Q = solution.root
    print(f"Calculated Heat Transfer Rate (Q): {Q:.2f} W")

    # --- Part 2: Solve the ODEs with the calculated Q and plot ---
    # Use 'solve_ivp' (the modern replacement for ode45) to solve the system of ODEs.
    t_span = [0, tend]
    y0 = [Two, Tpo]
    args_for_RadTemp = (Q, P1, P2, P3)
    
    # Solve the ODE system for the full duration
    sol = solve_ivp(
        RadTemp, 
        t_span, 
        y0, 
        args=args_for_RadTemp, 
        dense_output=True,
        t_eval=np.linspace(0, tend, 500) # Points to evaluate for a smooth plot
    )
    
    t = sol.t
    T_wall = sol.y[0]
    T_plate = sol.y[1]

    # --- Part 3: Plot the results ---
    plt.figure(figsize=(10, 6))
    plt.plot(t, T_wall, 'k-', label='Wall temperature')
    plt.plot(t, T_plate, 'k--', label='Plate temperature')
    
    # Add horizontal and vertical lines for target temperature and time
    plt.axhline(y=Te, color='k', linestyle=':', linewidth=1)
    plt.axvline(x=th, color='k', linestyle=':', linewidth=1)
    
    # Add labels and a legend
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.title('Transient Radiation Heating of a Plate')
    plt.legend(loc='upper left')
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Add a text box showing the calculated Q value
    plt.text(0.05, 0.85, f'Q = {Q:.2f} W', transform=plt.gca().transAxes,
             bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.7))

    plt.show()


def RadTemp(t, T, Q, P1, P2, P3):
    """
    Defines the system of ordinary differential equations.
    
    Args:
        t (float): Time.
        T (list/array): A list [T_wall, T_plate] of current temperatures.
        Q (float): Heat transfer rate.
        P1, P2, P3 (float): Physical constants.
        
    Returns:
        list: The rate of change [dT_wall/dt, dT_plate/dt].
    """
    Tw, Tp = T[0], T[1]  # Unpack wall and plate temperatures
    
    dTdt_wall = P1 * Q - P2 * (Tw**4 - Tp**4)
    dTdt_plate = -P3 * (Tp**4 - Tw**4)
    
    return [dTdt_wall, dTdt_plate]


def QGen(Q, Te, th, Two, Tpo, P1, P2, P3):
    """
    Objective function for the root-finder ('root_scalar').
    It calculates the difference between the target plate temperature (Te) and
    the actual plate temperature at time 'th' for a given heat input Q.
    The root-finder adjusts Q until this difference is zero.
    """
    # We only need to integrate up to the target time 'th'
    t_span = [0, th]
    y0 = [Two, Tpo]
    args_for_RadTemp = (Q, P1, P2, P3)
    
    # Solve the ODE system
    sol = solve_ivp(
        RadTemp, 
        t_span, 
        y0, 
        args=args_for_RadTemp, 
        dense_output=True # Needed to accurately find temperature at 'th'
    )
    
    # Evaluate the solution at the precise target time 'th'
    # sol.sol(th) returns [T_wall(th), T_plate(th)]
    plate_temp_at_th = sol.sol(th)[1]
    
    # Return the deviation that the root-finder will drive to zero
    return Te - plate_temp_at_th

# --- Run the main function ---
if __name__ == "__main__":
    transient_radiation_heating()