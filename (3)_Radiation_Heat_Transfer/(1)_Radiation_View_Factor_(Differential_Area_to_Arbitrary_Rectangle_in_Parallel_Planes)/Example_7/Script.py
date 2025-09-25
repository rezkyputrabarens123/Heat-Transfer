import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import dblquad

def Fd1_2(x_2a, x_2b, y_2a, y_2b, dz):
    """
    Calculates the view factor between a differential area and a finite rectangle.
    
    Args:
        x_2a, x_2b: x-coordinates defining the rectangle.
        y_2a, y_2b: y-coordinates defining the rectangle.
        dz: Separation distance.
        
    Returns:
        The view factor F.
    """
    # Define the kernel function as a local function.
    def kernel2(x, y):
        # The 'dist' parameter is captured from the outer function's scope (dz).
        # This is a key difference from MATLAB's nested functions.
        dist = dz
        
        # In Python, we use NumPy for vector operations.
        # Ensure x and y are numpy arrays for element-wise operations.
        x_arr = np.asarray(x)
        y_arr = np.asarray(y)
        
        # Calculate the position vectors S.
        S = np.vstack((x_arr, y_arr, dist * np.ones_like(x_arr)))
        
        # Define the normal vector n.
        n = np.array([0, 0, 1]).reshape(-1, 1)
        
        # Calculate the dot products.
        n_dot_S = np.dot(n.T, S)
        S_dot_S = np.sum(S * S, axis=0) # Equivalent to dot(S, S) in MATLAB
        
        # Calculate the kernel function.
        f = (n_dot_S**2) / (S_dot_S**2)
        
        # dblquad expects a single value, not a vector.
        # This reshaping and averaging ensures compatibility.
        return f.flatten()
        
    # Perform the double integration of the kernel over the rectangle.
    F, _ = dblquad(kernel2, x_2a, x_2b, y_2a, y_2b)
    
    # Return the final view factor divided by pi.
    return F / np.pi

if __name__ == '__main__':
    # Define the number of data points.
    N = 100
    
    # Create an array of separation distances.
    dz = np.linspace(0.1, 5, N)
    
    # Pre-allocate an array to store the calculated view factors.
    Fd12 = np.zeros(N)
    
    # Loop through each separation distance to calculate the view factor.
    for i in range(N):
        Fd12[i] = Fd1_2(0, 1, 0, 1, dz[i])
    
    # Plot the results.
    plt.plot(dz, Fd12, 'k-')
    plt.xlabel('Separation distance of surfaces')
    plt.ylabel('View factor')
    plt.title('View Factor for a Differential Area and Finite Rectangle')
    plt.grid(True)
    plt.show()