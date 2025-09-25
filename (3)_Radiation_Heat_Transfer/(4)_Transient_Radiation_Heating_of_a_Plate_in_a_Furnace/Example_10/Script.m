% Transient radiation heating of a plate in a furnace
% This script models the transient heating of a plate inside a furnace due to radiation.
% It determines the heat transfer rate (Q) required to bring the plate to a specific target temperature (Te) at a specific time (th).
% The problem is solved by using a numerical solver (fzero) to find the correct Q, and then a second solver (ode45) to plot the temperature profiles over time.

% Main function to run the simulation
function Transient_Radiation_Heating_of_a_Plate_in_a_Furnace

% Define physical constants related to heat transfer and system properties.
% These constants are specific to the problem setup (e.g., thermal mass, surface area, Stefan-Boltzmann constant).
P1 = 1.67e-5;  % Constant for wall temperature change rate, related to heat input Q.
P2 = 8.8e-14;   % Constant for wall temperature change rate, related to radiation from the plate.
P3 = 6.3e-13;   % Constant for plate temperature change rate, related to radiation from the wall.

% Define problem parameters and initial conditions.
Qguess = 1e5;   % Initial guess for the total heat transfer rate (W) to the furnace wall.
Te = 1100;      % Target plate temperature (K).
th = 600;       % Target time (s) to reach the target temperature.
tend = 660;     % Total duration of the simulation (s).
Two = 300;      % Initial wall temperature (K).
Tpo = 300;      % Initial plate temperature (K).

% Use 'fzero' to find the value of Q that makes the plate temperature at time 'th' equal to 'Te'.
% The 'QGen' function serves as the objective function for 'fzero'.
% It finds the root of the equation T_plate(th) - T_e = 0, where T_plate(th) is the plate temperature at time th.
Q = fzero(@QGen, Qguess, [], Te, th, Two, Tpo, tend, P1, P2, P3);

% Use 'ode45' to solve the system of ordinary differential equations (ODEs) for the given Q.
% 'RadTemp' is the function that defines the system of ODEs.
% The output 't' is a vector of time points, and 'T' is a matrix where T(:,1) is the wall temperature and T(:,2) is the plate temperature.
[t, T] = ode45(@RadTemp, [0, tend], [Two; Tpo], [], Q, Te, th, Two, Tpo, tend, P1, P2, P3);

% Plot the temperature profiles over time.
plot(t, T(:,1), 'k-', t, T(:,2), 'k--')

% Adjust the plot axes for better visualization and add markers.
z = axis; % Get the current axis limits.
hold on % Keep the current plot and add new plots on top.
plot([0, z(2)], [Te, Te], 'k.:', [th, th], [z(3), z(4)], 'k.:') % Plot dashed lines for target temperature and time.

% Add labels and a legend to the plot.
xlabel('Time (s)')
% Add a text label showing the calculated Q value.
text(0.05*z(2), 0.85*z(4), ['Q = ' num2str(Q,6) ' W'])
ylabel('Temperature (K)')
legend('Wall temperature', 'Plate temperature', 'Location', 'NorthWest')

% --- Nested Functions ---
% These functions are defined within the main function and have access to its variables.

% RadTemp: Defines the system of ordinary differential equations.
% It calculates the rate of change of temperature for the wall (dT(1)/dt) and the plate (dT(2)/dt).
% The equations are based on a simple energy balance considering heat input (Q) and radiative exchange between the wall and the plate.
function dTdt = RadTemp(t, T, Q, Te, th, Two, Tpo, tend, P1, P2, P3)
    % T(1) is the wall temperature, T(2) is the plate temperature.
    % dTdt(1) is the rate of change of wall temperature.
    % dTdt(2) is the rate of change of plate temperature.
    dTdt = [P1*Q-P2*(T(1)^4-T(2)^4); -P3*(T(2)^4-T(1)^4)];

% QGen: Objective function for 'fzero'.
% It solves the ODEs for a given Q and returns the difference between the target temperature and the actual plate temperature at time 'th'.
% 'fzero' will iteratively call this function, adjusting Q until this difference is approximately zero.
function PlateTempDev = QGen(Q, Te, th, Two, Tpo, tend, P1, P2, P3)
    % Solve the ODEs for the current guess of Q.
    [t,T] = ode45(@RadTemp, [0, tend], [Two;Tpo], [], Q,Te, th,Two,Tpo, tend, P1, P2, P3);
    % Interpolate to find the plate temperature at the exact time 'th'.
    % The 'spline' method is used for smooth interpolation.
    % PlateTempDev is the error term that 'fzero' tries to minimize.
    PlateTempDev = Te-interp1(t, T(:,2), th, 'spline');