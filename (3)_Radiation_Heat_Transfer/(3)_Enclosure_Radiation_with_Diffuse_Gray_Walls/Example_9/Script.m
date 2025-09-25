% Total heat transfer rate of a rectangular enclosure
% This script calculates the net radiative heat transfer rate for each surface of a rectangular enclosure.
% It assumes the surfaces are gray, diffuse emitters and reflectors, and a non-participating medium.

% Stefan-Boltzmann constant (W/m^2*K^4)
sigma = 5.6693e-8;

% Surface properties of the four sides of the rectangular enclosure.
% A: Area of each surface (m^2). The dimensions suggest a 3m x 5m enclosure, with two sides of 3m^2 and two of 5m^2.
% epsilon: Emissivity of each surface.
% T: Absolute temperature of each surface (K).
A = [3, 5, 3, 5];
epsilon = [0.7, 0.3, 0.85, 0.45];
T = [550, 700, 650, 600];

% View factors matrix (F_ij) for the enclosure.
% This is a 4x4 matrix where F_ij is the fraction of radiation leaving surface i that is incident on surface j.
% The negative signs are a convention used in this specific matrix formulation.
F = -[0, 0.3615, 0.277, 0.3615;...
0.2169, 0, 0.2169, 0.5662;...
0.277, 0.3615, 0, 0.3615;...
0.2169, 0.5662, 0.2169, 0];

% Initialize vectors for net heat transfer rate and a control vector.
% Q: Net heat transfer rate for each surface (W).
% c: A control vector, initialized to zeros. It's often used to handle different boundary conditions (e.g., prescribed temperature vs. prescribed heat flux).
% In this case, since it's all zeros, it assumes all surfaces have a prescribed temperature.
Q = zeros(1, length(A));
c = zeros(1, length(A));

% Calculate vector 'b' for the system of linear equations.
% This vector represents the radiative power from each surface.
% The expression is derived from the radiosity method for gray, diffuse surfaces with known temperatures.
b = sigma*epsilon./(1-epsilon).*(1-c).*T.^4+c.*Q./A;

% Calculate vector 'd' for the system of linear equations.
% This vector modifies the diagonal elements of the view factor matrix.
% It accounts for surface properties (emissivity) and boundary conditions.
d = (1-c).*1./(1-epsilon)+c;

% Modify the view factor matrix 'F' by adding 'd' to the diagonal.
% This step constructs the coefficient matrix for the linear system.
F = F+diag(d);

% Solve the linear system to find the radiosities (q0).
% q0: Radiosity of each surface (W/m^2). Radiosity is the total radiative energy leaving a surface per unit area per unit time.
% The equation is of the form F * q0 = b, so q0 = F^-1 * b.
q0 = F\b';

% Calculate the net radiative heat transfer rate for each surface (Q).
% Q is the difference between the emitted and incident radiation.
% This is derived from the radiosity equation, where Q = (epsilon*A / (1-epsilon)) * (sigma*T^4 - q0).
Q = A.*epsilon./(1-epsilon).*(1-c).*(sigma*T.^4-q0')

% Calculate the net heat flux for each surface (q).
% q: Net heat flux (W/m^2), which is the heat transfer rate per unit area.
q = Q./A