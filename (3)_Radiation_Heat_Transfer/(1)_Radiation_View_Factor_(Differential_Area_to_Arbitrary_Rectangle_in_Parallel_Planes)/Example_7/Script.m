% Differential_Area_and_a_Finite_Rectangle_in_Parallel_Planes
%
% This script calculates and plots the view factor between a differential
% area and a finite rectangle, both of which are located on parallel planes.
%
% The script uses two main functions:
% - Fd1_2: This function calculates the view factor F_{d1-2} between a
%          differential area dA1 and a finite rectangle A2.
% - kernel2: This is a sub-function used by Fd1_2 to compute the integrand
%            for the double integral.
%
% The view factor (F_{d1-2}) is a dimensionless quantity that represents
% the fraction of radiation leaving surface dA1 that strikes surface A2.
% For parallel planes, the view factor depends on the relative positions
% and separation distance of the surfaces.
%
% The script performs the following steps:
% 1. Calculates the view factor for two specific rectangle configurations
%    and prints the results to the command window.
% 2. Calculates the view factor for a range of separation distances and
%    plots the results to show how the view factor changes with distance.
%
% References:
% - Incropera, F. P., Dewitt, D. P., Bergman, T. L., & Lavine, A. S. (2007).
%   Fundamentals of Heat and Mass Transfer. John Wiley & Sons.
% - A_{i}-F_{i-j} View Factors for Common Geometries. (n.d.).
%   Retrieved from .

function Differential_Area_and_a_Finite_Rectangle_in_Parallel_Planes

% --- Step 1: Calculate view factors for specific cases ---
%
% Case 1: Rectangle located from x = -1 to 0 and y = -1 to 0,
%         at a separation distance of 5.
Set1 = Fd1_2(-1, 0, -1, 0, 5);
fprintf('View factor for case 1: Fd1-2 = %f\n', Set1);

% Case 2: Rectangle located from x = -1 to 1 and y = -1 to 1,
%         at a separation distance of 1.
Set2 = Fd1_2(-1, 1, -1, 1, 1);
fprintf('View factor for case 2: Fd1-2 = %f\n', Set2);

% --- Step 2: Plot view factor vs. separation distance ---
%
% Calculate the view factor for a rectangle (x=0 to 1, y=0 to 1) as the
% separation distance (dz) varies from 0.1 to 5.
N = 100;
dz = linspace(0.1, 5, N);
Fd12 = zeros(N,1);
for i = 1:N
    Fd12(i) = Fd1_2(0, 1, 0, 1, dz(i));
end

% Plot the results.
plot(dz, Fd12, 'k-')
xlabel('Separation distance of surfaces')
ylabel('View factor')
title('View Factor between a Differential Area and a Rectangle')

end % End of main function

% =========================================================================
% Function Definitions
% =========================================================================

function F = Fd1_2(x_2a, x_2b, y_2a, y_2b, dz)
% Fd1_2 calculates the view factor F_d1-2 from a differential area dA1 to a
% finite rectangle A2 on a parallel plane.
%
% The differential area dA1 is assumed to be at the origin (0,0) of a
% coordinate system. The finite rectangle A2 is defined by the coordinates
% (x_2a, x_2b) and (y_2a, y_2b) on a parallel plane at a separation distance dz.
%
% Inputs:
%   x_2a: Starting x-coordinate of the rectangle
%   x_2b: Ending x-coordinate of the rectangle
%   y_2a: Starting y-coordinate of the rectangle
%   y_2b: Ending y-coordinate of the rectangle
%   dz:   Separation distance between the planes
%
% Output:
%   F:    The view factor F_d1-2
%
% The view factor is calculated using the formula:
% F_{d1-2} = (1/pi) * integral_A2 [ (cos(theta_1) * cos(theta_2)) / r^2 ] dA2
% where dA2 = dx dy.
% Since dA1 and A2 are on parallel planes, the formula simplifies.
% The integrand is calculated by the sub-function 'kernel2'.
% dblquad is used for the numerical double integration.

% Integrate the kernel function over the rectangle's area.
F = dblquad(@kernel2, x_2a, x_2b, y_2a, y_2b, [], [], dz) / pi;

end % End of Fd1_2 function

% =========================================================================

function f = kernel2(x, y, dist)
% kernel2 is a sub-function that computes the integrand for the view factor
% calculation.
%
% The integrand represents the geometric relationship between the
% differential area and each point on the finite rectangle.
%
% The formula for the integrand is:
% (cos(theta_1) * cos(theta_2)) / r^2 = (dz^2) / (x^2 + y^2 + dz^2)^2
% This is derived from the geometric relationships between the two surfaces.
%
% A more general form is used here which applies the dot product between
% vectors. This form is more robust for general cases but simplifies for
% parallel planes.
%
% Inputs:
%   x:    x-coordinate vector of points on the rectangle
%   y:    y-coordinate vector of points on the rectangle
%   dist: Separation distance (dz)
%
% Output:
%   f:    The value of the integrand at each point (x, y)

% L is the number of points at which the integrand is calculated.
L = length(x);

% S is the vector from the differential area to a point on the rectangle.
% S = [x, y, dz]
S = [x; repmat(y, 1, L); dist*ones(1, L)];

% n is the normal vector of the differential area (pointing along the z-axis).
% n = [0, 0, 1]
n = repmat([0, 0, 1]', 1, L);

% Calculate the integrand value.
% The formula used here is dot(n, S)^2 / dot(S, S)^2, which simplifies to
% (dz^2) / (x^2 + y^2 + dz^2)^2 for parallel planes.
f = dot(n, S).^2 ./ dot(S, S).^2;

end % End of kernel2 function