% Transient Heat Conduction in an Infinite Solid Cylinder with Convection

Bi = 0.5;                     % Biot number (dimensionless)
Nroot = 15;                  % Number of roots to find

% Define the transcendental equation for root finding
CylinderRoots = @(x, Bi) x .* besselj(1, x) - Bi * besselj(0, x);

% Find the first Nroot roots of the equation using a helper function
r = FindZeros(CylinderRoots, Nroot, linspace(0, 50, 200), Bi);  % r is [Nroot x 1]

% Time domain (dimensionless time ?)
tau = linspace(0, 1.5, 20);  % tau is [1 x Nt]

% Create meshgrid for time and roots
[t, rt] = meshgrid(tau, r);  % t and rt are [Nroot x Nt]

% Transient exponential decay term
Fn = exp(-t .* rt.^2);       % Fn is [Nroot x Nt]

% Coefficients for each term in the series solution
cn = 2 * besselj(1, r) ./ (r .* (besselj(0, r).^2 + besselj(1, r).^2));  % [Nroot x 1]

% Expand coefficients across time domain
ccn = cn .* ones(1, length(tau));  % Broadcast cn across columns to match Fn

% Multiply coefficients with exponential decay
pro = ccn .* Fn;  % Element-wise multiplication, both [Nroot x Nt]

% Radial position (dimensionless radius ?)
rstar = linspace(0, 1, 20);  % [1 x Nr]

% Create meshgrid for radial position and roots
[R, rx] = meshgrid(rstar, r);  % R and rx are [Nroot x Nr]

% Bessel function evaluated at scaled radial positions
Jo = besselj(0, rx .* R);  % [Nroot x Nr]

% Compute temperature distribution ?(?, ?)
the = Jo' * pro;  % Jo' is [Nr x Nroot], pro is [Nroot x Nt] ? the is [Nr x Nt]

% Prepare meshgrid for plotting
[rr, tt] = meshgrid(rstar, tau);  % rr and tt are [Nt x Nr]

% Plot the temperature distribution
mesh(rr, tt, the')  % the' is [Nt x Nr] to match rr and tt
xlabel('\xi')       % Dimensionless radial position
ylabel('\tau')      % Dimensionless time
zlabel('\theta')    % Dimensionless temperature
view(49.5, -34)     % Set 3D view angle