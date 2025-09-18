function Example_3
    % Example_3: Solves transient heat conduction in a 1D slab
    % using MATLAB's PDE solver (pdepe) and compares with steady-state solution (bvp4c).

    % ---------------------------
    % PARAMETERS
    % ---------------------------
    Bi    = 0.1;   % Biot number (dimensionless convective heat transfer coefficient)
    Tr    = 0.55;  % Reference temperature ratio at the right boundary
    Sigma = 1;     % Source term for heat generation (dimensionless)

    % ---------------------------
    % SPATIAL AND TEMPORAL GRID
    % ---------------------------
    xi  = linspace(0, 1, 41);   % Spatial discretization: 41 points from 0 to 1
    tau = linspace(0, 1, 101);  % Time discretization: 101 points from 0 to 1

    % ---------------------------
    % SOLVE TRANSIENT PDE
    % ---------------------------
    % pdepe solves parabolic PDEs of the form:
    %   c*u_t = x^-m * d/dx(x^m * f) + s
    theta = pdepe(0, @pde1D, @pdeIC, @pdeBC, xi, tau, [], Bi, Tr, Sigma);

    % ============================================================
    % FIGURE 1: Transient temperature at selected spatial points
    % ============================================================
    z = 0:0.25:1;   % Selected positions along the slab (xi)
    figure(1)
    hold on
    for k = 1:length(z)
        % Find index of xi closest to z(k)
        [~, kk] = min(abs(xi - z(k)));

        % Plot temperature vs time at this position
        plot(tau, theta(:, kk), 'k-')

        % Add curve labels at specific points
        if k == 1
            text(0.5, 1.02*theta(end, kk), '\xi = 0.0 and 0.25')
        elseif k > 2
            text(0.5, theta(end, kk) + 0.02, ['\xi = ' num2str(xi(kk))])
        end
    end
    hold off
    axis([0 1 0.5 1])
    xlabel('\tau')
    ylabel('\theta')
    title('Transient Response at Selected Positions')

    % ============================================================
    % FIGURE 2: Spatial temperature profile at selected times
    % ============================================================
    figure(2)
    % Find time index when temperature at xi=0 is minimum
    [thmin, imin] = min(theta(:, 1));

    % Plot initial condition
    plot(xi, theta(1,:), 'k-', 'LineWidth', 2)
    hold on
    % Plot at tau = tau(2) and at minimum temperature time
    plot(xi, theta(2,:), 'k', xi, theta(imin,:), 'k:')

    % ---------------------------
    % STEADY-STATE SOLUTION
    % ---------------------------
    % Define initial guess for bvp4c
    solinit = bvpinit(linspace(0, 1, 20), [1 1]); 
    sol     = bvp4c(@barode, @barbc, solinit, [], Bi, Sigma);

    % Evaluate steady-state solution
    x = linspace(0, 1, 100);
    y = deval(sol, x);

    % Plot steady-state profile
    plot(x, y(1,:), 'k--');

    % Formatting
    xlabel('\xi')
    ylabel('\theta')
    legend(['\tau = 0 (Initial condition)'], ...
           ['\tau = ' num2str(tau(2))], ...
           ['\tau = ' num2str(tau(imin)) ' (Minimum at \xi = 0)'], ...
           '\tau > 2 (Steady state)', ...
           'Location', 'SouthWest')
    title('Spatial Profiles at Different Times')
    hold off
end

% ============================================================
% FUNCTIONS DEFINING STEADY-STATE ODE PROBLEM
% ============================================================
function dydx = barode(x, y, Bi, Sigma)
    % barode: Defines ODEs for steady-state heat conduction
    % y(1) = temperature, y(2) = dT/dx
    dydx = [ y(2);        % dT/dx
            -Sigma ];    % d^2T/dx^2 = -Sigma
end

function res = barbc(ya, yb, Bi, Sigma)
    % barbc: Boundary conditions for steady-state problem
    % Left BC: -Bi*T + dT/dx = 0
    % Right BC: T = Tr (0.55)
    res = [ ya(2) - Bi * ya(1);
            yb(1) - 0.55 ];
end

% ============================================================
% FUNCTIONS DEFINING TRANSIENT PDE
% ============================================================
function [c, f, s] = pde1D(x, t, u, DuDx, Bi, Tr, Sigma)
    % pde1D: PDE coefficients for pdepe
    c = 1;        % Coefficient of time derivative (u_t)
    f = DuDx;     % Flux term (heat conduction)
    s = Sigma;    % Source term (internal heat generation)
end

function T0 = pdeIC(x, Bi, Tr, Sigma)
    % pdeIC: Initial temperature profile along the slab
    T0 = 1 - 0.45 * x;  % Linear initial temperature profile
end

function [pl, ql, pr, qr] = pdeBC(xl, ul, xr, ur, t, Bi, Tr, Sigma)
    % pdeBC: Boundary conditions for transient PDE
    % Left boundary (convection): -Bi*ul + dT/dx = 0
    % Right boundary (Dirichlet): T = Tr
    pr = ur - Tr;
    qr = 0;
    pl = -Bi * ul;
    ql = 1;
end