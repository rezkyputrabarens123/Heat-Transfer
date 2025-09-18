function Natural_Convection_along_a_Heated_Plate
% Natural_Convection_along_a_Heated_Plate:
% Solves the boundary layer equations for natural convection along a
% vertical heated plate using the Blasius-type formulation with BVP solver.

Pr = [.07 .7 7];          % Prandtl numbers to consider
etaMax = [11, 8, 8];      % Maximum eta (similarity variable) for each Pr
xm = [10, 5, 5];          % x-axis limits for plotting
ym = [2, 0.8, 0.5];       % y-axis limits for stream function/velocity/shear
guess = [0 0 0 0 0];      % Initial guess for the BVP solver

for k = 1:3
    figure(k)
    
    % Initialize the boundary value problem with initial guess
    solinit = bvpinit(linspace(0, etaMax(k), 5), guess);
    
    % Solve the BVP using bvp4c
    sol = bvp4c(@NatConv, @NatConvBC, solinit, [], Pr(k));
    
    % Evaluate solution at fine grid points for plotting
    eta = linspace(0, etaMax(k), 300);
    y = deval(sol, eta);
    
    % --------------------------
    % Plot stream function, velocity, and shear
    % --------------------------
    subplot(2, 1, 1)
    plot(eta, y(1,:), '-.k', eta, y(2,:), '-k', eta, y(3,:), '--k')
    legend('Stream function, f = y_1', 'Velocity, df/d\eta = y_2',...
        'Shear, d^2f/d\eta^2 = y_3')
    axis([0 xm(k) -0.2 ym(k)])
    xlabel('\eta')
    ylabel('y_1, y_2, y_3')
    
    % --------------------------
    % Plot temperature and heat flux
    % --------------------------
    subplot(2, 1, 2)
    plot(eta, y(4,:), '-k', eta, y(5,:), '--k')
    legend('Temperature, T^* = y_4', 'Heat flux, dT^*/d\eta = y_5')
    axis([0 xm(k) -1.2 1])
    xlabel('\eta')
    ylabel('y_4, y_5')
end

% --------------------------
% BVP system for natural convection
% --------------------------
function ff = NatConv(eta, y, Pr)
    ff = [y(2);                 % dy1/deta = y2 (stream function derivative)
          y(3);                 % dy2/deta = y3 (velocity derivative)
          -3*y(1)*y(3)+2*y(2)^2-y(4); % dy3/deta = nonlinear term (shear eq.)
          y(5);                 % dy4/deta = y5 (temperature derivative)
          -3*Pr*y(1)*y(5)];    % dy5/deta = -3*Pr*f*dT/deta (heat eq.)

% --------------------------
% Boundary conditions for BVP
% --------------------------
function res = NatConvBC(ya, yb, Pr)
    res = [ya(1);      % f(0) = 0 (no slip at plate)
           ya(2);      % df/deta(0) = 0 (zero velocity at wall)
           ya(4)-1;    % T*(0) = 1 (plate temperature)
           yb(2);      % df/deta(inf) = 0 (velocity goes to 0)
           yb(4)];     % T*(inf) = 0 (ambient temperature)