function Blasius_Formulation

Pr = [0.07, 0.7, 7.0];        % Prandtl numbers
etaMax = [15, 8, 8];           % Maximum eta for each Pr
xm = [15, 5, 5];               % x-axis limit for plotting

for k = 1:3
    figure(k)
    
    % Initial guess for BVP solution
    solinit = bvpinit(linspace(0, etaMax(k), 8), [0, 0, 0, 0, 0]);
    
    % Solve boundary value problem for momentum and temperature
    sol = bvp4c(@(eta, y) BlasiusT(eta, y, Pr(k)), ...
                @(ya, yb) BlasiusTbc(ya, yb), solinit);
    
    % Evaluate solution on a fine eta grid
    eta = linspace(0, etaMax(k), 100);
    y = deval(sol, eta);
    
    % --- Plot momentum boundary layer results ---
    subplot(2, 1, 1)
    plot(eta, y(1,:), '-.k', eta, y(2,:), '-k', eta, y(3,:), '--k')
    xlabel('\eta')
    ylabel('y_1, y_2, y_3')
    legend('Stream function, f = y_1', 'Velocity, df/d\eta = y_2', ...
           'Shear, d^2f/d\eta^2= y_3')
    axis([0 xm(k) 0 2])
    
    % --- Plot thermal boundary layer results ---
    subplot(2, 1, 2)
    plot(eta, y(4,:), '-k', eta, y(5,:), '--k')
    axis([0 xm(k) 0 2])
    legend('Temperature, T^* = y_4', 'Heat flux, dT^*/d\eta = y_5')
    xlabel('\eta')
    ylabel('y_4, y_5')
end

% --- Nested ODE system ---
function F = BlasiusT(eta, y, Pr)
    F = [y(2); 
         y(3); 
         -0.5*y(1)*y(3); 
         y(5); 
         -Pr*0.5*y(1)*y(5)];
end

% --- Nested boundary conditions ---
function res = BlasiusTbc(ya, yb)
    res = [ya(1);     % f(0) = 0
           ya(2);     % f'(0) = 0
           ya(4);     % T*(0) = 0
           yb(2)-1;   % f'(?) = 1
           yb(4)-1];  % T*(?) = 1
end

end % main function