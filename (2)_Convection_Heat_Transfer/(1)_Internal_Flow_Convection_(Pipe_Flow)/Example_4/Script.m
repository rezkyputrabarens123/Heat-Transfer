% Heat Transfer Coefficient for Laminar Flow in a Pipe
% ----------------------------------------------------
% This program calculates the Nusselt number (Nu) and temperature profiles
% for fully developed laminar flow inside a pipe under two boundary
% conditions:
%   1. Constant wall temperature (isothermal wall)
%   2. Constant wall heat flux
% using MATLAB's PDE solver (pdepe).

% Main function
function Heat_Transfer_Coefficient_for_Laminar_Flow_in_a_Pipe

% -------------------------------
% Parameters
% -------------------------------
Tw = 40;          % Wall temperature [°C] (for constant Tw condition)
qw = 10;          % Wall heat flux [W/m^2] (for constant q condition)
Re = 40;          % Reynolds number (laminar flow, Re < 2300)
Pr = 5;           % Prandtl number
R  = 0.01;        % Pipe radius [m]
L  = 0.5;         % Pipe length [m]
k  = 0.6;         % Thermal conductivity of fluid [W/m-K]

Rt = 401;         % Number of radial grid points
zt = 50;          % Number of axial grid points
dxi = 1/(Rt-1);   % Radial step size (normalized radius)
xi = linspace(0, 1, Rt);   % Dimensionless radial coordinate (0=center, 1=wall)
zeta = linspace(0, 1, zt); % Dimensionless axial coordinate (0=inlet, 1=outlet)

% -------------------------------
% PDE Solution for Both Conditions
% -------------------------------
% solT -> solution for constant wall temperature
% solF -> solution for constant wall heat flux
solT = pdepe(1, @pdepde, @pdeic, @pdebcT, xi, zeta, [], Tw, qw, Re, Pr, R, L, k);
solF = pdepe(1, @pdepde, @pdeic, @pdebcF, xi, zeta, [], Tw, qw, Re, Pr, R, L, k);

% Arrays to store Nusselt numbers along zeta
NuT = zeros(zt,1);
NuF = zeros(zt,1);

% -------------------------------
% Chart 1: Temperature Distribution (Outlet Profile)
% -------------------------------
figure(1)
for i = 1:zt
    % -------- Constant Wall Temperature --------
    % Bulk mean temperature
    TmT = 4*trapz(xi, xi.*(1-xi.^2).*solT(i,:));
    % Temperature gradient near the wall
    dThdxiT = (solT(i,Rt)-solT(i,Rt-1)) / (dxi*(TmT-solT(i,Rt)));
    % Nusselt number definition
    NuT(i) = -2*dThdxiT;

    % -------- Constant Wall Heat Flux --------
    TmF = 4*trapz(xi, xi.*(1-xi.^2).*solF(i,:));
    dThdxiF = (solF(i,Rt)-solF(i,Rt-1)) / (dxi*(TmF-solF(i,Rt)));
    NuF(i) = -2*dThdxiF;
end

% Dimensionless temperature profile at outlet (zeta = 1)
ThT = (solT(end,:)-ones(1,Rt)*Tw)/(TmT-Tw);
ThF = (solF(end,:)-ones(1,Rt)*solF(end,Rt))/(TmF-solF(end,Rt));

% Plot dimensionless temperature distribution
plot(xi, ThT, 'k-', xi, ThF, 'k--')
xlabel('\xi (r/R)')
ylabel('\theta (dimensionless T)')
legend('Constant wall temperature', 'Constant wall heat flux')

% -------------------------------
% Chart 2: Nusselt Number along Pipe
% -------------------------------
figure(2)
plot(zeta, NuT, 'k-', zeta, NuF, 'k--')
xlabel('\zeta (z/L)')
ylabel('Nu (Nusselt number)')
ylim([0 6])
legend('Constant wall temperature', 'Constant wall heat flux')

% =====================================================
% Nested Functions for PDE Solver
% =====================================================

% PDE definition
% Equation form: c*dT/dz = d/dxi(f*dT/dxi) + s
function [c, f, s] = pdepde(xi, zeta, T, DTDxi, Tw, qw, Re, Pr, R, L, k)
    c = Re*Pr*R/L*(1-xi^2);  % coefficient (includes velocity profile)
    f = DTDxi;               % diffusion term
    s = 0;                   % no source term
end

% Initial condition (inlet temperature profile)
function T0 = pdeic(xi, Tw, qw, Re, Pr, R, L, k)
    T0 = 20; % Uniform initial fluid temperature [°C]
end

% Boundary condition: Constant wall temperature
function [pl, ql, pr, qr] = pdebcT(xil, Tl, xir, Tr, zeta, Tw, qw, Re, Pr, R, L, k)
    pl = 0;     % symmetry at centerline
    ql = 1;
    pr = Tr-Tw; % wall temperature fixed at Tw
    qr = 0;
end

% Boundary condition: Constant wall heat flux
function [pl, ql, pr, qr] = pdebcF(xil, Tl, xir, Tr, zeta, Tw, qw, Re, Pr, R, L, k)
    pl = 0;     % symmetry at centerline
    ql = 1;
    pr = -qw;   % impose constant heat flux
    qr = k/R;   % thermal conductivity / radius
end

end % End of main function