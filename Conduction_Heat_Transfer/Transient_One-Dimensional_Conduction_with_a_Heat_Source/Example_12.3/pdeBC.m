function [pl, ql, pr, qr] = pdeBC(xl, ul, xr, ur, t, Bi, Tr, Sigma)
    % pdeBC: Defines boundary conditions for the PDE solver (pdepe)
    %
    % General PDE form in pdepe requires:
    %   p(x,t,u) + q(x,t) * f(x,t,u,DuDx) = 0
    %
    % Inputs:
    %   xl, ul - Left boundary location and solution value
    %   xr, ur - Right boundary location and solution value
    %   t      - Time (?)
    %   Bi     - Biot number
    %   Tr     - Reference temperature ratio (applied at right boundary)
    %   Sigma  - Source/sink parameter (not used here, but kept for consistency)
    %
    % Outputs:
    %   pl, ql - Left boundary condition coefficients
    %   pr, qr - Right boundary condition coefficients
    %
    % Boundary conditions implemented:
    %   Left (x = 0):  ?'(0,t) = Bi * ?(0,t)
    %       ? pl = -Bi*ul, ql = 1
    %
    %   Right (x = 1): ?(1,t) = Tr (fixed temperature)
    %       ? pr = ur - Tr,  qr = 0

    % Right boundary (Dirichlet condition: ? = Tr)
    pr = ur - Tr;   % Enforces ?(x=1,t) = Tr
    qr = 0;

    % Left boundary (Robin/convective condition: ?' = Bi*?)
    pl = -Bi * ul;  % Residual from convective BC
    ql = 1;         % Coupled with flux term f = DuDx
end