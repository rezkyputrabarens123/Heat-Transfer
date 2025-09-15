function [c, f, s] = pde1D(x, t, u, DuDx, Bi, Tr, Sigma)
    % pde1D: Defines the PDE for transient 1D heat conduction
    %
    % PDE form used by pdepe:
    %   c(x,t,u,DuDx) * ?u/?t = ?/?x [ f(x,t,u,DuDx) ] + s(x,t,u,DuDx)
    %
    % Inputs:
    %   x     - Spatial coordinate (?)
    %   t     - Time (?)
    %   u     - Solution at (x,t), here u = ? (dimensionless temperature)
    %   DuDx  - Spatial derivative of u, i.e. d?/d?
    %   Bi    - Biot number (not used here, but kept for consistency)
    %   Tr    - Reference temperature ratio (not used here, included for consistency)
    %   Sigma - Source/sink parameter
    %
    % Outputs:
    %   c     - Coefficient for time derivative term
    %   f     - Flux term
    %   s     - Source term

    % Transient diffusion equation with source term:
    %   ??/?t = ?²?/??² + Sigma
    %
    % ? In pdepe form:
    %   c = 1                (time coefficient)
    %   f = DuDx             (flux = ??/??)
    %   s = Sigma            (constant source term)

    c = 1;        % Time derivative coefficient
    f = DuDx;     % Flux term (heat conduction)
    s = Sigma;    % Source term
end