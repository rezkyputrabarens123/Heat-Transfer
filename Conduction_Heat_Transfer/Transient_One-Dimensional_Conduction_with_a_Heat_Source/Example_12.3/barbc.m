function res = barbc(ya, yb, Bi, Sigma)
    % barbc: Defines the boundary conditions for the steady-state problem
    %
    % Inputs:
    %   ya    - Solution vector at the left boundary (x = 0):
    %           ya(1) = ?(0)   (dimensionless temperature at x=0)
    %           ya(2) = ?'(0)  (first derivative at x=0)
    %   yb    - Solution vector at the right boundary (x = 1):
    %           yb(1) = ?(1)   (dimensionless temperature at x=1)
    %           yb(2) = ?'(1)  (first derivative at x=1)
    %   Bi    - Biot number
    %   Sigma - Source/sink parameter (not used here, but included for consistency)
    %
    % Output:
    %   res   - Residual vector for boundary conditions:
    %           res(1) = ?'(0) - Bi * ?(0) = 0   (convective BC at x=0)
    %           res(2) = ?(1) - 0.55 = 0        (fixed temperature at x=1)

    res = [ ya(2) - Bi * ya(1);   % Left boundary: convective heat transfer
            yb(1) - 0.55 ];       % Right boundary: fixed temperature
end