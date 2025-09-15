function dydx = barode(x, y, Bi, Sigma)
    % barode: Defines the system of ODEs for the steady-state problem
    %
    % Inputs:
    %   x     - Independent variable (spatial coordinate ?)
    %   y     - Vector of dependent variables:
    %           y(1) = ? (dimensionless temperature)
    %           y(2) = d?/d? (first derivative)
    %   Bi    - Biot number (not directly used here, but kept for consistency)
    %   Sigma - Source/sink term parameter
    %
    % Output:
    %   dydx  - Column vector containing the derivatives:
    %           dydx(1) = d?/d?
    %           dydx(2) = d²?/d?²

    % System of first-order ODEs:
    %   d?/d?     = y(2)
    %   d²?/d?²   = -Sigma
    dydx = [ y(2);          % First equation: ?' = y2
            -Sigma ];       % Second equation: ?'' = -Sigma
end