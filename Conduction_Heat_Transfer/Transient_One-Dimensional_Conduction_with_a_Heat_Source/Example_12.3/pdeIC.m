function T0 = pdeIC(x, Bi, Tr, Sigma)
    % pdeIC: Defines the initial condition for the PDE solver (pdepe)
    %
    % Inputs:
    %   x     - Spatial coordinate (?), in the range [0,1]
    %   Bi    - Biot number (not used here, kept for consistency)
    %   Tr    - Reference temperature ratio (not used here)
    %   Sigma - Source/sink parameter (not used here)
    %
    % Output:
    %   T0    - Initial temperature distribution ?(x,0)
    %
    % Example here:
    %   ?(x,0) = 1 - 0.45·x
    %   ? Linear profile: starts at ?(0,0) = 1 and decreases with slope -0.45

    T0 = 1 - 0.45 * x;
end