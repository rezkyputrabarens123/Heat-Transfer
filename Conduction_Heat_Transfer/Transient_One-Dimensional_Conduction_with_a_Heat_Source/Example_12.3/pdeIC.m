function T0 = pdeIC(x, Bi, Tr, Sigma)
    % pdeIC: Defines the initial condition for the PDE solver (pdepe)
    T0 = 1 - 0.45 * x;
end