function res = barbc(ya, yb, Bi, Sigma)
    % barbc: Defines the boundary conditions for the steady-state problem
    res = [ ya(2) - Bi * ya(1);
            yb(1) - 0.55 ];
end