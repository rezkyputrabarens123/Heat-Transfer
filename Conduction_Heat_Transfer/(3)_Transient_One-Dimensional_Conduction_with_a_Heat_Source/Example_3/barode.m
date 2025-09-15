function dydx = barode(x, y, Bi, Sigma)
    % barode: Defines the system of ODEs for the steady-state problem
    dydx = [ y(2);
            -Sigma ];
end