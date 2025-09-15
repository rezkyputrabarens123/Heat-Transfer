function [pl, ql, pr, qr] = pdeBC(xl, ul, xr, ur, t, Bi, Tr, Sigma)
    % pdeBC: Defines boundary conditions for the PDE solver (pdepe)
    pr = ur - Tr;
    qr = 0;

    pl = -Bi * ul;
    ql = 1;
end