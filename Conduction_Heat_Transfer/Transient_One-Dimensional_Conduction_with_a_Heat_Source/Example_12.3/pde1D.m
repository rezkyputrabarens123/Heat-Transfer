function [c, f, s] = pde1D(x, t, u, DuDx, Bi, Tr, Sigma)
    % pde1D: Defines the PDE for transient 1D heat conduction
    c = 1;        % Time derivative coefficient
    f = DuDx;     % Flux term (heat conduction)
    s = Sigma;    % Source term
end