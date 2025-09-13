function Rt = FindZeros(FunName, Nroot, x, w)
% FindZeros locates up to Nroot zeros of a function specified by FunName
% within the interval defined by vector x, using weights w if needed.
% Inputs:
%   FunName - function handle or name of the function to evaluate
%   Nroot   - maximum number of roots to find
%   x       - vector of x values defining the search interval
%   w       - additional parameter passed to FunName
% Output:
%   Rt      - vector of root locations

% Evaluate the function at all x values
f = feval(FunName, x, w);

% Find indices where sign change occurs (indicating a zero crossing)
indx = find(f(1:end-1) .* f(2:end) < 0);

% Adjust Nroot if fewer sign changes are found
L = length(indx);
if L < Nroot
    warning('Requested %d roots, but only %d sign changes detected.', Nroot, L);
    Nroot = L;
end

% Preallocate output vector
Rt = zeros(Nroot, 1);

% Find roots using fzero within brackets where sign change occurs
for k = 1:Nroot
    bracket = [x(indx(k)), x(indx(k)+1)];
    Rt(k) = fzero(@(z) feval(FunName, z, w), bracket);
end
end