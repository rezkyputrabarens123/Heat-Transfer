function Example12_8
% This script calculates the radiation view factor (F_12) between two 
% parallel rectangular surfaces for two different geometric configurations.
% The view factor represents the fraction of radiation energy leaving
% surface 1 that directly strikes surface 2.
%
% The calculation involves solving a quadruple integral, which is handled
% here by nesting calls to MATLAB's 'integral2' function.

% --- Configuration 1 ---
% Two identical 2x2 squares, aligned and centered on the z-axis,
% separated by a distance of dz=2.
fprintf('Calculating view factor for Set 1...\n');
Set1 = F1_2(-1, 1, -1, 1, -1, 1, -1, 1, 2);
fprintf('View factor for aligned squares (Set 1): F_12 = %.4f\n\n', Set1);

% --- Configuration 2 ---
% Two adjacent 2x2 squares in different quadrants, separated by dz=2.
% NOTE: The original call had inverted integration limits for surface 2
% (e.g., x from 2 to 0). This has been corrected to be from 0 to 2 for clarity.
fprintf('Calculating view factor for Set 2...\n');
Set2 = F1_2(-2, 0, -2, 0, 0, 2, 0, 2, 2);
fprintf('View factor for adjacent squares (Set 2): F_12 = %.4f\n', Set2);
end

function F12 = F1_2(x1a, x1b, y1a, y1b, x2a, x2b, y2a, y2b, dz)
% Calculates the view factor F_12 between two parallel rectangles.
%
% Surface 1 is defined by: x in [x1a, x1b], y in [y1a, y1b] at z=0.
% Surface 2 is defined by: x in [x2a, x2b], y in [y2a, y2b] at z=dz.
%
% The view factor is given by the formula:
% F_12 = (1 / (A1 * pi)) * ?_A1 ?_A2 (cos?1 * cos?2 / r^2) dA2 dA1

% Calculate the area of the emitting surface (surface 1).
% The original code confusingly named this 'A2'. Renamed to 'A1'.
A1 = abs(x1b - x1a) * abs(y1b - y1a);

% Check for a non-physical surface.
if A1 == 0
    error('The area of surface 1 cannot be zero.');
end

% Perform the outer integration over the area of surface 1.
% The integrand is the 'OuterKernel' function, which itself computes the
% inner integral over surface 2 for each point on surface 1.
integral_value = integral2(@OuterKernel, x1a, x1b, y1a, y1b);

% Calculate the final view factor.
F12 = integral_value / (A1 * pi);

    % --- Nested Function: Outer Kernel ---
    function f_outer = OuterKernel(x1, y1)
        % This function is the integrand for the outer integral over surface 1.
        % The 'integral2' function calls this with matrices of x1 and y1 points.
        % Since the inner integral must be calculated for each individual point
        % (x1, y1), we use 'arrayfun' to loop over the input matrices.
        
        % For each point (x1_pt, y1_pt) from the input matrices, 'arrayfun'
        % calls the anonymous function, which performs the inner integration.
        f_outer = arrayfun(@(x1_pt, y1_pt) ...
            integral2(@(x2, y2) InnerKernel(x1_pt, y1_pt, x2, y2), x2a, x2b, y2a, y2b), ...
            x1, y1);
    end

    % --- Nested Function: Inner Kernel ---
    function f_inner = InnerKernel(x1, y1, x2, y2)
        % This function is the innermost integrand of the view factor equation.
        % It represents: (cos?1 * cos?2) / r^2 = dz^2 / ((x1-x2)^2 + (y1-y2)^2 + dz^2)^2
        %
        % Args:
        %   (x1, y1): A single point on surface 1.
        %   (x2, y2): Matrices of points on surface 2 (integration variables).
        
        % Calculate the squared distance (r^2) between the points.
        % This is vectorized to work efficiently with 'integral2'.
        r_squared = (x1 - x2).^2 + (y1 - y2).^2 + dz^2;
        
        % Calculate the integrand value.
        f_inner = dz^2 ./ (r_squared.^2);
    end

end