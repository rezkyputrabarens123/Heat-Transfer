function Example12_3
    % Example12_3: Transient heat conduction in a 1D slab 
    % using PDE solver (pdepe) and steady-state solution (bvp4c).

    % Parameters
    Bi    = 0.1;   % Biot number
    Tr    = 0.55;  % Reference temperature ratio
    Sigma = 1;     % Some physical parameter (thermal-related)

    % Discretization in space (xi) and time (tau)
    xi  = linspace(0, 1, 41);   % 41 spatial nodes between 0 and 1
    tau = linspace(0, 1, 101);  % 101 time points between 0 and 1

    % Solve PDE using pdepe
    theta = pdepe(0, @pde1D, @pdeIC, @pdeBC, xi, tau, [], Bi, Tr, Sigma);

    % ============================================================
    % FIGURE 1: Transient response at selected spatial positions
    % ============================================================
    z = 0:0.25:1;   % Selected positions (xi = 0, 0.25, 0.5, 0.75, 1)
    figure(1)
    hold on
    for k = 1:length(z)
        % Find index in xi closest to z(k)
        [~, kk] = min(abs(xi - z(k)));

        % Plot temperature ratio vs time
        plot(tau, theta(:, kk), 'k-')

        % Add text labels for curves
        if k == 1
            text(0.5, 1.02*theta(end, kk), '\xi = 0.0 and 0.25')
        elseif k > 2
            text(0.5, theta(end, kk) + 0.02, ['\xi = ' num2str(xi(kk))])
        end
    end
    hold off
    axis([0 1 0.5 1])
    xlabel('\tau')
    ylabel('\theta')
    title('Transient Response at Selected Positions')

    % ============================================================
    % FIGURE 2: Spatial profile at different times
    % ============================================================
    figure(2)
    % Find minimum value of theta at xi = 0 (first column)
    [thmin, imin] = min(theta(:, 1));

    % Plot at tau = 0 (initial condition)
    plot(xi, theta(1,:), 'k-', 'LineWidth', 2)
    hold on
    % Plot at tau = tau(2) and at the time when theta is minimum at xi=0
    plot(xi, theta(2,:), 'k', xi, theta(imin,:), 'k:')

    % Steady-state solution using boundary value solver (bvp4c)
    solinit = bvpinit(linspace(0, 1, 20), [1 1]); 
    sol     = bvp4c(@barode, @barbc, solinit, [], Bi, Sigma);

    % Evaluate steady-state solution
    x = linspace(0, 1, 100);
    y = deval(sol, x);

    % Plot steady-state profile
    plot(x, y(1,:), 'k--');

    % Formatting
    xlabel('\xi')
    ylabel('\theta')
    legend(['\tau = 0 (Initial condition)'], ...
           ['\tau = ' num2str(tau(2))], ...
           ['\tau = ' num2str(tau(imin)) ' (Minimum at \xi = 0)'], ...
           '\tau > 2 (Steady state)', ...
           'Location', 'SouthWest')
    title('Spatial Profiles at Different Times')
    hold off
end