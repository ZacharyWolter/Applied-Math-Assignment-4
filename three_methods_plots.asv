close all
% Plot true solution against approximate solutions for forward Euler,
% explicit midpoint, and Kutta's third-order method.
% Two subplots: x(t) on top, y(t) on bottom.
 
%% Initialize
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;
x0 = 8;
y0 = 0;
dxdt0 = 0;
dydt0 = 1.5;
V0 = [x0; y0; dxdt0; dydt0];

t_end = 40;
tspan = [0, t_end];
t_range = linspace(0, t_end, 1000);
h_ref = 0.1;

anon_grav_rate_func = @(t, V) gravity_rate_func(t, V, orbit_params);
 
% Analytical
V_list = compute_planetary_motion(t_range, V0, orbit_params);
x_analytical = V_list(:, 1);
y_analytical = V_list(:, 2);
 
%% Numerical 
methods = struct('name', {}, 'BT', {}, 'style', {}, 'width', {});
 
BT = struct('A', [0, 0, 0; 1/2, 0, 0; 0, 3/4, 0], 'B', [2/9, 1/3, 4/9], 'C', [0, 1/2, 3/4]);
methods(1) = struct('name', 'Ralston''s Third-Order Method', 'BT', BT, 'style', '--', 'width', 2.5);
 
BT = struct('A', [0, 0, 0; 0.5, 0, 0; 0, 3/4 0], 'B', [2/9; 1/3; 4/9], 'C', [0; 1/2; 3/4]);
methods(2) = struct('name', '3/8-Rule Fourth-Order Method', 'BT', BT, 'style', '-.', 'width', 2.5);
 
BT = struct('A', [0, 0, 0; 0.5, 0, 0; -1, 2, 0], 'B', [1/6; 2/3; 1/6], 'C', [0; 0.5; 1]);
methods(3) = struct('name', 'Kutta RK3', 'BT', BT, 'style', ':', 'width', 3);
 
results = cell(1, numel(methods));
for i = 1:numel(methods)
    [t_list, X_list, ~, ~] = explicit_RK_fixed_step_integration(anon_grav_rate_func, tspan, V0, h_ref, methods(i).BT);
    results{i} = struct('t', t_list, 'X', X_list);
end
 
% Plot
fig = figure('Position', [100 100 800 800]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, {'Comparison of Forward Euler, Midpoint, and RK3 Methods', ...
    'for a Celestial System'}, 'Interpreter', 'latex', 'FontSize', 16);

comp_labels = {'$x(t)\ (-)$', '$y(t)\ (-)$'};
comp_titles = {'$x$-Component of Planetary Motion', '$y$-Component of Planetary Motion'};
axs = gobjects(1, 2);

for i = 1:2
    ax = nexttile(tl);
    axs(i) = ax;
    hold(ax, 'on');

    if i == 1
        y_exact = x_analytical;
    else
        y_exact = y_analytical;
    end

    plot(ax, t_range, y_exact, 'k', 'LineWidth', 4, 'DisplayName', 'Analytical Solution');
    for J = 1:numel(methods)
        plot(ax, results{J}.t, results{J}.X(:, i), 'LineStyle', methods(J).style, 'LineWidth', methods(J).width, 'DisplayName', methods(J).name);
    end

    xlabel(ax, '$t\ (-)$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel(ax, comp_labels{i}, 'Interpreter', 'latex', 'FontSize', 16);
    title(ax, comp_titles{i}, 'Interpreter', 'latex', 'FontSize', 16);
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', 14);
    grid(ax, 'on');
    xlim(ax, [0, t_end]);
    if i == 1
        ylim(ax, [-4, 9]);
    end
end

% Single shared legend in its own tile at the top, under the main title
lgd = legend(axs(1), 'Orientation', 'horizontal', 'NumColumns', 2, 'Interpreter', 'latex', 'FontSize', 13);
lgd.Layout.Tile = 'north';

exportgraphics(fig, 'Runga-Kutta Compariso.png', 'Resolution', 600);
