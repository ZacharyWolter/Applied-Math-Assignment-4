close all

% Plot true solution against approximate solution for forward euler method,
% explicit midpoint method, and Kutta's third-order method

% Initialize
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;
x0 = 8;
y0 = 0;
dxdt0 = 0;
dydt0 = 1.5;
V0 = [x0;y0;dxdt0;dydt0];
t_range = linspace(0,30,100);
tspan = [t_range(1), t_range(end)];


h_ref = 0.2;


anon_grav_rate_func = @ (t, V) gravity_rate_func(t, V, orbit_params);

figure;

%%

% Analytical

V_list = compute_planetary_motion(t_range,V0,orbit_params);

axis equal; axis square;
axis([-20,20,-20,20])
hold on
plot(0,0,'ro','markerfacecolor','r','markersize', 9, 'DisplayName', 'Sun');
plot(V_list(:,1),V_list(:,2),'k', 'LineWidth', 4, 'DisplayName', 'Analytical Solution');


%%

% Euler
BT_struct = struct();
BT_struct.A = 0;
BT_struct.B = 1;
BT_struct.C = 0;

[t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(anon_grav_rate_func,tspan,V0,h_ref,BT_struct);

plot(X_list(:, 1), X_list(:, 2), 'LineStyle', '--', 'LineWidth', 2.5, 'DisplayName', 'Forward Euler');


%%

% Explicit midpoint method
BT_struct = struct();
BT_struct.A = [0 , 0; 0.5, 0];
BT_struct.B = [0; 1];
BT_struct.C = [0; 0.5];

[t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(anon_grav_rate_func,tspan,V0,h_ref,BT_struct);

plot(X_list(:, 1), X_list(:, 2), 'LineStyle', '-.', 'LineWidth', 2.5, 'DisplayName', 'Explicit Midpoint');




%%

% Kutta's third-order method
BT_struct = struct();
BT_struct.A = [0, 0, 0; 0.5, 0, 0; -1, 2, 0];
BT_struct.B = [1/6; 2/3; 1/6];
BT_struct.C = [0 ; 0.5; 1];


[t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(anon_grav_rate_func,tspan,V0,h_ref,BT_struct);

plot(X_list(:, 1), X_list(:, 2), 'LineStyle', ':', 'LineWidth', 2.5, 'DisplayName', 'RK3');
%% 

xlabel('$x\ (-)$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$y\ (-)$', 'Interpreter', 'latex', 'FontSize', 17);
title({'Comparison of Runge-Kutta (RK3), Forward Euler,', 'and Midpoint Methods for a Celestial System'}, 'Interpreter', 'latex', 'FontSize', 13)
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10.5);

xlim([-6 12.5])
ylim([-8 12])
axis square
grid on;




myfig = gcf();
exportgraphics(myfig, 'Runga-Kutta Comparison.png', 'Resolution', 600);