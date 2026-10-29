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


h_ref = 0.15;


anon_grav_rate_func = @ (t, V) gravity_rate_func(t, V, orbit_params);

figure;

%%

% Analytical

V_list = compute_planetary_motion(t_range,V0,orbit_params);

axis equal; axis square;
axis([-20,20,-20,20])
hold on
plot(0,0,'ro','markerfacecolor','r','markersize',5);
plot(V_list(:,1),V_list(:,2),'k');


%%

% Euler
BT_struct = struct();
BT_struct.A = 0;
BT_struct.B = 1;
BT_struct.C = 0;

[t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration ...
(anon_grav_rate_func,tspan,V0,h_ref,BT_struct);

plot(X_list(:, 1), X_list(:, 2));


%%

% Explicit midpoint method
BT_struct = struct();
BT_struct.A = [0 , 0; 0.5, 0];
BT_struct.B = [0; 1];
BT_struct.C = [0; 0.5];

[t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration ...
(anon_grav_rate_func,tspan,V0,h_ref,BT_struct);

plot(X_list(:, 1), X_list(:, 2));




%%

% Kutta's third-order method
BT_struct = struct();
BT_struct.A = [0, 0, 0; 0.5, 0, 0; -1, 2, 0];
BT_struct.B = [1/6; 2/3; 1/6];
BT_struct.C = [0 ; 0.5; 1];


[t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration ...
(anon_grav_rate_func,tspan,V0,h_ref,BT_struct);

plot(X_list(:, 1), X_list(:, 2));



%%
% Embedded Heun-Euler tableau
BT_struct = struct();
BT_struct.C = [0, 1];
BT_struct.B = [1/2, 1/2; 1, 0];
BT_struct.A = [0, 0; 1, 0];

error_desired = 0.05;


[t_list,X_list,h_avg, num_evals] = explicit_RK_variable_step_integration ...
(anon_grav_rate_func,tspan,V0,h_ref,BT_struct,p,error_desired);

plot(X_list(:, 1), X_list(:, 2), 'r--');