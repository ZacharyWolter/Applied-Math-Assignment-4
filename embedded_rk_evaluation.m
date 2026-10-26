function embedded_rk_evaluation
    close all;
    % Embedded RK method evaluation (Heun-Euler)

    % Embedded Heun-Euler tableau
    BT_struct = struct();
    BT_struct.C = [0, 1];
    BT_struct.B = [1/2, 1/2; 1, 0];
    BT_struct.A = [0, 0; 1, 0];

    % Function handles (BT_struct is defined above, so it is captured correctly)
    anon_rate = @(t,X) rate_func01(t,X);
    anon_sol  = @(t) solution01(t);
    G_func    = @(t,XA,h) RK_step_embedded(anon_rate,t,XA,h,BT_struct);

    % Local truncation error sweep
    h_vals = logspace(-5, 1, 100);
    t_ref  = 0;
    e_list1 = zeros(size(h_vals));
    e_list2 = zeros(size(h_vals));
    f_diff  = zeros(size(h_vals));
    for i = 1:numel(h_vals)
        [e_list1(i), e_list2(i)] = embedded_local_truncation_error(G_func, anon_sol, t_ref, h_vals(i));
        f_diff(i) = abs(anon_sol(t_ref + h_vals(i)) - anon_sol(t_ref));
    end

    plot_embedded_error(h_vals, e_list1, e_list2, f_diff);
end

function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
    X = cos(t);
end

function plot_embedded_error(h_vals, e_list1, e_list2, f_diff)
    errors1 = abs(e_list1);
    errors2 = abs(e_list2);
    errors3 = abs(e_list1 - e_list2);
    errors4 = f_diff(:)';

    colors = {[0 0.447 0.741], [0.85 0.325 0.098], [0.466 0.674 0.188]};

    figure;
    hold on;
        
    plot(h_vals, errors1, 'o-', 'Color', colors{1}, 'MarkerFaceColor', colors{1}, 'MarkerSize',8, 'LineWidth',2, 'DisplayName', '$X_{B1}$ (Heun, 2nd order) $\quad$');
    plot(h_vals, errors2, 'o-', 'Color', colors{2}, 'MarkerFaceColor', colors{2}, 'MarkerSize',6, 'LineWidth',2, 'DisplayName', '$X_{B2}$ (Euler, 1st order) $\quad$');
    plot(h_vals, errors3, 'o-', 'Color', colors{3}, 'MarkerFaceColor', colors{3}, 'MarkerSize',3, 'LineWidth',2, 'DisplayName', '$|X_{B1} - X_{B2}|$ (error est.) $\quad$');
    plot(h_vals, errors4, 'k--', 'LineWidth', 3, 'DisplayName', '$|f(t_{\mathrm{refm}}+h) - f(t_{\mathrm{ref}})|$ $\quad$');

    ax = gca;
    set(ax,'XScale','log','YScale','log', 'FontSize',16, 'TickLabelInterpreter','latex', 'LineWidth',1);
    xlim([(1e-5)/2 (1e2)/2])
    ylim([1e-11 1e3])
    grid on;
    box on;

    title('Local Truncation Error vs. Step Size for Embedded Heun--Euler', 'Interpreter','latex','FontSize',13);
    xlabel('Step Size (-)','Interpreter','latex','FontSize',18);
    ylabel('$\epsilon_{\mathrm{local}}$ (-)', 'Interpreter','latex','FontSize',20);
    legend('Location','southeast', 'Interpreter','latex','FontSize',10);

    exportgraphics(gcf, 'Local Error - embedded Heun Euler vs h.png', 'Resolution',600);
end