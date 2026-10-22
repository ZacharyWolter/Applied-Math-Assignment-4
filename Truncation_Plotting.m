function Truncation_Plotting

    rate_func = @rate_func01;
    solution  = @solution01;
    X0 = 10;

    names = {'Forward Euler', 'Backward Euler', 'Explicit Midpoint', 'Implicit Midpoint'};
    steps = {@forward_euler_step, @backward_euler_step, @explicit_midpoint_step, @implicit_midpoint_step};
    colors  = {[0 0.447 0.741], [0.85 0.325 0.098], [0.466 0.674 0.188], [0.494 0.184 0.556]};
    explicit = [1 3];          
    nM = 4;

    % local experiment 
    t_ref     = 1;
    h_list    = logspace(-5, 1, 100);
    fit_range = [.001, .01];

    % global experiment 
    tspan  = [0 5];
    h_glob = logspace(-3, -1, 20);   

    %  Local truncation error 
    local_err = zeros(nM, numel(h_list));
    for m = 1:nM
        [local_err(m,:), diff_err] = local_truncation_error(rate_func, solution, steps{m}, t_ref, h_list);
    end

    i_start = find(h_list >= fit_range(1), 1, 'first');
    i_end   = find(h_list <= fit_range(2), 1, 'last');
    idx = i_start:i_end;

    c_local = zeros(nM, 2);
    for m = 1:nM
        c_local(m,:) = polyfit(log10(h_list(idx)), log10(local_err(m,idx)), 1);
    end
    c_diff = polyfit(log10(h_list(idx)), log10(diff_err(idx)), 1);

    % Global truncation error 
    glob_err = zeros(nM, numel(h_glob));
    h_used   = zeros(nM, numel(h_glob));
    n_calls  = zeros(nM, numel(h_glob));
    X_final_true = solution(tspan(2));

    for m = 1:nM
        for k = 1:numel(h_glob)
            [~, X_list, h_avg, num_evals] = fixed_step_integration(rate_func, steps{m}, tspan, X0, h_glob(k));
            glob_err(m,k) = norm(X_list(end,:)' - X_final_true);
            h_used(m,k)   = h_avg;
            n_calls(m,k)  = num_evals;
        end
    end

    c_gh = zeros(nM, 2);   % global error vs h
    c_gn = zeros(nM, 2);   % global error vs number of calls
    for m = 1:nM
        c_gh(m,:) = polyfit(log10(h_used(m,(1:7))),  log10(glob_err(m,(1:7))), 1);
        c_gn(m,:) = polyfit(log10(n_calls(m,(1:7))), log10(glob_err(m,(1:7))), 1);
    end


    flag = 3;

    if flag == 1   % Plot 1 - Local Error: explicit methods vs h: with fits
        figure; hold on;
        fit_h = gobjects(nM,1);

        for m = explicit
            
            plot(h_list, local_err(m,:), 'o', 'Color', colors{m}, 'MarkerFaceColor', colors{m}, 'MarkerSize', 5, 'DisplayName', names{m});
            fit_h(m) = plot(h_list(idx), 10^c_local(m,2) * h_list(idx).^c_local(m,1), '-', 'Color', 'k', 'LineWidth', 2, 'DisplayName', sprintf('%s fit: $p = %.2f$', names{m}, c_local(m,1)));
        end

        plot(h_list, diff_err, 'o', 'color', colors{3}, 'MarkerFaceColor', colors{3}, 'MarkerSize', 7, 'DisplayName', '$|X(t+h)-X(t)|$');
        fit_diff_h = plot(h_list(idx), 10^c_diff(2) * h_list(idx).^c_diff(1), 'k-', 'LineWidth', 2, 'DisplayName', sprintf('$|X(t+h)-X(t)|$ fit: $p = %.2f$', c_diff(1)), 'Color', 'k');
        ax = gca;
        set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'LineWidth', 1);
        grid on; box on;
        pad_frac = 0.08;
        line_h = findobj(ax, 'Type', 'Line');
        xall = []; yall = [];

        for i = 1:numel(line_h)
            xall = [xall, line_h(i).XData(:)'];
            yall = [yall, line_h(i).YData(:)'];
        end

        xall = xall(xall > 0 & isfinite(xall));
        yall = yall(yall > 0 & isfinite(yall));
        xl_log = log10([min(xall), max(xall)]);
        yl_log = log10([min(yall), max(yall)]);
        xlim(10.^(xl_log + [-1 1]*pad_frac*diff(xl_log)));
        ylim(10.^(yl_log + [-1 1]*pad_frac*diff(yl_log)));


        xf = xlim(ax);
        for m = explicit
            set(fit_h(m), 'XData', xf, 'YData', 10^c_local(m,2) * xf.^c_local(m,1));
        end

        set(fit_diff_h, 'XData', xf, 'YData', 10^c_diff(2) * xf.^c_diff(1));

        title('Local Truncation Error vs. Step Size: Explicit Methods', 'Interpreter', 'latex', 'FontSize', 16);
        xlabel('$h\ (-)$', 'Interpreter', 'latex', 'FontSize', 17);
        ylabel('$\epsilon_{\mathrm{local}}\ (-)$', 'Interpreter', 'latex', 'FontSize', 19);
        legend('Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 10.5);

        myfig = gcf();
        exportgraphics(myfig, 'Local Error - explicit vs h - with fits.png', 'Resolution', 600);

    elseif flag == 2   % Plot 2 - Local Error: all methods vs h: no fits 
        figure; hold on;
        for m = 1:nM
            
            if m == 1 || m == 3
                plot(h_list, local_err(m,:), 'o', 'Color', colors{m}, 'MarkerFaceColor', colors{m}, 'MarkerSize', 6, 'DisplayName', names{m});
            else
                plot(h_list, local_err(m,:), 'o', 'Color', colors{m}, 'MarkerFaceColor', colors{m}, 'MarkerSize', 3, 'DisplayName', names{m});
            end
        end
    
        ax = gca;
        set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'LineWidth', 1);
        grid on; box on;
        pad_frac = 0.08;
        line_h = findobj(ax, 'Type', 'Line');
        xall = []; yall = [];
        for i = 1:numel(line_h)
            xall = [xall, line_h(i).XData(:)'];
            yall = [yall, line_h(i).YData(:)'];
        end
        xall = xall(xall > 0 & isfinite(xall));
        yall = yall(yall > 0 & isfinite(yall));
        xl_log = log10([min(xall), max(xall)]);
        yl_log = log10([min(yall), max(yall)]);
        xlim(10.^(xl_log + [-1 1]*pad_frac*diff(xl_log)));
        ylim(10.^(yl_log + [-1 1]*pad_frac*diff(yl_log)));
        title('Local Truncation Error vs. Step Size: All Methods', 'Interpreter', 'latex', 'FontSize', 16);
        xlabel('$h\ (-)$', 'Interpreter', 'latex', 'FontSize', 17);
        ylabel('$\epsilon_{\mathrm{local}}\ (-)$', 'Interpreter', 'latex', 'FontSize', 19);
        legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10.5);

        myfig = gcf();
        exportgraphics(myfig, 'Local Error - all vs h - no fits.png', 'Resolution', 600);
    
    elseif flag == 3   % Plot 3 - Global Error vs h: explicit methods: with fits
        figure; hold on;
        fit_h = gobjects(nM,1);
        for m = explicit
            plot(h_used(m,:), glob_err(m,:), 'o', 'Color', colors{m}, 'MarkerFaceColor', colors{m}, 'MarkerSize', 7, 'DisplayName', names{m});
            fit_h(m) = plot(h_used(m,:), 10^c_gh(m,2) * h_used(m,:).^c_gh(m,1), '-', 'Color', 'k', 'LineWidth', 2, 'DisplayName', sprintf('%s fit: $p = %.2f$', names{m}, c_gh(m,1)));
        end
    
        ax = gca;
        set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'LineWidth', 1);
        grid on; box on;
        pad_frac = 0.08;
        line_h = findobj(ax, 'Type', 'Line');
        xall = []; yall = [];
        for i = 1:numel(line_h)
            xall = [xall, line_h(i).XData(:)'];
            yall = [yall, line_h(i).YData(:)'];
        end
        xall = xall(xall > 0 & isfinite(xall));
        yall = yall(yall > 0 & isfinite(yall));
        xl_log = log10([min(xall), max(xall)]);
        yl_log = log10([min(yall), max(yall)]);
        xlim(10.^(xl_log + [-1 1]*pad_frac*diff(xl_log)));
        ylim(10.^(yl_log + [-1 1]*pad_frac*diff(yl_log)));

        % extend all fit lines across the full x-range of the axes
        xf = xlim(ax);
        for m = explicit
            set(fit_h(m), 'XData', xf, 'YData', 10^c_gh(m,2) * xf.^c_gh(m,1));
        end

        title('Global Truncation Error vs. Step Size: Explicit Methods', 'Interpreter', 'latex', 'FontSize', 16);
        xlabel('$h\ (-)$', 'Interpreter', 'latex', 'FontSize', 17);
        ylabel('$\epsilon_{\mathrm{global}}\ (-)$', 'Interpreter', 'latex', 'FontSize', 19);
        legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10.5);

        myfig = gcf();
        exportgraphics(myfig, 'Global Error - explicit vs h - with fits.png', 'Resolution', 600);
    
    elseif flag == 4   % Plot 4 - Global Error vs calls: all methods: with fits
    
        figure; hold on;
        fit_h = gobjects(nM,1);
        for m = explicit
            plot(n_calls(m,:), glob_err(m,:), 'o', 'Color', colors{m}, 'MarkerFaceColor', colors{m}, 'MarkerSize', 7, 'DisplayName', names{m});
            fit_h(m) = plot(n_calls(m,:), 10^c_gn(m,2) * n_calls(m,:).^c_gn(m,1), '-', 'Color', 'k', 'LineWidth', 2, 'DisplayName', sprintf('%s fit: $p = %.2f$', names{m}, c_gn(m,1)));
        end
    
        ax = gca;
        set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'LineWidth', 1);
        grid on; box on;
        pad_frac = 0.08;
        line_h = findobj(ax, 'Type', 'Line');
        xall = []; yall = [];
        for i = 1:numel(line_h)
            xall = [xall, line_h(i).XData(:)'];
            yall = [yall, line_h(i).YData(:)'];
        end
        xall = xall(xall > 0 & isfinite(xall));
        yall = yall(yall > 0 & isfinite(yall));
        xl_log = log10([min(xall), max(xall)]);
        yl_log = log10([min(yall), max(yall)]);
        xlim(10.^(xl_log + [-1 1]*pad_frac*diff(xl_log)));
        ylim(10.^(yl_log + [-1 1]*pad_frac*diff(yl_log)));

        % extend all fit lines across the full x-range of the axes
        xf = xlim(ax);
        for m = explicit
            set(fit_h(m), 'XData', xf, 'YData', 10^c_gn(m,2) * xf.^c_gn(m,1));
        end

        title('Global Truncation Error vs. Rate Function Calls: Explicit Methods', 'Interpreter', 'latex', 'FontSize', 13);
        xlabel('$\mathrm{Number\ of\ Rate\ Function\ Calls}\ (-)$', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('$\epsilon_{\mathrm{global}}\ (-)$', 'Interpreter', 'latex', 'FontSize', 16);
        legend('Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 9.75);

        myfig = gcf();
        exportgraphics(myfig, 'Global Error - explicit vs calls - with fits.png', 'Resolution', 600);
    
    elseif flag == 5   % Plot 5 - Global Error vs h: all methods: no fits
    
        figure; hold on;
        for m = 1:nM
            plot(h_used(m,:), glob_err(m,:), 'o-', 'Color', colors{m}, 'MarkerFaceColor', colors{m}, 'MarkerSize', 7, 'LineWidth', 1, 'DisplayName', names{m});
        end
    
        ax = gca;
        set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'LineWidth', 1);
        grid on; box on;
        pad_frac = 0.08;
        line_h = findobj(ax, 'Type', 'Line');
        xall = []; yall = [];
        for i = 1:numel(line_h)
            xall = [xall, line_h(i).XData(:)'];
            yall = [yall, line_h(i).YData(:)'];
        end
        xall = xall(xall > 0 & isfinite(xall));
        yall = yall(yall > 0 & isfinite(yall));
        xl_log = log10([min(xall), max(xall)]);
        yl_log = log10([min(yall), max(yall)]);
        xlim(10.^(xl_log + [-1 1]*pad_frac*diff(xl_log)));
        ylim(10.^(yl_log + [-1 1]*pad_frac*diff(yl_log)));
        title('Global Truncation Error vs. Step Size: All Methods', 'Interpreter', 'latex', 'FontSize', 15);
        xlabel('$h\ (-)$', 'Interpreter', 'latex', 'FontSize', 18);
        ylabel('$\epsilon_{\mathrm{global}}\ (-)$', 'Interpreter', 'latex', 'FontSize', 20);
        legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10);

        myfig = gcf();
        exportgraphics(myfig, 'Global Error - explicit vs h - no fits.png', 'Resolution', 600);
    
    elseif flag == 6 % Plot 6 - Global Error vs h: all methods: no fits
    
        figure; hold on;
        for m = 1:nM
            plot(n_calls(m,:), glob_err(m,:), 'o-', 'Color', colors{m}, 'MarkerFaceColor', colors{m}, 'MarkerSize', 7, 'LineWidth', 1, 'DisplayName', names{m});
        end
    
        ax = gca;
        set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 16, 'TickLabelInterpreter', 'latex', 'LineWidth', 1);
        grid on; box on;
        pad_frac = 0.08;
        line_h = findobj(ax, 'Type', 'Line');
        xall = []; yall = [];
        for i = 1:numel(line_h)
            xall = [xall, line_h(i).XData(:)'];
            yall = [yall, line_h(i).YData(:)'];
        end
        xall = xall(xall > 0 & isfinite(xall));
        yall = yall(yall > 0 & isfinite(yall));
        xl_log = log10([min(xall), max(xall)]);
        yl_log = log10([min(yall), max(yall)]);
        xlim(10.^(xl_log + [-1 1]*pad_frac*diff(xl_log)));
        ylim(10.^(yl_log + [-1 1]*pad_frac*diff(yl_log)));
        title('Global Truncation Error vs. Rate Function Calls: All Methods', 'Interpreter', 'latex', 'FontSize', 14);
        xlabel('$\mathrm{Number\ of\ Rate\ Function\ Calls}\ (-)$', 'Interpreter', 'latex', 'FontSize', 17);
        ylabel('$\epsilon_{\mathrm{global}}\ (-)$', 'Interpreter', 'latex', 'FontSize', 19);
        legend('Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 9.5);

        myfig = gcf();
        exportgraphics(myfig, 'Global Error - all vs calls - no fits.png', 'Resolution', 600);
    
    elseif flag == 7   % print tables
        fprintf('\n%-20s %12s %14s %16s\n', 'Method', 'p (local)', 'p (global,h)', 'p (global,calls)');
        fprintf('%-20s %12.3f\n', '|X(t+h)-X(t)|', c_diff(1));
        for m = 1:nM
            fprintf('%-20s %12.3f %14.3f %16.3f\n', names{m}, c_local(m,1), c_gh(m,1), c_gn(m,1));
        end
    end
end

function [local_err, diff_err] = local_truncation_error(rate_func, solution, step_func, t_ref, h_list)
    local_err = zeros(size(h_list));
    diff_err  = zeros(size(h_list));
    X_exact_now = solution(t_ref);

    for k = 1:length(h_list)
        h = h_list(k);
        X_exact_next = solution(t_ref + h);
        X_step = step_func(rate_func, t_ref, X_exact_now, h);

        local_err(k) = norm(X_step - X_exact_next);
        diff_err(k)  = norm(X_exact_next - X_exact_now);
    end
end

% Fixed step intergration consolidated from all solvers 
function [t_list,X_list,h_avg,num_evals] = fixed_step_integration(rate_func_in,step_func,tspan,X0,h_ref)
    N = ceil((tspan(2) - tspan(1))/h_ref);
    h = (tspan(2) - tspan(1))/N;

    t_list = linspace(tspan(1), tspan(2), N+1)';
    X_list = zeros(N+1, length(X0));
    X_list(1,:) = X0';

    num_evals = 0;
    for i = 1:N
        [XB, evals] = step_func(rate_func_in, t_list(i), X_list(i,:)', h);
        X_list(i+1,:) = XB';
        num_evals = num_evals + evals;
    end
    h_avg = h;
end



% Step functions (moved here for better naming convention) 
function [XB,num_evals] = forward_euler_step(rate_func_in,t,XA,h)
    XB = XA + h*rate_func_in(t,XA);
    num_evals = 1;
end

function [XB,num_evals] = explicit_midpoint_step(rate_func_in,t,XA,h)
    X_mid = XA + (h/2)*rate_func_in(t,XA);
    XB = XA + h*rate_func_in(t + h/2, X_mid);
    num_evals = 2;
end

function [XB,num_evals] = backward_euler_step(rate_func_in,t,XA,h)
    G = @(X) XA + h*rate_func_in(t+h, X) - X;
    dxtol = 1e-10;  ftol = 1e-10;  max_iter = 200;  dxmax = 1e7;
    [XB,~,num_evals] = multi_newton_solver(G,XA,dxtol,ftol,max_iter,dxmax);
end

function [XB,num_evals] = implicit_midpoint_step(rate_func_in,t,XA,h)
    G = @(X) XA + h*rate_func_in(t + h/2, (XA + X)/2) - X;
    dxtol = 1e-10;  ftol = 1e-10;  max_iter = 200;  dxmax = 1e7;
    [XB,~,num_evals] = multi_newton_solver(G,XA,dxtol,ftol,max_iter,dxmax);
end

function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
    X = cos(t);
end


% Multi Newton solver from assignment 1
function [x, exit_flag, num_evals] = multi_newton_solver(fun,x0,dxtol,ftol,max_iter,dxmax)
    x = x0(:);
    num_evals = 0;
    exit_flag = -1;

    for i = 1:max_iter
        fval = fun(x);
        num_evals = num_evals + 1;

        if norm(fval) < ftol
            exit_flag = 1;
            return
        end

        [J, n] = approximate_jacobian(fun, x, fval);
        num_evals = num_evals + n;

        dx = -J \ fval;
        x = x + dx;

        if norm(dx) > dxmax
            exit_flag = 0;
            return
        end
        if norm(dx) < dxtol
            exit_flag = 2;
            return
        end
    end
    disp("Newton solver did not converge within tolerance")
end

function [J, num_evals] = approximate_jacobian(fun, x, fval)
    n = length(x);
    J = zeros(length(fval), n);
    num_evals = 0;
    for j = 1:n
        delta = sqrt(eps) * max(1, abs(x(j)));
        xp = x;
        xp(j) = xp(j) + delta;
        J(:,j) = (fun(xp) - fval) / delta;
        num_evals = num_evals + 1;
    end
end
