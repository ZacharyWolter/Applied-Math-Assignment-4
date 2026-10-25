function Truncation_Plotting

    % Orbit setup
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;

    X0 = [8; 0; 0; 1.5];
    rate_func = @(t,V) gravity_rate_func(t,V,orbit_params);

    names = {'Forward Euler', 'Explicit Midpoint', 'Kutta Third Order'};
    colors = {[0 0.447 0.741], ...
              [0.85 0.325 0.098], ...
              [0.466 0.674 0.188]};
    nM = 3;

    % Butcher tableaux
    BT = cell(1,nM);

    BT{1} = struct('A', 0, ...
                   'B', 1, ...
                   'C', 0);

    BT{2} = struct('A', [0 0; 0.5 0], ...
                   'B', [0; 1], ...
                   'C', [0; 0.5]);

    BT{3} = struct('A', [0 0 0; 0.5 0 0; -1 2 0], ...
                   'B', [1/6; 2/3; 1/6], ...
                   'C', [0; 0.5; 1]);

    % Global experiment
    tspan = [0 30];
    h_glob = logspace(-4,-2,20);

    V_true = compute_planetary_motion(tspan,X0,orbit_params);
    X_final_true = V_true(end,:)';

    % Global truncation error
    glob_err = zeros(nM,numel(h_glob));
    h_used = zeros(nM,numel(h_glob));

    for m = 1:nM
        for k = 1:numel(h_glob)
            [~,X_list,h_avg,~] = ...
                explicit_RK_fixed_step_integration( ...
                rate_func,tspan,X0,h_glob(k),BT{m});

            % Position error at the final time
            glob_err(m,k) = norm( ...
                X_list(end,1:2)' - X_final_true(1:2));

            h_used(m,k) = h_avg;
        end
    end

    % Fit global error against step size
    c_gh = zeros(nM,2);
    fit_idx = 1:7;

    for m = 1:nM
        c_gh(m,:) = polyfit( ...
            log10(h_used(m,fit_idx)), ...
            log10(glob_err(m,fit_idx)),1);

        fprintf('%s: p = %.3f\n',names{m},c_gh(m,1));
    end


    % Print table of global convergence orders
    fprintf('\n%-25s %14s\n', 'Method', 'p (global,h)');
    for m = 1:nM
        fprintf('%-25s %14.3f\n', names{m}, c_gh(m,1));
    end

    flag = 5;

    if flag == 5
        figure;
        hold on;
        fit_h = gobjects(nM,1);

        for m = 1:nM
            plot(h_used(m,:),glob_err(m,:),'o-', ...
                'Color',colors{m}, ...
                'MarkerFaceColor',colors{m}, ...
                'MarkerSize',7, ...
                'LineWidth',1, ...
                'DisplayName',names{m});

            fit_h(m) = plot( ...
                h_used(m,fit_idx), ...
                10^c_gh(m,2) * h_used(m,fit_idx).^c_gh(m,1), ...
                'k--','LineWidth',2, ...
                'DisplayName',sprintf( ...
                '%s fit: $p = %.2f$',names{m},c_gh(m,1)));
        end

        ax = gca;
        set(ax,'XScale','log','YScale','log', ...
            'FontSize',16, ...
            'TickLabelInterpreter','latex', ...
            'LineWidth',1);

        grid on;
        box on;

        pad_frac = 0.08;
        xall = h_used(:);
        yall = glob_err(:);

        xall = xall(xall > 0 & isfinite(xall));
        yall = yall(yall > 0 & isfinite(yall));

        xl_log = log10([min(xall),max(xall)]);
        yl_log = log10([min(yall),max(yall)]);

        xlim(10.^(xl_log + [-1 1]*pad_frac*diff(xl_log)));
        ylim(10.^(yl_log + [-1 1]*pad_frac*diff(yl_log)));

        % Extend fit lines across the full x-range
        xf = xlim(ax);
        for m = 1:nM
            set(fit_h(m),'XData',xf, ...
                'YData',10^c_gh(m,2) * xf.^c_gh(m,1));
        end

        title('Global Position Error vs. Step Size', ...
            'Interpreter','latex','FontSize',15);
        xlabel('Step Size (-)','Interpreter','latex','FontSize',18);
        ylabel('$\epsilon_{\mathrm{global}}$ (-)', ...
            'Interpreter','latex','FontSize',20);
        legend('Location','southeast', ...
            'Interpreter','latex','FontSize',10);

        exportgraphics(gcf, ...
            'Global Error - explicit vs h - with fits.png', ...
            'Resolution',600);
    end
end