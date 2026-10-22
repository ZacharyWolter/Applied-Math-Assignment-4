function [epsilon, num_calls] = global_truncation_error(rate_func, step_func, X_func, tspan, h)
    X0 = X_func(tspan(1));
    [~, X_list, ~, num_calls] = fixed_step_integration(rate_func, step_func, tspan, X0, h);
    X_end = X_list(end,:).';
    epsilon = norm(X_end - X_func(tspan(2)));
end

