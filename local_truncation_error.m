function epsilon = embedded_local_truncation_error(G_func, X_func,t, h)
    X0 = X_func(t);
    G_func(t, X0, h)-X_func(t+h)
    epsilon = norm(G_func(t, X0, h)-X_func(t+h));
end


[XB1, XB2, num_evals]