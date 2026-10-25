function [epsilon1, epsilon2] = embedded_local_truncation_error(G_func, X_func,t, h)
    X0 = X_func(t);
    [XB1, XB2, ~] = G_func(t, X0, h);
    epsilon1 = norm(XB1-X_func(t+h));
    epsilon2 = norm(XB2-X_func(t+h));
end


