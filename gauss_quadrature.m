function [kab, fa] = gauss_quadrature(x_0, x_L, k, kappa, f, basis, dbasis, a, b, xi_q, w_q)

    % just do k=1 for now
    h_e = x_L - x_0;
    nq = k + 1; % use nq = k+1, advised from notes
    
    % initalize
    kab = 0;
    fa = 0;

    % loop over quadrature points
    for q = 1:nq
        xi = xi_q(q);
        w = w_q(q);
        x_phys = x_0 + (x_L - x_0)*(xi + 1)/2; % physical x location

        % evaluate k and f
        kab = kab + 2/h_e*w*kappa(x_phys)*dbasis{a}(xi)*dbasis{b}(xi);
        fa = fa + h_e/2*w*f(x_phys)*basis{a}(xi);
    end
end