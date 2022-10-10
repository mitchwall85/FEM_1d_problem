function [kab, fa] = gauss_quadrature(x_0, x_L, k, kappa, f, basis, dbasis)

    % just do k=1 for now
    h_e = abs(x_L - x_0);
    nq = k + 1; 
    xi_q = [-1/sqrt(3), 1/sqrt(3)];  % change for k~=1
    w_q = [1, 1]; % change for k~=1
    kappa_q = kappa(xi_q);
    f_q = f(xi_q);

    for q = 1:nq % use nq = k+1
        x = xi_q(q);
        kappa = kappa_q(q);
        f = f_q(q);
        w = w_q(q);

        kab = zeros(k+1, k+1);
        fa = zeros(k+1, 1);
        for a = 1:k+1
            for b = 1:k+1
                  kab(a,b) = kab(a,b) + 2/h_e*kappa*dbasis{a}(x)*dbasis{b}(x)*w;

            end
            fa(a) = fa(a) + 2/h_e*f*basis{a}(x)*w;

        end

    end

end