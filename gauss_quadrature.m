function [kab, fa] = gauss_quadrature(x_0, x_L, k, kappa, f, basis, dbasis, a, b)

    % just do k=1 for now
    h_e = abs(x_L - x_0);
    nq = k + 1; 
    xi_q = [-1/sqrt(3), 1/sqrt(3)];  % change for k~=1
    w_q = [1, 1]; % change for k~=1


    for q = 1:nq % use nq = k+1
        xi = xi_q(q);
        w = w_q(q);

        x_phys = x_0 + (x_L - x_0)*(xi + 1)/2; % physical x location
%         disp(x_phys)
        kappa_q = kappa(x_phys); % should be at physical locations
        f_q = f(x_phys);

        kab =  2/h_e*kappa_q*dbasis{a}(xi)*dbasis{b}(xi)*w;

        fa = 2/h_e*f_q*basis{a}(xi)*w;
    end

end