function [norm] = H1_norm(u, u_an, du_an, basis, dbasis, L, xi_q, w_q, IEN)
    
    [~, n_el] = size(IEN);
    k =  length(basis) - 1;
    dx = L/n_el;
    dof = n_el - 1 + n_el*(k - 1);
    
    prod = 0;
    for e = 1:n_el % for all elements in space
        x_0 = dx*(e-1);
        x_L = dx*e;
        for q = 1:k+1 % for each quadriture point in the elem
            xi = xi_q(q); % xi in paraent elem
            w = w_q(q); % weight in parent elem
            x_q = x_0 + (x_L - x_0)*(xi+1)/2; % transform to physical coords
            u_q = 0;
            du_q = 0;
            for d = 1:k+1 % output DOFs
                A = IEN(d,e); % forget about BCs right now. those should go to zero anyway
                u_q = u_q + u(A+1)*basis{d}(xi);     
                du_q = du_q + u(A+1)*dbasis{d}(xi)*2/dx;  
            end
%             prod = prod + w*((u_q - u_an(x_q))^2);
            prod = prod + w*((u_q - u_an(x_q))^2 + (du_q - du_an(x_q))^2);
        end
    end

    norm = sqrt(prod);

end