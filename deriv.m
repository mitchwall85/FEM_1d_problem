function [x, du] = deriv(x, u, dbasis, IEN, L)
    % evaluate derivative with shape functions
    [~, n_el] = size(IEN);
    k =  length(dbasis) - 1;
    dx = L/n_el;

    xi_elem_vect = linspace(-1,1,k+1);
    du = zeros(length(u),1);
    i = 0;
    for e = 1:n_el % for all elements in space
        i = i + 1;
        for xi_elem = 1:k % for each dof
            xi = xi_elem_vect(xi_elem); % xi location in parent elem
            for d = 1:k+1 % output DOFs
                A = IEN(d,e); % forget about BCs right now. those should go to zero anyway 
                du(i) = du(i) + u(A+1)*dbasis{d}(xi)*2/dx;  

            end

            i = i + 1;
            
            if e == n_el  && xi_elem == k % catch last elem
                for d = 1:k+1 % output DOFs
                    A = IEN(d,e); % forget about BCs right now. those should go to zero anyway 
                    du(i) = du(i) + u(A+1)*dbasis{d}(xi)*2/dx;  
    
                end
            end

        end
        i = i - 1;
    end

end