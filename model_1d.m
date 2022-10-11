function [x,u, norm] = model_1d(k, n_el, kappa, f, g_0, g_L, L, u_an, du_an)
%{
OUTPUTS:
x: location of nodes
u: finite element solutions at nodes

INPUTS:
k: polynomial degree to be used [1,2,3]
n_el: number of elements to use
kappa: function for k(x) at any point
f: function handle for computing f(x)
g_0: boundary condition at x=0
g_L: boundary condition at x=L
L: length of domain
u_an: analitcal soln  

%}

%% General Values
dof = n_el - 1 + n_el*(k - 1);
K = zeros(dof, dof);
F = zeros(dof, 1);
dx = L/n_el; % assumed evenly spaced
x = linspace(0, L, dof + 2);
n = dof; % subtract boundaries


%%
% hard code basis fcns
switch k
    case 1
        % basis
        basis{1} = @(x) 0.5 - 0.5*x;
        basis{2} = @(x) 0.5 + 0.5*x;
        % dbasis
        dbasis{1} = @(x) -0.5;
        dbasis{2} = @(x)  0.5;
        
        % quadrature
        xi_q = [-1/sqrt(3), 1/sqrt(3)];
        w_q  = [1, 1];

    case 2
        % basis
        basis{1} = @(x) -0.5*x + 0.5*x.^2;
        basis{2} = @(x) 1 - x.^2;
        basis{3} = @(x)  0.5*x + 0.5*x.^2;
        % dbasis
        dbasis{1} = @(x) -0.5 + x;
        dbasis{2} = @(x) -2*x;
        dbasis{3} = @(x)  0.5 + x;
%         x = linspace(-1,1,10); % test domain

        % quadrature
        xi_q = [-sqrt(3/5), 0, sqrt(3/5)];
        w_q  =  [ 5/9, 8/9, 5/9];

    case 3
        % basis
        basis{1} = @(x)   27/48*(-1/9 + 1/9*x + x.^2 - x.^3);
        basis{2} = @(x)   27/16*(1/3 - x - 1/3*x.^2 + x.^3);
        basis{3} = @(x)   27/16*(1/3 + x - 1/3*x.^2 - x.^3);
        basis{4} = @(x)  -9/16* (1/9 + x/9 - x.^2 - x.^3);
        % dbasis
        dbasis{1} = @(x)  27/48*(1/9 + 2*x - 3*x.^2);
        dbasis{2} = @(x)  27/16*( -1 - 2/3*x + 3*x.^2);
        dbasis{3} = @(x)  27/16*(1 - 2/3*x - 3*x.^2);
        dbasis{4} = @(x) -9/16* (1/9 - 2*x - 3*x.^2);
        % x = linspace(-1,1,10) % test domain

        % quadrature
        xi_q = [-sqrt((3/7) + 2/7*sqrt(6/5)), ...
                -sqrt((3/7) - 2/7*sqrt(6/5)), ...
                 sqrt((3/7) - 2/7*sqrt(6/5))...
                 sqrt((3/7) + 2/7*sqrt(6/5))];
        w_q  = [(18 - sqrt(30))/36, ...
                (18 + sqrt(30))/36, ...
                (18 + sqrt(30))/36, ...
                (18 - sqrt(30))/36];

end

%% Element Assembly
% Assemble IEN matrix
IEN = zeros(k+1,n_el);
for e = 1:n_el % loop over elements
    for a = 1:k+1 % loop over nodes
        IEN(a,e) = k*(e-1)+(a-1);
    end
end

% loop over elements 
for e = 1:n_el % loop over elements
    x_0 = dx*(e-1); % location of start of element
    x_L = dx*e; % location of end of element
    
    for a = 1:k+1 % loop over nodes
        A = IEN(a,e);
        if 1 <= A && A <= n % bounds might be 2 & n-1 here...
            for b = 1:k+1
                B = IEN(b,e);
%                 disp(["A: ",num2str(A)])
%                 disp(["B: ",num2str(B)])  
                if 1 <= B && B <= n
                    [kab,~] = gauss_quadrature(x_0, x_L, k, kappa, f, basis, dbasis, a, b, xi_q, w_q); % maybe I shouldnt be outputting F
                    K(A,B) = K(A,B) + kab;
                    stuff = 1;
                elseif B == 0
                    [kab,~] = gauss_quadrature(x_0, x_L, k, kappa, f, basis, dbasis, a, b, xi_q, w_q);
                    F(A) = F(A) - kab*g_0;
                    stuff = 1;
        
                elseif B == n+1
                    [kab,~] = gauss_quadrature(x_0, x_L, k, kappa, f, basis, dbasis, a, b, xi_q, w_q);
                    F(A) = F(A) - kab*g_L;
                    stuff = 1;
                
                end
            end
   
            [~,fa] = gauss_quadrature(x_0, x_L, k, kappa, f, basis, dbasis, a, b, xi_q, w_q);
            F(A) = F(A) + fa;
            stuff = 1;
        end
    end
end

d = K\F;
u = [g_0; d; g_L];

% need to transform basis fcns again
[norm] = H1_norm(u, u_an, du_an, basis, dbasis, L, xi_q, w_q, IEN);


end