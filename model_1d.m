function [x,u] = model_1d(k, n_el, kappa, f, g_0, g_L, L)
%{
OUTPUTS:
x: location of nodes
d: finite element solutions at nodes

INPUTS:
k: polynomial degree to be used [1,2,3]
n_el: number of elements to use
kappa: function for k(x) at any point
f: function handle for computing f(x)
g_0: boundary condition at x=0
g_L: boundary condition at x=L
L: length of domain

%}

%% General Values
n_nodes = n_el + 1;
K = zeros(n_nodes - 2, n_nodes - 2);
F = zeros(n_nodes - 2, 1);
dx = L/n_el; % assumed evenly spaced
x = linspace(0, L, n_nodes);


%%
% hard code basis fcns
% basis
basis{1} = @(x) 0.5 - 0.5*x;
basis{2} = @(x) 0.5 + 0.5*x;
% dbasis
dbasis{1} = @(x) -0.5;
dbasis{2} = @(x)  0.5;

%% Element Assembly
n = n_nodes - 2; % subtract boundaries

% Assemble IEN matrix
IEN = zeros(k+1,n_el);
for e = 1:n_el % loop over elements
    for a = 1:k+1 % loop over nodes
        IEN(a,e) = k*(e-1)+(a-1);
    end
end

% loop over elements 
for e = 1:n_el % loop over elements
    x_0 = n_el*(e-1); % location of start of element
    x_L = n_el*e; % location of end of element
    
    for a = 1:k+1 % loop over nodes
        A = IEN(a,e);
        if 1 <= A && A <= n % bounds might be 2 & n-1 here...
            for b = 1:k+1
                B = IEN(b,e);
                disp(A)
                disp(B)
                if 1 <= B && B <= n
                    [kab,fa] = gauss_quadrature(x_0, x_L, k, kappa, f, basis, dbasis);
                    K(A:A+k,B:B+k) = K(A:A+k,B:B+k) + kab;
                    stuff = 1;
                elseif B == 0
                    [kab,fa] = gauss_quadrature(x_0, x_L, k, kappa, f, basis, dbasis);
                    F(A) = F(A) - kab(2,2)*g_0; % jank fix to apply bc
        
                elseif b == n+1
                    [kab,fa] = gauss_quadrature(x_0, x_L, k, kappa, f, basis, dbasis);
                    F(A) = F(A) - kab(1,1)*g_L; % jank fix to apply bc
                
                end
            end
            [kab,fa] = gauss_quadrature(x_0, x_L, k, kappa, f, basis, dbasis);
            F(A:A+k) = F(A:A+k) + fa;
        end
    end
end


d = K\F;
u = zeros(n_el - 1, 1);
for i = 1:n_el - 1
    u = u + d(i)*basis{i}(x);
end



end