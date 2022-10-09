function [x,d] = model_1d(k, n_el, kappa, f, g_0, g_L, L)
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
K = zeros(n_el - 1, n_el - 1);
F = zeros(n_el, 1);
dx = L/n_el; % assumed evenly spaced


%% k = 1
% 3 element galerkin
x = linspace(0, L, n_el + 1);
basis = zeros(n_el - 1, n_el - 1);
for n = 1:length(basis)
    basis(n,n) = 1; 

end

%{


d = k\f;
u_gal1 = zeros(1,length(x_c));
for i = 1:2
    u_gal1 = u_gal1 + d(i)*basis(i,:);
end
u_gal1_err = 100*(u_gal1-u_an(1:14:end))./(u_an(1:14:end));
%}

%% Element Assembly
n = n_nodes - 1; % number of DOF

% Assemble IEN matrix
IEN_mat = zeros(n_el,k+1);
for e = 1:n_el % loop over elements
    for a = 1:k+1 % loop over nodes
        IEN_mat(e,a) = IEN(k,a,e);
    end
end

% loop over elements 
for e = 1:n_el % loop over elements
    for a = 1:k+1 % loop over nodes
        if 1 <= IEN_mat(e,a) && IEN_mat(e,a) <= n
            for b = 1:k+1
                B = IEN(b,e);
                if 1 <= B && B <= n
                    k = 2/dx*trapz(kappa(x0+xi)*dbasis(a)*db n asis(b));
                    
                    K(A,B) = K(A,B) + k_ab;
                
                elseif B == 0
                    F(A) = F(A) - k_ab*g_0;
        
                elseif b == n+1
                    F(A) = F(A) - k_ab*g_L;
                
                end
            end
            F(A) = F(A) + f = 2/dx*trapz(f(x0+xi)*basis(a))
        end
    end
end

a = 1;





end