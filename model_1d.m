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
x = [0, 1/3, 2/3, 1];   
x = linspace(0, L, n_el + 1);
basis = zeros(n_el - 1, n_el - 1);
for n = 1:length(basis)
    basis(n,n) = 1; 

end

k = zeros(n_nodes - 1, n_nodes - 1);
k = -[6, -3; -3, 6];
f = zeros(2,1);
for a = 1:2
    f(a) = trapz(x_c,x_c.*basis(a,:));
end

d = k\f;
u_gal1 = zeros(1,length(x_c));
for i = 1:2
    u_gal1 = u_gal1 + d(i)*basis(i,:);
end
u_gal1_err = 100*(u_gal1-u_an(1:14:end))./(u_an(1:14:end));




end