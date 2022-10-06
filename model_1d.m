function [x,d] = model_1d(k, n_el, kalla, f, g_0, g_L, L)
%{
OUTPUTS:
x: location of nodes
d: finite element solutions at nodes

INPUTS:
k: polynomial degree to be used
n_el: number of elements to use
kappa: function for k(x) at any point
f: function handle for computing f(x)
g_0: boundary condition at x=0
g_L: boundary condition at x=L
L: length of domain

%}


end