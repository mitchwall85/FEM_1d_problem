%% Main script for running models
% Mitch Wall
% ASEN 5007
%%
clear; clc; close all;

%% Problem 1:
% test
k = 1;
n_el = 5;
kappa = @(x) x;
f = @(x) x;
g_0 = 0;
g_L = 0;
L = 1;

[x_fem,u_fem] = model_1d(k, n_el, kappa, f, g_0, g_L, L);


%% Problem 2: Analitical Soln
% use model from HW1
% -u,xx = x

x = linspace(0,1,10);
u = @(x) x.^3 + sin(x);
f = @(x) 6*x - sin(x);

% plot
figure()
plot(x,u(x))
plot(x_fem, u_fem)
legend('Manufactured Solution','Galerkin Solution')
xlabel('x')
ylabel('u(x)')
title('P2: Analitical Solution Verification') % TODO is verification the right word? is this a manufactured soln?