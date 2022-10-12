%% Main script for running models
% Mitch Wall
% ASEN 5007
%%
clear; clc; close all;

%% Problem 2: Analitical Soln
% {
k_order = [1, 2, 3];
% n_el = [3]; % idk lets try these values
n_el = [3, 2^2, 2^3, 2^4, 2^5, 2^6, 2^7, 2^8];
kappa = @(x) 1;
f = @(x) -6*x + sin(x);
g_0 = 0;
g_L = 1.841470984807897;
L = 1;

norm_mat = zeros(length(n_el), length(k_order));
figure()
hold on
% analitical soln
u = @(x) x.^3 + sin(x);
du = @(x) 3*x.^2 + cos(x);
for k = 1:length(k_order)
    for el = 1:length(n_el) 
        % execute solution
        [x_fem, u_fem, du_fem, norm] = model_1d(k_order(k), n_el(el), kappa, f, g_0, g_L, L, u, du);
        % calculate norms
        norm_mat(el, k) = norm; 

        %{
        % block to plot indivdual solutions
        figure()
        hold on
        plot(x_fem, u(x_fem))
        scatter(x_fem, u_fem)
        legend('Manufactured Solution',['Galerkin Solution, k = ',num2str(k_order(k))]')
        grid on
        xlabel('x')
        ylabel('u(x)')
        title('P2: Analitical Solution Verification') % TODO is verification the right word? is this a manufactured soln?
        close(gcf)
        %}

    end
    plot(L./n_el, norm_mat(:,k), 'DisplayName', ['k = ',num2str(k_order(k))])

end

% plot convergence
xlabel('Element Spacing [m]')
ylabel('H1 Norm')
legend
set(gca,'XScale','log');
set(gca,'YScale','log');
grid on
axis equal
title('Grid Convergence for Manufactured Solution')


% return
%}
%% Problem 3: Rod Heating
% {
k = 1;
n_el = 100;
kappa = @(x) 385; % [w/m/c]
R = 0.025; % [m]
h_max = 1500; % [w/m^2]
L = 2;
f = @(x) 2/R*h_max*exp(-100*(x/L - 0.5).^2); % TODO change back
T_end = 30; % [c]
g_0 = T_end;
g_L = T_end;


% execute solution
[x_rod,u_rod, du_fem] = model_1d(k, n_el, kappa, f, g_0, g_L, L);

% plot
figure()
hold on
plot(x_rod, u_rod)
grid on
xlabel('x')
ylabel('T(x)')
title('P3: Rod Heating') % TODO is verification the right word? is this a manufactured soln?


% return
%}
%% Problem 4: Tensile Test
% {
k = 3;
n_el = 100;
E = 200e9;
A_max = 0.0001; % [m^2]
A_min = 0.00002; % [m^2]
L = 0.1;
kappa = @(x) E*(A_max - (A_max - A_min)*exp(-50*(x/L - 0.5)^2)); % [w/m/c]
f = @(x) 0; % TODO change back
u_end = 0.00002; % [c]
g_0 = 0;
g_L = u_end;


% execute solution
[x_rod, u_rod, du_fem] = model_1d(k, n_el, kappa, f, g_0, g_L, L);

% calculate stress
sigma = E*du_fem; % fix dimention issue

% plot
figure()
hold on
plot(x_rod, sigma/1e6) % convert to MPa
grid on
xlabel('x')
ylabel('\sigma (x) [MPa]')
title('P4: Tensile Test') % TODO is verification the right word? is this a manufactured soln?
%}
