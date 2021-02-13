% np3217 01333401
% Computational Mechanics Coursework

clear; clc; close all;
L = 1;
EA = 1e5;
EI = 1e2;

%% Plots
q = [1,2,5];
N = [50:5:200];
linearError    = zeros(11,3);
nonlinearError = linearError;
for i = 1:31
    for k = 1:3
        % linear solver
        displacements = linearFE(L, EI, q(k), N(i));
        vertical  = displacements(1:2:2*N(i)-2);
        rotations = displacements(2:2:2*N(i)-2);
        % non linear solver
        displacements_NL = nonLinearFE(L, EA, EI, q(k), N(i));
        vertical_NL  = displacements_NL(1:2:2*N(i)-2);
        rotations_NL = displacements_NL(2:2:2*N(i)-2);
        % error calculation
        error_x = 0:L/(N(i)-2):L;
        error_w = (q(k)/(24*EI)) * error_x.^2 .*(L^2 - 2*L*error_x + error_x.^2);
        linearError(i,k)    = max(abs(vertical - error_w'));
        nonlinearError(i,k) = max(abs(vertical_NL - error_w'));
    end
end

%% Question 2
figure; hold on; grid on
% q = 1kN
loglog(linearError(:,1), N, 'x-', 'linewidth', 2, 'color', [0 0.4470 0.7410])
% q = 2kN
loglog(linearError(:,2), N, 'o--', 'linewidth', 2, 'color', [0 0.4470 0.7410])
% q = 5kN
loglog(linearError(:,3), N, '+:', 'linewidth', 2, 'color', [0 0.4470 0.7410])
% prettiness
title('Linear Solver Convergence', 'interpreter', 'latex')
legend('Linear FE Solver, $q = 1 kN$', 'Linear FE Solver, $q = 2 kN$', ...
    'Linear FE Solver, $q = 5 kN$', 'interpreter', 'latex')
xlabel('Error', 'interpreter', 'latex')
ylabel('No. Elements', 'interpreter', 'latex')

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%% Question 5 figure 1
figure; hold on; grid on
% q = 1kN
loglog(nonlinearError(:,1), N, 'x-', 'linewidth', 2, 'color', [0.8500 0.3250 0.0980])
% q = 2kN
loglog(nonlinearError(:,2), N, 'o--', 'linewidth', 2, 'color', [0.8500 0.3250 0.0980])
% q = 5kN
loglog(nonlinearError(:,3), N, '+:', 'linewidth', 2, 'color', [0.8500 0.3250 0.0980])
% prettiness
title('Nonlinear Solver Convergence', 'interpreter', 'latex')
legend('Nonlinear FE Solver, $q = 1 kN$', 'Nonlinear FE Solver, $q = 2 kN$', ...
    'Nonlinear FE Solver, $q = 5 kN$', 'interpreter', 'latex')
xlabel('Error', 'interpreter', 'latex')
ylabel('No. Elements', 'interpreter', 'latex')

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%% Question 5 figure 2
figure; hold on; grid on
% q = 1kN
loglog(linearError(:,1), N, 'x-', 'linewidth', 2, 'color', [0 0.4470 0.7410])
loglog(nonlinearError(:,1), N, 'x-', 'linewidth', 2, 'color', [0.8500 0.3250 0.0980])
% q = 2kN
loglog(linearError(:,2), N, 'o--', 'linewidth', 2, 'color', [0 0.4470 0.7410])
loglog(nonlinearError(:,2), N, 'o--', 'linewidth', 2, 'color', [0.8500 0.3250 0.0980])
% q = 5kN
loglog(linearError(:,3), N, '+:', 'linewidth', 2, 'color', [0 0.4470 0.7410])
loglog(nonlinearError(:,3), N, '+:', 'linewidth', 2, 'color', [0.8500 0.3250 0.0980])
% prettiness
title('Comparison of Convergence', 'interpreter', 'latex')
legend('Linear FE Solver, $q = 1 kN$', 'Nonlinear FE Solver, $q = 1 kN$'...
    , 'Linear FE Solver, $q = 2 kN$', 'Nonlinear FE Solver, $q = 2 kN$' ...
    , 'Linear FE Solver, $q = 5 kN$', 'Nonlinear FE Solver, $q = 5 kN$' ...
    , 'interpreter', 'latex')
xlabel('Error', 'interpreter', 'latex')
ylabel('No. Elements', 'interpreter', 'latex')

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
%% Question 3
clear;
L = 1;
EA = 1e5;
EI = 1e2;
q = [1,2,5];
N = 60;

x = 0:1/58:L;

for i = 1:3
    % linear solver
    displacements = linearFE(L, EI, q(i), N);
    vertical(:,i)  = displacements(1:2:2*N-2);
    rotations(:,i) = displacements(2:2:2*N-2);
    % analytical
    w(:,i)   = (q(i)/(24*EI)) * x.^2 .*(L^2 - 2*L*x + x.^2);
end

figure; hold on; grid on
% analytical
plot(x, w(:,1), 'linewidth', 2, 'color', [0.4660 0.6740 0.1880])
plot(x, w(:,2), 'linewidth', 2, 'color', [0.4660 0.6740 0.1880])
plot(x, w(:,3), 'linewidth', 2, 'color', [0.4660 0.6740 0.1880])
% linear
plot(x, vertical(:,1), 'linewidth', 2, 'color', [0.3010 0.7450 0.9330])
plot(x, vertical(:,2), 'linewidth', 2, 'color', [0.3010 0.7450 0.9330])
plot(x, vertical(:,3), 'linewidth', 2, 'color', [0.3010 0.7450 0.9330])
% prettiness)
title('Linear FE Solver Displacement Comparison', 'interpreter', 'latex')
legend('Analytical Result, $q = 1 kN$', 'Analytical Result, $q = 2 kN$'...
    , 'Analytical Result, $q = 5 kN$', 'Linear FE Solver, $q = 1 kN$'...
    , 'Linear FE Solver, $q = 2 kN$', 'Linear FE Solver, $q = 5 kN$'...
    , 'interpreter', 'latex')
xlabel('Length of Beam, m', 'interpreter', 'latex')
ylabel('Vertical Displacement, m', 'interpreter', 'latex')
