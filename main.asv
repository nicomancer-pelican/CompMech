% np3217 01333401
% Computational Mechanics Coursework

% Main run file - assembly i.e. global strain matrix

% this calls element functions within a for loop (running through all
% elements) - uses index "e" for an arbitary number of N elemets
% this uses given functions:
% linear element stiffness: Kmat.m
% non linear element force: Fgeom.m
% first derivatives of Fgeom for Newton Raphson: Kgeom.m

% the bottom layer called from the above includes the shape functions (and
% material properties). shape.m, shapeder.m, shapder2.m to find matrix N 
% and its first and second derivatives.

%%
% inputs: [L, EA, EI q, N] = [beam length, axial stiffness, bending 
% stiffness, distibuted force, number of elements]
% Additionally: whether to use linear or non linear solver
% Are the BCs hard coded?

%% main driver
% ask user linear or non-linear
% get required inputs from user
% run appropriate solver
% display results

clear; clc; close all;
L = 1;
EA = 1e5;
EI = 1e2;
q = 1e3;
N = 200;

displacements = linearFE(L, EI, q, N);
displacements_NL = nonLinearFE(L, EA, EI, q, N);

% split rho into vertical displacements and moments
disc = 0:L/(N-2):1;
% linear
vertical  = displacements_NL(1:2:2*N-2);
rotations = displacements_NL(2:2:2*N-2);
moments   = gradient(displacements_NL(2:2:2*N-2))*EI;
% non linear
vertical_NL  = displacements_NL(1:2:2*N-2);
rotations_NL = displacements_NL(2:2:2*N-2);
moments_NL   = gradient(displacements_NL(2:2:2*N-2))*EI;

% analytical
x   = 0:0.01:L;
w   = (q/(24*EI)) * x.^2 .*(L^2 - 2*L*x + x.^2);
wp  = (q/(24*EI)) *(2*x*L^2 - 6*L*x.^2 + 4*x.^3);
wpp = (q/(24*EI)) *(2*L^2 - 12*L.*x + 12*x.^2);

% plots
figure; hold on; grid on
plot(x, w);
plot(disc, vertical);
plot(disc, vertical_NL);
legend('analytical', 'FE Linear', 'FE Non Linear')

figure; hold on; grid on;
plot(x, wp);
plot(disc, rotations);
plot(disc, rotations_NL);
legend('analytical', 'FE', 'FE Non Linear')

figure; hold on; grid on;
plot(x, wpp);
plot(disc, moments);
plot(disc, moments_NL);
legend('analytical', 'FE', 'FE Non Linear')

%% LINEAR FE SOLVER
function rho = linearFE(L, EI, q, N)
    % preallocate global stiffness matrix and global force vector 
    K   = zeros(2*N + 2, 2*N + 2);
    F   = zeros(2*N + 2, 1);
    
    % element properties
    L_e = L / N; %element length assuming equal length elements

    % stiffness matrix
    K = globalK(K, EI, N, L_e);
    
    % forces vector
    F = globalF(F, q, N, L_e);
    
    % apply BCs - fully fixed so remove 2 DoF at either end
    K = K(3:2*N, 3:2*N);
    F = F(3:2*N);
    
    % solve for unknown displacements and rotations
    rho = K\F;
end



%% NON-LINEAR FE SOLVER
function rho = nonLinearFE(L, EA, EI, q, N)
    % preallocate global stiffness matrices and global force vectors and
    % global displacement vector
    Km = zeros(2*N + 2, 2*N + 2);
    Kg = zeros(2*N + 2, 2*N + 2);
    Fm = zeros(2*N + 2, 1);
    Fg = zeros(2*N + 2, 1);
    rho_old = zeros(2*N + 2, 1);
    
    % element properties
    L_e = L / N; %element length assuming equal length elements
    
    % set arbitrary error to jumpstart loop
    eps = 1;
    
    % setup
    Km = globalK(Km, EI, N, L_e);
    Fm = globalF(Fm,  q, N, L_e);
    Kg = globalKgeom(rho_old, Kg, EI, N, L_e);
    Fg = globalFgeom(rho_old, Fg, EA, N, L_e);
    
    % apply BCs
    KmBC  = Km(3:2*N, 3:2*N);
    KgBC  = Kg(3:2*N, 3:2*N);
    FmBC  = Fm(3:2*N);
    FgBC  = Fg(3:2*N);
    rho_old_BC = rho_old(3:2*N);
    
    % functional
    G_old = KmBC*rho_old_BC - FmBC - FgBC;
    G_old_prime = KmBC + KgBC;
    
    while eps > 0.0001      
        % Newton Raphson iteration
        rho_new_BC = rho_old_BC - G_old_prime\G_old;
        
        % update rho (replace lost BC elements)
        rho_new = [0; 0; rho_new_BC; 0; 0];
        
        % calculate new Fg
        Fg = globalFgeom(rho_new, Fg, EA, N, L_e);
        Kg = globalKgeom(rho_new, Kg, EI, N, L_e);
        
        % apply BCs again
        KgBC = Kg(3:2*N, 3:2*N);
        FgBC = Fg(3:2*N);
        rho_new_BC = rho_new(3:2*N);
        
        % new functional
        G_new = KmBC*rho_new_BC - FmBC - FgBC;
        G_new_prime = KmBC + KgBC;
        
        % error calculation
        eps = max(abs(G_new - G_old));
        
        % update
        G_old = G_new;
        G_old_prime = G_new_prime;
        rho_old_BC = rho_new_BC;
    end
    rho = rho_new_BC;
end








