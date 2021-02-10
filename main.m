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
N = 100;

% displacements = linearFE(L, EA, EI, q, N);
displacements = nonLinearFE(L, EA, EI, q, N);

vertical = displacements(1:2:2*N-2);
moments  = gradient(displacements(2:2:2*N-2)*EI);


%% plots
x = 0:0.01:L;
w = 0.0001*((q/24*EI) * x.^2 .*(L^2 - 2*L*x + x.^2));
% wp = (q/12*EI) .*x .*(L^2 - 2*L*x + x.*2) + (q/24*EI) .* (-2*L + 2*x);
wpp = (q/EI)*(0.5 .*x.^2 - 0.5 * L .* x + (1/12)*L^2);
figure; hold on; grid on
plot(x, w);
plot([0:L/(N-2):1], vertical);
% plot(wp);
figure; hold on; grid on;
plot(x, wpp);
plot([0:L/(N-2):1], moments);

%% LINEAR FE SOLVER
function rho = linearFE(L, EA, EI, q, N)
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
    Km_new = zeros(2*N + 2, 2*N + 2);
    Kg_new = zeros(2*N + 2, 2*N + 2);
    Fm = zeros(2*N + 2, 1);
    Fg = zeros(2*N + 2, 1);
    rho = zeros(2*N + 2, 1);
    
    % element properties
    L_e = L / N; %element length assuming equal length elements
    
    % set arbitrary error to jumpstart loop
    eps = 1;
    
    % power counter for G_old
    i = 0;
    
    while eps > 0.01
        % linear part
        Km = globalK(Km_new, EI, N, L_e);
        Fm = globalF(Fm, q, N, L_e);
    
        % non linear part
        Kg = globalKgeom(rho, Kg_new, EA, N, L_e);
        Fg = globalFgeom(rho, Fg, EA, N, L_e);

        % apply BCs - fully fixed so remove 2 DoF at either end
        Km = Km(3:2*N, 3:2*N);
        Kg = Kg(3:2*N, 3:2*N);
        Fm = Fm(3:2*N);
        rho = rho(3:2*N);
        
        % initial functional G = K_M * rho - F_M
        G = Km*rho + ((Km + Kg)^i)*rho- Fm;
        
        % run Newton Raphson iteration
        rho_new = NewtonRaphson(rho, G, Km, Kg, Fm, EA, L);
        
        % update rho (replace lost BC elements)
        rho = [0; 0; rho_new; 0; 0];
        
        % calculate new Fg
        Fg = globalFgeom(rho, Fg, EA, N, L_e);
        
        % apply BC to Fg
        Fg = Fg(3:2*N);
        
        % new functional
        G_new = Km*rho_new - Fm - Fg;
        
        % error calculation
        eps = abs(G_new - G);
        
        % update other values
        Km_new = zeros(2*N + 2, 2*N + 2);
        Km_new(3:2*N, 3:2*N) = Km;
        Kg_new = zeros(2*N + 2, 2*N + 2);
        Kg_new(3:2*N, 3:2*N) = Kg;
        Fm = [0; 0; Fm; 0; 0];
        Fg = [0; 0; Fg; 0; 0];

        i = i + 1;
    end
end




%% UI sandbox
% 
% % setup window with grid, inputs panel and axes for results plot
% fig = uifigure('Name', 'AERO96015 - FE Solver');
% fig.Position = [200 200 1000 500];
% grid1 = uigridlayout(fig,[2 1]);
% grid1.RowHeight = {'1x',350};
% p = uipanel(grid1,'Title','FE Solver');
% ax = uiaxes(grid1);
% ax.XGrid = 'on'; ax.YGrid = 'on';
% 
% % grid in the inputs panel
% grid2 = uigridlayout(p,[4 3]);
% grid2.RowHeight = {22,22,22,22,22,22};
% grid2.ColumnWidth = {'1x','1x','1x','1x','1x','1x'};
% 
% % solver label
% solver = uilabel(grid2);
% solver.HorizontalAlignment = 'right';
% solver.Text = 'Solver';
% 
% % solver drop-down
% solver = uidropdown(grid2);
% solver.Items = {'Linear', 'Geometrically Non-Linear'};
% 
% % beam length label
% L = uilabel(grid2);
% L.HorizontalAlignment = 'right';
% L.Text = 'Beam length (m)';
% 
% % beam length field
% L = uieditfield(grid2, 'numeric');
% L.Value = 1;
% 
% % axial stiffness label
% EA = uilabel(grid2);
% EA.HorizontalAlignment = 'right';
% EA.Text = 'Axial stiffness, EA (N)';
% 
% % axial stiffness field
% EI = uieditfield(grid2, 'numeric');
% EI.Value = 10e5;
% 
% % bending stiffness label
% EI = uilabel(grid2);
% EI.HorizontalAlignment = 'right';
% EI.Text = 'Bending stiffness, EI (N m^2)';
% 
% % bending stiffness field
% EI = uieditfield(grid2, 'numeric');
% EI.Value = 10e2;
% 
% % distributed force label
% q = uilabel(grid2);
% q.HorizontalAlignment = 'right';
% q.Text = 'Distributed force (kN/m)';
% 
% % distributed force field
% q = uieditfield(grid2, 'numeric');
% q.Value = 1;
% 
% % no. elements label
% N = uilabel(grid2);
% N.HorizontalAlignment = 'right';
% N.Text = 'Number of elements';
% 
% % no. elements field
% N = uieditfield(grid2, 'numeric');
% N.Value = 1;
% 
% % Create a push button
% btn = uibutton(fig,'push',...
%                'ButtonPushedFcn', @(btn,event) linearFE(btn, L.Value, EA.Value, EI.Value, q.Value, N.value));
