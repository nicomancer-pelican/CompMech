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
EA = 10^5;
EI = 10^2;
q = 1000;
N = 2;

% displacements = linearFE(L, EA, EI, q, N);
displacements_NL = nonLinearFE(L, EA, EI, q, N);

vertical = displacements(1:2:2*N-2);
moments  = gradient(displacements(2:2:2*N-2)*EI);

figure
plot(vertical);

figure
plot(moments);

%% analytical
x = 0:0.01:L;
w = (q/24*EI) * x.^2 .*(L^2 - 2*L*x + x.^2);
wp = (q/12*EI) .*x .*(L^2 - 2*L*x + x.*2) + (q/24*EI) .* (-2*L + 2*x);
wpp = (q/EI)*(0.5 .*x.^2 - 0.5 * L .* x + (1/12)*L^2);
figure;
plot(w);
% plot(wp);
figure
plot(wpp);

%% LINEAR FE SOLVER
function rho = linearFE(L, EA, EI, q, N)
    % preallocate global stiffness matrix and global force vector 
    K   = zeros(2*N + 2, 2*N + 2);
    F   = zeros(2*N + 2, 1);
    
    % element properties
    L_e = L / N; %element length assuming equal length elements
    
    for e = 1:N
        % determine element stiffness matrix
        elementK = Kmat(EI, L_e);
        
        % place in global stiffness matrix
        K  = globalK(K, elementK, e);
        
        % element forces
        elementF = consistentNodalForces(q, L_e);
        
        % assemble global forces
        F = globalF(F, elementF, e);
    end
    
    % apply BCs - fully fixed so remove 2 DoF at either end
    K = K(3:2*N, 3:2*N);
    F = F(3:2*N);
    
    % solve for unknown displacements and rotations
    rho = K\F;
end



%% NON-LINEAR FE SOLVER
function rho = nonLinearFE(L, EA, EI, q, N)
    % preallocate global stiffness matrices and global force vectors
    Km = zeros(2*N + 2, 2*N + 2);
    Kg = zeros(2*N + 2, 2*N + 2);
    Fm = zeros(2*N + 2, 1);
    Fg = zeros(2*N + 2, 1);

    % element properties
    L_e = L / N; %element length assuming equal length elements
    
    for e = 1:N
        % determine element stiffness matrix
        elementKmat = Kmat(EI, L_e);
%         elementKgeo = NewtonRaphson(q, N, EA, L);
        
        % place in global stiffness matrix
        Km = globalK(Km, elementKmat, e);
%         Kg = globalK(Kg, elementKgeo, e);
        
        % element forces
        elementFmat = consistentNodalForces(q, L_e);
        elementFgeo = NewtonRaphson(q, N, EA, L);
        
        % assemble global forces
        Fm = globalF(Fm, elementFmat, e);
        Fg = globalF(Fg, elementFgeo, e);
    end
    
    % combine linear and nonlinear parts
    K = Km + Kg;
    F = Fm + Fg;
    
    % apply BCs - fully fixed so remove 2 DoF at either end
    K = K(3:2*N, 3:2*N);
    F = F(3:2*N);
    
    % solve for unknown displacements and rotations
    rho = K\F;
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
