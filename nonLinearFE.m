%NONLINEARFE linear FE solver
%   nonLinearFE (L, EI, q, N) performs a nonlinear finite element analysis
%   on a beam. The beam has 4 degrees of freedom, given by:
%     i,j=1: displacement at first node; 
%     i,j=2: rotation at first node;
%     i,j=3: displacement at second node;
%     i,j=4: rotation at second node. 
%   INPUTS: L  --> beam length
%           EA --> beam axial stiffness
%           EI --> beam bending stiffness
%           q  --> beam constant loading force
%           N  --> number of elements
%   OUTPUT: rho  --> nodal displacements and rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rho = nonLinearFE(L, EA, EI, q, N)
    % preallocate global stiffness matrices and global force vectors and
    % global displacement vector
    Km = zeros(2*N + 2);
    Kg = zeros(2*N + 2);
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
    
    while eps > 0.1      
        Kg = globalKgeom(rho_old, Kg, EA, N, L_e);
        Fg = globalFgeom(rho_old, Fg, EA, N, L_e);

        % functional
        G_old   = Km*rho_old - Fm - Fg;
        G_old_p = Km + Kg;

        % apply BCS
        rho_old_BC = rho_old(3:2*N);
        G_old_BC   = G_old(3:2*N);
        G_old_p_BC = G_old_p(3:2*N, 3:2*N);
        
        % Newton Raphson iteration with BCs
        rho_new = [0; 0; rho_old_BC - G_old_p_BC\G_old_BC; 0; 0];
        
        % updates using rho_new
        Fg_new = globalFgeom(rho_new, Fg, EA, N, L_e);
        G_new = Km*rho_new - Fm - Fg_new;
        rho_old = rho_new;
        
        % error calculation
        eps = max(abs(G_new - G_old));
    end
    rho = rho_new;
end

%eof