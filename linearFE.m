%LINEARFE linear FE solver
%   LinearFE (L, EI, q, N) performs a linear finite element analysis on a
%   beam. The beam has 4 degrees of freedom, given by:
%     i,j=1: displacement at first node; 
%     i,j=2: rotation at first node;
%     i,j=3: displacement at second node;
%     i,j=4: rotation at second node. 
%   INPUTS: L  --> beam length
%           EI --> beam bending stiffness
%           q  --> beam constant loading force
%           N  --> number of elements
%   OUTPUT: rho  --> nodal displacements and rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rho = linearFE(L, EI, q, N)
    % preallocate global stiffness matrix and global force vector 
    K   = zeros(2*N + 2);
    F   = zeros(2*N + 2, 1);
    
    % element properties
    L_e = L / N; %element length assuming equal length elements

    % stiffness matrix
    K = globalK(K, EI, N, L_e);
    
    % forces vector
    F = globalF(F, q, N, L_e);
    
    % apply BCs - fully fixed so remove 2 DoF at either end
    K(1:2,1:2) = 0; K(2*N:2*N+2, 2*N:2*N+2) = 0;
    F(1:2,1)   = 0; F(2*N:2*N+2)            = 0;
    
    % solve for unknown displacements and rotations
    rho = K\F;
end

%eof