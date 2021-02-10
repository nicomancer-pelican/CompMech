%GLOBALKGEOM global geometric stiffness matrix. 
%   globablKgeom (K, Ke, e) places the element stiffness matrix into the
%   appropriate position in the global stiffness matrix. Designed to be run
%   for every element iteratively. Once all elements have been iterated
%   over, the global stiffness matrix is complete.
%   INPUTS: Kg --> global stiffness matrix
%           Ke --> element stiffness matrix
%           e  --> element number
%   OUTPUT: Kg  --> updated global stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Kg = globalKgeom(rho, Kg, EA, N, L_e)
    for e = 1:N
        Ke = Kgeom(rho, EA, L_e);
        Kg(2*e-1 : 2*e+2, 2*e-1 : 2*e+2) = Kg(2*e-1 : 2*e+2, 2*e-1 : 2*e+2) + Ke;
    end
end

% eof