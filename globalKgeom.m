%GLOBALKGEOM global geometric stiffness matrix. 
%   globablKgeom (rho, Kg, EA, N, L_e) places the element stiffness matrix
%   into the appropriate position in the global stiffness matrix. Runs for
%   every element iteratively. Once all elements have been iterated over,
%   the global stiffness matrix is complete.
%   INPUTS: rho  --> displacement and rotation vector
%           Kg   --> global geometric stiffness matrix
%           EA   --> beam axial stiffness
%           N    --> number of elements
%           L_e  --> element length
%   OUTPUT: Kg  --> updated global geometric stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Kg = globalKgeom(rho, Kg, EA, N, L_e)
    for e = 1:N
        Ke = Kgeom(rho, EA, L_e);
        Kg(2*e-1 : 2*e+2, 2*e-1 : 2*e+2) = Kg(2*e-1 : 2*e+2, 2*e-1 : 2*e+2) + Ke;
    end
end

% eof