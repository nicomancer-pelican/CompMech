%GLOBALK global stiffness matrix. 
%   globablK (K, Ke, e) places the element stiffness matrix into the
%   appropriate position in the global stiffness matrix. Designed to be run
%   for every element iteratively. Once all elements have been iterated
%   over, the global stiffness matrix is complete.
%   INPUTS: K  --> global stiffness matrix
%           Ke --> element stiffness matrix
%           e  --> element number
%   OUTPUT: K  --> updated global stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = globalK(K, EI, N, L_e)
    for e = 1:N
        Ke = Kmat(EI, L_e);
        K(2*e-1 : 2*e+2, 2*e-1 : 2*e+2) = K(2*e-1 : 2*e+2, 2*e-1 : 2*e+2) + Ke;
    end
end

% eof