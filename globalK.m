%GLOBALK global stiffness matrix. 
%   globablK (K, EI, N, L_e) places the element stiffness matrix into the
%   appropriate position in the global stiffness matrix. Runs for every
%   element iteratively. Once all elements have been iterated
%   over, the global stiffness matrix is complete.
%   INPUTS: K   --> global stiffness matrix
%           EI  --> beam bending stiffness
%           N   --> number of elements
%           L_e --> element length
%   OUTPUT: K  --> updated global stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = globalK(K, EI, N, L_e)
    for e = 1:N
        Ke = Kmat(EI, L_e);
        K(2*e-1 : 2*e+2, 2*e-1 : 2*e+2) = K(2*e-1 : 2*e+2, 2*e-1 : 2*e+2) + Ke;
    end
end

% eof