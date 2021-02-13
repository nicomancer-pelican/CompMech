%GLOBALF global force vector
%   globablF (F, q, N, L_e) places the element force vector into the
%   appropriate position in the global force vector. Runs for every element
%   iteratively. Once all elements have been iterated
%   over, the global force vector is complete.
%   INPUTS: F   --> global force vector
%           q   --> applied force
%           N   --> number of elements
%           L_e --> element length
%   OUTPUT: F  --> updated global force vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = globalF(F, q, N, L_e)
    for e = 1:N
        Fe = consistentNodalForces(q, L_e);
        F(2*e-1: 2*e+2) = F(2*e-1: 2*e+2) + Fe;
    end
end

% eof