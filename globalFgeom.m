%GLOBALFGEOM global geometric force vector
%   globablFgeom (F, Fe, e) places the element force vector into the
%   appropriate position in the global force vector. Designed to be run
%   for every element iteratively. Once all elements have been iterated
%   over, the global force vector is complete.
%   INPUTS: Fg --> global force vector
%           Fe --> element force vector
%           e  --> element number
%   OUTPUT: Fg --> updated global force vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fg = globalFgeom(rho, Fg, EA, N, L_e)
    for e = 1:N
        Fe = Fgeom(rho, EA, L_e);
        Fg(2*e-1: 2*e+2) = Fg(2*e-1: 2*e+2) + Fe;
    end
end

% eof