%GLOBALF global force vector
%   globablF (F, Fe, e) places the element force vector into the
%   appropriate position in the global force vector. Designed to be run
%   for every element iteratively. Once all elements have been iterated
%   over, the global force vector is complete.
%   INPUTS: F  --> global force vector
%           Fe --> element force vector
%           e  --> element number
%   OUTPUT: K  --> updated global force vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = globalF(F, Fe, e)
    for i = 1:4
        F(e*2 - 2 + i) = F(e*2 - 2 + i) + Fe(i);
    end
end

% eof