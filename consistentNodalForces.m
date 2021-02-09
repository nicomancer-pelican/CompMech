%CONSISTENTNODALFORCES Equivalent nodal forces and moments
%   consistentNodalForces (q,L) finds the equivalent nodal forces for a
%   consistend loading q on an Euler beam element of length L. The beam has
%   4 degrees of freedom, given by:
%     i,j=1: displacement at first node; 
%     i,j=2: rotation at first node;
%     i,j=3: displacement at second node;
%     i,j=4: rotation at second node.
%   INPUTS: q --> consistent force per unit length [kN/m]
%           L --> element length
%   OUTPUT: nodalForces --> equivalent nodal forces and moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nodalForces = consistentNodalForces(q, L)
    nodalForces = zeros(4,1);
    for i=1:4
      Ffun = @(x) shape(x,L,i) * q;
      nodalForces(i) = quadl(Ffun,0,L);
    end
end

% eof