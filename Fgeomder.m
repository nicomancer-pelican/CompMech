%FGEOMDER 1st derivative of geometricallly-nonlinear structural function.
%   Fgeomder () 


function Fnlprime = Fgeomder(rhoe, EA, L)
%     for i=1:4
%         % Compute the derivative of the displacements in the element.
%         Wprime = @(x) shapeder(x,L,1)*rhoe(i) ...
%                +shapeder(x,L,2)*rhoe(i) ...
%                +shapeder(x,L,3)*rhoe(i) ...
%                +shapeder(x,L,4)*rhoe(i);
%         for j = 1:4
%             % Evaluate the elements of the geometric stiffness matrix.
%             Fnl = @(x) (-EA/2)*(Wprime(x).^3).*shapeder(x,L,j);
%             Fnlprime(j,i) = quadl(Fnl, 0, L);
%         end
%     end

%     Fnlprime = zeros(4,4);
%     for i = 1:4
%         Fun = @(x) (shapeder2(x,L,i)).^4;
%         Fnlprime(i,i) = (-3*EA/2) * rhoe(i)^2 * quadl(Fun,0,L);
%     end

    Fnl = Fgeom(rhoe, EA, L);
end