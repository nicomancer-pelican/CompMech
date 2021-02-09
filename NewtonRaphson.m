function r = NewtonRaphson(R_c, N, EA, L)
    % starting values
    R_0 = zeros(2*N, 1);
    r_0 = zeros(2*N, 1);
    
    % initial tangent stiffness
    K_t = Fgeomder(r_0, EA, L);
    
    % set arbitrary error to jumpstart loop
    eps = 1;
    while abs(eps) > 0.01
        % first estimate
        delta_R = R_c - R_0;
        delta_r = K_t\delta_R;

        % update displacement
        r_e = r_0 + delta_r;

        % use new displacement in function
        F_NL = Fgeom(r_e, EA, L);

        % check error
        eps = abs(F_NL - R_c);
        
        % new tangent stiffness
        K_t = Fgeomder(r_0, EA, L);
    end
    r = r_e;
end

%should output r as the displacements and K_t as the geometric stiffness

% eof