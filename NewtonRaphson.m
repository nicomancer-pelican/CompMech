function rho = NewtonRaphson(rho_old, G_old, Km, Kg, Fm, EA, L)
%     while eps > 0.01
        % tangent stiffness K_T = K_M + K_G
%         Kg = Kgeom(rho_old, EA, L);
        Kt = Km + Kg;
        
        % update rho
        rho_new = rho_old - Kt\G_old;
        
        % find delta rho
        delta_rho = rho_new - rho_old;
        
        rho = rho_new;
        
%         % update functional
%         G_new = Km*rho_new - Fm - Fgeom(rho_new, EA, L);
%         
%         % error calculation
%         eps = abs(G_new - G_old);
%         
%         % update values
%         G_old = G_new;
%         rho_old = rho_new;
%         
%         % power counter
%         i = i+1;
%     end
end

%should output r as the displacements and K_t as the geometric stiffness

% eof