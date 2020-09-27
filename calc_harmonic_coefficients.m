function output = calc_harmonic_coefficients(chainCode, nHarmonics)

    % This function returns the n-th set of four harmonic coefficients.
    % Input: 
    %   chain code (chainCode)
    %   number of harmonics (nHarmonics)
    % Output:
    %   coeffients of harmonics [an bn cn dn]

    k = size(chainCode, 2); % length of chain code
    
    l = calc_traversal_length(chainCode); % traversal length for each link
        
    L = l(k); % basic period = traveral length for entire chain (perimeter)
    
    % Store this value to make computation faster
    two_n_pi = 2 * nHarmonics * pi;
    
    %% Compute Harmonic cofficients: an, bn, cn, dn
    sigma_a = 0;
    sigma_b = 0;
    sigma_c = 0;
    sigma_d = 0;
        
    for p = 1 : k
        if (p > 1)
            lp_prev = l(p - 1);            
        else
            lp_prev = 0;
        end
        
        delta_d = calc_change_in_projections(chainCode(p));
        delta_x = delta_d(:,1);
        delta_y = delta_d(:,2);
        delta_l = calc_traversal_length(chainCode(p));
        
        q_x = delta_x / delta_l;
        q_y = delta_y / delta_l;
        
        sigma_a = sigma_a + q_x * (cos(two_n_pi * l(p) / L) - cos(two_n_pi * lp_prev / L));
        sigma_b = sigma_b + q_x * (sin(two_n_pi * l(p) / L) - sin(two_n_pi * lp_prev / L));
        sigma_c = sigma_c + q_y * (cos(two_n_pi * l(p) / L) - cos(two_n_pi * lp_prev / L));
        sigma_d = sigma_d + q_y * (sin(two_n_pi * l(p) / L) - sin(two_n_pi * lp_prev / L));   
    end
    
    r = L/(2*nHarmonics^2*pi^2);
    
    an = r * sigma_a;
    bn = r * sigma_b;
    cn = r * sigma_c;
    dn = r * sigma_d;
    
    %% Assign  to output
    output = [an bn cn dn];
end
