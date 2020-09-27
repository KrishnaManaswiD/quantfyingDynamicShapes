function [A0, C0] = calc_dc_components(chainCode)

    % Calculate DC components.
    % Input: 
    %   chain code (chain code)
    % Output: 
    %   A0 and C0 are bias coefficeis, corresponding to a frequency of zero.

    %% Maximum length of chain code
    k = size(chainCode, 2);
    
    %% Traversal length and change in projection length
    l = calc_traversal_length(chainCode);
    s = calc_change_in_projections(chainCode);
    
    %% Basic period of the chain code
    L = l(k);
    
    %% DC Components: A0, C0
    sum_a0 = 0;
    sum_c0 = 0;
    
    for p = 1 : k     

        delta_d = calc_change_in_projections(chainCode(p));
        delta_x = delta_d(:,1);
        delta_y = delta_d(:,2);
        delta_l = calc_traversal_length(chainCode(p));

        if (p > 1)       
                zeta = s(p - 1, 1) - delta_x / delta_l * l(p - 1);
                delta = s(p - 1, 2) - delta_y / delta_l * l(p - 1);
        else
                zeta = 0;
                delta = 0;
        end

        if (p > 1)
            sum_a0 = sum_a0 + delta_x / (2 * delta_l) * ((l(p))^2 - (l(p - 1))^2) + zeta * (l(p) - l(p-1));
            sum_c0 = sum_c0 + delta_y / (2 * delta_l) * ((l(p))^2 - (l(p - 1))^2) + delta * (l(p) - l(p-1));
        else
            sum_a0 = sum_a0 + delta_x / (2 * delta_l) * (l(p))^2 + zeta * l(p);
            sum_c0 = sum_c0 + delta_y / (2 * delta_l) * (l(p))^2 + delta * l(p);
        end
          
    end
    
    %% Assign  to output
    A0 = sum_a0 / L;
    C0 = sum_c0 / L;
end
