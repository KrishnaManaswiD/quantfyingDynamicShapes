function coefficients = fourier_approx(chainCode, nHarmonics, shouldNormalize)

    % This function generates coefficients of fourier approximation, given a chain code.
    % Input: 
    %   chain code (chainCode)
    %   number of harmonics (nHarmonics)
    %   whether to normalise or not (shouldNormalize)
    % Output:
    %   coeffients of harmonics n+1x4 matrix
    %   first row is [A0 0 C0 0]

    a = zeros(nHarmonics,1);
    b = zeros(nHarmonics,1);
    c = zeros(nHarmonics,1);
    d = zeros(nHarmonics,1);

    for i = 1 : nHarmonics       % loop over each harmonic
        harmonic_coeff = calc_harmonic_coefficients(chainCode, i);
        a(i,1) = harmonic_coeff(1, 1);
        b(i,1) = harmonic_coeff(1, 2);
        c(i,1) = harmonic_coeff(1, 3);
        d(i,1) = harmonic_coeff(1, 4);
    end

    [A0, C0] = calc_dc_components(chainCode); % bias components corresponding to zero frequency

    % Normalization procedure
    if (shouldNormalize == 1)   
        % Remove DC components
        A0 = 0;
        C0 = 0;
        
        % Compute theta1
        theta1 = 0.5 * atan(2 * (a(1) * b(1) + c(1) * d(1)) / ...
             	 (a(1)^2 + c(1)^2 - b(1)^2 - d(1)^2));
       
        costh1 = cos(theta1);
        sinth1 = sin(theta1);
             	 
        a_star_1 = costh1 * a(1) + sinth1 * b(1);
        b_star_1 = -sinth1 * a(1) + costh1 * b(1);
        c_star_1 = costh1 * c(1) + sinth1 * d(1);
        d_star_1 = -sinth1 * c(1) + costh1 * d(1);
       
        % Compute psi1 
        psi1 = atan(c_star_1 / a_star_1) ;
        
        % Compute E
        E = sqrt(a_star_1^2 + c_star_1^2);
        
        cospsi1 = cos(psi1);
        sinpsi1 = sin(psi1);
        
        for i = 1 : nHarmonics
            normalized = [cospsi1 sinpsi1; -sinpsi1 cospsi1] * [a(i) b(i); c(i) d(i)] * ... 
                        [cos(theta1 * i) -sin(theta1 * i); sin(theta1 * i) cos(theta1 * i)];
       
            a(i,1) = normalized(1,1) / E;
            b(i,1) = normalized(1,2) / E;
            c(i,1) = normalized(2,1) / E;
            d(i,1) = normalized(2,2) / E;
        end
        
    end  % end if normalized
    
    coefficients = [A0, 0, C0, 0];
    coefficients = [coefficients; [a, b, c, d]];
end