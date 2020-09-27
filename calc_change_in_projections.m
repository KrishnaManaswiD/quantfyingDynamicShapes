function output = calc_change_in_projections(chainCode)

% returns the changes in x and y projections of the chain, as a link in the
% chain code is traversed
%
% input: chain code (chainCode)
% output: 
%    if size(chainCode) is 1, [Dx, Dy]
%    if size(chainCode) > 1, ith row of ouput is [Sum_1toi(Dx), Sum_1toi(Dy)]

    sum_Dx = 0;
    sum_Dy = 0;
    
    p = zeros(size(chainCode,2),2);
    
    for i = 1 : size(chainCode, 2)        
        sum_Dx = sum_Dx + sign(6 - chainCode(i)) * sign(2 - chainCode(i));
        sum_Dy = sum_Dy + sign(4 - chainCode(i)) * sign(chainCode(i));
        p(i, 1) = sum_Dx;
        p(i, 2) = sum_Dy;
    end
    
    output = p;
    
end
