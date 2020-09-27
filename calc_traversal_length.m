function output = calc_traversal_length(chainCode)

    % The length of each link is either 1 or sqrt(2) depending on the 
    % orientation of the link.
    %
    % Input: chain code (chainCode)
    % Output:
    %   if size(chainCode) is 1, returns length of that link
    %   if size(chainCode) > 1, ith element is accumulated length up to link i. 

    sum_deltaLength = 0;
    l = zeros(size(chainCode,2),1);
    for i = 1 : size(chainCode, 2)
        sum_deltaLength = sum_deltaLength + 1 + ((sqrt(2)-1)/2).*(1-(-1).^chainCode(i));
        l(i,1) = sum_deltaLength;
    end
    
    output = l';
end
