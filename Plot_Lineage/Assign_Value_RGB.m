function [RGB] = Assign_Value_RGB( Value, RGB_range ) 
%
%Assign_Value_RGB = Give a value and a colormap, return the appropriate
%   color to use for plotting the inputted value
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

    % In rare exception, the mask created during analysis is empty, due to
    % some weirdly shaped cells being unable to generate a proper mask.
    % Value will be NaN, resulting in fatal error. In such cases, we will
    % assign a "empty" grey color for filling
    if isnan( Value )
        RGB(1) = 0.5;     RGB(2) = 0.5;     RGB(3) = 0.5;
        return
    end   

    % use abs() because half will be negative values: we search for the one
    % where the difference is closest to zero.
    diff = (abs(RGB_range(:,4) - Value));
    idx_min = find(diff == min(diff));
    if length(idx_min) >= 1           % in the rare exception when there are two ...    
        idx_min = idx_min(1);         % equal minima diff, take just the first
    end
    R = RGB_range(idx_min,1);      G = RGB_range(idx_min,2);    B = RGB_range(idx_min,3); 
    if Value > RGB_range(end,4) 
        R = RGB_range(end,1);      G = RGB_range(end,2);        B = RGB_range(end,3);        
    elseif Value < RGB_range(1,4)
        R = RGB_range(1,1);        G = RGB_range(1,2);          B = RGB_range(1,3);
    end
    RGB(1) = R/256;     RGB(2) = G/256;     RGB(3) = B/256;
    
end




