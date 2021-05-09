function [Val_A , Val_B] = ValueFinder_B3(clone, ChN, opt_C, algorithm)
% Extract from a clone entry the value to use for plotting 
% "Line type" = 3 - single value Ratio
%--------------------------------------------------------------------------   
    
Val_A = [];         % Numerator (of the Ratio)
Val_B = [];         % Denominator (of the Ratio)

% Evaluate the Value to apply inside the square
switch opt_C
    case 1      % Cytosol / Poles
        switch algorithm
         case 1
            Val_A = mean(clone.Fluor_Chan(ChN).IC( clone.Mask.Cytosol )) ;
            Val_B = mean(clone.Fluor_Chan(ChN).IC( logical(clone.Mask.MembPole_1 + clone.Mask.MembPole_2) )) ;
         case 2
            Val_A = mean(clone.Fluor_Chan(ChN).IC( clone.Fluor_Chan(ChN).Mask.Cytosol )) ;
            Val_B = mean(clone.Fluor_Chan(ChN).IC( logical(clone.Fluor_Chan(ChN).Mask.Pole_1 + clone.Fluor_Chan(ChN).Mask.Pole_2) )) ;
        end

    case 2      % Poles / Cytosol
        switch algorithm
         case 1
            Val_A = mean(clone.Fluor_Chan(ChN).IC( logical(clone.Mask.MembPole_1 + clone.Mask.MembPole_2) )) ;
            Val_B = mean(clone.Fluor_Chan(ChN).IC( clone.Mask.Cytosol )) ;                
         case 2
            Val_A = mean(clone.Fluor_Chan(ChN).IC( logical(clone.Fluor_Chan(ChN).Mask.Pole_1 + clone.Fluor_Chan(ChN).Mask.Pole_2) )) ;
            Val_B = mean(clone.Fluor_Chan(ChN).IC( clone.Fluor_Chan(ChN).Mask.Cytosol )) ;
        end

    case 3      % Old Pole / New Pole
        switch algorithm
         case 1
            Val_A = mean(clone.Fluor_Chan(ChN).IC( clone.Mask.MembPole_1 )) ;
            Val_B = mean(clone.Fluor_Chan(ChN).IC( clone.Mask.MembPole_2 )) ;                
         case 2
            Val_A = mean(clone.Fluor_Chan(ChN).IC( clone.Fluor_Chan(ChN).Mask.Pole_1 )) ;
            Val_B = mean(clone.Fluor_Chan(ChN).IC( clone.Fluor_Chan(ChN).Mask.Pole_2 )) ;
        end

    case 4      % New Pole / Old Pole
        switch algorithm
         case 1
            Val_A = mean(clone.Fluor_Chan(ChN).IC( clone.Mask.MembPole_2 )) ;
            Val_B = mean(clone.Fluor_Chan(ChN).IC( clone.Mask.MembPole_1 )) ;                
         case 2
            Val_A = mean(clone.Fluor_Chan(ChN).IC( clone.Fluor_Chan(ChN).Mask.Pole_2 )) ;
            Val_B = mean(clone.Fluor_Chan(ChN).IC( clone.Fluor_Chan(ChN).Mask.Pole_1 )) ;
        end

    case 5      % Cytosol / Membrane
        if algorithm == 1
            Val_A = mean(clone.Fluor_Chan(ChN).IC(  clone.Mask.Cytosol )) ;
            sumMask = logical(clone.Mask.MembLateral_1 + clone.Mask.MembLateral_2 ...
                            + clone.Mask.MembPole_1   + clone.Mask.MembPole_2);  
            Val_B = mean(clone.Fluor_Chan(ChN).IC( sumMask)) ;
        end

    case 6      % Membrane / Cytosol
        if algorithm == 1                
            sumMask = logical(clone.Mask.MembLateral_1 + clone.Mask.MembLateral_2 ...
                            + clone.Mask.MembPole_1   + clone.Mask.MembPole_2);  
            Val_A = mean(clone.Fluor_Chan(ChN).IC( sumMask)) ;
            Val_B = mean(clone.Fluor_Chan(ChN).IC(  clone.Mask.Cytosol )) ;
        end 

    case 7      % Lateral Membrane / Poles
        if algorithm == 1
            Val_A = mean(clone.Fluor_Chan(ChN).IC( logical(clone.Mask.MembLateral_1 + clone.Mask.MembLateral_2) )) ;
            Val_B = mean(clone.Fluor_Chan(ChN).IC( logical(clone.Mask.MembPole_1   + clone.Mask.MembPole_2)   )) ;
        end

    case 8      % Poles / Lateral Membrane
        if algorithm == 1
            Val_A = mean(clone.Fluor_Chan(ChN).IC( logical(clone.Mask.MembPole_1   + clone.Mask.MembPole_2)   )) ;
            Val_B = mean(clone.Fluor_Chan(ChN).IC( logical(clone.Mask.MembLateral_1 + clone.Mask.MembLateral_2) )) ;
        end
end %Opt_C

end