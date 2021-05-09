function [Value] = ValueFinder_B1(clone, ChN, opt_C, algorithm)
% Extract from a clone entry the value to use for plotting 
% "Line type" = 1 - single value
%--------------------------------------------------------------------------

Value = [];
% Evaluate the Value to apply inside the square    
switch opt_C
    case 1      % Whole Cell
        Value = mean(clone.Fluor_Chan(ChN).IC(  clone.Mask.Cell_body )) ; 

    case 2      % Cytosol
        if algorithm == 2
            Value = mean(clone.Fluor_Chan(ChN).IC(  clone.Fluor_Chan(ChN).Mask.Cytosol )) ; 
        else
            Value = mean(clone.Fluor_Chan(ChN).IC(  clone.Mask.Cytosol)) ;
        end

    case 3      % Old_PL                   
        switch algorithm
         case 1
            Value = mean(clone.Fluor_Chan(ChN).IC(  clone.Mask.MembPole_1 )) ; 
         case 2
            Value = mean(clone.Fluor_Chan(ChN).IC(  clone.Fluor_Chan(ChN).Mask.Pole_1 )) ;
        end

    case 4      % New_PL                    
        switch algorithm
         case 1
            Value = mean(clone.Fluor_Chan(ChN).IC(  clone.Mask.MembPole_2 )) ; 
         case 2
            Value = mean(clone.Fluor_Chan(ChN).IC(  clone.Fluor_Chan(ChN).Mask.Pole_2 )) ;
        end     

    case 5      % only algorithm 1 - Lateral Membrane  
        if algorithm == 1
            sumMask = logical(clone.Mask.MembLateral_1 + clone.Mask.MembLateral_2) ;
            Value = mean(clone.Fluor_Chan(ChN).IC( sumMask )) ; 
        end

    case 6      % only algorithm 1 - All Membrane                    
        if algorithm == 1
            sumMask = logical(clone.Mask.MembLateral_1 + clone.Mask.MembLateral_2 ...
                            + clone.Mask.MembPole_1   + clone.Mask.MembPole_2);
            Value = mean(clone.Fluor_Chan(ChN).IC( sumMask )) ; 
        end
end %Opt_C

end