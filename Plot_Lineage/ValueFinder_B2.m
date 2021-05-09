function [ Val_1, Val_2 ] = ValueFinder_B2(clone, ChN, opt_C, algorithm)
% Extract from a clone entry the value to use for plotting 
% "Line type" = 2 - bi-value
%--------------------------------------------------------------------------

Val_1 = [];         % upper value (in y-axis of lineage plot)
Val_2 = [];         % lower value (in y-axis of lineage plot)

% Evaluate the Value to apply inside the square    
switch opt_C
    case 1      % Old_PL / New_PL
        switch algorithm
         case 1
            Val_1 = mean(clone.Fluor_Chan(ChN).IC( clone.Mask.MembPole_1 )) ;
            Val_2 = mean(clone.Fluor_Chan(ChN).IC( clone.Mask.MembPole_2 )) ;
         case 2
            Val_1 = mean(clone.Fluor_Chan(ChN).IC( clone.Fluor_Chan(ChN).Mask.Pole_1 )) ;
            Val_2 = mean(clone.Fluor_Chan(ChN).IC( clone.Fluor_Chan(ChN).Mask.Pole_2 )) ;
        end 

    case 2      % New Pole / Cytosol
        switch algorithm
         case 1
            Val_1 = mean(clone.Fluor_Chan(ChN).IC( clone.Mask.MembPole_2 )) ;
            Val_2 = mean(clone.Fluor_Chan(ChN).IC( clone.Mask.Cytosol  ))  ;
         case 2
            Val_1 = mean(clone.Fluor_Chan(ChN).IC( clone.Fluor_Chan(ChN).Mask.Pole_2 )) ;
            Val_2 = mean(clone.Fluor_Chan(ChN).IC( clone.Fluor_Chan(ChN).Mask.Cytosol ))  ;
        end

    case 3      % Old_PL / Cytosol
        switch algorithm
         case 1
            Val_1 = mean(clone.Fluor_Chan(ChN).IC( clone.Mask.MembPole_1 )) ;
            Val_2 = mean(clone.Fluor_Chan(ChN).IC( clone.Mask.Cytosol  ))  ;
         case 2
            Val_1 = mean(clone.Fluor_Chan(ChN).IC( clone.Fluor_Chan(ChN).Mask.Pole_1 )) ;
            Val_2 = mean(clone.Fluor_Chan(ChN).IC( clone.Fluor_Chan(ChN).Mask.Cytosol ))  ;
        end

    case 4     % Poles / Cytosol
        switch algorithm
         case 1
            Val_1 = mean(clone.Fluor_Chan(ChN).IC( logical(clone.Mask.MembPole_1 + clone.Mask.MembPole_2) )) ;
            Val_2 = mean(clone.Fluor_Chan(ChN).IC(  clone.Mask.Cytosol ))  ;
         case 2
            Val_1 = mean(clone.Fluor_Chan(ChN).IC( logical(clone.Fluor_Chan(ChN).Mask.Pole_1 + clone.Fluor_Chan(ChN).Mask.Pole_2) )) ;
            Val_2 = mean(clone.Fluor_Chan(ChN).IC(  clone.Fluor_Chan(ChN).Mask.Cytosol ))  ;
        end            

    case 5     % All Membrane / Cytosol
        if algorithm == 1
            sumMask = logical(clone.Mask.MembLateral_1 + clone.Mask.MembLateral_2 ...
                            + clone.Mask.MembPole_1   + clone.Mask.MembPole_2);                  
            Val_1 = mean(clone.Fluor_Chan(ChN).IC( sumMask )) ;
            Val_2 = mean(clone.Fluor_Chan(ChN).IC(  clone.Mask.Cytosol )) ;          
        end

    case 6     % Lateral Membrane / Cytosol
        if algorithm == 1                
            Val_1 = mean(clone.Fluor_Chan(ChN).IC( logical(clone.Mask.MembLateral_1 + clone.Mask.MembLateral_2) )) ;
            Val_2 = mean(clone.Fluor_Chan(ChN).IC(  clone.Mask.Cytosol )) ;                  
        end
    case 7      % Lateral Membrane / Poles
        if algorithm == 1                 
            Val_1 = mean(clone.Fluor_Chan(ChN).IC( logical(clone.Mask.MembLateral_1 + clone.Mask.MembLateral_2) )) ;
            Val_2 = mean(clone.Fluor_Chan(ChN).IC( logical(clone.Mask.MembPole_1   + clone.Mask.MembPole_2) )) ; 

        end         
end %Opt_C
    
end