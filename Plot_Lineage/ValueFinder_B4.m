function [Value] = ValueFinder_B4(clone, ChN, opt_C)
% Extract from a clone entry the value to use for plotting 
% "Line type" = 4 - Segmentation
%--------------------------------------------------------------------------   
    
Value = [];         % Denominator (of the Ratio)

% Evaluate the Value to apply inside the square
if opt_C  <= 3          % Cell Segmentation
        Value = clone.Fluor_Chan(ChN).SegmentSig ;

elseif opt_C  >= 4      % Axial Segmentation
        Value = clone.Fluor_Chan(ChN).AxSig ;

end %Opt_C

end
