function [Xs, Ys] = Create_Segment_SquareBox( CellLength, ff, yrow, Dst_pl, stepL, PlotOpt_C )
% Create_Segment_SquareBox = Create a set of rectangular polygon that
% represent a cell segmented compartments and plot the value(s) at time
% point ff as a filled box. 
% It work for only 1 channel at a time
%
%  CellLength = total number of compartments/segments for the given cell
%  yrow   = store the y position for plotting the current cc-th clone
%  Dst_pl = distances between ploted cell lines 
%  stepL  = the segment size
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

% Centered at the "New pole" (Lower)   
if PlotOpt_C == 1  ||  PlotOpt_C == 4   
    for bb = 0 : CellLength
        Xs(bb+1,:) = [ ff-1, ff, ff, ff-1, ff-1 ];
        Ys(bb+1,:) = [(yrow*Dst_pl)+(stepL*bb)     , (yrow*Dst_pl)+(stepL*bb),    ... % last point = first point
                      (yrow*Dst_pl)+(stepL*(bb+1)) , (yrow*Dst_pl)+(stepL*(bb+1)),    (yrow*Dst_pl)+(stepL*bb) ];
    end
 
% Centered at the "old pole" (Upper)       
elseif PlotOpt_C == 3  ||  PlotOpt_C == 6
    for bb = 0 : CellLength
        Xs(bb+1,:) = [ ff-1, ff, ff, ff-1, ff-1 ];
        Ys(bb+1,:) = [(yrow*Dst_pl)-(stepL*bb)     , (yrow*Dst_pl)-(stepL*bb),    ... % last point = first point
                      (yrow*Dst_pl)-(stepL*(bb+1)) , (yrow*Dst_pl)-(stepL*(bb+1)),    (yrow*Dst_pl)-(stepL*bb) ];
    end
   
% Centered in the middle of the cell    
elseif PlotOpt_C == 2  ||  PlotOpt_C == 5    
    % Normalize position according to cell length
    normY = (yrow*Dst_pl) - round(CellLength/2)*stepL ;
    for bb = 0 : CellLength
        Xs(bb+1,:) = [ ff-1, ff, ff, ff-1, ff-1 ];
        Ys(bb+1,:) = [normY+(stepL*bb)     , normY+(stepL*bb),     ... 
                      normY+(stepL*(bb+1)) , normY+(stepL*(bb+1)), ...
                      normY+(stepL*bb) ];
    end
end

end

% s(bb+1,:) = [(yrow*Dst_pl)-round(Dst_pl/2)+(stepL*bb)     , (yrow*Dst_pl)-round(Dst_pl/2)+(stepL*bb),     ... 
%                   (yrow*Dst_pl)-round(Dst_pl/2)+(stepL*(bb+1)) , (yrow*Dst_pl)-round(Dst_pl/2)+(stepL*(bb+1)), ...
%                   (yrow*Dst_pl)-round(Dst_pl/2)+(stepL*bb) ];