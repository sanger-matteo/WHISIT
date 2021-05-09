function [Xs, Ys] = Create_1_SquareBox( subDiv, ff, yrow, Dst_pl, hg_pl )
% Create_1_SquareBox = Create a rectangular polygon for plotting the
% value(s) at time point ff as a filled box.
% Each channel is represented by one value at each time point
%
%  Dst_pl = distances between ploted cell lines 
%  hg_pl  = heights of the plotted cell lines
%  yrow   = store the y position for plotting the current cc-th clone
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


if length(subDiv) == 1
    % Create the square for plotting the value at time point Frm
    Xs{1} = [ ff-1, ff, ff, ff-1, ff-1 ];
    Ys{1} = [ yrow*Dst_pl,  yrow*Dst_pl, (yrow*Dst_pl)+hg_pl, (yrow*Dst_pl)+hg_pl, yrow*Dst_pl ];
  
elseif length(subDiv) == 2
    % Divide the square-"cell" into 2 sub-areas
    Xs{1} = [ ff-1, ff, ff, ff-1, ff-1 ] ;
    Xs{2} = [ ff-1, ff, ff, ff-1, ff-1 ] ;
    
    % upper half
    Ys{1} = [ (yrow*Dst_pl)+(hg_pl*1/2), (yrow*Dst_pl)+(hg_pl*1/2) , ... % last point = first point
              (yrow*Dst_pl)+(hg_pl*1)  , (yrow*Dst_pl)+(hg_pl*1),    (yrow*Dst_pl)+(hg_pl*1/2) ];
    % lower half
    Ys{2} = [ (yrow*Dst_pl)+(hg_pl*0)  , (yrow*Dst_pl)+(hg_pl*0) ,   ... % last point = first point
              (yrow*Dst_pl)+(hg_pl*1/2), (yrow*Dst_pl)+(hg_pl*1/2),  (yrow*Dst_pl)+(hg_pl*0) ];
       
end

end
