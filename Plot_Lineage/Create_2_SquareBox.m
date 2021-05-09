function [Xs, Ys] = Create_2_SquareBox( subDiv, ff, yrow, Dst_pl, hg_pl )
% Create_2_SquareBox = Create two rectangular polygons for plotting the
% values at time point ff as a filled box.
% Each channel is represented by two value at each time point
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
    % Divide the square-"cell" into 2 sub-areas
    Xs = [ ff-1, ff, ff, ff-1, ff-1 ] ;
    
    % upper half
    Ys{1,1} = [ (yrow*Dst_pl)+(hg_pl*1/2), (yrow*Dst_pl)+(hg_pl*1/2) , ... % last point = first point
              (yrow*Dst_pl)+(hg_pl*1)  , (yrow*Dst_pl)+(hg_pl*1),    (yrow*Dst_pl)+(hg_pl*1/2) ];
    % lower half
    Ys{1,2} = [ (yrow*Dst_pl)+(hg_pl*0)  , (yrow*Dst_pl)+(hg_pl*0) ,   ... % last point = first point
              (yrow*Dst_pl)+(hg_pl*1/2), (yrow*Dst_pl)+(hg_pl*1/2),  (yrow*Dst_pl)+(hg_pl*0) ];
       
elseif length(subDiv) == 2
    % Divide the square-"cell" into 2 sub-areas
    Xs = [ ff-1, ff, ff, ff-1, ff-1 ] ;
    
    % TOP upper half
    Ys{1,1} = [ (yrow*Dst_pl)+(hg_pl*3/4) , (yrow*Dst_pl)+(hg_pl*3/4) , ... % last point = first point
              (yrow*Dst_pl)+(hg_pl*1)   , (yrow*Dst_pl)+(hg_pl*1),    (yrow*Dst_pl)+(hg_pl*3/4) ];
    % TOP lower half 
    Ys{1,2} = [ (yrow*Dst_pl)+(hg_pl*2/4) , (yrow*Dst_pl)+(hg_pl*2/4) , ... % last point = first point
              (yrow*Dst_pl)+(hg_pl*3/4) , (yrow*Dst_pl)+(hg_pl*3/4) , (yrow*Dst_pl)+(hg_pl*2/4) ];
    
    % BOTTOM upper half      
    Ys{2,1} = [ (yrow*Dst_pl)+(hg_pl*1/4) , (yrow*Dst_pl)+(hg_pl*1/4) ,  ... % last point = first point
              (yrow*Dst_pl)+(hg_pl*2/4) , (yrow*Dst_pl)+(hg_pl*2/4) , (yrow*Dst_pl)+(hg_pl*1/4) ];
    % BOTTOM lower half
    Ys{2,2} = [ (yrow*Dst_pl)+(hg_pl*0)   , (yrow*Dst_pl)+(hg_pl*0)   ,  ... % last point = first point
              (yrow*Dst_pl)+(hg_pl*1/4) , (yrow*Dst_pl)+(hg_pl*1/4) , (yrow*Dst_pl)+(hg_pl*0) ];          
    
end

end
