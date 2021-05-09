function [t_DB] = Retrive_mesh_microTrack(cellList)
%%-------------------- RETRIVE DATA from cellList --------------------------
% Retrive and reorganize data from Track.mat files created by MICROBRTRACKER 
%
% Therefore the INPUT passed to the function will be either :
%   - cellList      or      cellList.meshData
%
% Originally cellList has frames at first level, and in second you access 
% all cells' data present at that frame.
%%--------------------------------------------------------------------------

    for f = 1:length(cellList)                      % go thourgh each frame
    for c = 1:length(cellList{f})                   % access all cells
        if  isempty(cellList{1,f}{1,c}) == 0  &&  length(cellList{1,f}{1,c}.mesh) > 1    % if it is NOT empty

            % NB: save in my_cell(...,c) and not f because we retrive all 
            % cells data and save in our struct frame by frame
            t_DB(f).cell(c).mesh = cellList{1,f}{1,c}.mesh; 
            % first colums are all xs, and second are all ys 
            t_DB(f).cell(c).coord = round( [ [cellList{1,f}{1,c}.mesh(:,1) ; flipud(cellList{1,f}{1,c}.mesh(:,3))] ...
                                           , [cellList{1,f}{1,c}.mesh(:,2) ; flipud(cellList{1,f}{1,c}.mesh(:,4))] ]);
        
            xs =  [cellList{1,f}{1,c}.mesh(:,1) ; flipud( cellList{1,f}{1,c}.mesh(:,3)) ];           % we retrive the x-coordinates, putting 1st and 3rd mesh column (right and left part of the mesh) in same vector
            ys =  [cellList{1,f}{1,c}.mesh(:,2) ; flipud( cellList{1,f}{1,c}.mesh(:,4)) ];           % the same for y-coordinates

            % left and right side of the cell outline coordinates
            xr =  [cellList{1,f}{1,c}.mesh(:,1)] ; 
            xl =  [cellList{1,f}{1,c}.mesh(:,3)] ; 
            yr =  [cellList{1,f}{1,c}.mesh(:,2)] ; 
            yl =  [cellList{1,f}{1,c}.mesh(:,4)] ; 
            mx = [xr(1) ; (xr(2:end-1) + xl(2:end-1)) /2 ; xr(end)] ;
            my = [yr(1) ; (yr(2:end-1) + yl(2:end-1)) /2 ; yr(end)]  ;
            
            % Cell axis coordinates, following middle cell line: 
            % first colums are all xs, and second are all ys
            t_DB(f).cell(c).geom.axis = [mx , my] ;
            t_DB(f).cell(c).geom.length = sum( sqrt( (mx(2:end)-mx(1:end-1)).^2 + (my(2:end)-my(1:end-1)).^2 )) ;            
            t_DB(f).cell(c).geom.area = polyarea(xs,ys) ;
            t_DB(f).cell(c).geom.polarity =  cellList{1,f}{1,c}.polarity ;
               
        end
    end
    end
    
end


%--------------------------------------------------------------------------
%    Copyright (c) 2016; WHISIT, Matteo Sangermani, All rights reserved
%--------------------------------------------------------------------------

