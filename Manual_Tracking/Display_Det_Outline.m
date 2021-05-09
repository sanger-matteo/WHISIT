function Display_Det_Outline
%
%Display_Det_Outline - Show the cells' outline detected using Oufticolor as 
%   visual aid during tracking
%
%   ---> Notes about variables used:
%   CellTracks = a global variable where that stores the information of
%   every cell track. 
%   scc = stores the ID-number of the currently selected cell-track
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

global APP_opt;	    
global CellTracks;     
global scc;       
global oDet;

mk_sz = APP_opt.t5_display_Marker  ;    % Chosen marker size   
switch APP_opt.t5_display_Marker        % according to mk_sz value ...
    case 10 ;    lw = 1 ;               % choose line width for cells' outline
    case 14 ;    lw = 1.4 ;
    case 18 ;    lw = 1.7 ;
end 
basic_clr = [0 .8 0] ;                  % Starndard track RGB color (Green)  
selected_clr = [.9 .4 0] ;              % Use special RGB color for selected scc-th track 
ff = APP_opt.t5_ff;                     % presently displayed frame

dist = [];                              % store the distances of each cell outline from the 
                                        % track point of scc-th track in current frame

% Find the column position for the scc-th cell in CellTracks cell array       
idx_scc = find(cell2mat(CellTracks( 1,: )) == scc) ;
% Store XY_t point for selected cell scc-th in present ff-th frame
X_t = CellTracks{2,idx_scc}(APP_opt.t5_ff, 1) ;
Y_t = CellTracks{2,idx_scc}(APP_opt.t5_ff, 2) ;
% remove '0' coordinates - not-tracked points will not be considered
X_t(X_t==0) = NaN;     
Y_t(Y_t==0) = NaN;
  
% If Cell Outline is NOT "none", all cells will have otline and with standard color
if APP_opt.t5_display_Outline ~= 0
 
  % Go through all cell's outline detected in ff-th frame
  for kk = 1 : length( oDet.cellList.meshData{ff} )
    if   ~isempty(oDet.cellList.meshData{ff}{kk})  && ...
         ~isempty( oDet.cellList.meshData{ff}{kk}.mesh )  && ...
         size(oDet.cellList.meshData{ff}{kk}.mesh,2) == 4         
     
       % Store XY_o points for kk_th outline and plot it
       X_o = [oDet.cellList.meshData{ff}{kk}.mesh(:,1) ; flipud(  oDet.cellList.meshData{ff}{kk}.mesh(:,3)) ];
       Y_o = [oDet.cellList.meshData{ff}{kk}.mesh(:,2) ; flipud(  oDet.cellList.meshData{ff}{kk}.mesh(:,4)) ];
       plot( X_o , Y_o , '-', 'Color', basic_clr , 'LineWidth', lw);
       
       % Whe we have POLE MARKING mode, we then draw the cell polarity as well
       if APP_opt.t5_display_Pole == 1
           display_PoleMarking( oDet.cellList.meshData{ff}{kk}.mesh, lw, basic_clr);
       end      
        
       % Calculate R2 between each XY_t and all XY_o points of kk-th cell 
       % ontour outline and store the minumum value
       dist(kk) = min(double( sqrt( abs(X_t - X_o).^2 + abs(Y_t - Y_o).^2 ) ));
       
       % If option to display outline of selected track is true, ...
       if APP_opt.t5_display_Outline == 2                       % check if XY_t is inside a cell contour
           if inpolygon(X_t,Y_t, X_o,Y_o)                       % and lot cell outline in selected_clr if XY_t is inside it             
               plot( X_o , Y_o , '-', 'Color', selected_clr , 'LineWidth', lw);
               if APP_opt.t5_display_Pole == 1
                   display_PoleMarking( oDet.cellList.meshData{ff}{kk}.mesh, lw, selected_clr);
               end
               isinc = 1 ;
           else
               isinc = 0 ; 
           end
       end        
    end    
  end %/for
  
  % If XY_t was not inside any cell outline, find the closest cell outline 
  idx_closest =  find(dist == min(dist));     
  if APP_opt.t5_display_Outline == 2  &&  isinc == 0  &&  ~isempty(idx_closest)...  % if not inside cell ...      
          && ~isnan(X_t)
      X_o = [oDet.cellList.meshData{ff}{idx_closest}.mesh(:,1) ; flipud(  oDet.cellList.meshData{ff}{idx_closest}.mesh(:,3)) ];
      Y_o = [oDet.cellList.meshData{ff}{idx_closest}.mesh(:,2) ; flipud(  oDet.cellList.meshData{ff}{idx_closest}.mesh(:,4)) ];
      plot( X_o , Y_o , '-', 'Color', selected_clr , 'LineWidth', lw);
      if APP_opt.t5_display_Pole == 1
        display_PoleMarking( oDet.cellList.meshData{ff}{idx_closest}.mesh, lw, selected_clr);
      end
  end
  
end %/if display_Outline ~= 0 

end



function display_PoleMarking( cMesh, LineW, Colour)
    % Show the "Marked" first pole of the cell outline
    X1 = cMesh(1,1);  
    Y1 = cMesh(1,2);  
    X2 = cMesh(end,1);  
    Y2 = cMesh(end,2);  
    d_shift = - LineW*10 ;                               % shift the point of lw pixels
    m = (Y2-Y1)/(X2-X1);                                 % m = slope of the pole-to-pole line
    t = d_shift / sqrt((X2-X1)^2  + (Y2-Y1)^2);          % ratio between the distance shift and pole-to-pole distance
    X_new = (1-t)*X1 + t*X2 ;                            % find the location of the displaced point from pole 1
    Y_new = (1-t)*Y1 + t*Y2 ;
    % plot the pole indicator line
    plot( [X1, X_new], [Y1, Y_new] , '.-', 'Color', Colour , 'LineWidth', LineW*2, 'MarkerSize', 10);
end




