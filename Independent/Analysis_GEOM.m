function [cData, R_Xs, R_Ys]  =  Analysis_GEOM(cData, extra_border, Img_BF, Img_CH1, Img_CH2, Img_CH3)
% Analysis_GEOM = Evaluate the geometric parameters (cell meshes, contour
%   and cell body axis), create and store cropped images of the specific
%   cell in the different channels
%
% cData.meshes = main variable generated by cell detection in Oufto.
%    This variable carries the cell's contour coordinates and is the most
%    important to define the cell body area. 
%    Pole 1 is the first row of meshes. Pole 2 is the last row of meshes.
%    (NOTE: those Pole points are repeated for the left and right side, to
%    close the cell perimeter).
%
% New sets of coordinates and mesh are calculated to be relative to the 
% cropped image(s) created, for example:
% .R_mesh
% .R_model
% .geom.axis
% .geom.R_axis
% .geom.length
% .geom.area
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


global APP_opt ;                        % Variable storing WHISIT options

% Ensure to have no duplicates and only unique points are in meshes:
% having duplicates can cause errors and crash at later stages of analysis.
cData.mesh = unique(cData.mesh,'rows');

% In rare exception, .mesh given by OUFTI have two common error that can
% cause error at some stages of analysis (especially Segmentation and
% Perimeter)
% *** 1 *** Meshes' points are not always "ordered" progressively from pole
% to pole, althought the meshes describe a correct perimeter.
% First and Last rows are the same point for left (1:2) and right (3:4)
% side of cell's meshes. If it is not so, there is a possible misordering.
check = cData.mesh(:,1:2) == cData.mesh(:,3:4);
check = unique(check','rows')' ;
% To fix this issue, we reorder all the points to go from pole 1 to 2.
cData.mesh = checkMeshesOrder( cData.mesh , check );

% *** 2 *** .mesh define the cell perimeter and is organized in a left and 
% right side of the cell, as follow: [leftX, leftY, rightX, rightY]
% (left/right refere to direction Pole 1 ---> Pole2)
% However, sometimes they are inverted with right first and then left (due
% to some undetected error in Oufti probably). So we must exchange columns
% of mesh to ensure the correct order of columns
cData.mesh = checkMeshesSide( cData.mesh );


% Gather cell meshes as XY coordinates with handier variable name
Xs = [cData.mesh(:,1) ; flipud(cData.mesh(:,3))] ;
Ys = [cData.mesh(:,2) ; flipud(cData.mesh(:,4))] ;  

% Store in cData a cropped image(s) of the cell in the Bright Field ( + extra-border) 
if APP_opt.t1_choose_BrightField == 1 
    cData.Bright_Field.IC  = imcrop( Img_BF  , [(min(Xs)-extra_border), (min(Ys)-extra_border), ((max(Xs)-min(Xs))+extra_border*2), ((max(Ys)-min(Ys))+extra_border*2) ]);
end

% Store in cData a cropped fluorescent image(s) with additional extra-border 
cData.Fluor_Chan(1).IC = imcrop( Img_CH1 , [(min(Xs)-extra_border), (min(Ys)-extra_border), ((max(Xs)-min(Xs))+extra_border*2), ((max(Ys)-min(Ys))+extra_border*2) ]);
% if there is a Channel 2 selected
if APP_opt.t1_choose_Chan_2 == 1          
    cData.Fluor_Chan(2).IC = imcrop( Img_CH2 , [(min(Xs)-extra_border), (min(Ys)-extra_border), ((max(Xs)-min(Xs))+extra_border*2), ((max(Ys)-min(Ys))+extra_border*2) ]);
end
% if there is a Channel 3 selected
if APP_opt.t1_choose_Chan_3 == 1        
    cData.Fluor_Chan(3).IC = imcrop( Img_CH3 , [(min(Xs)-extra_border), (min(Ys)-extra_border), ((max(Xs)-min(Xs))+extra_border*2), ((max(Ys)-min(Ys))+extra_border*2) ]);
end

sf = 1;      % <<<----- correction factor for the shift of coordinates: origin at [0;0] or at [1;1]
% Meshes and Models Relative to cropped image of the cell (.CHX.IC)
cData.R_mesh(:,1) = double( cData.mesh(:,1) - (min(Xs)-sf - extra_border) );
cData.R_mesh(:,2) = double( cData.mesh(:,2) - (min(Ys)-sf - extra_border) );
cData.R_mesh(:,3) = double( cData.mesh(:,3) - (min(Xs)-sf - extra_border) );
cData.R_mesh(:,4) = double( cData.mesh(:,4) - (min(Ys)-sf - extra_border) );

cData.R_model(:,1) = double( cData.model(:,1) - (min(Xs)-sf - extra_border) );
cData.R_model(:,2) = double( cData.model(:,2) - (min(Ys)-sf - extra_border) );
    
% For simplicity in coding we give shorter names to coordinates, all,
% left and right side of the cell outline coordinates
R_Xs = [cData.R_mesh(:,1) ; flipud(cData.R_mesh(:,3))] ;
R_Ys = [cData.R_mesh(:,2) ; flipud(cData.R_mesh(:,4))] ;
xr = double( cData.mesh(:,1) ); 
xl = double( cData.mesh(:,3) ); 
yr = double( cData.mesh(:,2) ); 
yl = double( cData.mesh(:,4) ); 
% Create the axis coordinates as the midpoint between the two meshes'
% sides. NOTE: if one pole was marked, the cData.mesh have already been
% reoriented correctly. Therefore, the axis will have same orientation
% as the meshes (first axes point == to first meshes' point)
mx = [xr(1) ; (xr(2:end-1) + xl(2:end-1)) /2 ; xr(end)] ;
my = [yr(1) ; (yr(2:end-1) + yl(2:end-1)) /2 ; yr(end)]  ;

% Cell axis coordinates, following middle cell line: 
% first colums are all xs, and second are all ys
cData.geom.axis = [mx , my] ;
cData.geom.length = sum( sqrt( (mx(2:end)-mx(1:end-1)).^2 + (my(2:end)-my(1:end-1)).^2 )) ;            
cData.geom.area = polyarea(R_Xs,R_Ys) ;

cData.geom.R_axis(:,1) = cData.geom.axis(:,1) - (min(Xs)-sf - extra_border) ;
cData.geom.R_axis(:,2) = cData.geom.axis(:,2) - (min(Ys)-sf - extra_border) ;


%     % --- ELLIPSE FITTING
%     xs = [cData.R_mesh(:,1) , cData.R_mesh(:,3)] ;
%     ys = [cData.R_mesh(:,2) , cData.R_mesh(:,4)] ;
%     % Ellipse fitting of the cell mesh outline (Relative coordinates)
%     [A, B, X0, Y0, Phi] = Geom_My_fit_ellipse__v2(xs,ys);   
%     cData.geom.A_Axis = A ;
%     cData.geom.B_Axis = B ;
%     cData.geom.X0 = X0 ;
%     cData.geom.Y0 =  Y0 ;
%     cData.geom.phi = Phi ;
%     cData.geom.Epsilon = (1 - B/A); 

   
end %main Func





%% --- Script-related functions --------------------------------------------
function reorderedMesh = checkMeshesOrder( oriMesh , check )   
% Reorganize the points in meshes to go progressively from pole to pole,
%
% We can find the poles because they are row in meshes where left (1:2)
% and right (3:4) side are exactly the same point. Starting from one
% pole we calculate the distance to all other points, find and store
% the index of the closest. Then we remove this points from the xy
% list, to ensure that in next iterations we do not go backward in our
% search. We repeate the cicle with the latest point and so on until we
% have no points left

% Initialize variables to store indexes and meshes coordinates 
idxA = [];
idxB = [];
xyA = [oriMesh(:,1), oriMesh(:,2), [1:length(oriMesh(:,2))]'] ; 
xyB = [oriMesh(:,3), oriMesh(:,4), [1:length(oriMesh(:,4))]'] ; 

% %%% ---> Check Plot 1  <----------------------- %%%

% !!! The Reordering we always start from Pole 1 !!!
% In this way we ensure to keep same orientation of the meshes as in the
% oriXY, and avoid flipping the poles.
if check(1) == 1                            % Take and start from Pole 1                     
    idxA(1) = 1 ;                           idxB(1) = 1 ;
    last_A = xyA(idxA(1),:);                last_B = xyB(idxB(1),:);
    xyA(idxA(1),:) = [] ;                   xyB(idxB(1),:) = [] ;  
    
% In any other case, find the nearest to the Pole 1 (the one with smaller index)
% Assign it a sPole 1 and start ordering from it
else   
    tempidx = find(check(1:end)==1);    
    idxA(1) = min(tempidx) ;                idxB(1) = min(tempidx) ;
    last_A = xyA(idxA(1),:);                last_B = xyB(idxB(1),:); 
    xyA(idxA(1),:) = [] ;                   xyB(idxB(1),:) = [] ; 
end


% Iteratively find the next nearest point, update the list of indexes
% and remove the point from the current list of remaining points.
while length(idxA) < size(oriMesh(:,1),1)
    dA = sqrt( (last_A(1) -xyA(1:end,1)).^2 + (last_A(2) -xyA(1:end,2)).^2 );
    [~,ii] = sort(dA);
    idxA(end+1) = xyA(ii(1), 3) ;
    last_A = xyA(ii(1), :); 
    xyA(ii(1),:) = [] ;    

    dB = sqrt( (last_B(1) -xyB(1:end,1)).^2 + (last_B(2) -xyB(1:end,2)).^2 );
    [~,ii] = sort(dB);
    idxB(end+1) = xyB(ii(1), 3) ;
    last_B = xyB(ii(1), :); 
    xyB(ii(1),:) = [] ; 
end

reorderedMesh = oriMesh(idxA,:);

% %%% ---> Check Plot 2  <----------------------- %%%

end %func checkMeshesOrder




function updatedMesh = checkMeshesSide( oriMesh )
% The function checkMeshesSide evaluate that the meshes are in the
% correct order, which is: [leftX, leftY, rightX, rightY]. 
% It evaluate the cross product between Pole1-to-Pole2 vector and
% Pole1-to-Side1 point vector to gain the signed angle between them:
% --- If angle theta > 0 then the first two column of mesh is correctly the 
%     "left" side of the cell.
% --- If angle theta < 0 then the first two column of mesh is the "right"
%     side of the cell. We then invert the order to restore the standard one.

% 3D vector of the Poles points
P1 = [ oriMesh( 1, 1) , oriMesh( 1, 2) , 0 ];
P2 = [ oriMesh(end, 1), oriMesh(end, 2), 0 ];
% 3D vector of the next sides points
side_1 = [oriMesh( 4, 1) , oriMesh( 4, 2) , 0];
side_2 = [oriMesh( 4, 3) , oriMesh( 4, 4) , 0];

% Pole vector (normalized to P1) and norm
uP = (P2 -P1);
normU = norm(uP);
% Side 1 vector (normalized to P1) and its norm
vS1 = (side_1 -P1);
normV = norm(vS1);  
% Calculate the signed angle: uP x vSide1
dotUV = dot(uP, vS1);
c = cross( uP , vS1);
theta_1 = sign(c(3)) *atan2d(norm(cross(uP,vS1)),dot(uP,vS1));

if theta_1 < 0          
    % %%% ---> Check Plot 3  <----------------------- %%%  
    % Side 1 is at the right of pole 1, invert  columns order to maintain
    % it as intended: [leftX, leftY, rightX, rightY]
    temp = [oriMesh(:,3), oriMesh(:,4), oriMesh(:,1), oriMesh(:,2)] ;
    updatedMesh = temp ;    
    % %%% ---> Check Plot 4  <----------------------- %%%
    
else
    % Mesh is already in correct order, return them unchanged.
    updatedMesh = oriMesh ;    
    
end

end % func checkMeshesSide




% % %%% ---> Check Plot 1 <---------------------- %%%
% figure(1)
% clf(1)
% hold on
% axis equal
% plot(xyA(:,1), xyA(:,2), '.-c')
% plot(xyB(:,1), xyB(:,2), '.-m')

% % %%% ---> Check Plot 2 <----------------------- %%%
% AA = [oriMesh(idxA,1), oriMesh(idxA,2)];
% BB = [oriMesh(idxA,3), oriMesh(idxA,4)];
% plot(AA(:,1), AA(:,2), '.-c');
% plot(BB(:,1), BB(:,2), '.-m');
% pause(1)
% for jj = 1 : length(AA)
%     plot(AA(jj,1), AA(jj,2), 'ob')
%     plot(BB(jj,1), BB(jj,2), 'or')
%     pause(0.1)
% end
% pause(0.1);

% %%% ---> Check Plot 3  <----------------------- %%%
% figure(1)
% clf;
% hold on;        
% axis equal; 
% plot(oriMesh(:,1), oriMesh(:,2) , 'm.' );
% plot(oriMesh(:,3), oriMesh(:,4) , 'c.' );
% plot(oriMesh(1,1), oriMesh(1,2) , 'o' , 'Color', [.0 .0 .0]);
% plot([oriMesh(1,1), oriMesh(end,1)], [oriMesh(1,2) , oriMesh(end,2)] , '.-k');
% plot([oriMesh(1,1), oriMesh( 4, 1)], [oriMesh(1,2), oriMesh( 4, 2)] , '.-r' );
% plot([oriMesh(1,1), oriMesh( 4, 3)], [oriMesh(1,2), oriMesh( 4, 4)] , '.-g' );

% %%% ---> Check Plot 4  <----------------------- %%%
% figure(2)
% clf;
% hold on;        
% axis equal; 
% plot(updatedMesh(:,1), updatedMesh(:,2) , 'm.' );
% plot(updatedMesh(:,3), updatedMesh(:,4) , 'c.' );
% plot(updatedMesh(1,1), updatedMesh(1,2) , 'o' , 'Color', [.0 .0 .0]);
% plot([updatedMesh(1,1), updatedMesh(end,1)], [updatedMesh(1,2) , updatedMesh(end,2)] , '.-k');
% plot([updatedMesh(1,1), updatedMesh( 4, 1)], [updatedMesh(1,2), updatedMesh( 4, 2)] , '.-r' );
% plot([updatedMesh(1,1), updatedMesh( 4, 3)], [updatedMesh(1,2), updatedMesh( 4, 4)] , '.-g' );
% pause(0.5)
