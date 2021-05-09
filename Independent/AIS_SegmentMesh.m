function [cData] = AIS_SegmentMesh(cData , S_step)
%SegmentMesh = Segment the cell via segmentation lines orthogonal to the
%	cell body axis with spacing between them of S_step [in pixels].
%
% INPUTS ------------------------------------------------------------------
% cData = Stores the data concerning a specific cell cc at frame ff
%         (normally is provided as cellList.meshData{ff}{cc})
%
% S_step = the segments length in [pixels]. to establish the segment length
%          the algorith calculate the length of the axis and divide it by
%          S_step to define equally spaced segments
%
% OUTPUT ------------------------------------------------------------------
% cData = updated cData is returned updated:
% ---.FrameSegmentLines : Store segementation lines that span the entire  
%          frame of image cData.Fluor_Chan(1).IC 
% ---.BodySegmentLines : Store segementation lines that cut the cell body
%
%     (columns order for segmentation lines: [ x1, y1, x2 ,y2])
%
% ---.Mask.Segment : stores all the segmentation masks as a stack of images,
%          where z-dimention represent the ordered sequence of segments.
%          The first segment is POLE 1, while the last is the POLE 2.
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

global APP_opt ;                        % Variable storing WHISIT options


%%%--->>> insert here CHECK PLOT 1 <<<---%%% 

% --- STEP 1 --- Create equally spaced cell body axis curve -----------------
% Store axes coordinates in new variable and smoothen the curve (using
% default moving average filter  = 5)
Ax = [smooth(cData.geom.R_axis(:,1)) , smooth(cData.geom.R_axis(:,2))];

% Generate a new axis with a denser number of points: equal to the original
% number of point times the multiplier factor (m_Fact).
m_Fact = 5 ;       
DensAx = evenspaced( Ax, size(Ax,1).*m_Fact);
DensAx = [smooth(DensAx(:,1)) , smooth(DensAx(:,2))];     % Smoothen the curve (default moving average filter = 5)

% Calculate the piece wise length distances along the axis DensAx
d = sqrt( (DensAx(1:end-1,1) -DensAx(2:end,1)).^2 + (DensAx(1:end-1,2) -DensAx(2:end,2)).^2 );

% Here we create a new axis curve with equally spaced N_seg points. Then,
% find closest DensAx point for each equally spaced axes point (EqSpaceAx)
% and store its index
N_seg = round(sum(d)/S_step) +1;    
EqSpaceAx = evenspaced( DensAx , N_seg );    
for ii = 2 : size(EqSpaceAx,1)-1
    d = sqrt( (EqSpaceAx(ii,1) -DensAx(:,1)).^2 + (EqSpaceAx(ii,2) -DensAx(:,2)).^2 );
    t_idx = find(d == min(d));
    idx(ii-1) = t_idx(1);
end
% Add the first and last point (the poles) as well
idx = unique( [1, idx, length(d)] );



% --- STEP 2 --- Define orthogonal segmentation lines -----------------------
% Using DensAx points, we create orthogonal lines passing throught the
% points defined by EqSpaceAx. These segmentation lines span the entire
% image .Fluor_Chan(1).IC

% Find max X and Y coordinates possible for plotting within the IC figure
[yR, xC] = size(cData.Fluor_Chan(1).IC);
% Temporary variables to store x and y points defining segmentation lines
xx = [] ; 
yy = [] ; 

% Skip the poles (see WHY after end of loop)
for jj = 2 : (length(idx)-1)

    % Tak correct indeces that do not get beyond size limits
    low_idx = idx(jj)-1 ;
    high_idx = idx(jj)+1 ;
    if low_idx < 2
       low_idx = 1;
    end
    if high_idx > length(DensAx)
       high_idx = length(DensAx);
    end
    % Define the two points and gather respective xs and ys in variables
    pA = [DensAx( low_idx ,1) , DensAx( low_idx ,2)];
    pB = [DensAx( high_idx ,1), DensAx( high_idx ,2)];
    xs = [pA(1), pB(1)];
    ys = [pA(2), pB(2)];

    % Calculate the normal line between the two points and store the new
    % two points that define the normal line
    normal = [mean(xs),mean(ys)] + null(pA-pB)';
    nxs = [mean(xs), normal(1)];
    nys = [mean(ys), normal(2)];

    % Calculate Parameter Vector for the line function y = mx + b that pass
    % through the two normal XY points
    coeff = [[1;1]  nxs(:)] \ nys(:);                    
    % coeff(2) == slope, m
    % coeff(1) == intercept, b

    % Now find the segmentation line with same slope and normal line and
    % passing through the point DensAx(idx)
    zz = size(xx,1)+1;
    [xx(zz,:) ,yy(zz,:)] = line_2Points( [1,xC], [1,yR], 0, coeff(1), coeff(2));

end % jj

% In previous for-loop we skipped the poles: it usually result in weirdly
% oriented lines because we only use 2 points to define the normal.
% Moreover, axes' ending the  can have strange "twist" at the very last
% point.
% Instead, for each pole we define a segmentation lines parallel line to
% the previous one and passing through the pole point (last or first point
% in the axis)

% Pole 1 segmentation line
slope = (yy(1,2)-yy(1,1))/(xx(1,2)-xx(1,1)) ;       
[xP1 ,yP1] = line_2Points([ 1,xC], [1,yR], DensAx(1,1), DensAx(1,2), slope);

% Pol 2 segmentation line
slope = (yy(end,2)-yy(end,1))/(xx(end,2)-xx(end,1)) ; 
[xP2 ,yP2] = line_2Points( [1,xC], [1,yR], DensAx(end,1), DensAx(end,2), slope);    

% Update the segementation lines variables
xx = double([ xP1 ; xx ; xP2 ]);
yy = double([ yP1 ; yy ; yP2 ]);            
% Store segementation lines in one matrix, with columns order: [ x1, y1, x2 ,y2]
cData.FrameSegmentLines = [ xx(:,1), yy(:,1), xx(:,2), yy(:,2) ];



% --- STEP 3 --- Define Cell body segmentation lines ------------------------
% Find segmentation line of the cell body,: simply find intercepts between
% the segmentation linesand the cell meshes (excluding Poles)
x_curve = [cData.R_mesh(:,1) ; flipud(cData.R_mesh(:,3))];
y_curve = [cData.R_mesh(:,2) ; flipud(cData.R_mesh(:,4))] ;
xxB = [];
yyB = [];
for hh = 2 : (size(xx,1)-1) 
    % Find all intersect points between two curves
    [ix,iy] = curveintersect(xx(hh,:),yy(hh,:) , x_curve,y_curve) ;

    if length(ix) >= 3
        % In rare cases, weirdly shaped cell's meshes can yield >= 3
        % intercept points. We must have only two points: so take the
        % couple that have longest segmentation line
        for mm = 1 : length(ix)
            dd(:,mm) = sqrt( (ix(mm) -ix).^2 + (iy(mm) - iy).^2 );
        end
        % Get lower triangular matrix of dd to avoid duplicated results
        dd = tril(dd);
        [p1, p2] = find(dd == max(max(dd)));
        ix = [ix(p1); ix(p2)] ;
        iy = [iy(p1); iy(p2)] ;
    end
    xxB = [xxB; ix'];
    yyB = [yyB; iy']; 

end
% Update the segementation lines variables
xxB = double([ xP1 ; xxB ; xP2 ]);
yyB = double([ yP1 ; yyB ; yP2 ]);            
% Store cell body segementation lines, with columns order: [ x1, y1, x2 ,y2]
cData.BodySegmentLines = [ xxB(:,1), yyB(:,1), xxB(:,2), yyB(:,2) ];



% --- STEP 4 --- Create segmentation Masks ----------------------------------
for kk = 1 : (size(xxB,1) -1)
    % Take the points defining the two segmentation lines
    pnts = [ xxB(kk,1)   , yyB(kk,1)     ; ...
               xxB(kk+1,1) , yyB(kk+1,1) ; ...
               xxB(kk+1,2) , yyB(kk+1,2) ; ...
               xxB(kk,2)   , yyB(kk,2)   ] ;
    % Define a bounding polygon around these points and take its coordinates       
    k = boundary(pnts(:,1), pnts(:,2));                
    polygo = pnts(k,:);
    % Define a temporary mask of the area between the two segmentation lines
    tempBW = roipoly( cData.Fluor_Chan(1).IC, polygo(:,1), polygo(:,2) );

    % We store all the segmentation masks as a 3D stack.
    cData.Mask.Segment(:,:,kk) = tempBW ;  

    % Now we can evaluate the average signal for each segment and store
    % it in variable ---.SegmentSig(xx) for all channels available
    cData.Fluor_Chan(1).SegmentSig(kk) = mean( cData.Fluor_Chan(1).IC( cData.Mask.Segment(:,:,kk) ));
    if APP_opt.t1_choose_Chan_2 == 1 
        cData.Fluor_Chan(2).SegmentSig(kk) = mean( cData.Fluor_Chan(2).IC( cData.Mask.Segment(:,:,kk) ));

    end             
    if APP_opt.t1_choose_Chan_3 == 1  &&  APP_opt.t1_CH3_Marker ~= 1      
        cData.Fluor_Chan(3).SegmentSig(kk) = mean( cData.Fluor_Chan(3).IC( cData.Mask.Segment(:,:,kk) ));
    end

    %%%--->>> insert here CHECK PLOT 2 <<<---%%%   

end

%%%--->>> insert here CHECK PLOT 3 <<<---%%%
   
end %main Func




function [qx ,qy] = line_2Points( limX, limY, x3, y3 , slope) 
% line_2Points = given an equation in the form: 
%   y = slope * (qx - x3) + y3
% 	the function find two points laying on the line that are within the 
%   limits limX and limY and passes through the point (x3,y3)
%
% The aim is to find two points that define a segmentation line, points
% whose coordinates must be within the field of view of the image frame.
% In our case limX and limY are the frame size, i.e.: [1, 33] and [1, 56]
% 
% If the line passes through a specific point, x3 and y3 defines it.
% If we simply know the slope and incercept:
%   - y3 act as the "intercept" point on axis y
%   - x3 is simply set to be = 0
  

% Define the ys points to be within the limY and find respective xs points
% that lies on the line at those y-limits (rows of image frame).
ys(1) = limY(1);
ys(2) = limY(2);        
xs(1) = (ys(1) - y3)/slope + x3;                
xs(2) = (ys(2) - y3)/slope + x3;

% Check that xs are within the x-limits (columns of image frame).
% If so they are outside we set xs and find the correct ys.
if xs(1) < limX(1)
    xs(1) = 1;
    ys(1) = slope* (xs(1) - x3)+ y3;
elseif xs(1) > limX(2)         
    xs(1) = limX(2) ;
    ys(1) = slope* (xs(1) - x3)+ y3;
end

if xs(2) < limX(1)
    xs(2) = 1;
    ys(2) = slope* (xs(2) - x3)+ y3;
elseif xs(2) > limX(2)        
    xs(2) = limX(2) ;
    ys(2) = slope* (xs(2) - x3)+ y3;
end   

% Return the xs and ys in separate variables
qx = [xs(1), xs(2)];
qy = [ys(1), ys(2)];
    
end





%% ------- CHECK PLOTs ----------------------------------------------------

% % ----------- CHECK PLOT 1 ----------------------------------------------
% hold off ;
% imshow(cData.Fluor_Chan(1).IC, [min(min(cData.Fluor_Chan(1).IC)), max(max(cData.Fluor_Chan(1).IC))] );
% hold on;        axis equal; 


% % ----------- CHECK PLOT 2 ----------------------------------------------
% imshow(tempBW);   axis equal;
% plot([cData.R_mesh(:,1) ; flipud(cData.R_mesh(:,3))] , ...
%      [cData.R_mesh(:,2) ; flipud(cData.R_mesh(:,4))] , '.' , 'Color',[.9 .8 0]);
% pause(0.1)
% 
% tempBW = cData.Mask.Cell_body + tempBW ;
% tempBW( tempBW == 1 ) = 0;
% 
% imshow(tempBW);
% hold on;        axis equal;
% plot([cData.R_mesh(:,1) ; flipud(cData.R_mesh(:,3))] , ...
%      [cData.R_mesh(:,2) ; flipud(cData.R_mesh(:,4))] , '.' , 'Color',[.9 .8 0]);
% plot(cData.R_mesh(1,1) , cData.R_mesh(1,2) , '*', 'MarkerSize',8, 'Color',[.2 .6 1]); 
% plot(xxB(kk:kk+1,:)' , yyB(kk:kk+1,:)' , '.-r')
% pause(0.1)


% % ----------- CHECK PLOT 3 ----------------------------------------------
% % Plot perimeter
% plot([cData.R_mesh(:,1) ; flipud(cData.R_mesh(:,3))] , ...
%      [cData.R_mesh(:,2) ; flipud(cData.R_mesh(:,4))] , '.' , 'Color',[.9 .8 0]);
% % Plot the cell body axis and points
% plot(DensAx(:,1), DensAx(:,2), '.-', 'Color', [.2 .6 1] );    
% plot(EqSpaceAx(:,1), EqSpaceAx(:,2), '*', 'Color', [.2 1 .5] );  
% % Plot all segmentation lines   
% px = [cData.BodySegmentLines(:,1) , cData.BodySegmentLines(:,3)] ;
% py = [cData.BodySegmentLines(:,2) , cData.BodySegmentLines(:,4)] ; 
% plot(px' , py' , '.-r')
% 
% hold off ;        
% pause(0.1)



