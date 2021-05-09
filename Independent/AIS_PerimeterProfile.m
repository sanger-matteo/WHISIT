function [cData] = AIS_PerimeterProfile(cData , in_shift)
%SegmentMesh = The function measures the average signal along the perimeter
% of a cell. 
% 
% From the cells' meshes we create an outer perimeter, whose points are
% evenly spaced. We then create a second perimeter, translating the outer
% points inward of a specific distance (in_shift). 
% Each matching outer-inner couple of points is used to draw a line segment
% and evaluate the average pixel values below such line. This is done for
% every point: as result we have an array that represent an open, linearize
% perimeter, with the first element being Pole 1. Pole 2 is roughly the 
% middle point of such line perimeter.
%
% IF in_shift is set to 0, the algorothm simply measures the pixel values
% below the cell mesh perimeter (outer).
%
% INPUTS ------------------------------------------------------------------
% cData = Stores the data concerning a specific cell cc at frame ff
%         (normally is provided as cellList.meshData{ff}{cc})
%
% in_shift = [pixel] width of the space between outer and inner perimeter.
%         Concentric profile line, going from outer to inner, will then be
%         drawn to measure the pixel value in this space at regual
%         intervals.
%
% density_pts = will determine the number of equally spaced points (Npts)
%         on the outer perimeters. The aim is to have the sampling profile
%         line evenly spaced to have correct and not overlapping sampling.
% 
%
% OUTPUT ------------------------------------------------------------------
% cData = updated cData is returned updated:
%   ---.Fluor_Chan(X).PerimSig = the array contains average signal values
%           of the "linear" perimeter at each sampling point 
%
%   ---.R_OutPeri = the relative coordinates of the outer perimeter
%   ---.R_InnPeri = the relative coordinates of the inner perimeter
%
% [The first pole is Pole 1 (or marked). Pole 2 is roughly located in the
%  midddle of the array (NOTE: because we create evenly spaced perimeters,
%  pole positions are not a perfect match with the ones in the meshes) ]
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


global APP_opt ;                        % Variable storing WHISIT options


% --- STEP 1 --- Create equally spaced cell body axis curve -----------------
% Store axes coordinates in new variable and smoothen the curve (using
% default moving average filter  = 5)
outXY = [ [cData.R_mesh(:,1) ; flipud(cData.R_mesh(:,3))] , ...
          [cData.R_mesh(:,2) ; flipud(cData.R_mesh(:,4))] ];
d = sqrt( (outXY(1:end-1,1)-outXY(2:end,1)).^2 + (outXY(1:end-1,2)-outXY(2:end,2)).^2 );

% Create a new cell mesh, where all points are evenly spaced , which is
% never the case for the original ones, since there is a left and right
% side. After calculating the cell perimeter length (in pixel) set what
% fraction [0.25-1] of points to create for a new evenly spaced perimeter.
% [70-75% is best to roughly have 1 point avery 1-1.25 pixels].
density_pts = APP_opt.AIS_PerimSpacing;     
Npts = round(sum(d) *density_pts);
outXY = double( evenspaced( outXY, Npts ) );


if in_shift == 0
    % --- 0 --- No "thickness" of the membrane --------------------------------------
    % There is no inner perimeter, we simply make it equal to outer 
    cData.R_OutPeri = outXY ; 
    cData.R_InnPeri = outXY ;
    
    % We simply measure the signal below the outer perimeter line
    profi = [];
    for kk = 1 : length(outXY)  
        profi(kk) = improfile( cData.Fluor_Chan(1).IC, outXY(kk,1), outXY(kk,2) );   
    end
    cData.Fluor_Chan(1).PerimSig = profi;        

    if APP_opt.t1_choose_Chan_2 == 1          
        profi = [];
        for kk = 1 : length(outXY)  
            profi(kk) = improfile( cData.Fluor_Chan(2).IC, outXY(kk,1), outXY(kk,2) );   
        end
        cData.Fluor_Chan(2).PerimSig = profi;  
    end

    if APP_opt.t1_choose_Chan_3 == 1  &&  APP_opt.t1_CH3_Marker ~= 1          
        profi = [];
        for kk = 1 : length(outXY)  
            profi(kk) = improfile( cData.Fluor_Chan(3).IC, outXY(kk,1), outXY(kk,2) );   
        end
        cData.Fluor_Chan(3).PerimSig = profi;  
    end

    
    
elseif in_shift >= 1    
    % --- 1 --- Define perimeter sampling lines --------------------------------
    % Define the inner perimeter, tranlating each outer point inward, of a
    % distance = inshift, andin the normal direction, in respect to the
    % "tangent" line at each specific point. 
    % Then we use each inner-outer points couple to draw sampling lines to
    % assess the average pixel value at each point of the perimeter

    % Temporary variables to store x and y points defining segmentation lines
    xx = [] ; 
    yy = [] ;
    inXY_A = [];
    inXY_B = [];
    
    for ii = 1 : size(outXY,1)   
        % Take correct indexes 
        if ii == 1
            low_idx = (size(outXY,1)-1) ;
        else
            low_idx = ii -1;
        end
        if ii == size(outXY,1)
            high_idx = 1;
        else
            high_idx = ii+1 ;
        end

        % Define the two points before and after the ii-th point of interest,
        % From pA_pB segment we extract the normal vector. 
        pA = [outXY( low_idx ,1) , outXY( low_idx ,2) ];
        pB = [outXY( high_idx ,1), outXY( high_idx ,2)];
        xs = [pA(1), pB(1)];
        ys = [pA(2), pB(2)];

        % Calculate the second point that, together with the point middle
        % to pA_pB segment, will form the vector normal.
        normal = [mean(xs),mean(ys)] - null(pA-pB)';
        nxs = [mean(xs), normal(1)];
        nys = [mean(ys), normal(2)];   

        % Find the angle between x-axis and normal vector: acos(adjecent/hypothenus)
        alp = acos( (nxs(2)-nxs(1)) / sqrt((nxs(2)-nxs(1))^2+(nys(2)-nys(1))^2) );
        % ---> CHECK angle alp ----- 
        % Calculate the shift in the normal direction to place new point
        dx = cos(alp) .* in_shift;
        dy = sin(alp) .* in_shift;    

        % Cols [1,2] of meshes are the "left" side of the cell, while [3,4]
        % are the "right" side. However, sometimes they are inverted.
        % (some inconsistency in Oufti, and perhaps when we flip poles).
        % Then, though the cross product does not change (AxB), one case
        % calculates the "small/inner" angle between two vectors and the
        % other the "bigger/outer" angle. Therefore, the direction of
        % normal vector sometime inward, but occasionally outward.
        % To overcome this problem we can use two tricks:
        % - when x of the A-point < than B-point, change direction of translation
        %   (This ensure we have a "consistent" direction, translating all
        %   points in the same normal direction)
        % - we calculate both the inward and outward translation and then
        %   find and choose the one with the smallest perimeter as the
        %   inner one.
        
        if pA(1) < pB(1)
            inXY_A(ii, 1) = outXY(ii ,1) + dx ; 
            inXY_A(ii, 2) = outXY(ii ,2) - dy ;
            inXY_B(ii, 1) = outXY(ii ,1) - dx ;
            inXY_B(ii, 2) = outXY(ii ,2) + dy ;        
        else
            inXY_A(ii, 1) = outXY(ii ,1) - dx ;
            inXY_A(ii, 2) = outXY(ii ,2) + dy ;
            inXY_B(ii, 1) = outXY(ii ,1) + dx ;
            inXY_B(ii, 2) = outXY(ii ,2) - dy ;        
        end

    end % ii

    % ----- OPTIONAL --------
    % Create a new evenly spaced inner perimeter (a simple translation as
    % we did in for-loop above can cluster some points or misposition them
    % a bit when they are too densily packed).
    % However, in this way the profile lines will not be "concentric", in
    % other words pointing to the cell center/axis; instead they will be
    % slightly tilted.
    % inXY = double( evenspaced(inXY, Npts) );
    % ------------------------

    dA = sum( sqrt( (inXY_A(1:end-1,1) -inXY_A(2:end,1)).^2 + (inXY_A(1:end-1,2) -inXY_A(2:end,2)).^2 ));
    dB = sum( sqrt( (inXY_B(1:end-1,1) -inXY_B(2:end,1)).^2 + (inXY_B(1:end-1,2) -inXY_B(2:end,2)).^2 ));
    if dA < dB
        inXY = inXY_A ;
    elseif dB < dA
        inXY = inXY_B ;
    end 
    
    % Store the outer and inner perimeters used for analysis. (see function
    % introductory comment for clarifications)
    cData.R_OutPeri = outXY ; 
    cData.R_InnPeri = inXY ;

    
    % Now we can calculate the profile lines along all the perimeter, going
    % from outer (original) to the new inner perimeter
    profi = [];
    for jj = 1 : length(inXY)  
        profi(jj) = mean(improfile( cData.Fluor_Chan(1).IC, ...
                    [inXY(jj,1), outXY(jj,1)]', [inXY(jj,2), outXY(jj,2)]' ));          
    end 
    cData.Fluor_Chan(1).PerimSig = profi;        

    if APP_opt.t1_choose_Chan_2 == 1          
        profi = [];
        for jj = 1 : length(inXY)  
            profi(jj) = mean(improfile( cData.Fluor_Chan(2).IC, ...
                        [inXY(jj,1), outXY(jj,1)]', [inXY(jj,2), outXY(jj,2)]' ));          
        end 
        cData.Fluor_Chan(2).PerimSig = profi;  
    end

    if APP_opt.t1_choose_Chan_3 == 1  &&  APP_opt.t1_CH3_Marker ~= 1          
        profi = [];
        for jj = 1 : length(inXY)  
            profi(jj) = mean(improfile( cData.Fluor_Chan(3).IC, ...
                        [inXY(jj,1), outXY(jj,1)]', [inXY(jj,2), outXY(jj,2)]' ));          
        end 
        cData.Fluor_Chan(3).PerimSig = profi; 
    end

    
end %in_shift


end %main Func




% ---> CHECK angle alp ------------------------------------------------------
% NN = null(pA-pB)';
% disp([rad2deg(alp), acosd(-NN(1)), asind(NN(2))])




% % ----------- CHECK PLOT 1 ------------------------------------------------
% hold off;
% imshow(cData.Fluor_Chan.IC, [min(min(cData.Fluor_Chan(1).IC)), max(max(cData.Fluor_Chan(1).IC))] );
% hold on;        axis equal; 
%  
% % % % % in the for loop
% % % %  plot(nxs, nys, '.-' , 'Color',[.0 .5 1], 'LineWidth', 1.5);
%          
% plot(inXY(:,1), inXY(:,2), '.-' , 'Color',[.9 .0 .0]);
% plot(inXY(1,1), inXY(1,2), '*' , 'Color',[.9 .0 .0]);
% 
% plot(outXY(:,1), outXY(:,2) , '.' , 'Color',[.9 .8 .0]);
% plot(outXY(1,1), outXY(1,2) , '*' , 'Color',[.9 .8 .0]);
%
% plot( [inXY(:,1), outXY(:,1)]', [inXY(:,2), outXY(:,2)]' , '-' , 'Color',[.0 .8 .0]);
