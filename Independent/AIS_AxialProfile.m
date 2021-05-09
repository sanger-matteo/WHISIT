function [cData] = AIS_AxialProfile(cData , Width)
%SegmentMesh = measure the signal profile along cells' body axis
%
% IF Width == 0, then we simply evaluate the signal below the axial line.
%
% IF Width >= 1, to measure an axial profile with a specific Width, we
% need to create two "axis" from the original axis curve. Each curve is
% translated from the original axis of a distance AxDstep ( = Width /2 ),
% in the direction of a vector normal to the original axis. 
% [each new axis moves in one of the two possible directions (+1 and -1)]
%
% For any given couple of point in the two new axis, we can draw a line
% segment and evaluave the average or maximum pixel value below the line.
% The resulting array of value is oriented with Pole 1 as first value, and
% Pole 2 as last. Moreover, the array length is equal to the original axis
% (which in turn is equal to the number of points for one side of the
% cell's meshes)
%
%
% INPUTS ------------------------------------------------------------------
% cData = Stores the data concerning a specific cell cc at frame ff
%         (normally is provided as cellList.meshData{ff}{cc})
%
% width = [pixels]. to establish the segment length
%
%
% OUTPUT ------------------------------------------------------------------
% cData = updated cData is returned updated: ---.AxSig : A single line
%         array storing the pixel values (avg or max) of the signal along 
%         the cell axis for the entire cell length
% 
% Store the two new axis lines used to evaluate the axial profile as:
% ---.geom.R_Naxis_1 = [ x, y ] ;
% ---.geom.R_Naxis_2 = [ x, y ] ;
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


global APP_opt ;                        % Variable storing WHISIT options

% Step distance for translating the axis. Because we create a left and
% right axis, we move translate each only half distance from the
% middle-original axis.
AxDstep = Width / 2 ;

ax = cData.geom.R_axis(:,1);
ay = cData.geom.R_axis(:,2);

if Width == 0
    % The two axis are simply equal to the original
    cData.geom.R_Naxis_1(:,1) = ax ;
    cData.geom.R_Naxis_1(:,2) = ay ;
    cData.geom.R_Naxis_2(:,1) = ax ;
    cData.geom.R_Naxis_2(:,2) = ay ;
    
    % Simply evaluate the profile line below original axis line
    profi = [];
    for jj = 1 : length(nax1)  
        if APP_opt.AIS_ProfileType == 1
           profi(jj) = mean(improfile( cData.Fluor_Chan(1).IC, ax, ay )); 
        elseif APP_opt.AIS_ProfileType == 2
           profi(jj) = max(improfile( cData.Fluor_Chan(1).IC, ax, ay )); 
        end          
    end 
    cData.Fluor_Chan(1).AxSig = profi;        

    if APP_opt.t1_choose_Chan_2 == 1          
        profi = [];
        for jj = 1 : length(nax1)
          if APP_opt.AIS_ProfileType == 1
              profi(jj) = mean(improfile( cData.Fluor_Chan(2).IC, ax, ay )); 
          elseif APP_opt.AIS_ProfileType == 2
              profi(jj) = max(improfile( cData.Fluor_Chan(2).IC, ax, ay )); 
          end          
        end    
        cData.Fluor_Chan(2).AxSig = profi;
    end

    if APP_opt.t1_choose_Chan_3 == 1  &&  APP_opt.t1_CH3_Marker ~= 1          
        profi = [];
        for jj = 1 : length(nax1)
          if APP_opt.AIS_ProfileType == 1
               profi(jj) = mean(improfile( cData.Fluor_Chan(3).IC, ax, ay )); 
          elseif APP_opt.AIS_ProfileType == 2
               profi(jj) = max(improfile( cData.Fluor_Chan(3).IC, ax, ay )); 
          end           
        end    
        cData.Fluor_Chan(3).AxSig = profi;
    end

    
    
elseif Width >= 0
    
    % Find the gradient in X and Y for the entire cell body axis curve
    d_ay = gradient( ax );
    d_ax = gradient( ay );

    % Take the derivative and create normal vector
    d_ax = (ax(2:end)-ax(1:end-1))/2;
    d_ay = (ay(2:end)-ay(1:end-1))/2;
    nvec = [-d_ay,d_ax];
% % % %     % If we get a gradient of zero (no change means we have 2 or more
% % % %     % consecutive points that are exactely the same), then we when
% % % %     % calculating the normal vector we get NaN as results. 
% % % %     % This cause an error and crash of analysis.
% % % %     % In theory in Analysis_AIS.m we ensure to avoid the error by taking 
% % % %     % unique coordinates. This second stap is set as second safeguard.
% % % % %     nvec(nvec == 0) = 0.0001;

    % Normalize the normal vector
    nvecN = bsxfun(@rdivide,nvec,sqrt(sum(nvec.^2,2)));

    % Take the arc-cosine to get the angle in radian of the projection of
    % the normal vector onto the x-axis. From this we can calculate the Dx
    % and Dy movement for the new axis in both directions
    angle = acos(nvecN(:,1));        
    moveDx = AxDstep .* cos(mean(angle));
    moveDy = AxDstep .* sin(mean(angle));
    nax1 = double( ax + moveDx );
    nay1 = double( ay + moveDy );        
    nax2 = double( ax - moveDx );
    nay2 = double( ay - moveDy );  
    % Store the two new "profile" axis; to allow easier plotting
    % (as two variables with columns [x1,y1] and [x2,y2])
    cData.geom.R_Naxis_1(:,1) = nax1 ;
    cData.geom.R_Naxis_1(:,2) = nay1 ;
    cData.geom.R_Naxis_2(:,1) = nax2 ;
    cData.geom.R_Naxis_2(:,2) = nay2 ;

    %%%--->>> insert here CHECK PLOT 1 <<<---%%%

    %%%--->>> insert here CHECK PLOT 2 <<<---%%%

    % Now we can calculate the profile line along all the normal line
    % going from one new axis to the other.
    profi = [];
    for jj = 1 : length(nax1)  
        if APP_opt.AIS_ProfileType == 1
           profi(jj) = mean(improfile( cData.Fluor_Chan(1).IC, ...
                            [nax1(jj),nax2(jj)]', [nay1(jj),nay2(jj)]' )); 
        elseif APP_opt.AIS_ProfileType == 2
           profi(jj) = max(improfile( cData.Fluor_Chan(1).IC, ...
                            [nax1(jj),nax2(jj)]', [nay1(jj),nay2(jj)]' )); 
        end          
    end 
    cData.Fluor_Chan(1).AxSig = profi;        

    if APP_opt.t1_choose_Chan_2 == 1          
        profi = [];
        for jj = 1 : length(nax1)
          if APP_opt.AIS_ProfileType == 1
              profi(jj) = mean(improfile( cData.Fluor_Chan(2).IC, ...
                                [nax1(jj),nax2(jj)]', [nay1(jj),nay2(jj)]' )); 
          elseif APP_opt.AIS_ProfileType == 2
              profi(jj) = max(improfile( cData.Fluor_Chan(2).IC, ...
                                [nax1(jj),nax2(jj)]', [nay1(jj),nay2(jj)]' )); 
          end          
        end    
        cData.Fluor_Chan(2).AxSig = profi;
    end

    if APP_opt.t1_choose_Chan_3 == 1  &&  APP_opt.t1_CH3_Marker ~= 1          
        profi = [];
        for jj = 1 : length(nax1)
          if APP_opt.AIS_ProfileType == 1
               profi(jj) = mean(improfile( cData.Fluor_Chan(3).IC, ...
                                [nax1(jj),nax2(jj)]', [nay1(jj),nay2(jj)]' )); 
          elseif APP_opt.AIS_ProfileType == 2
               profi(jj) = max(improfile( cData.Fluor_Chan(3).IC, ...
                                [nax1(jj),nax2(jj)]', [nay1(jj),nay2(jj)]' )); 
          end           
        end    
        cData.Fluor_Chan(3).AxSig = profi;
    end

end %main Func




% % ----------- CHECK PLOT 1 ------------------------------------------------
% % insert at indicated points of the code to test and debug whether cell's
% % outline and axis are identified properly
%
% axis equal
% hold on
% % Absolute mesh - whole frame
% plot(Xs,Ys,'b')
% plot(cData.geom.axis(:,1),cData.geom.axis(:,2),'r')
% % Relative mesh - cropped frame
% plot(cData.geom.R_axis(:,1),cData.geom.R_axis(:,2))
% plot(R_Xs,R_Ys, 'r') 
% plot(cData.R_model(:,1),cData.R_model(:,2),'b')
% 
% 
% % ----------- CHECK PLOT 2 ------------------------------------------------
% imshow(cData.Fluor_Chan(1).IC, [min(min(cData.Fluor_Chan(1).IC)), max(max(cData.Fluor_Chan(1).IC))] );
% hold on;        axis equal;        
% plot(ax, ay, '.-y');       
% plot(nax1, nay1, '.-r'); 
% plot(nax2, nay2, '.-r'); 
% plot(cData.R_model(:,1),cData.R_model(:,2),'.-c')
% plot( [nax1,nax2]',  [nay1,nay2]', '-g');        
% quiver(x,y,-dy, dx, 'Color', [0 .8 .5])
% quiver(x,y, dy,-dx, 'Color', [0 .5 .8])  
% hold off;  

