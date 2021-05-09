function cData = InF_Fluor_Analysis_PL_v2(cData, ff, cc, I_CH1, I_CH2)
%Analysis_PL = analyse the signal inside cells using algorithm PL
% (Polar Signal). The cell area is analysed to identify two compartments: 
% the two, distinct, poles of a cell and the remaining cytosoic area 
%
%
% INPUTS ------------------------------------------------------------------
% cData = Stores the data concerning a specific cell cc at frame ff
%         (normally in the provided as cellList.meshData{ff}{cc})
% 
%
% OUTPUT ------------------------------------------------------------------
% cData = update cData is return back
% 
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-

% --- Code INFO -----------------------------------------------------------
% .CHX.IC = The function create a small cropped image from each
%           fluorescence image provided with the specific cc-th cell in the
%           center and stores it. Thus, at the end the Results.mat file
%           is a library of "cells" containing the raw images as well,
%           allowing to perform further analysis without the need of the
%           stack of images.
%
% Masks are generated to define specific cell areas. 
% They are logical matrix (0 or 1) with value of 1 if the pixel belong to
% the area of interest. Masks allow to easily analyse fluorescence in
% specific area by, since the (logic) mask can be used as a list of indexes
% to select specific areas in the fluorescence image. For example the
% average signal in a cell is simply  >>  mean(...CHX.IC( ...CellMask ));
%
% .Mask_wCell = this is a mask of the whole cell body . This is defined using
%           outline coordinates as deternined in Oufti. 
%
% .Mask_PL_1 and .Mask_PL_1 = the two masks define the corresponding pole areas.
%           The algorithm first search for a "bright" spot or cluster of
%           signal around the pole area. Once found using criteria given by
%           the user, it creates an area to around it that will define the
%           mask(s)
%
% .Mask_pCyto = This mask is defined by subtracting Mask_wCell - Mask_Memb
%
% (Although specific to Analysis_M2C, the algorithm need to define .Mask_Memb
% because it will be the starting point for defining the masks Mask_PL_X) 
%
% .Mask_Memb = this mask, at its minimum, correspond exactely with the cell
%           perimeter. Therefore is just one pixel width. According to the
%           parameter used by the user, the membrane area can be extended
%           inside the cell to comprise a width of X pixel. 
%
% .Mask_mCyto = This mask is defined by subtracting Mask_wCell - Mask_Memb
%
% The algorithm analyse several important geometric parameter of the cell
% and its outline identified during detection in Oufti. All coorinates and
% mesh are recalculated to be relative to the newly cropped image
% .R_mesh
% .R_model
% .geom.axis
% .geom.length            
% .geom.area
% .geom.R_axis%
%



%% --- Initialize ---------------------------------------------------------
global APP_opt ;                        % Variable storing WHISIT options

extra_border = APP_opt.BorderBox ;      % extra border area to add to a cell cropped image
Memb_sz = APP_opt.Memb_sz ;             % number of circle to repeat to extend the membrane area: roughly each 
                                        % cicle add a pixel width from perimeter inward to the cell to membrane area
                                        
Search_Cr_CH1 = APP_opt.CH1_Search_Circle ;     % Cell pole coordinates to chose for circle fit
sc_fc_CH1 = APP_opt.CH1_ScaleFact ;             % [%], scale foactor for pole circle
PL_Rad_CH1 = APP_opt.CH1_Peak_R ;               % Signal at pole radius

if APP_opt.t1_choose_Chan == 2                      % if there is a Chan 2 
    Search_Cr_CH2 = APP_opt.CH2_Search_Circle ;        % Thres_Cr = 10 ;
    sc_fc_CH2 = APP_opt.CH2_ScaleFact ;                % sc_fc = 0.7 ;   
    PL_Rad_CH2 = APP_opt.CH2_Peak_R ;                  % PL_Rad = 2.5 ;
end



%% -- Body of the Function --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                      

if ~isempty(cData.model)  &&  size(cData.mesh,2) == 4    
    % Save the option parammeters used for analysis
    cData.info.Original_Frame   = ff ;
    cData.info.Original_cellID  = cc ;    
    cData.info.Analysis_Opt.extra_border = APP_opt.BorderBox ;
    cData.info.Analysis_Opt.Memb_sz = APP_opt.Memb_sz ;
            
    % Find Min and Max in X-Y axis for the c-th-cell modelinates
    Xs = [cData.mesh(:,1) ; flipud(cData.mesh(:,3))]  ;
    Ys = [cData.mesh(:,2) ; flipud(cData.mesh(:,4))] ;   
    
    % Crop cell fluorescent signal of with additional (ex)-border
    cData.CH1.IC = imcrop( I_CH1 , [(min(Xs)-extra_border), (min(Ys)-extra_border), ((max(Xs)-min(Xs))+extra_border*2), ((max(Ys)-min(Ys))+extra_border*2) ]);
    % if there is a Chan 2 selected
    if APP_opt.t1_choose_Chan == 2 
        cData.CH2.IC = imcrop( I_CH2 , [(min(Xs)-extra_border), (min(Ys)-extra_border), ((max(Xs)-min(Xs))+extra_border*2), ((max(Ys)-min(Ys))+extra_border*2) ]);
    end    
    sf = 1;      % <<<----- correction factor for the shift between BF and Fluor. due
                 % to difference in modelinates system: origin at [0;0] or at [1;1]
    % modelinates Relative to cropped image of the cell (.CHX.IC)
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
%          
% --- Create the Masks that define cytosol and membrane area --------------
%--------------------------------------------------------------------------
    
    % ---> 1 Cell Mask 
    % Whether I use CH1.IC or CH2.IC does not change anything here: 
    % I just need IC for as blank mask to create BW, and both CH has same sizes
    BW = roipoly( cData.CH1.IC, R_Xs, R_Ys );
%     BW = roipoly( cData.CH1.IC, cData.R_model(:,1), cData.R_model(:,2) );

    % poly2mask = Convert region of interest (ROI) polygon to region mask
%     BW = poly2mask( [cData.R_mesh(:,1); cData.R_mesh(:,3)], [cData.R_mesh(:,2); cData.R_mesh(:,4)] ,...
%                     size(cData.CH1.IC,1), size(cData.CH1.IC,2)) ;

    %---------R_P_model------Relative_Perimeter_model----------------------
    % Using round(.mesh) to define .model (see Retrive_BF_mesh.m), 
    % does not give a continuos line of pixel around the cell. But
    % with the mask we can can extract the perimeter and re-define cell 
    % coordinates so that they describe a continuous line.
    
    tempPx = regionprops(bwperim(BW),'PixelList');
    cData.R_P_model = tempPx.PixelList ;                   
    
    % NOTE: coordinates are ordered in respect to x-axis. Is not a problem
    % for .model, we never need them ordered (compared to .mesh)
    %----------------------------------------------------------------------
        
    % ---> 2 Temporary Masks
    % define contour perimeter
    t_BWm = bwperim(BW);
    for i = 1 : Memb_sz -1
        t_BWc = imfill(t_BWm,'holes') - t_BWm;  % create temporay cytosol mask:
        in_per =  bwperim(t_BWc) ;              % subtract 'holes' to perimeter to have new Cytosol
        t_BWm = t_BWm + in_per ;                % inner hole perimiter ...  
    end                                         % ... add it to Membrane area
    
    % ---> 3 Final Cytosol and Membrane Masks and Save them in my_DB     
    %BWm = t_BWm + in_per;
    BWm = logical( t_BWm );
    BWc = logical( imfill(t_BWm,'holes') - BWm );      
    cData.Mask_wCell = BW;
    cData.Mask_mCyto = BWc;
    cData.Mask_Memb = BWm;

    
% --- FLUORO CHANNEL 1 --- 111111111111111111111111111111111111111111111111
% --- POLE: detection and measurements ------------------------------------
%--------------------------------------------------------------------------
    % We know that foci signal can only be at poles. Thus, we measure first
    % fit poles to circle and define an area where to search for the signal
        
% ---> POLE 1 (OLD) search area
    if size(cData.R_mesh, 1) <= Search_Cr_CH1
        xs = [ cData.R_mesh(1:end, 1); cData.R_mesh(1:end, 3) ];
        ys = [ cData.R_mesh(1:end, 2); cData.R_mesh(1:end, 4) ];
    else
        xs = [ cData.R_mesh(1:Search_Cr_CH1, 1); cData.R_mesh(1:Search_Cr_CH1, 3) ];
        ys = [ cData.R_mesh(1:Search_Cr_CH1, 2); cData.R_mesh(1:Search_Cr_CH1, 4) ];
    end
    [Cnt, R] = fit_circle( [xs,ys] ) ;
    ang = 0:0.1:2*pi;      xp = (R*sc_fc_CH1) *cos(ang);      yp = (R*sc_fc_CH1) *sin(ang);
    x_p = cData.R_mesh(1,1);         D_x = (x_p + Cnt(1)) / 2.08 ;
    y_p = cData.R_mesh(1,2);         D_y = (y_p + Cnt(2)) / 2.08 ;
    % create temp_mask of pole circle search area
    BW_CH1_sc_1  = roipoly( cData.CH1.IC, D_x+xp , D_y+yp ) ;
    
% ---> POLE 2 (NEW) search area
    if size(cData.R_mesh, 1) <= Search_Cr_CH1
        xs = [ cData.R_mesh(end-1:end, 1); cData.R_mesh(end-1:end, 3) ];
        ys = [ cData.R_mesh(end-1:end, 2); cData.R_mesh(end-1:end, 4) ];
    else
        xs = [ cData.R_mesh(end-Search_Cr_CH1:end, 1); cData.R_mesh(end-Search_Cr_CH1:end, 3) ];
        ys = [ cData.R_mesh(end-Search_Cr_CH1:end, 2); cData.R_mesh(end-Search_Cr_CH1:end, 4) ];
    end
    [Cnt, R] = fit_circle( [xs,ys] ) ;
    ang = 0:0.1:2*pi;      xp = (R*sc_fc_CH1) *cos(ang);      yp = (R*sc_fc_CH1) *sin(ang);
    x_p = cData.R_mesh(end,1);       D_x = (x_p + Cnt(1)) / 1.92  ;
    y_p = cData.R_mesh(end,2);       D_y = (y_p + Cnt(2)) / 1.92 ;
    % create temp_mask of pole circle search area
    BW_CH1_sc_2  = roipoly( cData.CH1.IC, D_x+xp , D_y+yp ) ;
    
    % We save those mask, just for plotting
    cData.CH1.M_sc_PL_1 = BW_CH1_sc_1;
    cData.CH1.M_sc_PL_2 = BW_CH1_sc_2;
    
% --->  Search highest value pixel. 
    PL1_CH1 = cData.CH1.IC .* BW_CH1_sc_1 ;                 % Apply mask to isolate raw-value signal at the poles area
    PL2_CH1 = cData.CH1.IC .* BW_CH1_sc_2 ;
    [iy_1, ix_1] = find(PL1_CH1==(max(max(PL1_CH1)))) ;     % Now at each pole find the max(pixel_value)
    [iy_2, ix_2] = find(PL2_CH1==(max(max(PL2_CH1)))) ;
    iy_1 = iy_1(1);    iy_2 = iy_2(1);    
    ix_1 = ix_1(1);    ix_2 = ix_2(1);
    % Then create a smaller area of interest with the higherst pixel at the
    % center. This is now the detected polar signal   
    xp = ix_1+ (PL_Rad_CH1*cos(ang));      yp = iy_1+ (PL_Rad_CH1*sin(ang));
    BW_PL1_CH1 = roipoly( cData.CH1.IC, xp , yp );
    xp = ix_2+ (PL_Rad_CH1*cos(ang));     yp = iy_2+ (PL_Rad_CH1*sin(ang));
    BW_PL2_CH1 = roipoly( cData.CH1.IC, xp , yp );
   
    
    % In case that there are no polar foci and signal is cytosolic, the pole
    % are tend to be much inside the cell. In order to avoid that we need
    % restrict the area of search and look for highest value pixel, but
    % closer to the membrane

    % Create perimeters and test if the pole area cross the cell membrane
    t_BWm = bwperim(BW);
    t_PL1 = bwperim(BW_PL1_CH1);
    t_PL2 = bwperim(BW_PL2_CH1);
    BW_CH1_sc_1 = ( cData.CH1.IC .*  (BW_CH1_sc_1 & BWm) ) ;
    if numel( find((t_BWm & t_PL1)==1)) <= 2
        [iy_1, ix_1] = find(BW_CH1_sc_1==(max(max(BW_CH1_sc_1)))) ;
        iy_1 = iy_1(1);    ix_1 = ix_1(1); 
        xp = ix_1+ (PL_Rad_CH1*cos(ang));      yp = iy_1+ (PL_Rad_CH1*sin(ang));
        BW_PL1_CH1 = roipoly( cData.CH1.IC, xp , yp );
    end

    BW_CH1_sc_2 = ( cData.CH1.IC .*  (BW_CH1_sc_2 & BWm) ) ;
    if numel( find((t_BWm & t_PL1)==1)) <= 2
        [iy_2, ix_2] = find(BW_CH1_sc_2==(max(max(BW_CH1_sc_2)))) ;
        iy_2 = iy_2(1);    ix_2 = ix_2(1); 
        xp = ix_2+ (PL_Rad_CH1*cos(ang));      yp = iy_2+ (PL_Rad_CH1*sin(ang));
        BW_PL2_CH1 = roipoly( cData.CH1.IC, xp , yp );
    end

    % We save the masks of the pole signal and define cell area without the poles
    cData.CH1.Mask_PL_1 = BW_PL1_CH1;        % OLD pole
    cData.CH1.Mask_PL_2 = BW_PL2_CH1;        % NEW polw
    cData.CH1.Mask_pCyto = logical(BW - ((BW & BW_PL1_CH1) + (BW & BW_PL2_CH1))) ;

% --- FLUORO CHANNEL 2 --- 222222222222222222222222222222222222222222222222
% --- POLE: detection and measurements ------------------------------------
%--------------------------------------------------------------------------    
    if APP_opt.t1_choose_Chan == 2         % ANALYSE Chan 2
        % We know that foci signal can only be at poles. Thus, we measure first
        % fit poles to circle and define an area where to search for the signal

    % ---> POLE 1 (OLD) search area
        xs = [ cData.R_mesh(1:Search_Cr_CH2,1); cData.R_mesh(1:Search_Cr_CH2,3) ];
        ys = [ cData.R_mesh(1:Search_Cr_CH2,2); cData.R_mesh(1:Search_Cr_CH2,4) ];
        [Cnt, R] = fit_circle( [xs,ys] ) ;
        ang = 0:0.1:2*pi;      xp = (R*sc_fc_CH2) *cos(ang);      yp = (R*sc_fc_CH2) *sin(ang);
        x_p = cData.R_mesh(1,1);         D_x = (x_p + Cnt(1)) / 2.08 ;
        y_p = cData.R_mesh(1,2);         D_y = (y_p + Cnt(2)) / 2.08 ;
        % create temp_mask of pole circle search area
        BW_CH2_sc_1  = roipoly( cData.CH2.IC, D_x+xp , D_y+yp ) ;

    % ---> POLE 2 (NEW) search area
        xs = [ cData.R_mesh(end-Search_Cr_CH2:end,1); cData.R_mesh(end-Search_Cr_CH2:end,3) ];
        ys = [ cData.R_mesh(end-Search_Cr_CH2:end,2); cData.R_mesh(end-Search_Cr_CH2:end,4) ];
        [Cnt, R] = fit_circle( [xs,ys] ) ;
        ang = 0:0.1:2*pi;      xp = (R*sc_fc_CH2) *cos(ang);      yp = (R*sc_fc_CH2) *sin(ang);
        x_p = cData.R_mesh(end,1);       D_x = (x_p + Cnt(1)) / 1.92  ;
        y_p = cData.R_mesh(end,2);       D_y = (y_p + Cnt(2)) / 1.92 ;
        % create temp_mask of pole circle search area
        BW_CH2_sc_2  = roipoly( cData.CH2.IC, D_x+xp , D_y+yp ) ;

        % We save those mask, just for plotting
        cData.CH2.M_sc_PL_1 = BW_CH2_sc_1;
        cData.CH2.M_sc_PL_2 = BW_CH2_sc_2;

    % --->  Search highest value pixel. 
        PL1_CH2 = cData.CH2.IC .* BW_CH2_sc_1 ;                  % Apply mask to isolate raw-value signal at the poles area
        PL2_CH2 = cData.CH2.IC .* BW_CH2_sc_2 ;
        [iy_1, ix_1] = find(PL1_CH2==(max(max(PL1_CH2)))) ;      % Now at each pole find the max(pixel_value)
        [iy_2, ix_2] = find(PL2_CH2==(max(max(PL2_CH2)))) ;
        iy_1 = iy_1(1);    iy_2 = iy_2(1);    
        ix_1 = ix_1(1);    ix_2 = ix_2(1);
        % Then create a smaller area of interest with the higherst pixel at the
        % center. This is now the detected polar signal   
        xp = ix_1+ (PL_Rad_CH2*cos(ang));      yp = iy_1+ (PL_Rad_CH2*sin(ang));
        BW_PL1_CH2 = roipoly( cData.CH2.IC, xp , yp );
        xp = ix_2+ (PL_Rad_CH2*cos(ang));     yp = iy_2+ (PL_Rad_CH2*sin(ang));
        BW_PL2_CH2 = roipoly( cData.CH2.IC, xp , yp );


        % In case that there are no polar foci and signal is cytosolic, the pole
        % are tend to be much inside the cell. In order to avoid that we need
        % restrict the area of search and look for highest value pixel, but
        % closer to the membrane

        % Create perimeters and test if the pole area cross the cell membrane
        t_BWm = bwperim(BW);
        t_PL1 = bwperim(BW_PL1_CH2);
        t_PL2 = bwperim(BW_PL2_CH2);
        BW_CH2_sc_1 = ( cData.CH2.IC .*  (BW_CH2_sc_1 & BWm) ) ;
        if numel( find((t_BWm & t_PL1)==1)) <= 2
            [iy_1, ix_1] = find(BW_CH2_sc_1==(max(max(BW_CH2_sc_1)))) ;
            iy_1 = iy_1(1);    ix_1 = ix_1(1); 
            xp = ix_1+ (PL_Rad_CH2*cos(ang));      yp = iy_1+ (PL_Rad_CH2*sin(ang));
            BW_PL1_CH2 = roipoly( cData.CH2.IC, xp , yp );
        end

        BW_CH2_sc_2 = ( cData.CH2.IC .*  (BW_CH2_sc_2 & BWm) ) ;
        if numel( find((t_BWm & t_PL1)==1)) <= 2
            [iy_2, ix_2] = find(BW_CH2_sc_2==(max(max(BW_CH2_sc_2)))) ;
            iy_2 = iy_2(1);    ix_2 = ix_2(1); 
            xp = ix_2+ (PL_Rad_CH2*cos(ang));      yp = iy_2+ (PL_Rad_CH2*sin(ang));
            BW_PL2_CH2 = roipoly( cData.CH2.IC, xp , yp );
        end

        % We save the masks of the pole signal and define cell area without the poles
        cData.CH2.Mask_PL_1 = BW_PL1_CH2;       % OLD pole
        cData.CH2.Mask_PL_2 = BW_PL2_CH2;       % NEW polw
        cData.CH2.Mask_pCyto = logical(BW - ((BW & BW_PL1_CH2) + (BW & BW_PL2_CH2))) ;   
    
    end
    
end

end


%--------------------------------------------------------------------------
%       Copyright (c) 2016, Matteo Sangermani, All rights reserved
%--------------------------------------------------------------------------

