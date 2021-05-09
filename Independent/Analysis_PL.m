function cData = Analysis_PL(cData, ff, cc, Img_BF, Img_CH1, Img_CH2, Img_CH3)
%Analysis_PL = analyse the signal inside cells using algorithm PL
% (Polar Signal). The cell area is analysed to identify two compartments: 
% the two, distinct, poles of a cell and the remaining cytosoic area 
%
% Called primarily by t1_MAIN_Analysis
%
% INPUTS ------------------------------------------------------------------
% cData = Stores the data concerning a specific cell cc at frame ff
%         (normally is provided as cellList.meshData{ff}{cc})
% I_CH1 = the image of the Channel 1 
% I_CH2 = the image of the Channel 2
% I_CH3 = the image of the Channel 3
% ff = frame number
% cc = cell number at ff-th frame 
%
% OUTPUT ------------------------------------------------------------------
% cData = updated cData is returned
% 
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------



% --- INFO on important Variables ----------------------------------------- 
% .Fluor_Chan(X).IC = The function create a small cropped image the
%       specific cc-th cell in the each fluorescence image provided and
%       stores it. In the end the Results.mat file is a library of
%       containing the raw images as well, allowing to perform further
%       analysis without the need of the original images.
%
% Masks are generated to define specific cell areas. 
% They are logical matrix (0 or 1) with value of 1 if the pixel belong to
% the area of interest. Masks allow to easily analyse fluorescence because
% the (logic) mask can be used as a list of indexes to select specific
% areas in the fluorescence image. For example the average signal in a cell
% is simply  >>  mean(...Fluor_Chan(X).IC( ...Mask ));
%
% .Mask.Cell_body = this is a mask of the cell body area, defined using 
%       cell outline (meshes) as deternined in Oufti.
%
% (Although used in Analysis_M2P, the algorithm must define .Mask.Memb_All
% because it will be the starting point for defining the masks Mask.Pole_X) 
% .Mask.Memb_All = this mask, at its minimum, correspond exactely with the
%       cell perimeter (thichness 1 pixel). According to the parameter used
%       by the user, the membrane area can be extended inside the cell to
%       comprise a width of X pixel.
%
% .Mask.Cytosol = This mask is defined as:  .Cell_body - .Memb_All
%
% .Fluor_Chan(X).Mask.Pole_X = two masks define the two pole areas. Those
%       masks are defined using the information from fluorescence,
%       therefore they are specific to each channel. The algorithm first
%       search for a "bright" spot or cluster of signal around the pole
%       area. Once found it creates an area to around it that will define
%       the mask(s)
%
% .Fluor_Chan(X).Mask.Cytosol = This mask is channel specific Cytosol. 
%       Defined as : .Cell_body - (.Memb_All + .Pole_1 + .Pole_2 )
%
%
% The algorithm store important geometric parameter of the cell in .geom.
% New sets of coordinates and mesh are calculated to be relative to the 
% cropped image(s) created, for example:
% .R_mesh
% .R_model
% .geom.axis
% .geom.R_axis
% .geom.length
% .geom.area



%% --- Initialize ---------------------------------------------------------
global APP_opt ;                                    % Variable storing WHISIT options

extra_border = APP_opt.BorderBox ;                  % extra border area to add to a cell cropped image

% "Memb" parameters define the membrane areas' thicknesses in [pixel].
% Number of circle to repeat to extend the membrane area for Mask.Memb_All 
% roughly during analysis it create a perimeters and add,  inward to the
% cell, XYW pixel in width to membrane area.
Memb_Thick = APP_opt.PL_AllMemb ;
% Thickness of polar membrane, for creating polar search areas, but not Mask.Memb_All                                        
PoleMemb_Thick = APP_opt.PL_PoleMemb ; 
% [0-1], scale factor to enlarge of dimish the diameter of the poles' search circle
Scale_Search_Circle = APP_opt.PL_ScaleFact ;  
% Define the pole Masks arount the brightest pixel, with a radius of ...
Rad_PoleArea = APP_opt.PL_PoleRad ;                    



%% -- Body of the Function --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                      

if ~isempty(cData.model)  &&  size(cData.mesh,2) == 4   
    
    % Save the option parammeters used for analysis  
    cData.info.Analysis_Opt.extra_border = APP_opt.BorderBox ;
    cData.info.Analysis_Opt.Memb_sz = APP_opt.PL_AllMemb ;
            
    % Evaluate the geometric parameters (cell meshes, contour and cell body 
    % axis); and create and store cropped images of the specific cell in
    % the different channels
    [cData, R_Xs, R_Ys] =  Analysis_GEOM(cData, extra_border, Img_BF, Img_CH1, Img_CH2, Img_CH3) ;
    
    
    % ---> Create Cell body Mask ------------------------------------------
    % We create .Mask.Cell_body to define the area of entire cell
    BW = roipoly( cData.Fluor_Chan(1).IC, R_Xs, R_Ys );
    
    % Relative Pixel-wise Perimeter of the entire cell
    tempPx = regionprops(bwperim(BW),'PixelList');
    cData.R_P_model = tempPx.PixelList ;
        
    % ---> Create Cytosol and Membrane Masks ------------------------------ 
    % Iterate below process to create a membrane masks of desired thickness 
    t_BWm = bwperim(BW);                        % define cell body perimeter
    for i = 1 : Memb_Thick -1
        t_BWc = imfill(t_BWm,'holes') - t_BWm;  % create temporay cytosol mask:
        in_per =  bwperim(t_BWc) ;              % subtract 'holes' to perimeter to have new Cytosol
        t_BWm = t_BWm + in_per ;                % inner hole perimiter ...  
    end                                         % ... add it to Membrane area  
    BWm = logical( t_BWm );
    BWc = logical( imfill(t_BWm,'holes') - BWm );   
    % Store the masks
    cData.Mask.Cell_body = BW;
    cData.Mask.Cytosol = BWc;
    cData.Mask.Memb_All = BWm;
    
    
% --- POLE: detection and measurements ------------------------------------
    % First we create a new membrane area with the thickness specified by
    % PoleMemb_Thick rather than Memb_Thick
    t_BWm = bwperim(BW);                        % define cell body perimeter
    for i = 1 : PoleMemb_Thick -1
        t_BWc = imfill(t_BWm,'holes') - t_BWm;  % create temporay cytosol mask:
        in_per =  bwperim(t_BWc) ;              % subtract 'holes' to perimeter to have new Cytosol
        t_BWm = t_BWm + in_per ;                % inner hole perimiter ...  
    end                                         % ... add it to Membrane area  
    Mask_PoleMembr = logical( t_BWm );
    

    % ---> Define search areas at the poles
    % We know that foci signal can only be at poles. Thus, we first define
    % areas where to search for the signal at each pole. These area are
    % defined by intersercting the Mask_PoleMembr with a search circle of
    % diamter equal to the average cell width
    % We place the center of the seracg cicle a bit off from the "pole" point
    % at about the second/3 point of the axis away from the pole
          
    % Calculate the distances between each couple of point (left-right) 
    % defining the cell meshes. Then calculate the median
    cell_w = sqrt( (cData.R_mesh(:,1)-cData.R_mesh(:,3)).^2 + (cData.R_mesh(:,2)-cData.R_mesh(:,4)).^2 );
    m_width = median(cell_w);
    % Create a circle coordinates with center (0;0) and radius that is half
    % the median cell width and scaled by scale_factor(sc_fc)
    ang = 0:0.1:2*pi;     
    x_cir = ((m_width/2)*Scale_Search_Circle) *cos(ang);      
    y_cir = ((m_width/2)*Scale_Search_Circle) *sin(ang);
    
    % We place the center of the cell cicle a bit off from the "pole",
    % roughly at the 3-rd point on the axis away from each resp. pole
    x_Center_1 = cData.geom.R_axis( 3, 1);
    y_Center_1 = cData.geom.R_axis( 3, 2);
    x_Center_2 = cData.geom.R_axis( end-3, 1);
    y_Center_2 = cData.geom.R_axis( end-3, 2);
    
    % Create final search Masks for pole 1 and 2
    BW_search_1 = roipoly( cData.Fluor_Chan(1).IC, x_Center_1+x_cir , y_Center_1+y_cir ) ;
    BW_search_1 = (BW_search_1 + Mask_PoleMembr );      
    BW_search_1((BW_search_1 == 1)) = 0;
    BW_search_1((BW_search_1 == 2)) = 1;    
    BW_search_2 = roipoly( cData.Fluor_Chan(1).IC, x_Center_2+x_cir , y_Center_2+y_cir ) ;
    BW_search_2 = (BW_search_2 + Mask_PoleMembr );      
    BW_search_2((BW_search_2 == 1)) = 0;
    BW_search_2((BW_search_2 == 2)) = 1;
    
    
% --- FLUORO CHANNEL 1 --- 111111111111111111111111111111111111111111111111   
    % Store search Masks for pole 1 and 2
    cData.Fluor_Chan(1).Mask.searchAPole_1 = logical(BW_search_1) ;  
    cData.Fluor_Chan(1).Mask.searchAPole_2 = logical(BW_search_2) ;

  % --->  Search highest value pixel.  
    % Apply mask to isolate raw-value signal at the poles area  
    PL1_CH1 = cData.Fluor_Chan(1).IC .* BW_search_1;  
    PL2_CH1 = cData.Fluor_Chan(1).IC .* BW_search_2;
    % Now at each pole find the max(pixel_value)
    [iy_1, ix_1] = find(PL1_CH1==(max(max(PL1_CH1)))) ;       
    [iy_2, ix_2] = find(PL2_CH1==(max(max(PL2_CH1)))) ;
    iy_1 = iy_1(1);    iy_2 = iy_2(1);    
    ix_1 = ix_1(1);    ix_2 = ix_2(1);
    
    % Then create a smaller area of interest with the higherst pixel at the
    % center. This is now the detected polar signal   
    x_cir = ix_1+ (Rad_PoleArea*cos(ang));      
    y_cir = iy_1+ (Rad_PoleArea*sin(ang));
    BW_PL1_CH1 = roipoly( cData.Fluor_Chan(1).IC, x_cir , y_cir );
    x_cir = ix_2+ (Rad_PoleArea*cos(ang));     
    y_cir = iy_2+ (Rad_PoleArea*sin(ang));
    BW_PL2_CH1 = roipoly( cData.Fluor_Chan(1).IC, x_cir , y_cir );
       
    % In case that there are no polar foci and signal is cytosolic, the
    % pole tend to be inside the cell. In order to avoid this we need to
    % restrict the area of search and look for highest value pixel, but
    % closer to the membrane

    % Create perimeters and test if the pole area cross the cell membrane
    t_BWm = bwperim(BW);
    t_PL1 = bwperim(BW_PL1_CH1);
    t_PL2 = bwperim(BW_PL2_CH1);
    BW_search_1 = ( cData.Fluor_Chan(1).IC .*  (BW_search_1 & BWm) ) ;
    if numel( find((t_BWm & t_PL1)==1)) <= 2
        [iy_1, ix_1] = find(BW_search_1==(max(max(BW_search_1)))) ;
        iy_1 = iy_1(1);    ix_1 = ix_1(1); 
        x_cir = ix_1+ (Rad_PoleArea*cos(ang));     
        y_cir = iy_1+ (Rad_PoleArea*sin(ang));
        BW_PL1_CH1 = roipoly( cData.Fluor_Chan(1).IC, x_cir , y_cir );
    end

    BW_search_2 = ( cData.Fluor_Chan(1).IC .*  (BW_search_2 & BWm) ) ;
    if numel( find((t_BWm & t_PL1)==1)) <= 2
        [iy_2, ix_2] = find(BW_search_2==(max(max(BW_search_2)))) ;
        iy_2 = iy_2(1);    ix_2 = ix_2(1); 
        x_cir = ix_2+ (Rad_PoleArea*cos(ang));      
        y_cir = iy_2+ (Rad_PoleArea*sin(ang));
        BW_PL2_CH1 = roipoly( cData.Fluor_Chan(1).IC, x_cir , y_cir );
    end

    % We save the masks of the pole signal and define cell area without the poles
    cData.Fluor_Chan(1).Mask.Pole_1 = BW_PL1_CH1;        % OLD pole
    cData.Fluor_Chan(1).Mask.Pole_2 = BW_PL2_CH1;        % NEW polw
    cData.Fluor_Chan(1).Mask.Cytosol = logical(BW - ((BW & BW_PL1_CH1) + (BW & BW_PL2_CH1))) ;




% --- FLUORO CHANNEL 2 --- 222222222222222222222222222222222222222222222222
    if APP_opt.t1_choose_Chan_2 == 1         % evaluate Channel 2
        % Store search Masks for pole 1 and 2
        cData.Fluor_Chan(2).Mask.searchAPole_1 = logical(BW_search_1) ;
        cData.Fluor_Chan(2).Mask.searchAPole_2 = logical(BW_search_2) ;  

      % --->  Search highest value pixel. 
        % Apply mask to isolate raw-value signal at the poles area
        PL1_CH2 = cData.Fluor_Chan(2).IC .* BW_search_1 ;         
        PL2_CH2 = cData.Fluor_Chan(2).IC .* BW_search_2 ;
        % Now at each pole find the max(pixel_value)
        [iy_1, ix_1] = find(PL1_CH2==(max(max(PL1_CH2)))) ;       
        [iy_2, ix_2] = find(PL2_CH2==(max(max(PL2_CH2)))) ;
        iy_1 = iy_1(1);    iy_2 = iy_2(1);    
        ix_1 = ix_1(1);    ix_2 = ix_2(1);
        % Then create a smaller area of interest with the higherst pixel at the
        % center. This is now the detected polar signal   
        x_cir = ix_1+ (Rad_PoleArea*cos(ang));     
        y_cir = iy_1+ (Rad_PoleArea*sin(ang));
        BW_PL1_CH2 = roipoly( cData.Fluor_Chan(2).IC, x_cir , y_cir );
        x_cir = ix_2+ (Rad_PoleArea*cos(ang));     
        y_cir = iy_2+ (Rad_PoleArea*sin(ang));
        BW_PL2_CH2 = roipoly( cData.Fluor_Chan(2).IC, x_cir , y_cir );

        % In case that there are no polar foci and signal is cytosolic, the
        % pole tend to be inside the cell. In order to avoid this we need
        % to restrict the area of search and look for highest value pixel,
        % but closer to the membrane

        % Create perimeters and test if the pole area cross the cell membrane
        t_BWm = bwperim(BW);
        t_PL1 = bwperim(BW_PL1_CH2);
        t_PL2 = bwperim(BW_PL2_CH2);
        BW_search_1 = ( cData.Fluor_Chan(2).IC .*  (BW_search_1 & BWm) ) ;
        if numel( find((t_BWm & t_PL1)==1)) <= 2
            [iy_1, ix_1] = find(BW_search_1==(max(max(BW_search_1)))) ;
            iy_1 = iy_1(1);    ix_1 = ix_1(1); 
            x_cir = ix_1+ (Rad_PoleArea*cos(ang));      
            y_cir = iy_1+ (Rad_PoleArea*sin(ang));
            BW_PL1_CH2 = roipoly( cData.Fluor_Chan(2).IC, x_cir , y_cir );
        end

        BW_search_2 = ( cData.Fluor_Chan(2).IC .*  (BW_search_2 & BWm) ) ;
        if numel( find((t_BWm & t_PL1)==1)) <= 2
            [iy_2, ix_2] = find(BW_search_2==(max(max(BW_search_2)))) ;
            iy_2 = iy_2(1);    ix_2 = ix_2(1); 
            x_cir = ix_2+ (Rad_PoleArea*cos(ang));     
            y_cir = iy_2+ (Rad_PoleArea*sin(ang));
            BW_PL2_CH2 = roipoly( cData.Fluor_Chan(2).IC, x_cir , y_cir );
        end

        % We save the masks of the pole signal and define cell area without the poles
        cData.Fluor_Chan(2).Mask.Pole_1 = BW_PL1_CH2;       % OLD pole
        cData.Fluor_Chan(2).Mask.Pole_2 = BW_PL2_CH2;       % NEW polw
        cData.Fluor_Chan(2).Mask.Cytosol = logical(BW - ((BW & BW_PL1_CH2) + (BW & BW_PL2_CH2))) ;   
    
    end % if Chan 2




% --- FLUORO CHANNEL 3 ---3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-% 
    if APP_opt.t1_choose_Chan_3 == 1  &&  APP_opt.t1_CH3_Marker ~= 1
        % Store search Masks for pole 1 and 2
        cData.Fluor_Chan(3).Mask.searchAPole_1 = logical(BW_search_1) ;  
        cData.Fluor_Chan(3).Mask.searchAPole_2 = logical(BW_search_2) ;  
        
        % --->  Search highest value pixel. 
        % Apply mask to isolate raw-value signal at the poles area
        PL1_CH3 = cData.Fluor_Chan(3).IC .* BW_search_1 ;         
        PL2_CH3 = cData.Fluor_Chan(3).IC .* BW_search_2 ;
        % Now at each pole find the max(pixel_value)
        [iy_1, ix_1] = find(PL1_CH3==(max(max(PL1_CH3)))) ;       
        [iy_2, ix_2] = find(PL2_CH3==(max(max(PL2_CH3)))) ;
        iy_1 = iy_1(1);    iy_2 = iy_2(1);    
        ix_1 = ix_1(1);    ix_2 = ix_2(1);
        % Then create a smaller area of interest with the higherst pixel at the
        % center. This is now the detected polar signal   
        x_cir = ix_1+ (Rad_PoleArea*cos(ang));      
        y_cir = iy_1+ (Rad_PoleArea*sin(ang));
        BW_PL1_CH3 = roipoly( cData.Fluor_Chan(3).IC, x_cir , y_cir );
        x_cir = ix_2+ (Rad_PoleArea*cos(ang));    
        y_cir = iy_2+ (Rad_PoleArea*sin(ang));
        BW_PL2_CH3 = roipoly( cData.Fluor_Chan(3).IC, x_cir , y_cir );

        % In case that there are no polar foci and signal is cytosolic, the
        % pole tend to be inside the cell. In order to avoid this we need
        % to restrict the area of search and look for highest value pixel,
        % but closer to the membrane

        % Create perimeters and test if the pole area cross the cell membrane
        t_BWm = bwperim(BW);
        t_PL1 = bwperim(BW_PL1_CH3);
        t_PL2 = bwperim(BW_PL2_CH3);
        BW_search_1 = ( cData.Fluor_Chan(3).IC .*  (BW_search_1 & BWm) ) ;
        if numel( find((t_BWm & t_PL1)==1)) <= 2
            [iy_1, ix_1] = find(BW_search_1==(max(max(BW_search_1)))) ;
            iy_1 = iy_1(1);    ix_1 = ix_1(1); 
            x_cir = ix_1+ (Rad_PoleArea*cos(ang));      
            y_cir = iy_1+ (Rad_PoleArea*sin(ang));
            BW_PL1_CH3 = roipoly( cData.Fluor_Chan(3).IC, x_cir , y_cir );
        end

        BW_search_2 = ( cData.Fluor_Chan(3).IC .*  (BW_search_2 & BWm) ) ;
        if numel( find((t_BWm & t_PL1)==1)) <= 2
            [iy_2, ix_2] = find(BW_search_2==(max(max(BW_search_2)))) ;
            iy_2 = iy_2(1);    ix_2 = ix_2(1); 
            x_cir = ix_2+ (Rad_PoleArea*cos(ang));      
            y_cir = iy_2+ (Rad_PoleArea*sin(ang));
            BW_PL2_CH3 = roipoly( cData.Fluor_Chan(3).IC, x_cir , y_cir );
        end

        % We save the masks of the pole signal and define cell area without the poles
        cData.Fluor_Chan(3).Mask.Pole_1 = BW_PL1_CH3;       % OLD pole
        cData.Fluor_Chan(3).Mask.Pole_2 = BW_PL2_CH3;       % NEW polw
        cData.Fluor_Chan(3).Mask.Cytosol = logical(BW - ((BW & BW_PL1_CH3) + (BW & BW_PL2_CH3))) ;   
    
    end % if Chan 3
    

end % if cData ~ empty


end % func






