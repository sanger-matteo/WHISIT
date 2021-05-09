function cData = Analysis_M2P(cData, ff, cc, Img_BF, Img_CH1, Img_CH2, Img_CH3)
%Analysis_PL = analyse the signal inside cells using algorithm M2P
% (Membrane 2 Poles). The cell area is analysed to identify several
% compartments: the membrane area, extending from the cell perimeter into
% the cell for X pixels, and the cytosol area, everything else that is not
% membrane. In turn, the membrane area is divided in 4 parts: the two polar
% membrane areas and the two "lateral" sides along the cell body.
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
% .Mask.Cell_body = this is a mask of the cell body area, defined using 
%       cell outline (meshes) as deternined in Oufti.
%
% .Mask.Memb_All = this mask, at its minimum, correspond exactely with the
%       cell perimeter (thichness 1 pixel). According to the parameter used
%       by the user, the membrane area can be extended inside the cell to
%       comprise a width of X pixel.
%
% .Mask.Cytosol = This mask is defined as:  .Cell_body - .Memb_All
%
% According to the parameter used by the user, the membrane areas extend
% inside the cell to comprise a width of X pixel. 
%
% .Mask.MembPole_1 and .Mask.MembPole_2 = the two masks define the pole 
%       areas. The algorithm intersect the membrane mask with thickness (d)
%       and a search circle (with specific radius) to determine the extend
%       of the areas.
%
% .Mask.MembLateral_1 and .Mask.MembLateral_2 = the two masks define the 
%       two lateral membrane area, found by subtracting .Mask.MembPole_X 
%       and .Cytosol. The thickness is determined by the parameter LateralMemb_Thick
%
% .Mask.Memb_All = this mask, at its minimum, correspond exactely with the
%       cell perimeter (thichness 1 pixel). According to the parameter used
%       by the user, the membrane area can be extended inside the cell to
%       comprise a width of X pixel.
%
% .Mask.Memb_All = this algorithat, unlike the others, define the mask as
%       the sum of .Mask.MembPole_1+...2 + .Mask.MembLateral_1 +...2
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
%
% Note: we use cData.Fluor_Chan(1).IC to define roipoly. Since, channel 1 
%       is mandatory to run the analysis, so we no not need to check that
%       it is provided




%% --- Initialize ---------------------------------------------------------
global APP_opt ;                         % Variable storing WHISIT options

extra_border = APP_opt.BorderBox ;                  % extra border area to add to a cell cropped image

% "Memb" parameters define the membrane areas' thicknesses in [pixel].
% Number of circle to repeat to extend the membrane area for Mask.Memb_All 
% roughly during analysis it create a perimeters and add,  inward to the
% cell, XYW pixel in width to membrane area.
LateralMemb_Thick = APP_opt.M2P_LatMemb ;   % [px], Thickness lateral areas 
PoleMemb_Thick = APP_opt.M2P_PoleMemb  ;    % [px], Thickness polar areas                                         
Rad_PoleArea = APP_opt.M2P_PoleRad ;        % [px], polar circle radius



%% -- Body of the Function --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                      

if ~isempty(cData.model)  &&  size(cData.mesh,2) == 4   
    
    % Save the option parammeters used for analysis
    cData.info.Analysis_Opt.extra_border = APP_opt.BorderBox ;
    cData.info.Analysis_Opt.Memb_sz = APP_opt.M2P_LatMemb ;
            
    % Evaluate the geometric parameters (cell meshes, contour and cell body 
    % axis); and create and store cropped images of the specific cell in
    % the different channels
    [cData, R_Xs, R_Ys] =  Analysis_GEOM(cData, extra_border, Img_BF, Img_CH1, Img_CH2, Img_CH3) ;
    
    
    % ---> Cell body Mask ------------------------------------------------
    % We create .Mask.Cell_body to define the area of entire cell
    BW = roipoly( cData.Fluor_Chan(1).IC, R_Xs, R_Ys );
    cData.Mask.Cell_body = BW;
    
    % Relative Pixel-wise Perimeter of the entire cell
    tempPx = regionprops(bwperim(BW),'PixelList');
    cData.R_P_model = tempPx.PixelList ;
     

    % First, create two temporary masks for the membrane:
    % - one whole membrane with thickness LateralMemb_Thick: 
    %       this is used to create the two .Mask.MembPole_X
    % - one whole membrane with thickness PL_Rad_CH1: 
    %       this is used to create the two .Mask.MembLateral_X
    
    % ---> Temporary polar Membrane Mask ----------------------------------    
    % Iterate below process to create a membrane masks of desired thickness 
    t_BWm = bwperim(BW);                        % define cell body perimeter
    for i = 1 : PoleMemb_Thick -1
        t_BWc = imfill(t_BWm,'holes') - t_BWm;  % create temporary cytosol mask:
        in_per =  bwperim(t_BWc) ;              % subtract 'holes' to perimeter to have new Cytosol
        t_BWm = t_BWm + in_per ;                % inner hole perimiter add it to Membrane area
    end
    pole_BWm = logical( t_BWm );

    % ---> Temporary lateral Membrane Mask -------------------------------- 
    % First, we create a perimeter mask for each side of the cell
    % (right/left), using the respective side of coordinates and the axis. 
    % Then we subtract each mask from the entire perimeter areas 
    % (built using both meshes sides).
    % This approach is faster and guarantee 100% that we always create two
    % and only two, separate lateral membrane area
    
    t_BWm = bwperim(BW);                        % define whole perimeter
     
    regionXs = [cData.R_mesh(:,1) ; flipud(cData.geom.R_axis(:,1))];
    regionYs = [cData.R_mesh(:,2) ; flipud(cData.geom.R_axis(:,2))];
    % define side 1-2 of perimeter
    t_BWm_1 = bwperim( roipoly( cData.Fluor_Chan(1).IC, regionXs, regionYs ) );      
    
    regionXs = [cData.R_mesh(:,3) ; flipud(cData.geom.R_axis(:,1))];
    regionYs = [cData.R_mesh(:,4) ; flipud(cData.geom.R_axis(:,2))]; 
    % define side 3-4 of perimeter
    t_BWm_2 = bwperim( roipoly( cData.Fluor_Chan(1).IC, regionXs, regionYs ) );
    
    % iterate below process to create desired thickness  
    for i = 1 : LateralMemb_Thick -1                          
        % both sides (whole perimeter)
        t_BWc = imfill(t_BWm,'holes') - t_BWm;      % create temporay cytosol mask:
        in_per =  bwperim(t_BWc) ;                  % subtract 'holes' to perimeter to have new Cytosol
        t_BWm = t_BWm + in_per ;                    % inner hole perimiter add it to Membrane area
        
        % side 1-2
        t_BWc = imfill(t_BWm_1,'holes') - t_BWm_1;  % create temporay cytosol mask:
        in_per =  bwperim(t_BWc) ;                  % subtract 'holes' to perimeter to have new Cytosol
        t_BWm_1 = t_BWm_1 + in_per ;                % inner hole perimiter add it to Membrane area
        
        % side 3-4
        t_BWc = imfill(t_BWm_2,'holes') - t_BWm_2;  % create temporay cytosol mask:
        in_per =  bwperim(t_BWc) ;                  % subtract 'holes' to perimeter to have new Cytosol
        t_BWm_2 = t_BWm_2 + in_per ;                % inner hole perimiter add it to Membrane area
    end
    
    % Convert the masks to logical ( where they "coincide": == 2 )
    lat_wBW = logical( t_BWm );
    lat_BWm_1 = logical( (t_BWm_1 + lat_wBW) == 2 );
    lat_BWm_2 = logical( (t_BWm_2 + lat_wBW) == 2 );

    
% --- DEFINE the final MASKs ----------------------------------------------
    % Using Var_Rad radius, draw two circles around the two poles,
    % necessary to define the extent of the pole areas

    % Pole 1 (OLD) search area
    ang = 0:0.1:2*pi;      x_cir = Rad_PoleArea *cos(ang);      y_cir = Rad_PoleArea *sin(ang);
    x_sO = cData.R_mesh(1,1);
    y_sO = cData.R_mesh(1,2);
    % create temp_mask of pole circle search area
    BW_search_PL_1 = roipoly( cData.Fluor_Chan(1).IC, x_sO+x_cir , y_sO+y_cir ) ;

    % Pole 2 (NEW) search area
    ang = 0:0.1:2*pi;      x_cir = Rad_PoleArea *cos(ang);      y_cir = Rad_PoleArea *sin(ang);
    x_sO = cData.R_mesh(end,1);
    y_sO = cData.R_mesh(end,2);
    % create temp_mask of pole circle search area
    BW_search_PL_2 = roipoly( cData.Fluor_Chan(1).IC, x_sO+x_cir , y_sO+y_cir ) ;
    
    
    % Now that we have polar membranes mask and search areas, we can
    % intersect them to find create final polar mask
    cData.Mask.MembPole_1= logical( (BW_search_PL_1 + pole_BWm) == 2) ;
    cData.Mask.MembPole_2 = logical( (BW_search_PL_2 + pole_BWm) == 2) ;

    % We need to create "fake" pole mask with thickness equal to that used
    % for the lateral membrane masks This can then be subtracted from the
    % temporary unified lateral Memb mask to create the final ones
    fake_PL_1 = logical( (BW_search_PL_1 + lat_wBW) == 2) ;
    fake_PL_2 = logical( (BW_search_PL_2 + lat_wBW) == 2) ;
    unified_lateral_Memb = logical((fake_PL_1 + fake_PL_2 + lat_wBW) == 1) ;
    
    % Finaly, create the two lateral membrane areas by intersecting with
    % the unified lateral Memb, where we already subtracted the poles    
    cData.Mask.MembLateral_1 = logical((unified_lateral_Memb + lat_BWm_1) == 2);
    cData.Mask.MembLateral_2 = logical((unified_lateral_Memb + lat_BWm_2) == 2);
    
    cData.Mask.Memb_All = logical(cData.Mask.MembPole_1+ cData.Mask.MembPole_2 +...
                               cData.Mask.MembLateral_1 + cData.Mask.MembLateral_2 );
    cData.Mask.Cytosol = logical( imfill(cData.Mask.Memb_All,'holes') - cData.Mask.Memb_All );

end % if cData ~ empty

end % func






