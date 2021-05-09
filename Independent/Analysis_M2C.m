function cData = Analysis_M2C(cData, ff, cc, Img_BF, Img_CH1, Img_CH2, Img_CH3)
%Analysis_M2C = analyse the signal inside cells using algorithm M2C
% (Membrane 2 Cytosol). The cell area is analysed to identify two
% compartments: the membrane area, extending from the cell perimeter into
% the cell for X pixels, and the cytosol area, everything else that is not
% membrane.
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
% .Mask.Memb_All = this mask, at its minimum, correspond exactely with the
%       cell perimeter (thichness 1 pixel). According to the parameter used
%       by the user, the membrane area can be extended inside the cell to
%       comprise a width of X pixel.
%
% .Mask.Cytosol = This mask is defined as:  .Cell_body - .Memb_All
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

global APP_opt ;                        % Variable storing WHISIT options

extra_border = APP_opt.BorderBox ;      % extra border area to add to a cell cropped image
Memb_sz = APP_opt.Memb_sz ;             % number of circle to repeat to extend the membrane area: roughly each 
                                        % cicle add a pixel width from perimeter inward to the cell to membrane area
                                       
%% -- Body of the Function --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                      

if ~isempty(cData.model)  &&  size(cData.mesh,2) == 4 
    
    % Save the option parammeters used for analysis
    cData.info.Analysis_Opt.extra_border = APP_opt.BorderBox ;
    cData.info.Analysis_Opt.Memb_sz = APP_opt.Memb_sz ;
            
    % Gather cell meshes as XY coordinates with handier variable name
    Xs = [cData.mesh(:,1) ; flipud(cData.mesh(:,3))] ;
    Ys = [cData.mesh(:,2) ; flipud(cData.mesh(:,4))] ; 
  
    % Store in cData a cropped image(s) of the cell in the Bright Field ( + extra-border) 
    if APP_opt.t1_choose_BrightField == 1 
        cData.Bright_Field.IC  = imcrop( Img_BF  , [(min(Xs)-extra_border), (min(Ys)-extra_border), ((max(Xs)-min(Xs))+extra_border*2), ((max(Ys)-min(Ys))+extra_border*2) ]);
    end

    % Storein cData a cropped fluorescent image(s) with additional extra-border 
    cData.Fluor_Chan(1).IC = imcrop( Img_CH1 , [(min(Xs)-extra_border), (min(Ys)-extra_border), ((max(Xs)-min(Xs))+extra_border*2), ((max(Ys)-min(Ys))+extra_border*2) ]);
    % if there is a Channel 2 selected
    if APP_opt.t1_choose_Chan_2 == 1          
        cData.Fluor_Chan(2).IC = imcrop( Img_CH2 , [(min(Xs)-extra_border), (min(Ys)-extra_border), ((max(Xs)-min(Xs))+extra_border*2), ((max(Ys)-min(Ys))+extra_border*2) ]);
    end
    % if there is a Channel 3 selected
    if APP_opt.t1_choose_Chan_3 == 1          
        cData.Fluor_Chan(3).IC = imcrop( Img_CH3 , [(min(Xs)-extra_border), (min(Ys)-extra_border), ((max(Xs)-min(Xs))+extra_border*2), ((max(Ys)-min(Ys))+extra_border*2) ]);
    end
   
    sf = 1;      % <<<----- correction factor for the shift between BF and Fluor. due
                 % to difference in modelinates system: origin at [0;0] or at [1;1]
    % modelinates Relative to cropped image of the cell (.CHX.IC)
    cData.R_mesh(:,1) = cData.mesh(:,1) - (min(Xs)-sf - extra_border) ;
    cData.R_mesh(:,2) = cData.mesh(:,2) - (min(Ys)-sf - extra_border) ;
    cData.R_mesh(:,3) = cData.mesh(:,3) - (min(Xs)-sf - extra_border) ;
    cData.R_mesh(:,4) = cData.mesh(:,4) - (min(Ys)-sf - extra_border) ;
    
    cData.R_model(:,1) = cData.model(:,1) - (min(Xs)-sf - extra_border) ;
    cData.R_model(:,2) = cData.model(:,2) - (min(Ys)-sf - extra_border) ;
    
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
%          
% --- Create the Masks that define cytosol and membrane area --------------
%--------------------------------------------------------------------------

    % ---> Whole Cell Mask ------------------------------------------------
    % We create .Mask.Cell_body to define the area of entire cell
    BW = roipoly( cData.Fluor_Chan(1).IC, R_Xs, R_Ys );
    
    % Relative Pixel-wise Perimeter of the entire cell
    tempPx = regionprops(bwperim(BW),'PixelList');
    cData.R_P_model = tempPx.PixelList ;
        
    % ---> 2 Temporary Masks
    % define contour perimeter
    t_BWm = bwperim(BW);
    for i = 1 : Memb_sz -1
        t_BWc = imfill(t_BWm,'holes') - t_BWm;  % create temporay cytosol mask:
        in_per =  bwperim(t_BWc) ;              % subtract 'holes' to perimeter to have new Cytosol
        t_BWm = t_BWm + in_per ;                % inner hole perimiter ...  
    end                                         % ... add it to Membrane area
    
    % ---> 3 Final Cytosol and Membrane Masks and Save them in cData     
    %BWm = t_BWm + in_per;
    BWm = logical( t_BWm );
    BWc = logical( imfill(t_BWm,'holes') - BWm );      
    cData.Mask.Cell_body = BW;
    cData.Mask.Cytosol = BWc;
    cData.Mask.Memb_All = BWm;

    
end % if cData not empty

end


%--------------------------------------------------------------------------
%       Copyright (c) 2016, Matteo Sangermani, All rights reserved
%--------------------------------------------------------------------------

