function cData = Analysis_M2C(cData, ff, cc, I_CH1, I_CH2)
%Analysis_M2C = analyse the signal inside cells using algorithm M2C
% (Membrane 2 Cytosol). The cell area is analysed to identify two
% compartments: the membrane area, extending from the cell perimeter into
% the cell for X pixels, and the cytosol area, everything else that is not
% membrane.
%
%
% INPUTS ------------------------------------------------------------------
% cData = Stores the data concerning a specific cell cc at frame ff
%         (normally in the provided as cellList.meshData{ff}{cc})
% 
%
% OUTPUT ------------------------------------------------------------------
% cData = update cData is return back
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
% .M_wCell = this is a mask of the whole cell body . This is defined using
%           outline coordinates as deternined in Oufti. 
%
% .M_Memb = this mask, at its minimum, correspond exactely with the cell
%           perimeter. Therefore is just one pixel width. According to the
%           parameter used by the user, the membrane area can be extended
%           inside the cell to comprise a width of X pixel. 
%
% .M_Cyto = This mask is defined by subtracting M_wCell - M_Memb
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
%                     size(cData.CH1.IC,1), size(cData.CH1.IC,2))

    %---------R_P_model------Relative_Perimeter_model----------------------
    % Using round(.mesh) to define .model (see Retrive_BF_mesh.m), 
    % does not give a continuos line of pixel around the cell. But
    % with the mask we can can extract the perimeter and re-define cell 
    % coordinates so that they describe a continuous line.
    %
    tempPx = regionprops(bwperim(BW),'PixelList');
    cData.R_P_model = tempPx.PixelList ;   
    %
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
    
    % ---> 3 Final Cytosol and Membrane Masks and Save them in cData     
    %BWm = t_BWm + in_per;
    BWm = logical( t_BWm );
    BWc = logical( imfill(t_BWm,'holes') - BWm );      
    cData.Mask_wCell = BW;
    cData.Mask_mCyto = BWc;
    cData.Mask_Memb = BWm;

    
end % if cData not empty

end


%--------------------------------------------------------------------------
%       Copyright (c) 2016, Matteo Sangermani, All rights reserved
%--------------------------------------------------------------------------

