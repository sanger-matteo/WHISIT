function cData = Analysis_AIS(cData, ff, cc, I_CH1, I_CH2)
%Analysis_AIS = analyse the signal inside cells using algorithm AIS
% (Average Signal). This simply involves finding the average intensity
% signal inside the cells for channel 1 (and 2).
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
% .CellMask = this is a mask of the cell body area. This is defined using 
%           outline coordinates as deternined in Oufti. It is simply a 
%           logic matrix (0 or 1 values) with size equal to .CHX.IC.
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

%% -- Body of the Function --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                      

if ~isempty(cData.model)  &&  size(cData.mesh,2) == 4
    % Save the option parammeters used for analysis
    cData.info.Original_Frame   = ff ;
    cData.info.Original_cellID  = cc ;    
    cData.info.Analysis_Opt.extra_border = APP_opt.BorderBox ;
    cData.info.Analysis_Opt.Memb_sz = -1 ;          % not used for this algorithm
            
    % Find Min and Max in X-Y axis for the c-th-cell modelinates
    Xs = [cData.mesh(:,1) ; flipud(cData.mesh(:,3))] ;
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

    
%     % -------- OPTIONAL ------------
%     % Ellipse fitting of the cell mesh outline (Relative modelinates)
%     [A, B, X0, Y0, Phi] = Geom_My_fit_ellipse__v2(xs,ys);   
%     cData.geom.A_Axis = A ;
%     cData.geom.B_Axis = B ;
%     cData.geom.X0 = X0 ;
%     cData.geom.Y0 =  Y0 ;
%     cData.geom.phi = Phi ;
%     cData.geom.Epsilon = (1 - B/A); 
    

% --- Create the Masks that define the whole cytosol area -----------------
%--------------------------------------------------------------------------

    % ---> 1 Cell Mask 
    % Whether I use CH1.IC or CH2.IC does not change anything here: 
    % We create a Mask BW, which can be used to extract the average signal
    % from both channels.
    BW = roipoly( cData.CH1.IC, R_Xs, R_Ys );
    cData.Mask_wCell = BW;    % The mask of the whole cell area
    
    %---------R_P_model------Relative_Perimeter_model----------------------
    tempPx = regionprops(bwperim(BW),'PixelList');
    cData.R_P_model = tempPx.PixelList ; 

    
%%%--->>> insert here CHECK PLOT 1 <<<---%%%

end

end


%--------------------------------------------------------------------------
%       Copyright (c) 2016, Matteo Sangermani, All rights reserved
%--------------------------------------------------------------------------



%% ----------- CHECK PLOTs ------------------------------------------------
% % insert at indicated points of the code to test and debug whether cell's
% outline and axis are identified properly
% axis equal
% hold on
% % Absolute mesh - whole frame
% plot(Xs,Ys,'b')
% plot(cData.geom.axis(:,1),cData.geom.axis(:,2),'r')
% % Relative mesh - cropped frame
% plot(cData.geom.R_axis(:,1),cData.geom.R_axis(:,2))
% plot(R_Xs,R_Ys, 'r') 
% plot(cData.R_model(:,1),cData.R_model(:,2),'b')
