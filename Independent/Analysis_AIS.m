function cData = Analysis_AIS(cData, ff, cc, Img_BF, Img_CH1, Img_CH2, Img_CH3)
%Analysis_AIS = analyse the signal inside cells using algorithm AIS
% (Average Signal). This is primarily the average intensity signal inside
% cells. If the options are selected, it will also Find Spots and extract
% cell's signal profile along the main axis.
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
% .Mask.Cell_body = this is a mask of the cell body area, defined using cell
%       outline (meshes) as deternined in Oufti. It is a logic matrix
%       (0 or 1 values) with size equal to .IC.





%% --- Initialize ---------------------------------------------------------

global APP_opt ;                        % Variable storing WHISIT options

extra_border = APP_opt.BorderBox ;      % extra border area to add to a cell cropped image

%% -- Body of the Function --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                      

if ~isempty(cData.model)  &&  size(cData.mesh,2) == 4
    
    % Save the option parammeters used for analysis
    cData.info.Analysis_Opt.extra_border = APP_opt.BorderBox ;
    cData.info.Analysis_Opt.Memb_sz = -1 ;          % not used for this algorithm
         
    % Evaluate the geometric parameters (cell meshes, contour and cell body 
    % axis); and create and store cropped images of the specific cell in
    % the different channels
    [cData, R_Xs, R_Ys] =  Analysis_GEOM(cData, extra_border, Img_BF, Img_CH1, Img_CH2, Img_CH3) ;
    
    
    
    % ---> Cell Body Mask -------------------------------------------------
    % We create .Mask.Cell_body to define the area of entire cell
    BW = roipoly( cData.Fluor_Chan(1).IC, R_Xs, R_Ys );
    cData.Mask.Cell_body = BW;
    
    % Relative Pixel-wise Perimeter of the entire cell
    tempPx = regionprops(bwperim(BW),'PixelList');
    cData.R_P_model = tempPx.PixelList ; 
    

    % ---> Find Spots -----------------------------------------------------
    % If spots was detected in Oufti, their data is extracted and stored in 
    % cData. We also create a normalize position of each spot, calculated
    % as ratio between position and cell length 
    % (0 and 1 being the two poles, if pole Marking was done this results 
    % will be oriented accordingly)
    % --- !!! ONLY POSSIBLE FOR CHANNEL 3, when used as Pole MARKER !!! ---
    if APP_opt.t1_choose_Chan_3 == 1  &&  APP_opt.t1_CH3_Marker ~= 1  &&  isfield(cData,'spots')  
        for ii = 1 : length(cData.spots)
            cData.spots(ii).Rx = double( cData.spots(ii).x - (min(Xs)-sf - extra_border) );
            cData.spots(ii).Ry = double( cData.spots(ii).y - (min(Ys)-sf - extra_border) );
            cData.spots(ii).l_ratio = cData.spots(ii).l/cData.geom.length;
        end
    end
    
            
    % --- Extract cell's sectional signal ---------------------------------
    % Analyse the singal distribution along the cell body axis
    if APP_opt.AIS_EvalAxisProfile == 1     
        cData = AIS_AxialProfile(cData , APP_opt.AIS_AxisWidth);
        % Store the cell body axis length of any cell analysed
        % (we need to normalize the results when saving in Save_txt_AIS.m)
        APP_opt.Len_AxSig(end+1) = length(cData.Fluor_Chan(1).AxSig);
    end 

      
    % --- SEGMENTATION ---------------------------------------------------- 
    % Perform cell segmentation
    if APP_opt.AIS_EvalSegmentation == 1 
        cData = AIS_SegmentMesh(cData , APP_opt.AIS_SegmentLength ) ;
        % Store the total number of sgements of any cell analysed
        % (we need to normalize the results when saving in Save_txt_AIS.m)        
        APP_opt.Len_SegmSig(end+1) = size(cData.Fluor_Chan(1).SegmentSig ,2);
    end
    
    
    % --- Cell Border profile ---------------------------------------------
    % Perform signal analysis of cell's perimeter
    if APP_opt.AIS_EvalPerim == 1 
        cData = AIS_PerimeterProfile(cData , APP_opt.AIS_PerimWidth ) ;
        % Store the perimeter length of any cell analysed
        % (we need to normalize the results when saving in Save_txt_AIS.m)
        APP_opt.Len_PerimSig(end+1) = length(cData.Fluor_Chan(1).PerimSig);
    end

  
end % if cData ~ empty

end % func






%% -------- OPTIONAL ------------

% % Ellipse fitting of the cell mesh outline (Relative modelinates)
%
% [A, B, X0, Y0, Phi] = Geom_My_fit_ellipse__v2(Xs,Ys);   
% cData.geom.A_Axis = A ;
% cData.geom.B_Axis = B ;
% cData.geom.X0 = X0 ;
% cData.geom.Y0 =  Y0 ;
% cData.geom.phi = Phi ;
% cData.geom.Epsilon = (1 - B/A); 
    

