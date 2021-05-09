function GUI_StartUp(app)
%GUI_StartUp = Initiate all the major variable of the Graphic User
%   Interface as well as the APP_opt, a variable that stores the global
%   options and parameters to run WHISIT
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

clc;
            
% Variable APP_opt store all the major information on the status of WHISIT
% and the choises, parameters chosen by the user in the GUI
global APP_opt ;

% Initialize an empty cell_Track for manual tracking
global CellTracks;          CellTracks = {};

% Determine where your m-file's folder is and add that folder plus all the
% subfolders to the path.
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

% my_DataBase variable, store and save analysis
APP_opt.my_DB = '' ; 

% Buttons Options
APP_opt.START_t1 = 0 ;             
APP_opt.PLOT_t3 = 0 ;
APP_opt.START_t5 = 0 ;
APP_opt.BREAK = 1 ;            
APP_opt.ERROR = 0 ;

APP_opt.name_delimiters = {'_', '.'};

%% --- TAB 1 --- Paths and Folders -------------------------------------------------

% Path, folder and prefix used for the Bright Field images 
APP_opt.t1_path_BF = ''; 
APP_opt.t1_foldName_BF = '' ;
APP_opt.t1_Prefix_BF = '';

% Path, folder and prefix used for the Fluorescence channel 1 images
APP_opt.t1_path_CH1 = '' ;
APP_opt.t1_foldName_CH1 = '' ; 
APP_opt.t1_Prefix_CH1 = '';          

% Path, folder and prefix used for the Fluorescence channel 2 images
APP_opt.t1_path_CH2 = '';         
APP_opt.t1_foldName_CH2 = '' ;
APP_opt.t1_Prefix_CH2 = '';     

% Path, folder and prefix used for the Fluorescence channel 2 images
APP_opt.t1_path_CH3 = '';         
APP_opt.t1_foldName_CH3 = '' ;
APP_opt.t1_Prefix_CH3 = '';     

% Path and filename used for the cell Detection.mat file
APP_opt.t1_path_Det = '' ;
APP_opt.t1_file_Det = '' ;            

APP_opt.t1_exp_name = '' ;      % User defined experiment name

% Path and filename used for the Detection_2_Track file
APP_opt.t1_path_Res_D2T = '' ;
APP_opt.t1_file_Res_D2T = '' ;

APP_opt.BorderBox = app.t1_Edit_BorderBox.Value ;                % extra border area to add to a cell cropped image
APP_opt.plot_pause = app.t1_Edit_WaitPlot.Value ;                % 2 [sec], display analysis plot for X sec,


%% --- TAB 1 --- Variables and input parameters for analysis --------------------

APP_opt.AIS_EvalAxisProfile  = app.t1_Check_AxisProfile.Value ;       % Chooice to analyse the profile line
APP_opt.AIS_AxisWidth        = app.t1_Edit_AxisWidth.Value;           % [px] width of the axial profile line analysed
APP_opt.AIS_ProfileType      = 1 ;                                    %  what type of profile line saved: 1 = mean; 2 = max value

APP_opt.AIS_EvalSegmentation = app.t1_Check_Segmentation.Value;       % Chooice to perform cell segmentation
APP_opt.AIS_SegmentLength    = app.t1_Edit_SegmentWidth.Value;        % length of each segment (in [pixels])

APP_opt.AIS_EvalPerim        = app.t1_Check_PerimeterSig.Value ;      % Chooice to analyse cell border signal
APP_opt.AIS_PerimWidth       = app.t1_Edit_PerimeterWidth.Value ;     % width in [pixel] to consider for the analysis
APP_opt.AIS_PerimSpacing     = app.t1_Edit_PerimeterSpacing.Value ;   % Fraction [0.25-1] of points to create for a new evenly spaced perimeter

APP_opt.Len_AxSig       = [];        % Length of Profile line of all cells         
APP_opt.Len_SegmSig     = [];        % Length of Segmentation line of all cells         
APP_opt.Len_BorderSig   = [];        % Length of Boder perimeter line of all cells         

% "Memb" parameters define membrane areas' thicknesses in [pixel]: to
% create a perimeters and add, inward to the cell, XYW pixel in width.
APP_opt.M2P_LatMemb     = app.t1_Edit_M2P_LatMemb.Value;         % [px], Thickness lateral areas         
APP_opt.M2P_PoleMemb    = app.t1_Edit_M2P_PoleMemb.Value ;       % [px], Thickness polar areas 
APP_opt.M2P_PoleRad     = app.t1_Edit_M2P_PoleRad.Value  ;       % [px], polar circle radius
        
APP_opt.PL_AllMemb      = app.t1_Edit_PL_Memb.Value ;            % [px], Thickness all membrane area              
APP_opt.PL_PoleMemb     = app.t1_Edit_PL_PoleMemb.Value ;        % [px], Thickness polar areas 
APP_opt.PL_PoleRad      = app.t1_Edit_PL_PoleRad.Value ;         % [px], Radius for final polar area  
APP_opt.PL_ScaleFact    = app.t1_Edit_PL_ScaleFact.Value  ;       % 0.7 , [%] Scale factor for search diameter ( = Avg cell width)  


APP_opt.choice_plot = 0 ;             % Display analysis process by plotting for each cell
APP_opt.choice_AIS = 1 ;              % Average Intensity Signal
APP_opt.choice_M2P = 0 ;              % Memb_2_Cyto
APP_opt.choice_PL = 0 ;               % PolarSign
% Analyse one or two fluorescent channels
APP_opt.t1_choose_BrightField = 0;    % Tab_1:  1 = Crop and store;     0 = do not analyse;    
APP_opt.t1_choose_Chan_2 = 0;         % Tab_1:  1 = analyse;            0 = do not analyse;    
APP_opt.t1_choose_Chan_3 = 0;         % Tab_1:  1 = analyse;            0 = do not analyse;    
% Choose whether to use channel 3 as a marker for a pole
APP_opt.t1_CH3_Marker = 0 ;
   
% If yes, perform analysis to create the cloneList file
APP_opt.t1_choose_CreateCloneList = app.t1_CheckBox_CreateCloneList.Value;




%% --- TAB 3 --- Paths and Variables for plotting ----------------------------------------------           
% Path and filename used for the clone cell's Detection.mat file
APP_opt.t3_path_cloneList = '' ;
APP_opt.t3_file_cloneList = '' ; 
APP_opt.t3_path_cloneFolder = '' ;
APP_opt.t3_fold_cloneFolder = '' ;

% Values and Parameters for plotting lineages (tab 3)       
% Inizialize the values for the choosing the plotting type 
value = app.t3_RadioBox_PlotType.SelectedObject;
Str_value = strsplit(value.Text, '-'); 
APP_opt.t3_PlotOpt_A = str2num(Str_value{1}(1));

value = app.t3_RadioBox_LineType.SelectedObject;
Str_value = strsplit(value.Text, '-');
APP_opt.t3_PlotOpt_B = str2num(Str_value{1}(1));  

Str_value = strsplit(app.t3_DropDownValue.Value, '-');
APP_opt.t3_PlotOpt_C = str2num(Str_value{1}(1));   

APP_opt.t3_PlotOpt_FounderN =  [];     % Cell Lineage to plot

APP_opt.t3_choose_ChannelMode = 1 ;     % Plot channel 1, 2 or both (Independetly) or plot the ratio of 1/2 or 2/1
APP_opt.t3_choose_NormColormap = 1 ;


APP_opt.t3_fpm_CH1 = app.t3_EditField_fpmCH1.Value;            
APP_opt.t3_display_PL  = app.t3_CheckBox_DisplayPoles.Value ;
APP_opt.t3_display_ExI = app.t3_CheckBox_DisplayExI.Value ;      
APP_opt.t3_display_cID = 0 ;


APP_opt.t3_height_PlotLine = app.t3_EditField_LineSize.Value;
APP_opt.t3_Dist_btw_PlotLines = app.t3_EditField_LinesDistance.Value; 
APP_opt.t3_BorderLineWidth = app.t3_EditField_Borderline.Value;
APP_opt.t3_B4_LenStep = app.t3_EditField_SegmentSize.Value;

APP_opt.t3_ColorMap_LUT_CH1 = (app.t3_DropDown_LUT_CH1.Value);            
APP_opt.t3_ColorMap_LUT_CH2 = (app.t3_DropDown_LUT_CH2.Value);
APP_opt.t3_Value_Range_CH1 = app.t3_ValueRange_CH1.Value ;
APP_opt.t3_Value_Range_CH2 = app.t3_ValueRange_CH2.Value ;

% Save plotted data in a .txt file
APP_opt.t3_choose_Save_txt = app.t3_Save_txt_CheckBox.Value;
APP_opt.t3_choose_AutoSave_Plot = 0 ;
APP_opt.t3_choose_Display_Plot = app.t3_Display_Plot.Value ;
APP_opt.t3_choose_Summary_Plot = app.t3_Summary_Plot.Value ;
% User defined experiment name
APP_opt.t3_exp_name = '';                               
% Analyse and plot a single file or all the .mat files in a folder        
APP_opt.t3_Generation_VS_Single = 1 ;       
% Store a list of the .mat file for intergeneration analysis        
APP_opt.t3_intergen_srcFiles = struct();            



%% --- TAB 5 --- Paths and Variables for Manual Tracking ------------------------------------
APP_opt.t5_ff = 1 ;        % keep track of frame number in manual tracking

% Path, folder and prefix used for the Bright Field images
APP_opt.t5_path_BF = '';              
APP_opt.t5_foldName_BF = '' ; 
APP_opt.t5_Prefix_BF = '';

% Path and filename used for the cell Detection.mat file           
APP_opt.t5_path_Det = '';             
APP_opt.t5_fileName_Det = '' ;

% Path and filename used for the manual tracking .txt file
APP_opt.t5_path_Track = '' ;          
APP_opt.t5_fileName_Track = '' ;

% Give a costum name when saving the track.txt file
APP_opt.t5_exp_name = '';  

% Variable store what to do after clicking on the image: Cell Tracking mode,
% (= 1) assign coordinates to cell track. Pole Switching mode (= 2) allows
% to change cell polarity
APP_opt.t5_Mode_CellTrack = 1;
APP_opt.t5_display_Pole = 1 ;  

APP_opt.t5_display_Det = app.t5_choose_Outline.Value ;
APP_opt.t5_display_Track = 0 ;            % Disply track type
APP_opt.t5_display_Color = [.9 .2 0];     % Display track with color RGB  
APP_opt.t5_display_Marker = 10 ;          % Displayed marker size  
APP_opt.t5_display_Outline = 0 ;          % cell outline style 

% Displayed colored frame border to aid tracking 
APP_opt.t5_display_TrackAid = app.t5_showAid_CheckBox.Value;

 

% The ColorOrder property of the axes contains the color order. To change
% the color order. Here we set different default values for the ColorOrder
% property.
co48 = [ 0.3922    0.7843    1.0000;     1.0000    0.0392         0;
         0.0784    0.8627         0;     1.0000    0.8431         0;
         0.1373    0.5490    1.0000;     1.0000    0.5137    0.9804;
         1.0000    0.5882         0;     0.1569    0.7451    0.4706;
         0.8627    0.0784    0.3137;     0.7451    0.9412    0.2549;
         1.0000    0.7059    0.1961;          0    0.7059    1.0000;
         0.8039    0.8039         0;     0.8627    0.3137    0.2353;
         1.0000    0.7137    0.7569;     0.4706    0.9412         0;
         0.3922    0.5882    0.9412;     1.0000    0.4157    0.4157;
         0.9608    0.2039    0.7059;     0.9412    0.5490    0.4706;
         0.2549    0.8824    0.8157;     0.9412    0.9412    0.7843;
         1.0000    0.5098    0.6667;     0.7451    0.2431    1.0000;
         0.9608    0.3137         0;     0.7451    0.9412    0.2549;
         0.2745    0.3922    0.9020;          0    0.9608    1.0000;
         0.8235    0.4314    0.3137;     0.5137    0.4353    1.0000;
         1.0000    0.2431    0.5882;     1.0000    0.9255    0.5490;
         0.6784    1.0000    0.0784;     0.1961    0.6275    0.7843;
              0    0.7843    0.5882;     0.7843    0.6784         0;
         0.7804    0.0824    0.5216;     0.9490    0.9490         0;
              0         0    0.9804;     0.8549    0.4392    0.8392;
         1.0000         0    1.0000;     0.5882    0.9608    1.0000;
         0.6275    1.0000    0.5882;          0    0.8078    0.8235;
         0.6078    0.1961    0.9333;     0.5294    0.8078    0.9804;
         0.8627    0.7843    0.3529;     0.6275    0.4706    0.8627;
         0.2353    1.0000    0.3922 ] ;
     
set(groot,'defaultAxesColorOrder',co48)

end






