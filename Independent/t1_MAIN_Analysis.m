function t1_MAIN_Analysis(app)
%t1_MAIN_Analysis = MAIN function that guide guide trhough all the steps of
%analysis of the fluorescence signal in tab 1 of WHISIT
%
% The function can analyse any detection file whether det.mat was combined
% or not with a manual tracking. Information about the tracking are used
% only after analysis of the fluorescence channel(s), when the cloneList is
% created.
%
% INPUT ------------------------------------------------------------------
% Are gathered from the selection and path provided using WHISIT interface.
% The major ones include:
% - BF stack of images
% - Channel 1 stack of images
% - Channel 2 stack of images
%   (all must be .tif and in the name format is AnyExpName_XXXX.tif, where 
%   XXXX defining individual frames. Must start at 1)
%
% - DETection.mat file: result from cell detection in the BF channel done 
%   using Outfi. cellList is the most important variable that will be used
%   for the analysis in this function
%
% - Several parameters, general or specific to an algorithm, can be chosen
%   by the user using WHISIT interface.
% 
% OUTPUT ------------------------------------------------------------------
% Unless stopped during analysis, several outputs are automatically 
% generated:
% - Res.mat
% - Data.txt
% - clone_List.mat
%
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-



% --- STEP 1 --------------------------------------------------------------
% Find total number of digits at ending number of filenames (BF and CH_X)
%
% --- STEP 2 --------------------------------------------------------------
% Check that stack and files have equal sizes
%
% --- STEP 3 --------------------------------------------------------------
% Perform analysis of fluorescence signal and of cell detection outline
%
% --- STEP 4 --------------------------------------------------------------
% Save files: the results.mat, the tab-separated.txt and, if necessary, the
% cloneList.mat.


global APP_opt;                        % Variable storing WHISIT options



%% --- STEP 1 --------------------------------------------------------------

% Find .tif files inside BF stack or CH_1/2 stack the filename and calcolate the  
% total number of digits at ending number (filename_xxxx.tif) 
% [filename must be is same as foldName]

%--- BF DIGITS -----------------------------------------------------------%  
tot_BFdig = TotDigits_in_Filename(APP_opt.t1_srcFiles_BF, APP_opt.name_delimiters);

if tot_BFdig == -1              % then we have not found any tif correctly named
    app.TextOUT.Value = sprintf('\n%s\n%s',  '!!! FL stack contain uncorrect .tif filename !!!',...
                                'Format should be in form:   filename_xxxx.tif');
    APP_opt.ERROR = 1;          
    return;
end 

%--- CH_1 DIGITS ---------------------------------------------------------%
tot_CH1dig = TotDigits_in_Filename(APP_opt.t1_srcFiles_CH1, APP_opt.name_delimiters);

if tot_CH1dig == -1                 % then we have not found any tif correctly named
    app.TextOUT.Value = sprintf('\n%s\n%s',  '!!! FL stack contain uncorrect .tif filename !!!',...
                                'Format should be in form:   filename_xxxx.tif');
    APP_opt.ERROR = 1;          
    return;
end

%--- CH_2 DIGITS ---------------------------------------------------------%
tot_CH2dig = -1;
if APP_opt.t1_choose_Chan == 2      % ANALYSE Chan 2
    tot_CH2dig = TotDigits_in_Filename(APP_opt.t1_srcFiles_CH2, APP_opt.name_delimiters);
    
    if tot_CH2dig == -1             % then we have not found any tif correctly named
        app.TextOUT.Value = sprintf('\n%s\n%s',  '!!! FL stack contain uncorrect .tif filename !!!',...
                                    'Format should be in form:   filename_xxxx.tif');
        APP_opt.ERROR = 1;          
        return;
    end
end



%% --- STEP 2 -------------------------------------------------------------
% Check that all images have same size and same number of images. Also make
% sure that the detection.mat file has length equal to the number of stack
% images provided
if length(APP_opt.t1_srcFiles_BF) ~= length(APP_opt.t1_srcFiles_CH1)
    app.TextOUT.Value = sprintf('\n%s\n%s',   '!!! Unequal number of frames between the BF and channel 1 !!!'); 
    APP_opt.ERROR = 1;          return;
end

if APP_opt.t1_choose_Chan == 2
    if length(APP_opt.t1_srcFiles_BF) ~= length(APP_opt.t1_srcFiles_CH2)
        app.TextOUT.Value = sprintf('\n%s\n%s',   '!!! Unequal number of frames between the BF and channel 2 !!!'); 
        APP_opt.ERROR = 1;          return;
    end
    if length(APP_opt.t1_srcFiles_CH1) ~= length(APP_opt.t1_srcFiles_CH2)
        app.TextOUT.Value = sprintf('\n%s\n%s',   '!!! Unequal number of frames between the two and channels !!!'); 
        APP_opt.ERROR = 1;          return;
    end
end

szCH(2,:) = size(imread([APP_opt.t1_srcFiles_CH1(1).folder '/' APP_opt.t1_srcFiles_CH1(1).name]));
szCH( ~any(szCH,2), : ) = [];           % remove rows where are present zeros
rows = all(szCH(:,2) == szCH(1,2)) ;
cols = all(szCH(:,2) == szCH(1,2)) ;
if rows ~= 1 && cols ~= 1
    app.TextOUT.Value = sprintf('\n%s\n%s',  '!!! Frames are not equal size !!!');
    APP_opt.ERROR = 1;          return;
end



%% --- STEP 3 --------------------------------------------------------------
% Main body of the function, analyse the fluorescence and cell's outline

% --- Load DET.mat file from Oufti
app.TextOUT.Value = sprintf('\n%s',  'Loading Det.mat file' );
load([APP_opt.t1_path_Det , APP_opt.t1_file_Det ]);
app.TextOUT.Value = sprintf('\n%s',  'Processing ...' );

% Count number total number of cell in entire cellList.meshData. Will
% be used as count-down and progression to display in GUI of WHISIT
tot_N_cells = sum(cellfun(@ length, cellList.meshData)) ;

I_CH1 = [];       % Initialize images ... 
I_CH2 = [];       % (so that we do not pass empty arguments) 
n_c = 0;          % number of cell being analyse (out of tot_N_cells)
tot_frames = length(cellList.meshData) ;

for ff = 1 : length(cellList.meshData)      % go through all frames

    % Assess number of digit and open next Fluorescence frame 
    N_dig = length(num2str(ff)); 
    N1_null = repmat('0', [1, tot_CH1dig - N_dig]);
    I_CH1 = double(imread([ APP_opt.t1_path_CH1, APP_opt.t1_foldName_CH1 ,'/' APP_opt.t1_Prefix_CH1 '_' N1_null num2str(ff), '.tif'])) ; 
    % for each frame calculate the average background signal, the minimum and maximum pixel value
    [BkGr_Ch1, min_px_Ch1, max_px_Ch1] = BckGrd_GaussFilt( I_CH1, 2, 4);

    if APP_opt.t1_choose_Chan == 2          % if we analyse two channels 
        N2_null = repmat('0', [1, tot_CH2dig - N_dig]);
        I_CH2 = double(imread([ APP_opt.t1_path_CH2, APP_opt.t1_foldName_CH2 ,'/' APP_opt.t1_Prefix_CH2 '_' N2_null num2str(ff), '.tif'])) ; 
        % for each frame calculate the average background signal, the minimum and maximum pixel value
        [BkGr_Ch2, min_px_Ch2, max_px_Ch2] = BckGrd_GaussFilt( I_CH2, 2, 4);
    end   

    for cc = 1 : length(cellList.meshData{ff})             % go through each cell in ff frame
        if ~isempty(cellList.meshData{ff}{cc})  &&  ~isempty(cellList.meshData{ff}{cc}.mesh) 
            cellList.meshData{ff}{cc}.CH1.BkGr_Fr = BkGr_Ch1;       % Avg backgroung pixel value
            cellList.meshData{ff}{cc}.CH1.Min_px_Fr = min_px_Ch1;   % Min pixel value in whole frame
            cellList.meshData{ff}{cc}.CH1.Max_px_Fr = max_px_Ch1;   % Max pixel value in whole frame

            if APP_opt.t1_choose_Chan == 2          % if we analyse two channels
                cellList.meshData{ff}{cc}.CH2.BkGr_Fr = BkGr_Ch2;       % Avg backgroung pixel value
                cellList.meshData{ff}{cc}.CH2.Min_px_Fr = min_px_Ch2;   % Min pixel value in whole frame
                cellList.meshData{ff}{cc}.CH2.Max_px_Fr = max_px_Ch2;   % Max pixel value in whole frame
            end

            % Analyse Fluorescence channel(s) according to the algorithm chosen
            if     APP_opt.choice_AIS == 1 
                cellList.meshData{1,ff}{1,cc} = Analysis_AIS(cellList.meshData{1,ff}{1,cc}, ff, cc, I_CH1, I_CH2) ;
            elseif APP_opt.choice_M2C == 1 
                cellList.meshData{1,ff}{1,cc} = Analysis_M2C(cellList.meshData{1,ff}{1,cc}, ff, cc, I_CH1, I_CH2) ;
            elseif APP_opt.choice_PL  == 1  
                cellList.meshData{1,ff}{1,cc} = Analysis_PL( cellList.meshData{1,ff}{1,cc}, ff, cc, I_CH1, I_CH2) ;
            end

            % --- Display the analisys in composite figure(s)
            if APP_opt.choice_plot == 1        
                N_dig = length(num2str(ff)); 
                N_null = repmat('0', [1, tot_BFdig - N_dig]);
                Img_BF = double(imread([ APP_opt.t1_path_BF, APP_opt.t1_foldName_BF ,'/' APP_opt.t1_Prefix_BF '_' N_null num2str(ff), '.tif'])) ;

                if     APP_opt.choice_AIS == 1 ;     Display_Cell_AIS(cellList.meshData{1,ff}{1,cc}, Img_BF ) ;
                elseif APP_opt.choice_M2C == 1 ;     Display_Cell_M2C(cellList.meshData{1,ff}{1,cc}, Img_BF ) ;
                elseif APP_opt.choice_PL  == 1 ;     Display_Cell_PL( cellList.meshData{1,ff}{1,cc}, Img_BF ) ;     
                end        
            end        
        end
        
        if APP_opt.BREAK == 1       % if STOP button is pressed,
            return                  % interrupt function
        end 

        % Update the analysis counter in GUI of WHISIT
        app.TextOUT.Value = sprintf('%s\n%s', [ 'Cell  ', num2str(cc), '  in Frame  ' num2str(ff) '  of  ' num2str(tot_frames)], ...
                             [ 'Total cells analysed:  ', num2str(n_c) ,'  of  ', num2str(tot_N_cells) ]);
        n_c = n_c + 1 ;

    end  %for cc
    
end  %for ff


%% --- STEP 4 --------------------------------------------------------------
% SAVE the analysis performed (Res.mat) and create the eithr CellData.txt 
% or cloneList.mat, depending whether the Manual Tracking was performed

% Save the analysis parameters used insid the Results or clone_List file
% (we will need them especially for Lineage_Plotting (tab_3))
WHISIT_parameters.choice_AIS = APP_opt.choice_AIS ;
WHISIT_parameters.choice_M2C = APP_opt.choice_M2C ;
WHISIT_parameters.choice_PL  = APP_opt.choice_PL ;

if APP_opt.t1_choose_ManualTrack ~= 1         % if NO manual tracking was performed
    % --- SAVE analysis Res.mat file 
    app.TextOUT.Value = sprintf('\n%s',   'Saving Results ... ');
    if isempty(APP_opt.t1_exp_name)
        save( [APP_opt.t1_path_CH1  ,'/' , 'Res'] ,...
               'cellList', 'cellListN', 'WHISIT_parameters' ,... 
               'coefPCA', 'mCell', 'p', 'paramString', 'rawPhaseFolder', 'shiftfluo', 'shiftframes', 'weights');
    else
        save( [ APP_opt.t1_path_CH1 ,'/' , APP_opt.t1_exp_name , '_', 'Res'],...
               'cellList', 'cellListN', 'WHISIT_parameters' ,... 
               'coefPCA', 'mCell', 'p', 'paramString', 'rawPhaseFolder', 'shiftfluo', 'shiftframes', 'weights');
    end

    % --- SAVE .txt file 
    if     APP_opt.choice_AIS == 1 ;       Save_txt_AIS(cellList) ;
    elseif APP_opt.choice_M2C == 1 ;       Save_txt_M2C(cellList) ;
    elseif APP_opt.choice_PL  == 1 ;       Save_txt_PL( cellList) ; 
    end
    
elseif APP_opt.t1_choose_ManualTrack == 1     % if manual tracking was done, we DO NOT save .txt file

    % --- CREATE and SAVE Res_2.mat
    app.TextOUT.Value = sprintf('\n%s',  'Creating clone_List ... ');
    
    [cellList.meshData, LG] = Lng_Analysis_Tracks( cellList.meshData , cellTrack );
    
    app.TextOUT.Value = sprintf('\n%s',   'Saving Results ... ');        
    % we also include cellTrack and LG in Res_Det2Track.mat. This file cell
    % detection file is still compatible with Oufti, where the user can go
    % back to check and change cell polarity if wished or needed

    if isempty(APP_opt.t1_exp_name)
        save( [APP_opt.t1_path_CH1  ,'/' , 'Res_Det2Track'] ,...
               'cellTrack', 'cellList', 'cellListN', 'WHISIT_parameters', 'LG',... 
               'coefPCA', 'mCell', 'p', 'paramString', 'rawPhaseFolder', 'shiftfluo', 'shiftframes', 'weights');
    else
        save( [ APP_opt.t1_path_CH1 ,'/' , APP_opt.t1_exp_name , '_', 'Res_Det2Track'] ,...
               'cellTrack', 'cellList', 'cellListN', 'WHISIT_parameters', 'LG',... 
               'coefPCA', 'mCell', 'p', 'paramString', 'rawPhaseFolder', 'shiftfluo', 'shiftframes', 'weights');
    end
    
end % if manual_tracking

end






