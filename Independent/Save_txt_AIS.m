function  Save_txt_AIS( cData, WHISIT_P, Path_Save, Exp_Name )  
%Save_txt_AIS = Save data analysed via algorithm AIS (Average Signal) in a
%   tab separated .txt file. Also, if selected in GUI, the function will
%   save two additional .txt files: Spot.txt and Profile.txt.
%
%
% INPUTS ------------------------------------------------------------------
% cData = (normally in the provided as cellList). Stores the analysis for 
%         all cells analysed  
%         cData.meshData - is a cell array which list all the frames
%         cData.meshData{ff} - each element contain the list of cells in
%                              the specific ff.th frame
%         cData.meshData{ff}{cc} - the data concerning a specific cell cc
%                                  at frame ff
%
% WHISIT_P = variable contains all the parameters used for analysis. This
%         will be used to decide which data to save and how.
%
% Path_Save = path to where to store the .txt file(s)
%
% Exp_Name = if provided, assign a costum name Prefix to the .txt file(s)
%
%
% OUTPUT ------------------------------------------------------------------
% Three tab separated .txt file can be generated as output:
% Cell_Data.txt = stores in each row the data relative to a specific cell.
%       This include the average signal value, areas, length...
% Spot.txt = contain the information of each spot detected in the analysis.
%       Each line is one detected spot
% Profile.txt = store the signal profile along the main axis of the cells
%       analysed. Each line is a single cell result
% Segment.txt = store the average signal from each the cell segmentation
%       process. Each line is a single cell result
%
%
% N.B.: In .txt file(s) empty position are filled with value -1. This allow 
%       to have only numeric results to easily export/import and handling
%       in i.e. Excell. It is possible because all values generated during 
%       analysis can never be negative numbers (although some can be 0.0 ).
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


%% --- Initialize ---------------------------------------------------------


% Create file Cell_Data.txt to store cell entries -----------------------------------------------------
if isempty(Exp_Name)
    filename_cell_txt = [ Path_Save ,'/', 'Cell_Data.txt'];
else
    filename_cell_txt = [ Path_Save ,'/' , Exp_Name , '_Cell_Data.txt'];
end
% Initialize first row with labels in each column 
file_C = fopen(filename_cell_txt, 'w+');
fprintf( file_C, 'Frame \tCell_Num \tCell_Length \tArea_Cell \t' );
fprintf( file_C, '\tCH1_Avg_Int_Cell \tCH1_Avg_BkGr_noise \tCH1_Min_px_value_Frame \tCH1_Max_px_value_Frame \t' );
fprintf( file_C, '\tCH2_Avg_Int_Cell \tCH2_Avg_BkGr_noise \tCH2_Min_px_value_Frame \tCH2_Max_px_value_Frame \t' );
fprintf( file_C, '\tCH3_Avg_Int_Cell \tCH3_Avg_BkGr_noise \tCH3_Min_px_value_Frame \tCH3_Max_px_value_Frame \t' );
fprintf( file_C, '\tSet_Polarity \t\n' );


% Create Profile.txt to store axial_Line entries -------------------------------------------------------
if WHISIT_P.AIS.eval_AxProfile == 1 
    if isempty(Exp_Name)
        filename_CH1_txt = [ Path_Save ,'/', 'Profile_CH1.txt'];
    else
        filename_CH1_txt = [ Path_Save ,'/' , Exp_Name , '_Profile_CH1.txt'];
    end
    % Initialize first row with labels for in column 
    file_AxSigCH1 = fopen(filename_CH1_txt, 'w+');
    fprintf( file_AxSigCH1, 'Frame \tCell Num \tCell Length \tArea Cell \tAxis_Width \t' );
    fprintf( file_AxSigCH1, '    \tProfile_Line_CH1 \t\n' );
    
    if WHISIT_P.eval_Channel_2 == 1 
        if isempty(Exp_Name)
            filename_CH2_txt = [ Path_Save ,'/', 'Profile_CH2.txt'];
        else
            filename_CH2_txt = [ Path_Save ,'/' , Exp_Name , '_Profile_CH2.txt'];
        end
        % Initialize first row with labels in each column 
        file_AxSigCH2 = fopen(filename_CH2_txt, 'w+');
        fprintf( file_AxSigCH2, 'Frame \tCell Num \tCell Length \tArea Cell \tAxis_Width \t' );
        fprintf( file_AxSigCH2, '    \tProfile_Line_CH2 \t\n' );
    end

    if WHISIT_P.eval_Channel_3 == 1  &&  WHISIT_P.MarkedPole ~= 1 
        if isempty(Exp_Name)
            filename_CH3_txt = [ Path_Save ,'/', 'Profile_CH3.txt'];
        else
            filename_CH3_txt = [ Path_Save ,'/' , Exp_Name , '_Profile_CH3.txt'];
        end
        % Initialize first row with labels in each column 
        file_AxSigCH3 = fopen(filename_CH3_txt, 'w+');
        fprintf( file_AxSigCH3, 'Frame \tCell Num \tCell Length \tArea Cell \tAxis_Width \t' );
        fprintf( file_AxSigCH3, '    \tProfile_Line_CH3 \t\n' );
    end
end


% Create Segmentation.txt to store cell segmentaton entries -------------------------------------------------------
if WHISIT_P.AIS.eval_Segment == 1 
    if isempty(Exp_Name)
        filename_CH1_txt = [ Path_Save ,'/', 'Segment_CH1.txt'];
    else
        filename_CH1_txt = [ Path_Save ,'/' , Exp_Name , '_Segment_CH1.txt'];
    end
    % Initialize first row with labels for in column 
    file_SegmentSigCH1 = fopen(filename_CH1_txt, 'w+');
    fprintf( file_SegmentSigCH1, 'Frame \tCell Num \tCell Length \tArea Cell \tSegment_Size_px\t' );
    fprintf( file_SegmentSigCH1, '    \tSegment_AvgSig_CH1 \t\n' );
    
    if WHISIT_P.eval_Channel_2 == 1 
        if isempty(Exp_Name)
            filename_CH2_txt = [ Path_Save ,'/', 'Segment_CH2.txt'];
        else
            filename_CH2_txt = [ Path_Save ,'/' , Exp_Name , '_Segment_CH2.txt'];
        end
        % Initialize first row with labels in each column 
        file_SegmentSigCH2 = fopen(filename_CH2_txt, 'w+');
        fprintf( file_SegmentSigCH2, 'Frame \tCell Num \tCell Length \tArea Cell \tSegment_Size_px\t' );
        fprintf( file_SegmentSigCH2, '    \tSegment_AvgSig_CH2 \t\n' );
    end

    if WHISIT_P.eval_Channel_3 == 1  &&  WHISIT_P.MarkedPole ~= 1 
        if isempty(Exp_Name)
            filename_CH3_txt = [ Path_Save ,'/', 'Segment_CH3.txt'];
        else
            filename_CH3_txt = [ Path_Save ,'/' , Exp_Name , '_Segment_CH3.txt'];
        end
        % Initialize first row with labels in each column 
        file_SegmentSigCH3 = fopen(filename_CH3_txt, 'w+');
        fprintf( file_SegmentSigCH3, 'Frame \tCell Num \tCell Length \tArea Cell \tSegment_Size_px\t' );
        fprintf( file_SegmentSigCH3, '    \tSegment_AvgSig_CH3 \t\n' );
    end
end


% Create Perimeter.txt to store cell's perimeter profile entries -------------------------------------------------------
if WHISIT_P.AIS.eval_Perim == 1 
    if isempty(Exp_Name)
        filename_CH1_txt = [ Path_Save ,'/', 'Perimeter_CH1.txt'];
    else
        filename_CH1_txt = [ Path_Save ,'/' , Exp_Name , '_Perimeter_CH1.txt'];
    end
    % Initialize first row with labels for in column 
    file_PerimeterSigCH1 = fopen(filename_CH1_txt, 'w+');
    fprintf( file_PerimeterSigCH1, 'Frame \tCell Num \tCell Length \tArea Cell \tPerimeter_Width_px \tSampling_points_distance_px \t' );
    fprintf( file_PerimeterSigCH1, '    \tPerimeter_AvgSig_CH1 \t\n' );
    
    if WHISIT_P.eval_Channel_2 == 1 
        if isempty(Exp_Name)
            filename_CH2_txt = [ Path_Save ,'/', 'Perimeter_CH2.txt'];
        else
            filename_CH2_txt = [ Path_Save ,'/' , Exp_Name , '_Perimeter_CH2.txt'];
        end
        % Initialize first row with labels in each column 
        file_PerimeterSigCH2 = fopen(filename_CH2_txt, 'w+');
        fprintf( file_PerimeterSigCH2, 'Frame \tCell Num \tCell Length \tArea Cell \tPerimeter_Width_px \tSampling_points_distance_px \t' );
        fprintf( file_PerimeterSigCH2, '    \tPerimeter_AvgSig_CH2 \t\n' );
    end

    if WHISIT_P.eval_Channel_3 == 1  &&  WHISIT_P.MarkedPole ~= 1 
        if isempty(Exp_Name)
            filename_CH3_txt = [ Path_Save ,'/', 'Perimeter_CH3.txt'];
        else
            filename_CH3_txt = [ Path_Save ,'/' , Exp_Name , '_Perimeter_CH3.txt'];
        end
        % Initialize first row with labels in each column 
        file_PerimeterSigCH3 = fopen(filename_CH3_txt, 'w+');
        fprintf( file_PerimeterSigCH3, 'Frame \tCell Num \tCell Length \tArea Cell \tPerimeter_Width_px \tSampling_points_distance_px \t' );
        fprintf( file_PerimeterSigCH3, '    \tPerimeter_AvgSig_CH3 \t\n' );
    end
end


% Create Spot.txt to store fluorescence spots entries -----------------------------------------------------
% --- !!! ONLY POSSIBLE FOR CHANNEL 3, when used as Pole MARKER !!! ---
if WHISIT_P.eval_Channel_3 == 1  &&  WHISIT_P.MarkedPole ~= 1  &&  isfield(cData,'spots') 
    % set dummy variable to enable faster check within the for-loops,
    % rather than doin a triple check (as above) every time
    SaveSpotDataCH3 = 1;            
    if isempty(Exp_Name)
        filename_Spt_txt = [ Path_Save ,'/', 'Spots_CH3.txt'];
    else
        filename_Spt_txt = [ Path_Save ,'/' , Exp_Name , '_Spots_CH3.txt'];
    end
    % Initialize first row with labels for in column 
    file_Spt = fopen(filename_Spt_txt, 'w+');
    fprintf( file_Spt, 'Frame \tCell Num \tCell Length \tArea Cell \t' );
    fprintf( file_Spt, '    \tSpot Number \tSpot_Axes_Pos \tSpot_Axes_D \tSpot_Axes_Ratio \tAbs_X \tAbs_Y \t\n' );
else
    SaveSpotDataCH3 = 0;
end


%% --- Data Extraction ----------------------------------------------------
for ff = 1 : length(cData.meshData)               % go through all frames
for cc = 1 : length(cData.meshData{ff})           % go through each cell
    
% temporary arrays to store data for each respective .txt file
t_C = [];     
t_axCH1 = [];
t_axCH2 = [];
t_axCH3 = [];
t_segCH1 = [];
t_segCH2 = [];
t_segCH3 = [];
t_periCH1 = [];
t_periCH2 = [];
t_periCH3 = [];
t_Spt = [];

if  ~isempty(cData.meshData{ff}{cc}) ...
        && ~isempty(cData.meshData{ff}{cc}.model) ...
        &&  size(cData.meshData{ff}{cc}.mesh, 2) == 4      % if it is not empty 
    
    t_C( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
    t_C( 2) =  cData.meshData{ff}{cc}.info.Oufti_cellID ;
    t_C( 3) =  cData.meshData{ff}{cc}.geom.length ;
        
    % Pixel Area Membrane and Cytosol and Cell
    t_C( 4) =  sum(sum( cData.meshData{ff}{cc}.Mask.Cell_body ));  
    
    t_C( 5) = -1;    
    % --- CHANNEL 1 -------1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1- % 
    % Average Pixel signal of Membrane and Cytosol and Cell: Absolute
    t_C( 6) =  mean(cData.meshData{ff}{cc}.Fluor_Chan(1).IC( cData.meshData{ff}{cc}.Mask.Cell_body ));    
    % Average background noise value in Frame
    t_C( 7) =  cData.meshData{ff}{cc}.Fluor_Chan(1).Avg_BkGr_value ;
    % Min and Max background pixel value in Frame
    t_C( 8) =  cData.meshData{ff}{cc}.Fluor_Chan(1).Min_px_value;
    t_C( 9) =  cData.meshData{ff}{cc}.Fluor_Chan(1).Max_px_value;
    
    t_C( 10) = -1;
    % --- CHANNEL 2 ---------2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2- %  
    if WHISIT_P.eval_Channel_2 == 1
        % Average Pixel signal of Membrane and Cytosol and Cell: Absolute
        t_C( 11) =  mean(cData.meshData{ff}{cc}.Fluor_Chan(2).IC( cData.meshData{ff}{cc}.Mask.Cell_body ));
        % Average background noise value in Frame
        t_C( 12) =  cData.meshData{ff}{cc}.Fluor_Chan(2).Avg_BkGr_value ;
        % Min and Max background pixel value in Frame
        t_C( 13) =  cData.meshData{ff}{cc}.Fluor_Chan(2).Min_px_value;
        t_C( 14) =  cData.meshData{ff}{cc}.Fluor_Chan(2).Max_px_value;
    else
        % skip all this columns 
        t_C( 11: 14) = -1;
    end

    t_C( 15) = -1;
    % --- CHANNEL 3 ---------3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3- %  
    if WHISIT_P.eval_Channel_3 == 1  &&  WHISIT_P.MarkedPole ~= 1
        % Average Pixel signal of Membrane and Cytosol and Cell: Absolute
        t_C( 16) =  mean(cData.meshData{ff}{cc}.Fluor_Chan(3).IC( cData.meshData{ff}{cc}.Mask.Cell_body ));
        % Average background noise value in Frame
        t_C( 17) =  cData.meshData{ff}{cc}.Fluor_Chan(3).Avg_BkGr_value ;
        % Min and Max background pixel value in Frame
        t_C( 18) =  cData.meshData{ff}{cc}.Fluor_Chan(3).Min_px_value;
        t_C( 19) =  cData.meshData{ff}{cc}.Fluor_Chan(3).Max_px_value;
    else
        % skip all this columns 
        t_C( 16: 19) = -1;
    end  

    t_C( 20) = -1;
    t_C( 21) = cData.meshData{ff}{cc}.polarity ;
    
    t_C( 22) = -1;   
    
    % --- SAVE file_C as single entry-row
    fprintf( file_C , '%f\t', t_C(:) );
    fprintf( file_C , '\n' );
    

    
    
    
%% --- Axial Profile Line --------------------------------------------------
    if WHISIT_P.AIS.eval_AxProfile == 1         % if user choose the option

        % --- CHANNEL 1 -------1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1- %
        t_axCH1( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
        t_axCH1( 2) =  cData.meshData{ff}{cc}.info.Oufti_cellID ;
        t_axCH1( 3) =  cData.meshData{ff}{cc}.geom.length ;

        % Pixel Area Membrane and Cytosol and Cell
        t_axCH1( 4) =  sum(sum( cData.meshData{ff}{cc}.Mask.Cell_body ));

        % Axis width used for analysis
        t_axCH1( 5) =  WHISIT_P.AIS.AxWidth ;        
        t_axCH1( 6) = -1;        
        
        % We add zeros left-right in order to center the profiles line and 
        % so all cells are aligned around midpoint. Emptpy values are = -1
        Max_length_Line = max( WHISIT_P.AIS.Len_AxSig );
        aLine = (cData.meshData{ff}{cc}.Fluor_Chan(1).AxSig) ;    
        Z = Max_length_Line - length(aLine);

        if mod(Z,2) == 0
            plus = ones(1, (Z/2))*(-1);
            aLine = [plus, aLine , plus];
        else
            plus = ones(1, fix(Z/2))*(-1);      % fix(), division returning integer
            aLine = [plus, aLine , plus, -1];
        end
        % Write in cell array the profile line of the foci
        for pp = 1 : length(aLine)
            t_axCH1( 7+pp ) = aLine(pp) ;
        end

        % --- SAVE file_C as single entry-row      
        fprintf( file_AxSigCH1 , '%f\t', t_axCH1(:) );
        fprintf( file_AxSigCH1 , '\n' );

        
        % --- CHANNEL 2 ---------2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2- %  
        if WHISIT_P.eval_Channel_2 == 1
            t_axCH2( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
            t_axCH2( 2) =  cData.meshData{ff}{cc}.info.Oufti_cellID ;
            t_axCH2( 3) =  cData.meshData{ff}{cc}.geom.length ;

            % Pixel Area Membrane and Cytosol and Cell
            t_axCH2( 4) =  sum(sum( cData.meshData{ff}{cc}.Mask.Cell_body ));

            % Axis width used for analysis
            t_axCH2( 5) =  WHISIT_P.AIS.AxWidth ;
            t_axCH2( 6) = -1;        

            % We add zeros left-right in order to center the profiles line and 
            % so all cells are aligned around midpoint. Emptpy values are = -1
            Max_length_Line = max( WHISIT_P.AIS.Len_AxSig );
            aLine = (cData.meshData{ff}{cc}.Fluor_Chan(2).AxSig) ;    
            Z = Max_length_Line - length(aLine);

            if mod(Z,2) == 0
                plus = ones(1, (Z/2))*(-1);
                aLine = [plus, aLine , plus];
            else
                plus = ones(1, fix(Z/2))*(-1);      % fix(), division returning integer
                aLine = [plus, aLine , plus, -1];
            end
            % Write in cell array the profile line of the foci
            for pp = 1 : length(aLine)
                t_axCH2( 7+pp ) = aLine(pp) ;
            end
            
            % --- SAVE file_C as single entry-row
            fprintf( file_AxSigCH2 , '%f\t', t_axCH2(:) );
            fprintf( file_AxSigCH2 , '\n' );
        end  
        

        % --- CHANNEL 3 ---------3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3- %  
        if WHISIT_P.eval_Channel_3 == 1  &&  WHISIT_P.MarkedPole ~= 1
            t_axCH3( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
            t_axCH3( 2) =  cData.meshData{ff}{cc}.info.Oufti_cellID ;
            t_axCH3( 3) =  cData.meshData{ff}{cc}.geom.length ;

            % Pixel Area Membrane and Cytosol and Cell
            t_axCH3( 4) =  sum(sum( cData.meshData{ff}{cc}.Mask.Cell_body ));

            % Axis width used for analysis
            t_axCH3( 5) =  WHISIT_P.AIS.AxWidth ;
            t_axCH3( 6) = -1;        

            % We add zeros left-right in order to center the profiles line and 
            % so all cells are aligned around midpoint. Emptpy values are = -1
            Max_length_Line = max( WHISIT_P.AIS.Len_AxSig );
            aLine = (cData.meshData{ff}{cc}.Fluor_Chan(3).AxSig) ;    
            Z = Max_length_Line - length(aLine);

            if mod(Z,2) == 0
                plus = ones(1, (Z/2))*(-1);
                aLine = [plus, aLine , plus];
            else
                plus = ones(1, fix(Z/2))*(-1);      % fix(), division returning integer
                aLine = [plus, aLine , plus, -1];
            end
            % Write in cell array the profile line of the foci
            for pp = 1 : length(aLine)
                t_axCH3( 7+pp ) = aLine(pp) ;
            end
            
            % --- SAVE file_C as single entry-row
            fprintf( file_AxSigCH3 , '%f\t', t_axCH3(:) );
            fprintf( file_AxSigCH3 , '\n' );
        end  
    
    end %end Profile_Line
  


    

%% --- Cell Segmentation average signal --------------------------------------------------
    if WHISIT_P.AIS.eval_Segment == 1         % if user choose the option

        % --- CHANNEL 1 -------1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1- %
        t_segCH1( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
        t_segCH1( 2) =  cData.meshData{ff}{cc}.info.Oufti_cellID ;
        t_segCH1( 3) =  cData.meshData{ff}{cc}.geom.length ;

        % Pixel Area Membrane and Cytosol and Cell
        t_segCH1( 4) =  sum(sum( cData.meshData{ff}{cc}.Mask.Cell_body ));

        % Axis width used for analysis
        t_segCH1( 5) =  WHISIT_P.AIS.SegmentLength ;        
        t_segCH1( 6) = -1;   
        
        % We add zeros left-right in order to center the segmentation array  
        % and all cells are aligned around midpoint. Emptpy values are = -1        
        Max_length_Line = max( WHISIT_P.AIS.Len_SegmSig);
        aLine = (cData.meshData{ff}{cc}.Fluor_Chan(1).SegmentSig ) ;    
        Z = Max_length_Line - length(aLine);

        if mod(Z,2) == 0
            plus = ones(1, (Z/2))*(-1);
            aLine = [plus, aLine , plus];
        else
            plus = ones(1, fix(Z/2))*(-1);      % fix(), division returning integer
            aLine = [plus, aLine , plus, -1];
        end
        % Write in cell array the profile line of the foci
        for pp = 1 : length(aLine)
            t_segCH1( 7+pp ) = aLine(pp) ;
        end

        % --- SAVE file_C as single entry-row      
        fprintf( file_SegmentSigCH1 , '%f\t', t_segCH1(:) );
        fprintf( file_SegmentSigCH1 , '\n' );

        
        % --- CHANNEL 2 ---------2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2- %  
        if WHISIT_P.eval_Channel_2 == 1
            t_segCH2( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
            t_segCH2( 2) =  cData.meshData{ff}{cc}.info.Oufti_cellID ;
            t_segCH2( 3) =  cData.meshData{ff}{cc}.geom.length ;

            % Pixel Area Membrane and Cytosol and Cell
            t_segCH2( 4) =  sum(sum( cData.meshData{ff}{cc}.Mask.Cell_body ));

            % Axis width used for analysis
            t_segCH1( 5) =  WHISIT_P.AIS.SegmentLength ;        
            t_segCH1( 6) = -1;   

            % We add zeros left-right in order to center the segmentation array  
            % and all cells are aligned around midpoint. Emptpy values are = -1        
            Max_length_Line = max( WHISIT_P.AIS.Len_SegmSig);
            aLine = (cData.meshData{ff}{cc}.Fluor_Chan(2).SegmentSig ) ;    
            Z = Max_length_Line - length(aLine);

            if mod(Z,2) == 0
                plus = ones(1, (Z/2))*(-1);
                aLine = [plus, aLine , plus];
            else
                plus = ones(1, fix(Z/2))*(-1);      % fix(), division returning integer
                aLine = [plus, aLine , plus, -1];
            end
            % Write in cell array the profile line of the foci
            for pp = 1 : length(aLine)
                t_segCH2( 7+pp ) = aLine(pp) ;
            end
            
            % --- SAVE file_C as single entry-row
            fprintf( file_SegmentSigCH2 , '%f\t', t_segCH2(:) );
            fprintf( file_SegmentSigCH2 , '\n' );
        end  
        

        % --- CHANNEL 3 ---------3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3- %  
        if WHISIT_P.eval_Channel_3 == 1  &&  WHISIT_P.MarkedPole ~= 1
            t_segCH3( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
            t_segCH3( 2) =  cData.meshData{ff}{cc}.info.Oufti_cellID ;
            t_segCH3( 3) =  cData.meshData{ff}{cc}.geom.length ;

            % Pixel Area Membrane and Cytosol and Cell
            t_segCH3( 4) =  sum(sum( cData.meshData{ff}{cc}.Mask.Cell_body ));

            % Axis width used for analysis
            t_segCH1( 5) =  WHISIT_P.AIS.SegmentLength ;        
            t_segCH1( 6) = -1;   
            
            % We add zeros left-right in order to center the segmentation array  
            % and all cells are aligned around midpoint. Emptpy values are = -1         
            Max_length_Line = max( WHISIT_P.AIS.Len_SegmSig);
            aLine = (cData.meshData{ff}{cc}.Fluor_Chan(3).SegmentSig ) ;    
            Z = Max_length_Line - length(aLine);

            if mod(Z,2) == 0
                plus = ones(1, (Z/2))*(-1);
                aLine = [plus, aLine , plus];
            else
                plus = ones(1, fix(Z/2))*(-1);      % fix(), division returning integer
                aLine = [plus, aLine , plus, -1];
            end
            % Write in cell array the profile line of the foci
            for pp = 1 : length(aLine)
                t_segCH3( 7+pp ) = aLine(pp) ;
            end
            
            % --- SAVE file_C as single entry-row
            fprintf( file_SegmentSigCH3 , '%f\t', t_segCH3(:) );
            fprintf( file_SegmentSigCH3 , '\n' );
        end  
    
    end %end Profile_Line
  

    

    
%% --- Perimeter Profile Line --------------------------------------------------
    if WHISIT_P.AIS.eval_Perim == 1         % if user choose the option

        % --- CHANNEL 1 -------1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1- %
        t_periCH1( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
        t_periCH1( 2) =  cData.meshData{ff}{cc}.info.Oufti_cellID ;
        t_periCH1( 3) =  cData.meshData{ff}{cc}.geom.length ;
        % Pixel Area of Cell body
        t_periCH1( 4) =  sum(sum( cData.meshData{ff}{cc}.Mask.Cell_body ));

        % Axis width used for analysis
        t_periCH1( 5) =  WHISIT_P.AIS.PerimWidth ;
        aLine = (cData.meshData{ff}{cc}.Fluor_Chan(1).PerimSig ) ;
        t_periCH1( 6) =  (length(aLine)/WHISIT_P.AIS.PerimSpacing) / length(aLine) ;
        t_periCH1( 7) = -1;   
        
        % We add zeros left-right in order to center the segmentation array  
        % and all cells are aligned around midpoint. Emptpy values are = -1        
        Max_length_Line = max( WHISIT_P.AIS.Len_PerimSig);
            
        Z = Max_length_Line - length(aLine);

        if mod(Z,2) == 0
            plus = ones(1, (Z/2))*(-1);
            aLine = [plus, aLine , plus];
        else
            plus = ones(1, fix(Z/2))*(-1);      % fix(), division returning integer
            aLine = [plus, aLine , plus, -1];
        end
        % Write in cell array the profile line of the foci
        for pp = 1 : length(aLine)
            t_periCH1( 8+pp ) = aLine(pp) ;
        end

        % --- SAVE file_C as single entry-row      
        fprintf( file_PerimeterSigCH1 , '%f\t', t_periCH1(:) );
        fprintf( file_PerimeterSigCH1 , '\n' );

        
        % --- CHANNEL 2 ---------2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2- %  
        if WHISIT_P.eval_Channel_2 == 1
            t_periCH2( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
            t_periCH2( 2) =  cData.meshData{ff}{cc}.info.Oufti_cellID ;
            t_periCH2( 3) =  cData.meshData{ff}{cc}.geom.length ;
            % Pixel Area of Cell body
            t_periCH2( 4) =  sum(sum( cData.meshData{ff}{cc}.Mask.Cell_body ));

            % Axis width used for analysis
            t_periCH2( 5) =  WHISIT_P.AIS.PerimWidth ;
            aLine = (cData.meshData{ff}{cc}.Fluor_Chan(2).PerimSig ) ;
            t_periCH2( 6) =  (length(aLine)/WHISIT_P.AIS.PerimSpacing) / length(aLine) ;
            t_periCH2( 7) = -1;   

            % We add zeros left-right in order to center the segmentation array  
            % and all cells are aligned around midpoint. Emptpy values are = -1        
            Max_length_Line = max( WHISIT_P.AIS.Len_PerimSig);

            Z = Max_length_Line - length(aLine);

            if mod(Z,2) == 0
                plus = ones(1, (Z/2))*(-1);
                aLine = [plus, aLine , plus];
            else
                plus = ones(1, fix(Z/2))*(-1);      % fix(), division returning integer
                aLine = [plus, aLine , plus, -1];
            end
            % Write in cell array the profile line of the foci
            for pp = 1 : length(aLine)
                t_periCH2( 8+pp ) = aLine(pp) ;
            end

            % --- SAVE file_C as single entry-row      
            fprintf( file_PerimeterSigCH2 , '%f\t', t_periCH2(:) );
            fprintf( file_PerimeterSigCH2 , '\n' );
        end  
        

        % --- CHANNEL 3 ---------3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3- %  
        if WHISIT_P.eval_Channel_3 == 1  &&  WHISIT_P.MarkedPole ~= 1
            t_periCH3( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
            t_periCH3( 2) =  cData.meshData{ff}{cc}.info.Oufti_cellID ;
            t_periCH3( 3) =  cData.meshData{ff}{cc}.geom.length ;
            % Pixel Area of Cell body
            t_periCH3( 4) =  sum(sum( cData.meshData{ff}{cc}.Mask.Cell_body ));

            % Axis width used for analysis
            t_periCH3( 5) =  WHISIT_P.AIS.PerimWidth ;
            aLine = (cData.meshData{ff}{cc}.Fluor_Chan(3).PerimSig ) ;
            t_periCH3( 6) =  (length(aLine)/WHISIT_P.AIS.PerimSpacing) / length(aLine) ;
            t_periCH3( 7) = -1;   

            % We add zeros left-right in order to center the segmentation array  
            % and all cells are aligned around midpoint. Emptpy values are = -1        
            Max_length_Line = max( WHISIT_P.AIS.Len_PerimSig);

            Z = Max_length_Line - length(aLine);

            if mod(Z,2) == 0
                plus = ones(1, (Z/2))*(-1);
                aLine = [plus, aLine , plus];
            else
                plus = ones(1, fix(Z/2))*(-1);      % fix(), division returning integer
                aLine = [plus, aLine , plus, -1];
            end
            % Write in cell array the profile line of the foci
            for pp = 1 : length(aLine)
                t_periCH3( 8+pp ) = aLine(pp) ;
            end

            % --- SAVE file_C as single entry-row      
            fprintf( file_PerimeterSigCH3 , '%f\t', t_periCH3(:) );
            fprintf( file_PerimeterSigCH3 , '\n' );
        end  
    
    end %end Profile_Line
  

    
    
    
%% --- Analysis of Spots --------------------------------------------------      
    % --- !!! ONLY POSSIBLE FOR CHANNEL 3, when used as Pole MARKER !!! ---
    if SaveSpotDataCH3 == 1

        for pp = 1 : length(cData.meshData{ff}{cc}.spots.l)
            t_Spt( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
            t_Spt( 2) =  cData.meshData{ff}{cc}.info.Oufti_cellID ;
            t_Spt( 3) =  cData.meshData{ff}{cc}.geom.length ;

            % Pixel Area Membrane and Cytosol and Cell
            t_Spt( 4) =  sum(sum( cData.meshData{ff}{cc}.Mask.Cell_body ));

            t_Spt( 5) = -1;

            t_Spt( 6) = length(cData.meshData{ff}{cc}.spots.l) ;        % number of spots
            t_Spt( 7) = cData.meshData{ff}{cc}.spots.l(pp) ;            % spot axial position
            t_Spt( 8) = cData.meshData{ff}{cc}.spots.d(pp) ;            % spot distance from cell's axes
            t_Spt( 9) = cData.meshData{ff}{cc}.spots.l_ratio(pp) ;      % spot axial position as ratio of cell length
            t_Spt(10) = cData.meshData{ff}{cc}.spots.x(pp) ;
            t_Spt(11) = cData.meshData{ff}{cc}.spots.y(pp) ;
                        
            % --- Save a single entry
            fprintf( file_Spt , '%f\t', t_Spt(:) );
            fprintf( file_Spt , '\n' );
        end
    end

end
end
end


% Close handles to all .txt files created
fclose(file_C) ;

if WHISIT_P.AIS.eval_AxProfile == 1 
    fclose(file_AxSigCH1) ;
    if WHISIT_P.eval_Channel_2 == 1 
        fclose(file_AxSigCH2) ;
    end
    if WHISIT_P.eval_Channel_3 == 1  &&  WHISIT_P.MarkedPole ~= 1 
        fclose(file_AxSigCH3) ;
    end
end

if WHISIT_P.AIS.eval_Segment == 1 
    fclose(file_SegmentSigCH1) ;
    if WHISIT_P.eval_Channel_2 == 1 
        fclose(file_SegmentSigCH2) ;
    end
    if WHISIT_P.eval_Channel_3 == 1  &&  WHISIT_P.MarkedPole ~= 1 
        fclose(file_SegmentSigCH3) ;
    end
end

% --- !!! ONLY POSSIBLE FOR CHANNEL 3, when used as Pole MARKER !!! ---
if SaveSpotDataCH3 == 1 
    fclose(file_Spt) ;
end


end






