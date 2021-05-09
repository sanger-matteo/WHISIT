function  Save_txt_PL( cData, WHISIT_P, Path_Save, Exp_Name )
%Save_txt_PL = Save data analysed via algorithm PL (Polar Signal) in a
% tab separated .txt file.
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
% Path_Save = path to where to store the .txt file
%
% Exp_Name = if provided, assign a costum name Prefix to the .txt file
%
%
% OUTPUT ------------------------------------------------------------------
% Cell_Data.txt = tab separated .txt file, where each row is the data
%                 relative to a specific cell
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

% Select appropriate filename for Cell_Data.txt
if isempty(Exp_Name)
    filename_cell_txt = [ Path_Save ,'/', 'Cell_Data.txt'];
else
    filename_cell_txt = [ Path_Save ,'/' , Exp_Name , '_', 'Cell_Data.txt'];
end

% Initialize first row with labels for each column necessary
file_C = fopen(filename_cell_txt, 'w+');
fprintf( file_C, 'Frame \tCell Num \tCell Length \tCell_Area \t' );
fprintf( file_C, '      \tCH1_Cyto_AvgSig \tCH1_Cyto_Area \tCH1_Cyto_Sig_sum \t' );
fprintf( file_C, 'CH1_P1_AvgSig \tCH1_P1_Area \tCH1_P1_Sig_sum \tCH1_P1-Cyto_ratio \t' );
fprintf( file_C, 'CH1_P2_AvgSig \tCH1_P2_Area \tCH1_P2_Sig_sum \tCH1_P2-Cyto_ratio \t' );
fprintf( file_C, '      \tCH1_Avg_BkGr_noise \tCH1_Min_px_value \tCH1_Max_px_value \t' );

fprintf( file_C, '\tCH2_Cyto_AvgSig \tCH2_Cyto_Area \tCH2_Cyto_Sig_sum \t' );
fprintf( file_C, 'CH2_P1_AvgSig \tCH2_P1_Area \tCH2_P1_Sig_sum \tCH2_P1-Cyto_ratio \t' );
fprintf( file_C, 'CH2_P2_AvgSig \tCH2_P2_Area \tCH2_P2_Sig_sum \tCH2_P2-Cyto_ratio \t' );
fprintf( file_C, '      \tCH2_Avg_BkGr_noise \tCH2_Min_px_value \tCH2_Max_px_value \t' );

fprintf( file_C, '\tCH3_Cyto_AvgSig \tCH3_Cyto_Area \tCH3_Cyto_Sig_sum \t' );
fprintf( file_C, 'CH3_P1_AvgSig \tCH3_P1_Area \tCH3_P1_Sig_sum \tCH3_P1-Cyto_ratio \t' );
fprintf( file_C, 'CH3_P2_AvgSig \tCH3_P2_Area \tCH3_P2_Sig_sum \tCH3_P2-Cyto_ratio \t' );
fprintf( file_C, '      \tCH3_Avg_BkGr_noise \tCH3_Min_px_value \tCH3_Max_px_value \t' );

fprintf( file_C, '      \tSet_Polarity \t\n' );



%% --- Data Extraction ----------------------------------------------------

for ff = 1 : length(cData.meshData)               % go through all frames

for cc = 1 : length(cData.meshData{ff})           % go through each cell
t_C = [];               % temporary matrix to handle and store data 

if  ~isempty(cData.meshData{ff}{cc}) ...
        && ~isempty(cData.meshData{ff}{cc}.model) ...
        &&  size(cData.meshData{ff}{cc}.mesh, 2) == 4      % if it is not empty 
  
    t_C( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
    t_C( 2) =  cData.meshData{ff}{cc}.info.Oufti_cellID ;
    t_C( 3) =  cData.meshData{ff}{cc}.geom.length ;       % [px]
    t_C( 4) =  cData.meshData{ff}{cc}.geom.area ;         % [px]
    
    % --- CHANNEL 1 ---------1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1- %  

    t_C( 5) = -1;        % --- SKIP the column
    % Average intensity inside Cytosol area
    t_C( 6) =  mean(cData.meshData{ff}{cc}.Fluor_Chan(1).IC( cData.meshData{ff}{cc}.Fluor_Chan(1).Mask.Cytosol )) ;
    % Area of Cytosol Mask
    t_C( 7) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(1).Mask.Cytosol ));   
    % Signal sum of Cytosol Mask
    t_C( 8) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(1).IC .* cData.meshData{ff}{cc}.Fluor_Chan(1).Mask.Cytosol ));

    % Average intensity, Area and Signal sum in Cytosol Mask
    avg_PL_1 = mean(cData.meshData{ff}{cc}.Fluor_Chan(1).IC( cData.meshData{ff}{cc}.Fluor_Chan(1).Mask.Pole_1  ));
    avg_PL_2 = mean(cData.meshData{ff}{cc}.Fluor_Chan(1).IC( cData.meshData{ff}{cc}.Fluor_Chan(1).Mask.Pole_2  ));
   
    t_C( 9 ) =  avg_PL_1;    
    t_C( 10) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(1).Mask.Pole_1 ));    
    t_C( 11) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(1).IC .* cData.meshData{ff}{cc}.Fluor_Chan(1).Mask.Pole_1 ));
    t_C( 12) =  avg_PL_1 / mean(cData.meshData{ff}{cc}.Fluor_Chan(1).IC( cData.meshData{ff}{cc}.Fluor_Chan(1).Mask.Cytosol ));
    
    t_C( 13) =  avg_PL_2;    
    t_C( 14) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(1).Mask.Pole_2 ));    
    t_C( 15) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(1).IC .* cData.meshData{ff}{cc}.Fluor_Chan(1).Mask.Pole_2 ));
    t_C( 16) =  avg_PL_2 / mean(cData.meshData{ff}{cc}.Fluor_Chan(1).IC( cData.meshData{ff}{cc}.Fluor_Chan(1).Mask.Cytosol ));

    t_C( 17) = -1;        % --- SKIP the column
    
    % Average background noise value in Frame
    t_C( 18)  =  cData.meshData{ff}{cc}.Fluor_Chan(1).Avg_BkGr_value ;
    % Min and Max background pixel value in Frame
    t_C( 19) =  cData.meshData{ff}{cc}.Fluor_Chan(1).Min_px_value;
    t_C( 20) =  cData.meshData{ff}{cc}.Fluor_Chan(1).Max_px_value; 
    
    
    % --- CHANNEL 2 -----------2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2- %  
    if WHISIT_P.eval_Channel_2 == 1       
        t_C( 21) = -1;        % --- SKIP the column

        % Average intensity inside Cytosol area
        t_C( 22) =  mean(cData.meshData{ff}{cc}.Fluor_Chan(2).IC( cData.meshData{ff}{cc}.Fluor_Chan(2).Mask.Cytosol )) ;
        % Area of Cytosol Mask
        t_C( 23) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(2).Mask.Cytosol ));   
        % Signal sum of Cytosol Mask
        t_C( 24) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(2).IC .* cData.meshData{ff}{cc}.Fluor_Chan(2).Mask.Cytosol ));

        % Average intensity, Area and Signal sum in Cytosol Mask
        avg_PL_1 = mean(cData.meshData{ff}{cc}.Fluor_Chan(2).IC( cData.meshData{ff}{cc}.Fluor_Chan(2).Mask.Pole_1  ));
        avg_PL_2 = mean(cData.meshData{ff}{cc}.Fluor_Chan(2).IC( cData.meshData{ff}{cc}.Fluor_Chan(2).Mask.Pole_2  ));

        t_C( 25) =  avg_PL_1;    
        t_C( 26) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(2).Mask.Pole_1 ));    
        t_C( 27) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(2).IC .* cData.meshData{ff}{cc}.Fluor_Chan(2).Mask.Pole_1 ));
        t_C( 28) =  avg_PL_1 / mean(cData.meshData{ff}{cc}.Fluor_Chan(2).IC( cData.meshData{ff}{cc}.Fluor_Chan(2).Mask.Cytosol ));
        
        t_C( 29) =  avg_PL_2; 
        t_C( 30) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(2).Mask.Pole_2 ));    
        t_C( 31) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(2).IC .* cData.meshData{ff}{cc}.Fluor_Chan(2).Mask.Pole_2 ));
        t_C( 32) =  avg_PL_2 / mean(cData.meshData{ff}{cc}.Fluor_Chan(2).IC( cData.meshData{ff}{cc}.Fluor_Chan(2).Mask.Cytosol ));

        t_C( 33) = -1;        % --- SKIP the column

        % Average background noise value in Frame
        t_C( 34)  =  cData.meshData{ff}{cc}.Fluor_Chan(2).Avg_BkGr_value ;
        % Min and Max background pixel value in Frame
        t_C( 35) =  cData.meshData{ff}{cc}.Fluor_Chan(2).Min_px_value;
        t_C( 36) =  cData.meshData{ff}{cc}.Fluor_Chan(2).Max_px_value;

    else
        t_C(21:36) = -1;        % --- SKIP the column
        
    end % end_if Channel 2


    % --- CHANNEL 3 -----------3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3- %
    if WHISIT_P.eval_Channel_3 == 1  &&  WHISIT_P.MarkedPole ~= 1
        t_C( 37) = -1;        % --- SKIP the column

        % Average intensity inside Cytosol area
        t_C( 38) =  mean(cData.meshData{ff}{cc}.Fluor_Chan(3).IC( cData.meshData{ff}{cc}.Fluor_Chan(3).Mask.Cytosol )) ;
        % Area of Cytosol Mask
        t_C( 39) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(3).Mask.Cytosol ));   
        % Signal sum of Cytosol Mask
        t_C( 40) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(3).IC .* cData.meshData{ff}{cc}.Fluor_Chan(3).Mask.Cytosol ));

        % Average intensity, Area and Signal sum in Cytosol Mask
        avg_PL_1 = mean(cData.meshData{ff}{cc}.Fluor_Chan(3).IC( cData.meshData{ff}{cc}.Fluor_Chan(3).Mask.Pole_1  ));
        avg_PL_2 = mean(cData.meshData{ff}{cc}.Fluor_Chan(3).IC( cData.meshData{ff}{cc}.Fluor_Chan(3).Mask.Pole_2  ));

        t_C( 41) =  avg_PL_1;    
        t_C( 42) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(3).Mask.Pole_1 ));    
        t_C( 43) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(3).IC .* cData.meshData{ff}{cc}.Fluor_Chan(3).Mask.Pole_1 ));
        t_C( 44) =  avg_PL_1 / mean(cData.meshData{ff}{cc}.Fluor_Chan(3).IC( cData.meshData{ff}{cc}.Fluor_Chan(3).Mask.Cytosol ));
       
        t_C( 45) =  avg_PL_2; 
        t_C( 46) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(3).Mask.Pole_2 ));    
        t_C( 47) =  sum(sum( cData.meshData{ff}{cc}.Fluor_Chan(3).IC .* cData.meshData{ff}{cc}.Fluor_Chan(3).Mask.Pole_2 ));
        t_C( 48) =  avg_PL_2 / mean(cData.meshData{ff}{cc}.Fluor_Chan(3).IC( cData.meshData{ff}{cc}.Fluor_Chan(3).Mask.Cytosol ));

        t_C( 49) = -1;        % --- SKIP the column

        % Average background noise value in Frame
        t_C( 50)  =  cData.meshData{ff}{cc}.Fluor_Chan(3).Avg_BkGr_value ;
        % Min and Max background pixel value in Frame
        t_C( 51) =  cData.meshData{ff}{cc}.Fluor_Chan(3).Min_px_value;
        t_C( 52) =  cData.meshData{ff}{cc}.Fluor_Chan(3).Max_px_value;

    else
        t_C(37:52) = -1;        % --- SKIP the column
        
    end % end_if Channel 3

   
    % Show polarity of the cell
    t_C(53) = -1;
    t_C(54) = cData.meshData{ff}{cc}.polarity ;
    

    % --- SAVE file_C as single entry-row  
    fprintf( file_C , '%f\t', t_C(:) );
    fprintf( file_C , '\n' );
    
end % if cData not empty
end % cc
end % ff

fclose(file_C) ;    % Close file handle

end






