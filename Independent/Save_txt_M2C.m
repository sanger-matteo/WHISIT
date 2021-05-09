function  Save_txt_M2C( cData , WHISIT_P, Path_Save, Exp_Name )
%Save_txt_M2C = Save data analysed via algorithm M2C (Membrane 2 Cytosol) 
% in a tab separated .txt file.
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
% OUTPUT ------------------------------------------------------------------
% Cell_Data.txt = tab separated .txt file, where each row is the data
%                 relative to a specific cell
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
fprintf( file_C, 'Frame \tCell_Num \tCell_Length \tArea_Membr \tArea_Cyto \tArea_Cell \t' );

fprintf( file_C, '     \tCH1_Avg_Int Membr \tCH1_Avg_Int_Cyto \tCH1_Avg_Int_Cell \t' );
fprintf( file_C, 'CH1_Avg_BkGr_noise \tCH1_Min_px_value_Frame \tCH1_Max_px_value_Frame \t' );

fprintf( file_C, '     \tCH2_Avg_Int_Membr \tCH2_Avg_Int_Cyto \tCH2_Avg_Int_Cell \t' );
fprintf( file_C, 'CH2_Avg_BkGr_noise \tCH2_Min_px_value_Frame \tCH2_Max_px_value_Frame \t' );

fprintf( file_C, '     \tCH3_Avg_Int_Membr \tCH3_Avg_Int_Cyto \tCH3_Avg_Int_Cell \t' );
fprintf( file_C, 'CH3_Avg_BkGr_noise \tCH3_Min_px_value_Frame \tCH3_Max_px_value_Frame \t\n' );



%% --- Data Extraction ----------------------------------------------------
for ff = 1 : length(cData.meshData)               % go through all frames
for cc = 1 : length(cData.meshData{ff})           % go through each cell    
t_C = [];               % temporary matricx to handle and store data 

if  ~isempty(cData.meshData{ff}{cc}) ...
        && ~isempty(cData.meshData{ff}{cc}.model) ...
        &&  size(cData.meshData{ff}{cc}.mesh, 2) == 4      % if it is not empty 
    
    t_C( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
    t_C( 2) =  cData.meshData{ff}{cc}.info.Oufti_cellID ;
    t_C( 3) =  cData.meshData{ff}{cc}.geom.length ;
        
    % Pixel Area Membrane and Cytosol and Cell
    t_C( 4) =  sum(sum( cData.meshData{ff}{cc}.Mask.Memb_All ));
    t_C( 5) =  sum(sum( cData.meshData{ff}{cc}.Mask.Cytosol ));
    t_C( 6) =  sum(sum( cData.meshData{ff}{cc}.Mask.Cell_body )); 
  
  
    % --- CHANNEL 1 ---------1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1- % 
    
    t_C( 7) = -1;        % --- SKIP the column
    % Average Pixel signal of Membrane and Cytosol and Cell: Absolute
    t_C( 8)  =  mean(cData.meshData{ff}{cc}.Fluor_Chan(1).IC( cData.meshData{ff}{cc}.Mask.Memb_All ));
    t_C( 9)  =  mean(cData.meshData{ff}{cc}.Fluor_Chan(1).IC( cData.meshData{ff}{cc}.Mask.Cytosol ));
    t_C( 10) =  mean(cData.meshData{ff}{cc}.Fluor_Chan(1).IC( cData.meshData{ff}{cc}.Mask.Cell_body ));
    
    % Average background noise value in Frame
    t_C( 11) =  cData.meshData{ff}{cc}.Fluor_Chan(1).Avg_BkGr_value ;
    % Min and Max background pixel value in Frame
    t_C( 12) =  cData.meshData{ff}{cc}.Fluor_Chan(1).Min_px_value;
    t_C( 13) =  cData.meshData{ff}{cc}.Fluor_Chan(1).Max_px_value;
    
   
    % --- CHANNEL 2 -----------2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2- %  
    if WHISIT_P.eval_Channel_2 == 1   

        t_C( 14) = -1;        % --- SKIP the column
        % Average Pixel signal of Membrane and Cytosol and Cell: Absolute
        t_C( 15) =  mean(cData.meshData{ff}{cc}.Fluor_Chan(2).IC( cData.meshData{ff}{cc}.Mask.Memb_All ));
        t_C( 16) =  mean(cData.meshData{ff}{cc}.Fluor_Chan(2).IC( cData.meshData{ff}{cc}.Mask.Cytosol ));
        t_C( 17) =  mean(cData.meshData{ff}{cc}.Fluor_Chan(2).IC( cData.meshData{ff}{cc}.Mask.Cell_body ));

        % Average background noise value in Frame
        t_C( 18) =  cData.meshData{ff}{cc}.Fluor_Chan(2).Avg_BkGr_value ;
        % Min and Max background pixel value in Frame
        t_C( 19) =  cData.meshData{ff}{cc}.Fluor_Chan(2).Min_px_value;
        t_C( 20) =  cData.meshData{ff}{cc}.Fluor_Chan(2).Max_px_value;
        
    else
        % --- SKIP all this columns 
        t_C( 14 : 20) = -1;

    end % end_if Channel 2


    % --- CHANNEL 3 -----------3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3-3- %  
    if WHISIT_P.eval_Channel_3 == 1  &&  WHISIT_P.MarkedPole ~= 1

        t_C( 21) = -1;        % --- SKIP the column
        % Average Pixel signal of Membrane and Cytosol and Cell: Absolute
        t_C( 22) =  mean(cData.meshData{ff}{cc}.Fluor_Chan(3).IC( cData.meshData{ff}{cc}.Mask.Memb_All ));
        t_C( 23) =  mean(cData.meshData{ff}{cc}.Fluor_Chan(3).IC( cData.meshData{ff}{cc}.Mask.Cytosol ));
        t_C( 24) =  mean(cData.meshData{ff}{cc}.Fluor_Chan(3).IC( cData.meshData{ff}{cc}.Mask.Cell_body ));

        % Average background noise value in Frame
        t_C( 25)  =  cData.meshData{ff}{cc}.Fluor_Chan(3).Avg_BkGr_value ;
        % Min and Max background pixel value in Frame
        t_C( 26) =  cData.meshData{ff}{cc}.Fluor_Chan(3).Min_px_value;
        t_C( 27) =  cData.meshData{ff}{cc}.Fluor_Chan(3).Max_px_value;
        
    else
        % --- SKIP all this columns 
        t_C( 21 : 27) = -1;

    end % end_if Channel 3

    
    % --- SAVE file_C as single entry-row  
    fprintf( file_C , '%f\t', t_C(:) );
    fprintf( file_C , '\n' );
    
end % if cData not empty
end % cc
end % ff

fclose(file_C) ;    % Close file handle

end






