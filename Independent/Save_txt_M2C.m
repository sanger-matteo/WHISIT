function  Save_txt_M2C( cData )
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
% OUTPUT ------------------------------------------------------------------
% Cell_Data.txt = tab separated .txt file, where each row is the data
%                 relative to a specific cell
% N.B.: Empty position are filled with value -1. This allow to keep have
%       only numeric results which allows to easily export/import and
%       handling of the results (i.e. Excell). 
%       This is possible because cell length, area, pixel value and all
%       other results generated during analysis can never be negative
%       numbers ( although some can have value 0.0 ).
%
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|- 


%% --- Initialize ---------------------------------------------------------

global APP_opt;                        % Variable storing WHISIT options

% Select appropriate filename for Cell_Data.txt
if isempty(APP_opt.t1_exp_name)
    filename_cell_txt = [ APP_opt.t1_path_Det ,'/', 'Cell_Data.txt'];
else
    filename_cell_txt = [ APP_opt.t1_path_Det ,'/' , APP_opt.t1_exp_name , '_', 'Cell_Data.txt'];
end

% Initialize first row with labels for each column necessary
file_C = fopen(filename_cell_txt, 'w+');
fprintf( file_C, 'Frame \tCell Num \tCell Length \tArea Membr \tArea Cyto \tArea Cell \t' );

fprintf( file_C, '     \tCH1_Avg Int Membr \tCH1_Avg Int Cyto \tCH1_Avg Int Cell \t' );
fprintf( file_C, 'CH1 Avg_BkGr_noise \tCH1_Min_px_value_Frame \tCH1_Max_px_value_Frame \t' );

fprintf( file_C, '     \tCH2_Avg Int Membr \tCH2_Avg Int Cyto \tCH2_Avg Int Cell \t' );
fprintf( file_C, 'CH2 Avg_BkGr_noise \tCH2_Min_px_value_Frame \tCH2_Max_px_value_Frame \t\n' );



%% --- Data Extraction ----------------------------------------------------
for ff = 1 : length(cData.meshData)               % go through all frames
for cc = 1 : length(cData.meshData{ff})           % go through each cell    
t_C = [];               % temporary matricx to handle and store data 

if  ~isempty(cData.meshData{ff}{cc}) ...
        && ~isempty(cData.meshData{ff}{cc}.model) ...
        &&  size(cData.meshData{ff}{cc}.mesh, 2) == 4      % if it is not empty 
    
    t_C( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
    t_C( 2) =  cData.meshData{ff}{cc}.info.Original_cellID ;
    t_C( 3) =  cData.meshData{ff}{cc}.geom.length ;
        
    % Pixel Area Membrane and Cytosol and Cell
    t_C( 4) =  sum(sum( cData.meshData{ff}{cc}.Mask_Memb ));
    t_C( 5) =  sum(sum( cData.meshData{ff}{cc}.Mask_mCyto ));
    t_C( 6) =  sum(sum( cData.meshData{ff}{cc}.Mask_wCell )); 
    
    % --- CHANNEL 1 ---------1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1- % 
    
    t_C( 7) = -1;        % --- SKIP the column

    % Average Pixel signal of Membrane and Cytosol and Cell: Absolute
    t_C( 8) =  mean(cData.meshData{ff}{cc}.CH1.IC( cData.meshData{ff}{cc}.Mask_Memb ));
    t_C( 9) =  mean(cData.meshData{ff}{cc}.CH1.IC( cData.meshData{ff}{cc}.Mask_mCyto ));
    t_C( 10) =  mean(cData.meshData{ff}{cc}.CH1.IC( cData.meshData{ff}{cc}.Mask_wCell ));
    
    % Average background noise value in Frame
    t_C( 11)  =  cData.meshData{ff}{cc}.CH1.BkGr_Fr ;
    % Min and Max background pixel value in Frame
    t_C( 12) =  cData.meshData{ff}{cc}.CH1.Min_px_Fr ;
    t_C( 13) =  cData.meshData{ff}{cc}.CH1.Max_px_Fr ;
    
   
    % --- CHANNEL 2 -----------2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2- %  
    if APP_opt.t1_choose_Chan == 2          % if we analyse two channels
        t_C( 14) = -1;        % --- SKIP the column

        % Average Pixel signal of Membrane and Cytosol and Cell: Absolute
        t_C( 15) =  mean(cData.meshData{ff}{cc}.CH2.IC( cData.meshData{ff}{cc}.Mask_Memb ));
        t_C( 16) =  mean(cData.meshData{ff}{cc}.CH2.IC( cData.meshData{ff}{cc}.Mask_mCyto ));
        t_C( 17) =  mean(cData.meshData{ff}{cc}.CH2.IC( cData.meshData{ff}{cc}.Mask_wCell ));

        % Average background noise value in Frame
        t_C( 18)  =  cData.meshData{ff}{cc}.CH2.BkGr_Fr ;
        % Min and Max background pixel value in Frame
        t_C( 19) =  cData.meshData{ff}{cc}.CH2.Min_px_Fr ;
        t_C( 20) =  cData.meshData{ff}{cc}.CH2.Max_px_Fr ;
        
    else
        % --- SKIP all this columns 
        t_C( 14 : 20) = -1;
    end % Channel 2

    
    % Save a single Cell entry
    fprintf( file_C , '%f\t', t_C(:) );
    fprintf( file_C , '\n' );
    
end % if cData not empty
end % cc
end % ff

fclose(file_C) ;


end






