function  Save_txt_PL( cData )
%Save_txt_PL = Save data analysed via algorithm PL (Polar Signal) in a
% tab separated .txt file.
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
fprintf( file_C, 'Frame \tCell Num \tCell Length \tCell_Area \t' );
fprintf( file_C, '      \tCH1_Cyto_AvgSig \tCH1_Cyto_Area \tCH1_Cyto_Sig_sum \t' );
fprintf( file_C, 'CH1_P1_AvgSig \tCH1_P1_Area \tCH1_P1_Sig_sum \tCH1_P1-Cyto_ratio \t' );
fprintf( file_C, 'CH1_P2_AvgSig \tCH1_P2_Area \tCH1_P2_Sig_sum \tCH1_P2-Cyto_ratio \t' );
fprintf( file_C, '      \tCH1_Avg_BkGr_noise \tCH1_Min_px_value \tCH1_Max_px_value \t' );

fprintf( file_C, '\tCH2_Cyto_AvgSig \tCH2_Cyto_Area \tCH2_Cyto_Sig_sum \t' );
fprintf( file_C, 'CH2_P1_AvgSig \tCH2_P1_Area \tCH2_P1_Sig_sum \tCH2_P1-Cyto_ratio \t' );
fprintf( file_C, 'CH2_P2_AvgSig \tCH2_P2_Area \tCH2_P2_Sig_sum \tCH2_P2-Cyto_ratio \t' );
fprintf( file_C, '      \tCH2_Avg_BkGr_noise \tCH2_Min_px_value \tCH2_Max_px_value \t\n' );



%% --- Data Extraction ----------------------------------------------------

for ff = 1 : length(cData.meshData)               % go through all frames

for cc = 1 : length(cData.meshData{ff})           % go through each cell
t_C = [];               % temporary matrix to handle and store data 

if  ~isempty(cData.meshData{ff}{cc}) ...
        && ~isempty(cData.meshData{ff}{cc}.model) ...
        &&  size(cData.meshData{ff}{cc}.mesh, 2) == 4      % if it is not empty 
  
    t_C( 1) =  cData.meshData{ff}{cc}.info.Original_Frame ;
    t_C( 2) =  cData.meshData{ff}{cc}.info.Original_cellID ;
    t_C( 3) =  cData.meshData{ff}{cc}.geom.length ;       % [px]
    t_C( 4) =  cData.meshData{ff}{cc}.geom.area ;         % [px]
    
    % --- CHANNEL 1 ---------1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1- %  

    t_C( 5) = -1;        % --- SKIP the column

    % Here we save for Cytosol, Pole 1 and Pole 2 the following:
    %---Average signal     
    %---Areas
    %---Signal sum (all-pixels in are)
    %---Ratio of signal PL/Cyto
    t_C( 6) =  mean(cData.meshData{ff}{cc}.CH1.IC( cData.meshData{ff}{cc}.CH1.Mask_pCyto )) ;
    t_C( 7) =  sum(sum( cData.meshData{ff}{cc}.CH1.Mask_pCyto ));    
    t_C( 8) =  sum(sum( cData.meshData{ff}{cc}.CH1.IC .* cData.meshData{ff}{cc}.CH1.Mask_pCyto ));
    
    avg_PL_1 = mean(cData.meshData{ff}{cc}.CH1.IC( cData.meshData{ff}{cc}.CH1.Mask_PL_1  ));
    avg_PL_2 = mean(cData.meshData{ff}{cc}.CH1.IC( cData.meshData{ff}{cc}.CH1.Mask_PL_2  ));
    if avg_PL_1 > avg_PL_2    
        t_C( 9 ) =  avg_PL_1;    
        t_C( 10) =  sum(sum( cData.meshData{ff}{cc}.CH1.Mask_PL_1 ));    
        t_C( 11) =  sum(sum( cData.meshData{ff}{cc}.CH1.IC .* cData.meshData{ff}{cc}.CH1.Mask_PL_1 ));
        t_C( 12) =  avg_PL_1 / mean(cData.meshData{ff}{cc}.CH1.IC( cData.meshData{ff}{cc}.CH1.Mask_pCyto ));

        t_C( 13) =  avg_PL_2;    
        t_C( 14) =  sum(sum( cData.meshData{ff}{cc}.CH1.Mask_PL_2 ));    
        t_C( 15) =  sum(sum( cData.meshData{ff}{cc}.CH1.IC .* cData.meshData{ff}{cc}.CH1.Mask_PL_2 ));
        t_C( 16) =  avg_PL_2 / mean(cData.meshData{ff}{cc}.CH1.IC( cData.meshData{ff}{cc}.CH1.Mask_pCyto ));

    elseif avg_PL_2 >= avg_PL_1    
        t_C( 9 ) =  avg_PL_2;    
        t_C( 10) =  sum(sum( cData.meshData{ff}{cc}.CH1.Mask_PL_2 ));    
        t_C( 11) =  sum(sum( cData.meshData{ff}{cc}.CH1.IC .* cData.meshData{ff}{cc}.CH1.Mask_PL_2 ));
        t_C( 12) =  avg_PL_2 / mean(cData.meshData{ff}{cc}.CH1.IC( cData.meshData{ff}{cc}.CH1.Mask_pCyto ));

        t_C( 13) =  avg_PL_1;    
        t_C( 14) =  sum(sum( cData.meshData{ff}{cc}.CH1.Mask_PL_1 ));    
        t_C( 15) =  sum(sum( cData.meshData{ff}{cc}.CH1.IC .* cData.meshData{ff}{cc}.CH1.Mask_PL_1 ));
        t_C( 16) =  avg_PL_1 / mean(cData.meshData{ff}{cc}.CH1.IC( cData.meshData{ff}{cc}.CH1.Mask_pCyto ));
    end
    
    t_C( 17) = -1;        % --- SKIP the column
    
    % Average background noise value in Frame
    t_C( 18)  =  cData.meshData{ff}{cc}.CH1.BkGr_Fr ;
    % Min and Max background pixel value in Frame
    t_C( 19) =  cData.meshData{ff}{cc}.CH1.Min_px_Fr ;
    t_C( 20) =  cData.meshData{ff}{cc}.CH1.Max_px_Fr ; 
    
    
    % --- CHANNEL 2 -----------2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2-2- %  
    if APP_opt.t1_choose_Chan == 2         % if we analyse two channels        
        t_C( 21) = -1;        % --- SKIP the column

        % As above, here we save for Cytosol, Pole 1 and Pole 2 the following:
        %---Average signal     
        %---Areas
        %---Signal sum (all-pixels in are)
        %---Ratio of signal PL/Cyto
        t_C( 22) =  mean(cData.meshData{ff}{cc}.CH2.IC( cData.meshData{ff}{cc}.CH2.Mask_pCyto )) ;
        t_C( 23) =  sum(sum( cData.meshData{ff}{cc}.CH2.Mask_pCyto ));    
        t_C( 24) =  sum(sum( cData.meshData{ff}{cc}.CH2.IC .* cData.meshData{ff}{cc}.CH2.Mask_pCyto ));

        avg_PL_1 = mean(cData.meshData{ff}{cc}.CH2.IC( cData.meshData{ff}{cc}.CH2.Mask_PL_1  ));
        avg_PL_2 = mean(cData.meshData{ff}{cc}.CH2.IC( cData.meshData{ff}{cc}.CH2.Mask_PL_2  ));
        if avg_PL_1 > avg_PL_2    
            t_C( 25) =  avg_PL_1;    
            t_C( 26) =  sum(sum( cData.meshData{ff}{cc}.CH2.Mask_PL_1 ));    
            t_C( 27) =  sum(sum( cData.meshData{ff}{cc}.CH2.IC .* cData.meshData{ff}{cc}.CH2.Mask_PL_1 ));
            t_C( 28) =  avg_PL_1 / mean(cData.meshData{ff}{cc}.CH2.IC( cData.meshData{ff}{cc}.CH2.Mask_pCyto ));

            t_C( 29) =  avg_PL_2;    
            t_C( 30) =  sum(sum( cData.meshData{ff}{cc}.CH2.Mask_PL_2 ));    
            t_C( 31) =  sum(sum( cData.meshData{ff}{cc}.CH2.IC .* cData.meshData{ff}{cc}.CH2.Mask_PL_2 ));
            t_C( 32) =  avg_PL_2 / mean(cData.meshData{ff}{cc}.CH2.IC( cData.meshData{ff}{cc}.CH2.Mask_pCyto ));

        elseif avg_PL_2 >= avg_PL_1    
            t_C( 25) =  avg_PL_2;    
            t_C( 26) =  sum(sum( cData.meshData{ff}{cc}.CH2.Mask_PL_2 ));    
            t_C( 27) =  sum(sum( cData.meshData{ff}{cc}.CH2.IC .* cData.meshData{ff}{cc}.CH2.Mask_PL_2 ));
            t_C( 28) =  avg_PL_2 / mean(cData.meshData{ff}{cc}.CH2.IC( cData.meshData{ff}{cc}.CH2.Mask_pCyto ));

            t_C( 29) =  avg_PL_1;    
            t_C( 30) =  sum(sum( cData.meshData{ff}{cc}.CH2.Mask_PL_1 ));    
            t_C( 31) =  sum(sum( cData.meshData{ff}{cc}.CH2.IC .* cData.meshData{ff}{cc}.CH2.Mask_PL_1 ));
            t_C( 32) =  avg_PL_1 / mean(cData.meshData{ff}{cc}.CH2.IC( cData.meshData{ff}{cc}.CH2.Mask_pCyto ));
        end

        t_C( 33) = -1;        % --- SKIP the column

        % Average background noise value in Frame
        t_C( 34)  =  cData.meshData{ff}{cc}.CH2.BkGr_Fr ;
        % Min and Max background pixel value in Frame
        t_C( 35) =  cData.meshData{ff}{cc}.CH2.Min_px_Fr ;
        t_C( 36) =  cData.meshData{ff}{cc}.CH2.Max_px_Fr ; 
        
    else
        t_C(21:36) = -1;        % --- SKIP the column
        
    end
    
    % Save a single Cell entry
    fprintf( file_C , '%f\t', t_C(:) );
    fprintf( file_C , '\n' );
    
end % if cData not empty
end % cc
end % ff

fclose(file_C) ;

end






