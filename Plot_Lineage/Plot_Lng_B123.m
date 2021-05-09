function Plot_Lng_B123
% Here we Plot the whole genalogy tree for a single lineage
%
%
%
%
%--------------------------------------------------------------------------
% cN_End_Line and cS_End_Line allow us to to decide in which order to plot
% the clones plotting order. Starting with the end-of-line cells
% (cN_End_Line ) we plot every coupled daughter cells and then find the
% update the list with the common ancestor. Doing this and progressing
% "back in time" we complete the lineage tree by exclusion, untill only the
% founder cell is 
% Those two list will keep changing as we procede in the tree. We also  keep  
% using ID_Raw_array as its c_ID have same index as in clone_List, making 
% search more easy.
%--------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

% Rembember:
% -> OLD_pole = Mask.Pole_1
% -> NEW_pole = Mask.Pole_2

  
%%
% ----- INITIALIZE --------------------------------------------------------
% Gather plotting parameters, the distribution of signal values and create
% a corresponding RGB colormap

global APP_opt ;

hg_pl = APP_opt.t3_height_PlotLine ;         % heights of the plot end lines
Dst_pl = APP_opt.t3_Dist_btw_PlotLines ;     % distances between plot end lines 
LinW =  APP_opt.t3_BorderLineWidth ;         % LineWidth of line connecting ancestor-to-daughters
Color_BorderLines = [.3 .3 .3]; 
LAST_Frm = 1;                                % highest value frame, use to find max x-axis length
fpmRate = 1 / APP_opt.t3_fpm_CH1 ;

% Load the clone_List
load([APP_opt.t3_path_cloneList, APP_opt.t3_file_cloneList]);   


% Consider option APP_opt.t3_choose_ChannelMode, what channel(s) we wish to
% plot: we can plot either individual channels, both in same plot or a
% ratio of the two.
% Each clone can be representated by one line containing the information
% for a single channel or a ratio of the two; or a clone can be
% representated by a line divided in two to plot the information for both
% independent channels.
% For this reason we use two variables to help discriminate all the cases:
% - ChNum   : carries the channel(s) we need to create the plot
% - SubLine : carries the number of "lines" to plot for each clone

if APP_opt.t3_choose_ChannelMode == 1
    ChNum = 1;
    SubLine = 1 ;
    LUT(1) = APP_opt.t3_ColorMap_LUT_CH1 ;
    vRange{1} = APP_opt.t3_Value_Range_CH1 ;
    
elseif APP_opt.t3_choose_ChannelMode == 2
    ChNum = 2 ;
    SubLine = 1 ;
    LUT(1) = APP_opt.t3_ColorMap_LUT_CH2 ;
    vRange{1} = APP_opt.t3_Value_Range_CH2 ;
    
elseif APP_opt.t3_choose_ChannelMode == 3
    ChNum = [1, 2] ;
    SubLine = [1, 2];
    LUT(1) = APP_opt.t3_ColorMap_LUT_CH1 ;
    LUT(2) = APP_opt.t3_ColorMap_LUT_CH2 ;
    vRange{1} = APP_opt.t3_Value_Range_CH1;
    vRange{2} = APP_opt.t3_Value_Range_CH2;
    
elseif APP_opt.t3_choose_ChannelMode == 4
    ChNum = [1, 2] ;
    SubLine = 1 ;
    LUT(1) = APP_opt.t3_ColorMap_LUT_CH1 ;
    vRange{1} = APP_opt.t3_Value_Range_CH1;
    
elseif APP_opt.t3_choose_ChannelMode == 5
    ChNum = [1, 2] ;
    SubLine = 1 ;
    LUT(1) = APP_opt.t3_ColorMap_LUT_CH1 ;
    vRange{1} = APP_opt.t3_Value_Range_CH1;
    
end

% Gather the distribution of all values (for each channel) needed to plot. 
% Note if APP_opt.t3_choose_ChannelMode >= 4, Hist_Values_cloneList
% return the ratio between the two channels.
for kk = ChNum
    [ Mu, Sigma, Distr_Vals ] = Hist_Values_cloneList(0, clone_List);
end
% For each channel create an RGB colormap, with corresponding range of
% values associated to each color.
for kk = SubLine
    % Load the color LUT. For ratios, we select LUT of channel 1.
    % [LUT folder should be in same folder of the scripts Plot_Lineage]
    pathParts_LUT = strsplit(mfilename('fullpath'), {'/','\'} ) ;        % mfilename = take path of currently running script.
    path_LUT = fullfile(pathParts_LUT{1,1:end-1},'\') ;                  % fullfile  = build full filename from string parts
    filename_LUT = [path_LUT '/LUT/' 'LUT_' LUT(kk) '.txt'] ;
    
    RGB_range{kk} = textread(filename_LUT) ;
    [ V_Min{kk}, V_Max{kk}, RGB_range{kk}] = Find_ColorValueRange( ...
                    Mu{kk}, Sigma{kk}, Distr_Vals{kk} , RGB_range{kk}, vRange{kk});                
end
                       
% Create ordered list of the clones (cs_All) and the end-of-line clone (cS_EndLine) 
% This will allow to plot in the correct order, from the "latest" to the founder cell
[cS_EndLine, cN_EndLine, cS_All] = fnc_Organize_cloneList(clone_List) ;


% --- SAVE_step_0 ------ Initialize -------------------------------------
if APP_opt.t3_choose_Save_txt == 1   
    % Create a subfolder to store .txt files, if it does not exsist    
    PathFolder = [APP_opt.t3_path_cloneList, '/Plotted_Values'];
    if exist(PathFolder) == 0
        mkdir(PathFolder) ;
    end    
    % Select appropriate file name for all Data.txt file to create
    if isempty(APP_opt.t3_exp_name)
        Path_txt = [ PathFolder ,'/'];
    else
        Path_txt = [ PathFolder ,'/', APP_opt.t3_exp_name, '_' ];
    end
    
    % Initialize variable to collect all values of data plotted
    % [cell array in which each element is a channel; for ratios we
    %  collect their values and not individual channels]
    for kk = SubLine       
        tree_Values{kk} = cell(4,length(clone_List));    
    end
    
    % Find the latest time point for entire clone list, necessary to save data
    % in .txt file by ensuring that each clone data has same length.
    last_timepoint = -1 ;
    for cc = 1 : size(clone_List,2)
        if last_timepoint <= clone_List{cc}{1}.frame_last
           last_timepoint = clone_List{cc}{1}.frame_last ;
        end
    end
end





if APP_opt.t3_PlotOpt_A == 1            % if GUI Plot Type: '1 - Lineage Tree'
%% 
% ***** PLOT LINEAGE TREE ************************************************************************

% Make figure visible of invisible according to user's choice    
if APP_opt.t3_choose_Display_Plot == 0
    h1f = figure('Position', [100 300 1200 600], 'Visible','off'); 
elseif APP_opt.t3_choose_Display_Plot == 1
    h1f = figure('Position', [100 300 1200 600]);    
end
hold on;
ax = gca; 
LAST_Frm = 1;                           % highest value frame, use to find max x-axis length
y_H_Row = [1:length(cN_EndLine)];   


%% --- STEP 1 --- Plot all End_Line clones ---------------------------------
% First we plot all end line clones, meaning all the clones that do not
% generate daughter cells at the end of their "life".

for cc = 1 : length(cN_EndLine)
    % Find where the position of the cc-th EndLine in the array
    idx_EndL = find(strcmp(cS_All , [APP_opt.t3_PlotOpt_NLineage cS_EndLine{cc}] ) == 1);
    % Check if there are two cells with identical ID name
    if length(idx_EndL)>2                           % if true, 
        fprintf('ERROR, two identical cell IDs');   % ERROR and return
        return
    end        

    % yrow = store the y position for plotting the current cc-th clone
    yrow = y_H_Row(cc); 
    
    % Plot all time points. We use two counters for "time":
    % - ff  = stores the array position of each element (/frame) of
    %         the clone_List{XXX}: simply goes from 1 to last frame
    % - Frm = x-axis "real" time variable, used to place the extracted
    %         Value at the correct time point of the lineage tree
    for ff = 1 : length(clone_List{idx_EndL})
        
        Frm = ff + clone_List{idx_EndL}{1}.frame_birth;
        if Frm > LAST_Frm;     LAST_Frm = Frm;     end

        % If choosen, plot HUD: cell poles orientation and cell_ID number
        if ff == 1      % at frame of birth;
           fnc_Plot_HUD( clone_List{idx_EndL}{1}, Frm, yrow, Dst_pl, hg_pl , LinW)
        end

        % Define the coordinates of the "box" square polygon and find the
        % signal value of the desired subcellular compartment to plot.
        [Xs, Ys, Value] = Find_Value_XYSquare( clone_List{idx_EndL}{ff}, Frm, yrow, Dst_pl, hg_pl, ChNum , SubLine, ...
                                               APP_opt.t3_choose_ChannelMode, APP_opt.t3_PlotOpt_B, APP_opt.t3_PlotOpt_C, APP_opt.algorithm );
        
        % Now, evaluate the RGB_Color for each "Line" that we wish to plot for each clone
        if APP_opt.t3_PlotOpt_B == 1  |  APP_opt.t3_PlotOpt_B == 3
            for kk = SubLine            
                cRGB  = Assign_Value_RGB( Value{kk} , RGB_range{kk} );
                % Fill squares with evaluated RGB_Color
                fill(Xs{kk}, Ys{kk}, [cRGB(1), cRGB(2), cRGB(3)], 'LineStyle', 'none');
            end
        elseif APP_opt.t3_PlotOpt_B == 2
            for kk = SubLine
                for ii = 1 : 2      % column 1 = Upper;   column 2 = Lower.
                    cRGB  = Assign_Value_RGB( Value{kk}(ii) , RGB_range{kk} );
                    % Fill squares with evaluated RGB_Color
                    fill(Xs, Ys{kk,ii}, [cRGB(1), cRGB(2), cRGB(3)], 'LineStyle', 'none');
                end
            end
        end
        
        % --- Extra Info ---- Place dots over cell_line to represent cell's
        % eccentricity value [greys scale: black = 0/circle; white = 1/ellipse)]
        if APP_opt.t3_display_ExI == 1
            eXs = [clone_List{idx_EndL}{ff}.mesh(:,1); flipud(clone_List{idx_EndL}{ff}.mesh(:,3))];
            eYs = [clone_List{idx_EndL}{ff}.mesh(:,2); flipud(clone_List{idx_EndL}{ff}.mesh(:,4))];       
            [la, sa] = Geom_My_fit_ellipse__v2( eXs, eYs );
            % Calculate ellipse eccentricity (sa = short semi-axis; la = long semi-axis)
            % Multipling by factor *0.9 restrict grey scale to a range 0:0.9 and 
            % it avoids having "white" marker (which are hardly visible)
            El_Ec = sqrt( 1 - (sa^2/la^2)) .*0.9 ;
            % The ratio cannot be bigger than 1
            if El_Ec > 0.9 ;      El_Ec = 0.9 ;     end
            El_Ec_Col = [El_Ec El_Ec El_Ec] ;            
            plot( Frm -0.5, (yrow*Dst_pl +hg_pl) +hg_pl/3 , '.', 'MarkerSize', 8, 'LineWidth', 0.2, 'Color', El_Ec_Col);
        end
                
        % --- SAVE_step_1 ---- Collect plotted value of all end-of-line clones
        if APP_opt.t3_choose_Save_txt == 1   
            for kk = SubLine
                % Save the value(s) plotted
                tree_Values{kk}{1, idx_EndL} = [tree_Values{kk}{1, idx_EndL} , Value{kk}'];     
                % [px] Cell axial-length       
                tree_Values{kk}{2, idx_EndL} = [tree_Values{kk}{2, idx_EndL}, clone_List{idx_EndL}{ff}.geom.length] ;    
                % [px] Cell area
                tree_Values{kk}{3, idx_EndL} = [tree_Values{kk}{3, idx_EndL}, clone_List{idx_EndL}{ff}.geom.area] ;
                % y-position in plot  
                tree_Values{kk}{4, idx_EndL} = yrow*Dst_pl ;
            end
        end

    end % ff

    % Draw two horizontal lines that bound the plotted data for the clone-Line
    plot( [Frm, Frm-length(clone_List{idx_EndL})], [yrow*Dst_pl,  yrow*Dst_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
    plot( [Frm, Frm-length(clone_List{idx_EndL})], [(yrow*Dst_pl)+hg_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
            
    % Draw vertical line at the end, to enclose all clone-Line
    plot( [Frm, Frm], [yrow*Dst_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);

    % If necessary draw a line to separate the two channels plot-lines
    if APP_opt.t3_choose_ChannelMode == 3
        plot( [Frm, Frm-length(clone_List{idx_EndL})], [(yrow*Dst_pl)+(hg_pl*1/2),  (yrow*Dst_pl)+(hg_pl*1/2)], 'LineWidth', LinW/2, 'Color', Color_BorderLines);
    end
    
end % cc



%% --- STEP 2 --- Plot all remaining clones ---------------------------------
% This step will plot all the remaining clones. The algorithm find clone
% couples using cS_EndLine: a couple is two clonse that share same ID root
% (.XYZ ). From this it is possible to identify the common ancestror and
% plot such cell (which would be .XY ). The cS_EndLine list is updated and
% the algorithm proceed with the next couple, and so on until we end up to
% the progenitor cell itself. 
% [This must use while-loop to ensure we plot all clones]

while length(cS_EndLine)>=2
    % cc is the counter of the cell position in arrays, such as cS_EndLine and clone_List 
    cc = 0 ; 
    
    while cc < length(cS_EndLine)
        cc = cc+1;
        
        % Find the two twins: ID name string as well as index positionin array
        % If a twin ID name end in 1, the other must be 2, and viceversa
        TWstr_A = cS_EndLine{cc}(2:end);
        if  TWstr_A(end) == '1'
            TWstr_B = [cS_EndLine{cc}(2:end-1) '2'];
        elseif  TWstr_A(end) == '2'
            TWstr_B = [cS_EndLine{cc}(2:end-1) '1'];   
        end
        idx_TA =  find(strcmp(cS_EndLine, ['.' TWstr_A] ) == 1) ;
        idx_TB =  find(strcmp(cS_EndLine, ['.' TWstr_B] ) == 1) ;
        
        % Plot the common ancestor
        if ~isempty(idx_TA) && ~isempty(idx_TB)            
            % Find ancestor ID name and postion in array
            str_Anc = [APP_opt.t3_PlotOpt_NLineage '.' TWstr_A(1:end-1)] ;
            idx_Anc = find(strcmp(cS_All , str_Anc ) == 1);            
            % The last time point of the ancestor is the birth of the two twins
            idx_birth = clone_List{idx_Anc}{1}.frame_last ;  
            
            % yrow = store the y position for plotting the current cc-th clone
            yrow = (y_H_Row(idx_TA) + y_H_Row(idx_TB)) /2 ;
                      
            % Plot vertical line at idx_birth (cell divisions), connecting 
            % two daughter cells (twins) in lineage tree
            Xs = [ idx_birth+1 , idx_birth+1 ];
            if y_H_Row(idx_TA) < y_H_Row(idx_TB)
                Ys = [ y_H_Row(idx_TA)*Dst_pl,  (y_H_Row(idx_TB)*Dst_pl)+hg_pl];
            elseif y_H_Row(idx_TA) > y_H_Row(idx_TB)
                Ys = [ y_H_Row(idx_TB)*Dst_pl,  (y_H_Row(idx_TA)*Dst_pl)+hg_pl];
            end
            plot(Xs, Ys, '-', 'LineWidth', LinW, 'Color', Color_BorderLines); 
    
            % Plot all time points. We use two counters for "time":
            % - ff  = stores the array position of each element (/frame) of  
            %         the clone_List{XXX}: simply goes from 1 to last frame
            % - Frm = x-axis "real" time variable, used to place the extracted
            %         Value at the correct time point of the lineage tree
            for ff = 1 : length(clone_List{idx_Anc})
                Frm = ff + clone_List{idx_Anc}{1}.frame_birth;
                if Frm > LAST_Frm;     LAST_Frm = Frm;     end

                % If choosen, plot HUD: cell poles orientation and cell_ID number
                if ff == 1      % at frame of birth;
                   fnc_Plot_HUD( clone_List{idx_Anc}{1}, Frm, yrow, Dst_pl, hg_pl , LinW)
                end 
                
                % Define the coordinates of the "box" square polygon and find the
                % signal value of the desired subcellular compartment to plot.
                [Xs, Ys, Value] = Find_Value_XYSquare( clone_List{idx_Anc}{ff}, Frm, yrow, Dst_pl, hg_pl, ChNum , SubLine, ...
                                                       APP_opt.t3_choose_ChannelMode, APP_opt.t3_PlotOpt_B, APP_opt.t3_PlotOpt_C, APP_opt.algorithm );

                % Now, evaluate the RGB_Color for each "Line" that we wish to plot for each clone
                if APP_opt.t3_PlotOpt_B == 1  |  APP_opt.t3_PlotOpt_B == 3
                    for kk = SubLine            
                        cRGB  = Assign_Value_RGB( Value{kk} , RGB_range{kk} );
                        % Fill squares with evaluated RGB_Color
                        fill(Xs{kk}, Ys{kk}, [cRGB(1), cRGB(2), cRGB(3)], 'LineStyle', 'none');
                    end
                elseif APP_opt.t3_PlotOpt_B == 2
                    for kk = SubLine
                        for ii = 1 : 2      % column 1 = Upper;   column 2 = Lower.
                            cRGB  = Assign_Value_RGB( Value{kk}(ii) , RGB_range{kk} );
                            % Fill squares with evaluated RGB_Color
                            fill(Xs, Ys{kk,ii}, [cRGB(1), cRGB(2), cRGB(3)], 'LineStyle', 'none');
                        end
                    end
                end
                                
                % --- Extra Info ---- Place dots over cell_line to represent cell's
                % eccentricity value [greys scale: black = 0/circle; white = 1/ellipse)]
                if APP_opt.t3_display_ExI == 1
                    eXs = [clone_List{idx_Anc}{ff}.mesh(:,1); flipud(clone_List{idx_Anc}{ff}.mesh(:,3))];
                    eYs = [clone_List{idx_Anc}{ff}.mesh(:,2); flipud(clone_List{idx_Anc}{ff}.mesh(:,4))];        
                    [la, sa] = Geom_My_fit_ellipse__v2( eXs, eYs );
                    % Calculate ellipse eccentricity (sa = short semi-axis; la = long semi-axis)
                    % Multipling by factor *0.9 restrict grey scale to a range 0:0.9 and 
                    % it avoids having "white" marker (which are hardly visible)
                    El_Ec = sqrt( 1 - (sa^2/la^2)) .*0.9 ;
                    % The ratio cannot be bigger than 1
                    if El_Ec > 0.9 ;      El_Ec = 0.9 ;     end
                    El_Ec_Col = [El_Ec El_Ec El_Ec] ;            
                    plot( Frm -0.5, (yrow*Dst_pl +hg_pl) +hg_pl/3 , '.', 'MarkerSize', 8, 'LineWidth', 0.2, 'Color', El_Ec_Col);
                end
                
                % --- SAVE_step_2 ---- Collect plotted value of all ancestors
                if APP_opt.t3_choose_Save_txt == 1  
                    for kk = SubLine
                        % Save the value(s) plotted
                        tree_Values{kk}{1, idx_Anc} = [tree_Values{kk}{1, idx_Anc} , Value{kk}'];                                                
                        % [px] Cell axial-length       
                        tree_Values{kk}{2, idx_Anc} = [tree_Values{kk}{2, idx_Anc}, clone_List{idx_Anc}{ff}.geom.length] ;    
                        % [px] Cell area
                        tree_Values{kk}{3, idx_Anc} = [tree_Values{kk}{3, idx_Anc}, clone_List{idx_Anc}{ff}.geom.area] ;      
                        % y-position in plot  
                        tree_Values{kk}{4, idx_Anc} = yrow*Dst_pl ; 
                    end
                end
                
            end % ff
                           
            % Draw two horizontal lines that bound the plotted data for the clone-Line
            plot( [Frm, Frm-length(clone_List{idx_Anc})], [yrow*Dst_pl,  yrow*Dst_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
            plot( [Frm, Frm-length(clone_List{idx_Anc})], [(yrow*Dst_pl)+hg_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
            
            % Draw vertical line at the end, to enclose all end_Lines
            plot( [Frm, Frm], [yrow*Dst_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);

            % If necessary draw a line to separate the two channels plot-lines
            if APP_opt.t3_choose_ChannelMode == 3
                plot( [Frm, Frm-length(clone_List{idx_Anc})], [(yrow*Dst_pl)+(hg_pl*1/2),  (yrow*Dst_pl)+(hg_pl*1/2)], 'LineWidth', LinW/2, 'Color', Color_BorderLines);
            end
    
            % UPDATE cS_EndLine ----------------------------------------
            % Remove the ID of the couple we just plotted and...
            % (we do not need to update cN_EndLine)
            cS_EndLine{idx_TA} = '';
            cS_EndLine{idx_TB} = '';
            y_H_Row(idx_TA) = yrow;            
            y_H_Row(idx_TB) = [];
            % ... put their common ancestor ID_anc in the list, since now it is a "end-line"                  
            [String] = strsplit(str_Anc, '.');
            cS_EndLine{idx_TA} = [ '.' String{2} ];
            % ... and we remove empty elements from the list
            cS_EndLine = cS_EndLine(~cellfun(@isempty, cS_EndLine)) ;

        end
    end  % while cc
end % while length

% Last idx_Anc is lineage progenitor: draw vertical line at the start, to enclose Plot_Line
plot( [clone_List{idx_Anc}{1}.frame_birth , clone_List{idx_Anc}{1}.frame_birth],...
      [yrow*Dst_pl , (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);


  
  
  
  

elseif APP_opt.t3_PlotOpt_A == 2            % if GUI Plot Type: '2 - Individual Clones'
%% 
% ***** PLOT CLONES INDEPENDENTLY ****************************************************************

% Make figure visible of invisible according to user choice   
if APP_opt.t3_choose_Display_Plot == 0
    h1f = figure('Position', [100 300 1200 600], 'Visible','off'); 
elseif APP_opt.t3_choose_Display_Plot == 1
    h1f = figure('Position', [100 300 1200 600]);    
end
hold on;       
ax = gca; 
% In independent clones plotting, there are as many as the number of clones
y_H_Row = [1:length(clone_List)];
   

% ----- Plot all clones independently
for cc = 1 : length(clone_List)    
    % yrow = store the y position where the current cc-th cell is plotted
    %        (in independent clone-plot this is just cc)
    yrow = cc ;
    
    % ff = only one counter for "time", since we plot each clone starting 
    % from x=0, array position and time (x-axis) coincides.
    for ff = 1 : length(clone_List{cc})        
        if ff > LAST_Frm;     LAST_Frm = ff;     end

        % If choosen, plot HUD: cell poles orientation and cell_ID number
        if ff == 1      % at frame of birth;
           fnc_Plot_HUD( clone_List{cc}{1}, ff, yrow, Dst_pl, hg_pl , LinW)
        end
        
        % Define the coordinates of the "box" square polygon and find the
        % signal value of the desired subcellular compartment to plot.
        [Xs, Ys, Value] = Find_Value_XYSquare( clone_List{cc}{ff}, ff, yrow, Dst_pl, hg_pl, ChNum , SubLine, ...
                                               APP_opt.t3_choose_ChannelMode, APP_opt.t3_PlotOpt_B, APP_opt.t3_PlotOpt_C, APP_opt.algorithm );

        % Now, evaluate the RGB_Color for each "Line" that we wish to plot for each clone
        if APP_opt.t3_PlotOpt_B == 1  |  APP_opt.t3_PlotOpt_B == 3
            for kk = SubLine            
                cRGB  = Assign_Value_RGB( Value{kk} , RGB_range{kk} );
                % Fill squares with evaluated RGB_Color
                fill(Xs{kk}, Ys{kk}, [cRGB(1), cRGB(2), cRGB(3)], 'LineStyle', 'none');
            end
        elseif APP_opt.t3_PlotOpt_B == 2
            for kk = SubLine
                for ii = 1 : 2      % column 1 = Upper;   column 2 = Lower.
                    cRGB  = Assign_Value_RGB( Value{kk}(ii) , RGB_range{kk} );
                    % Fill squares with evaluated RGB_Color
                    fill(Xs, Ys{kk,ii}, [cRGB(1), cRGB(2), cRGB(3)], 'LineStyle', 'none');
                end
            end
        end 

        % --- Extra Info ---- Place dots over cell_line to represent cell's
        % eccentricity value [greys scale: black = 0/circle; white = 1/ellipse)]
        if APP_opt.t3_display_ExI == 1
            eXs = [clone_List{cc}{ff}.mesh(:,1); flipud(clone_List{cc}{ff}.mesh(:,3))];
            eYs = [clone_List{cc}{ff}.mesh(:,2); flipud(clone_List{cc}{ff}.mesh(:,4))];        
            [la, sa] = Geom_My_fit_ellipse__v2( eXs, eYs );
            % Calculate ellipse eccentricity (sa = short semi-axis; la = long semi-axis)
            % Multipling by factor *0.9 restrict grey scale to a range 0:0.9 and 
            % it avoids having "white" marker (which are hardly visible)
            El_Ec = sqrt( 1 - (sa^2/la^2)) .*0.9 ;
            % The ratio cannot be bigger than 1
            if El_Ec > 0.9 ;      El_Ec = 0.9 ;     end
            El_Ec_Col = [El_Ec El_Ec El_Ec] ;            
            plot( ff -0.5, (yrow*Dst_pl +hg_pl) +hg_pl/3 , '.', 'MarkerSize', 8, 'LineWidth', 0.2, 'Color', El_Ec_Col);
        end
                
        % --- SAVE_step_1-2 ---- Collect plotted value of all clones
        if APP_opt.t3_choose_Save_txt == 1   
            for kk = SubLine
                % Save the value(s) plotted
                tree_Values{kk}{1, cc} = [tree_Values{kk}{1, cc} , Value{kk}'] ;                   
                % [px] Cell axial-length       
                tree_Values{kk}{2, cc} = [tree_Values{kk}{2, cc}, clone_List{cc}{ff}.geom.length] ;    
                % [px] Cell area
                tree_Values{kk}{3, cc} = [tree_Values{kk}{3, cc}, clone_List{cc}{ff}.geom.area] ;
                % y-position in plot  
                tree_Values{kk}{4, cc} = yrow*Dst_pl ; 
            end
        end
    end % ff
    
    % Draw two horizontal lines that bound the plotted data for the clone-Line
    plot( [0, length(clone_List{cc})], [yrow*Dst_pl,  yrow*Dst_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
    plot( [0, length(clone_List{cc})], [(yrow*Dst_pl)+hg_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
    
    % Draw vertical line at the STRAT and END, to enclose all end_Lines
    plot( [0, 0], [yrow*Dst_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
    plot( [ff, ff], [yrow*Dst_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);

    % If necessary draw a line to separate the two channels plot-lines
    if APP_opt.t3_choose_ChannelMode == 3
        plot( [ff, ff-length(clone_List{cc})], [(yrow*Dst_pl)+(hg_pl*1/2),  (yrow*Dst_pl)+(hg_pl*1/2)], 'LineWidth', LinW/2, 'Color', Color_BorderLines);
    end
    
end % cc

end % if Plot





%% *** SAVE the FIGURE and DATA *******************************************

% --- FINISH layout and SAVE Plot -----------------------------------------
h1f.Color = [1 1 1];
ax.TickDir = 'out';                 ax.YTick = [] ;                 
ax.XColor = Color_BorderLines ;     ax.YColor = 'none' ;
ax.LineWidth = LinW ;               ax.FontSize = 15 ;
% Create extra space at end of x axis, where place colorbar legends (simply extending of  MaxX+Delta_X)
XTicksPos = linspace(ax.XLim(1), ax.XLim(2), length(ax.XTick)) ./ fpmRate ;
XTicksPos = [XTicksPos , XTicksPos(end)+(XTicksPos(end)-XTicksPos(end-1))];
% Stretch the plot a bit in the x-axis direction
ax.XLim(2)     = ax.XTick(end) + ( ax.XTick(end)-ax.XTick(end-1)) ;
ax.XTickLabel  = XTicksPos ;


% ---> PLOT colorbars (for all channels)
for kk = SubLine
    % Use correct label for the colorbar legends
    if APP_opt.t3_choose_ChannelMode <= 3
        barLabel = ['Channel ' num2str(ChNum(kk))] ; 
    elseif APP_opt.t3_choose_ChannelMode == 4
        barLabel = 'Ratio 1-2' ;
    elseif APP_opt.t3_choose_ChannelMode == 5
        barLabel = 'Ratio 2-1' ;
    end
    
    % To show two separate colorbar legends, we need create second x-axis.
    % Colormap associate legend to one axis only, and if run it twice it
    % simply overwrite and not creating a second one.
    switch kk
        case 1
            barAX = ax ;            
        case 2
            ax1_pos = ax.Position;                  
            barAX = axes('Position', ax1_pos, 'XAxisLocation','top', 'YAxisLocation','right', 'Color','none');
            barAX.YTick = [] ;                barAX.YColor = 'none';
            barAX.XTick = [] ;                barAX.XColor = 'none';            
    end    
    TickPos = linspace(0,1,7);
    ValuePos = floor(size(RGB_range{kk},1) .* TickPos) ;
    ValuePos(ValuePos == 0) = 1;
    ValuePos = RGB_range{kk}(ValuePos,4) ;
    
    % Choose a readable number format for the legends
    StrBar = {};
    for vv = 1 : size(ValuePos,1)
        if ValuePos(vv) < 1000             
            StrBar{vv} = num2str(ValuePos(vv), '%.2f') ;
        else
            StrBar{vv} = num2str(ValuePos(vv), '%3.2e') ;
        end
    end
    % Plot the colorbar, adjust the position and size. Channel 1 is the
    % lower one (if there are twwo channels, 2 is positioned above the 1)
    m  = colormap(barAX, RGB_range{kk}(:,1:3) ./255);
    cp = colorbar(barAX, 'FontSize', 12 , 'Ticks',TickPos, 'TickLabels', StrBar );
    cp.Location     = 'eastoutside';
    cp.Label.String = barLabel;
    cp.Position(1)  = 0.85;
    cp.Position(2)  = ((kk-1)*0.4) + cp.Position(2) + cp.Position(3)*2;
    cp.Position(4)  = cp.Position(4)*0.4;    
    cp.LineWidth    = LinW;
    cp.Color        = Color_BorderLines ;
end
% Render the picture, improving quality of .tif images and vectorial .pdf
h1f.Renderer = 'Painters';


% Select appropriate filename ending according to analysis performed
terminator = [num2str(APP_opt.t3_PlotOpt_A) '_' num2str(APP_opt.t3_PlotOpt_B) '_' num2str(APP_opt.t3_PlotOpt_C)];
if APP_opt.t3_choose_AutoSave_Plot ~= 0
    if isempty(APP_opt.t3_exp_name)
        filename_plot = [ APP_opt.t3_path_cloneList ,'/', 'Plot_' terminator ];
    else
        filename_plot = [ APP_opt.t3_path_cloneList ,'/' , APP_opt.t3_exp_name , '_', 'Plot_' terminator ];
    end
    if     APP_opt.t3_choose_AutoSave_Plot == 1;       print(h1f, filename_plot ,'-dpdf');
    elseif APP_opt.t3_choose_AutoSave_Plot == 2;       print(h1f, filename_plot , '-r300', '-dtiffn');
    elseif APP_opt.t3_choose_AutoSave_Plot == 3;       print(h1f, filename_plot ,'-dsvg');
    end
end



% --- SAVE_step_3 ---- Save Plotted data in a .txt file -----------------
if APP_opt.t3_choose_Save_txt == 1 
    
    % Reorganize data to have all have same "length" and data at the correct "column" position.

    for kk = SubLine
        for cc = 1 : size(clone_List,2)
            for rr = [1,2,3]
                nRows = size(tree_Values{kk}{rr,cc},1);        % number of rows == values saved in variable
                if APP_opt.t3_PlotOpt_A == 1
                    tree_Values{kk}{rr,cc} = [ repmat(-1, nRows, clone_List{cc}{1}.frame_birth -1) , ...
                                               tree_Values{kk}{rr,cc}, ...
                                               repmat(-1, nRows, last_timepoint - length([repmat(-1, nRows, clone_List{cc}{1}.frame_birth-1), tree_Values{kk}{rr,cc}])) ] ;
                elseif APP_opt.t3_PlotOpt_A == 2
                    tree_Values{kk}{rr,cc} = [ tree_Values{kk}{rr,cc}, repmat(-1, nRows, last_timepoint - length( tree_Values{kk}{rr,cc})) ];
                end
            end %rr
        end %cc

        % Sort data and save it as .txt in the same vertical order as it was plotted
        [~,indx] = sort( cell2mat(tree_Values{kk}(4,:)) , 'descend');

        % ---- Fluorescent Signal Data ------------------------------------
        % We do inside kk for-loop, because we can have 2 channels, thus two
        % separate files must be generated
        % ValueRows: number of rows == values plotted per line
        ValueRows = size(tree_Values{1,1}{1,1},1);        
        file_S = fopen([Path_txt 'Signal_CH', num2str(ChNum(kk)) , '_S_' terminator '.txt' ], 'w+');     
        fprintf( file_S , 'Clone ID\tTime Birth\tFrame-to-Min\t' );            % First row is the time point numbers
        fprintf( file_S , '%f\t', [1 : last_timepoint]' );
        fprintf( file_S , '\n' );
        for ii = indx
            fprintf( file_S , '%s\t', clone_List{ii}{1}.ID_clone );
            fprintf( file_S , '%f\t', clone_List{ii}{1}.frame_birth );
            fprintf( file_S , '%f\t', 1/ fpmRate );
            for nn = 1 : ValueRows
                if nn == 1 
                    fprintf( file_S , '%f\t', tree_Values{kk}{1,ii}(nn,:) );
                else
                    fprintf( file_S , '\t \t \t');
                    fprintf( file_S , '%f\t', tree_Values{kk}{1,ii}(nn,:) );
                end
            fprintf( file_S , '\n' );
            end %nn
        end %ii
        fclose(file_S) ;         
        
    end %kk

    % tree_Values{1}: either channel has same infos, first element is always present
    % ---- Cell axial-length Data -------------------------------------
    file_L = fopen([Path_txt 'DataL_' terminator '.txt' ], 'w+');  
    fprintf( file_L , 'Clone ID\tTime Birth\tFrame-to-Min\t' );            % First row is the time point numbers
    fprintf( file_L , '%f\t', [1 : last_timepoint]' );
    fprintf( file_L , '\n' );
    for ii = indx
        fprintf( file_L , '%s\t', clone_List{ii}{1}.ID_clone );
        fprintf( file_L , '%f\t', clone_List{ii}{1}.frame_birth );
        fprintf( file_L , '%f\t', 1/ fpmRate );
        fprintf( file_L , '%f\t', tree_Values{1}{2,ii} );
        fprintf( file_L , '\n' );
    end
    fclose(file_L) ;

    % ---- Cell area Data --------------------------------------------
    file_A = fopen([Path_txt 'DataA_' terminator '.txt' ], 'w+');    
    fprintf( file_A , 'Clone ID\tTime Birth\tFrame-to-Min\t' );            % First row is the time point numbers
    fprintf( file_A , '%f\t', [1 : last_timepoint]' );
    fprintf( file_A , '\n' );
    for ii = indx
        fprintf( file_A , '%s\t', clone_List{ii}{1}.ID_clone );
        fprintf( file_A , '%f\t', clone_List{ii}{1}.frame_birth );
        fprintf( file_A , '%f\t', 1/ fpmRate );
        fprintf( file_A , '%f\t', tree_Values{1}{3,ii} );        
    end
    fclose(file_A) ;
    
end % Save .txt

end % MAIN fnc




