function Plot_Lng_M4
%  Here we Plot the whole genalogy tree for a single lineage
%
%  SEGMENTATION
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

% OLD_pole = Mask_PL_1
% NEW_pole = Mask_PL_2
  
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

% Convert the clone_List in a track_List
[ track_List , info_Track ] = Convert_Clone2Track(clone_List) ;

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
    % There cannot be two channels at the same time
    
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

% find the Max_CellLen, the number of segments the longest cell has
if APP_opt.t3_PlotOpt_C <= 3 
    [ Seg_Len, Max_CellLen ] = Find_NumberSegments(track_List, max(ChNum), 1) ;
elseif APP_opt.t3_PlotOpt_C >= 4 
    [ Seg_Len, Max_CellLen ] = Find_NumberSegments(track_List, max(ChNum), 2) ;
end
% In segmentation plot the "cell-line" is enlarged and can overlapping with
% other features of the plot. Thus, we take the 0.75% of Dst_pl , which it
% helps to avoid overlap, especially for additional Extra Info and Pole ID
stp_Seg = (Dst_pl*0.75) / Max_CellLen ;


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
    tree_Values= cell(4,length(track_List));    

    
    % Find the latest time point for entire clone list, necessary to save data
    % in .txt file by ensuring that each clone data has same length.
    last_timepoint = -1 ;
    for cc = 1 : size(track_List,2)
        if last_timepoint <= track_List{cc}{1}.frame_last
           last_timepoint = track_List{cc}{1}.frame_last ;
        end
    end
end


%%
% ***** PLOT ALL TRACKS ****************************************************************

% Make figure visible of invisible according to user choice   
if APP_opt.t3_choose_Display_Plot == 0
    h1f = figure('Position', [100 300 1200 600], 'Visible','off'); 
elseif APP_opt.t3_choose_Display_Plot == 1
    h1f = figure('Position', [100 300 1200 600]);    
end
hold on;       
ax = gca; 
% In independent clones plotting, there are as many as the number of clones
y_H_Row = [1:length(track_List)];
   

% ----- Plot all clones independently
for cc = 1 : length(track_List)    
    % yrow = store the y position where the current cc-th cell is plotted
    %        (in independent clone-plot this is just cc)
    yrow = cc ;
    
    % ff = only one counter for "time", since we plot each clone starting 
    % from x=0, array position and time (x-axis) coincides.
    for ff = 1 : length(track_List{cc})        
        if ff > LAST_Frm;     LAST_Frm = ff;     end

        % If choosen, plot HUD: cell poles orientation and cell_ID number
        if ff == 1      % at frame of birth;
           fnc_Plot_HUD( track_List{cc}{1}, ff, yrow, Dst_pl, hg_pl , LinW)
        end
        
        % Define the coordinates of the "box" square polygon and find the
        % signal value of the desired subcellular compartment to plot.
        [Xs, Ys, Value] = Find_Value_XYSquare( track_List{cc}{ff}, ff, yrow, Dst_pl, stp_Seg, ChNum , SubLine, ...
                                               APP_opt.t3_choose_ChannelMode, APP_opt.t3_PlotOpt_B, APP_opt.t3_PlotOpt_C, APP_opt.algorithm,  Seg_Len{cc}{ff} );

        % Transfor matrix of segments into a cell array, which is the format
        % compatible with the function Assign_Value_RGB  
        Value = mat2cell(Value{:}, 1, ones(1,length(Value{:})));

        % Now, evaluate the RGB_Color for each "Line" that we wish to plot for each clone
        for rr = 1 : Seg_Len{cc}{ff}            
            cRGB  = Assign_Value_RGB( Value{rr} , RGB_range{1} );
            % Fill squares with evaluated RGB_Color
            fill(Xs(rr,:), Ys(rr,:), [cRGB(1), cRGB(2), cRGB(3)], 'LineStyle', 'none');
        end

        % --- Extra Info ---- Place dots over cell_line to represent cell's
        % eccentricity value  [greys scale: black = 0/circle; white = 1/ellipse)]
        if APP_opt.t3_display_ExI == 1
            eXs = [track_List{cc}{ff}.mesh(:,1); flipud(track_List{cc}{ff}.mesh(:,3))];
            eYs = [track_List{cc}{ff}.mesh(:,2); flipud(track_List{cc}{ff}.mesh(:,4))];       
            [la, sa] = Geom_My_fit_ellipse__v2( eXs, eYs );
            % Calculate ellipse eccentricity (sa = short semi-axis; la = long semi-axis)
            % Multipling by factor *0.9 restrict grey scale to a range 0:0.9 and 
            % it avoids having "white" marker (which are hardly visible)
            El_Ec = sqrt( 1 - (sa^2/la^2)) .*0.9 ;
            % The ratio cannot be bigger than 1
            if El_Ec > 0.9 ;      El_Ec = 0.9 ;     end
            El_Ec_Col = [El_Ec El_Ec El_Ec] ;  
            % establish the Y-Position where to place the label
            if APP_opt.t3_PlotOpt_C  == 1  ||   APP_opt.t3_PlotOpt_C  == 4 
                Y_Position = (yrow*Dst_pl -hg_pl/3) ;
            elseif APP_opt.t3_PlotOpt_C  == 3  ||   APP_opt.t3_PlotOpt_C  == 6
                Y_Position = (yrow*Dst_pl +hg_pl/3) ;
            elseif APP_opt.t3_PlotOpt_C  == 2  ||   APP_opt.t3_PlotOpt_C  == 5
                Y_Position = Ys(end, 3) +hg_pl/3 ;
            end
            plot( ff -0.5, Y_Position , '.', 'MarkerSize', 8, 'LineWidth', 0.2, 'Color', El_Ec_Col);
        end
        
        % --- SAVE_step_1-2 ---- Collect plotted value of all clones        
        % Save the value(s) plotted
        if APP_opt.t3_choose_Save_txt == 1   
            % Save the value(s) plotted
            tree_Values{1, cc} = [tree_Values{1, cc} , {Value}];     
            % [px] Cell axial-length       
            tree_Values{2, cc} = [tree_Values{2, cc} , track_List{cc}{ff}.geom.length] ;    
            % [px] Cell area
            tree_Values{3, cc} = [tree_Values{3, cc} , track_List{cc}{ff}.geom.area] ;
            % y-position in plot  
            tree_Values{4, cc} = [tree_Values{4, cc} ,  ff ] ;
        end
                    
        % Draw two horizontal lines that bound the plotted data for the clone-Line
        plot( [ff-1, ff], Ys(1, 1:2), 'LineWidth', LinW, 'Color', Color_BorderLines);
        plot( [ff-1, ff], Ys(end-1, 3:4), 'LineWidth', LinW, 'Color', Color_BorderLines);
        if ff == 1
            % Draw vertical line at the first frame, to enclose cell
            plot( [ff-1, ff-1], [Ys(1, 1), Ys(end-1, 3)], 'LineWidth', LinW, 'Color', Color_BorderLines);
        end
        
    end % ff     

    % Draw vertical line at the END, to enclose all end_Lines
    plot( [ff, ff], [Ys(1, 1), Ys(end-1, 3)], 'LineWidth', LinW, 'Color', Color_BorderLines);
    
    % Plot a marker and track ID where the cc-th track generated new track(s)
    if ~isempty(info_Track{2,cc})  &&  APP_opt.t3_display_cID >= 1
      Dau_tracks = cellTrack{6,cc}( find(cellTrack{6,cc}~=0) );
      if ~isempty(Dau_tracks)
          % establish the Y-Position where to place the label
          if APP_opt.t3_PlotOpt_C  == 1  ||   APP_opt.t3_PlotOpt_C  == 4 
              Y_Position = (yrow*Dst_pl -hg_pl/3) ;
          elseif APP_opt.t3_PlotOpt_C  == 3  ||   APP_opt.t3_PlotOpt_C  == 6
              Y_Position = (yrow*Dst_pl +hg_pl/3) ;
          elseif APP_opt.t3_PlotOpt_C  == 2  ||   APP_opt.t3_PlotOpt_C  == 5
              Y_Position = Ys(end, 3) +hg_pl/3 ;
          end
          for dd = 1 : length(info_Track{2,cc})
              fdiv = info_Track{2,cc}(dd);
              plot( (fdiv - info_Track{1,cc}) -0.5 , Y_Position +hg_pl/1.5 , ...
                  'v', 'MarkerSize', 6, 'MarkerFaceColor', Color_BorderLines, 'MarkerEdgeColor', Color_BorderLines);
              if  APP_opt.t3_display_cID >= 1
                  text( (fdiv - info_Track{1,cc}) +0.5 , Y_Position +hg_pl/1.5 , ...
                        num2str(Dau_tracks(dd)), 'FontSize', hg_pl);
              end
          end
      end
    end
    
end % cc





%% *** SAVE the FIGURE and DATA *******************************************

% --- FINISH layout and SAVE Plot -----------------------------------------
h1f.Color = [1 1 1];
ax.TickDir = 'out';                 ax.YTick = [] ;                 
ax.XColor = Color_BorderLines ;     ax.YColor = 'none' ;
ax.LineWidth = LinW ;               ax.FontSize = 15 ;
% Create extra space at end of x axis, where place colorbar legends (simply extending of  MaxX+Delta_X)
XTicksPos = linspace(ax.XLim(1), ax.XLim(2), length(ax.XTick)) ./ fpmRate ;
XTicksPos = [XTicksPos , XTicksPos(end)+(XTicksPos(end)-XTicksPos(end-1))];
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
   
    % We save all data in a single file       
    if APP_opt.t3_choose_ChannelMode <= 2
        FileTag = ['CH' num2str(ChNum(kk))] ; 
    elseif APP_opt.t3_choose_ChannelMode == 4
        FileTag = 'R1-2' ;
    elseif APP_opt.t3_choose_ChannelMode == 5
        FileTag = 'R2-1' ;
    end
    file_S = fopen([Path_txt 'Data_Segment_', FileTag , '_' terminator '.txt' ], 'w+');     fprintf( file_S , 'Track ID\tTime Birth\tFrame-to-Min\tFrame\tCell_Legnth\tCellArea\t' );           
    fprintf( file_S , '%f\t', [1 : Max_CellLen]' );
    fprintf( file_S , '\n' );

    % Reorganize data to have all have same "length" and data at the correct "column" position.
    for cc = 1 : size(track_List,2)
        for ff = 1 : size(track_List{1,cc},2)
            fprintf( file_S , '%s\t', num2str(track_List{cc}{1}.ID_ManualTrack) );
            fprintf( file_S , '%f\t', track_List{cc}{1}.frame_birth );
            fprintf( file_S , '%f\t', 1/ fpmRate );
            fprintf( file_S , '%f\t', tree_Values{4,cc}(ff) );        
            fprintf( file_S , '%f\t', tree_Values{2,cc}(ff) );
            fprintf( file_S , '%f\t', tree_Values{3,cc}(ff) );   
            
            % Create the line with segmentation values, adding -1 before and
            % after to fill a Line and ensure that they all the same length
            Val =  cell2mat(tree_Values{1,cc}{1,ff}) ;                   
            Z = Max_CellLen - length(Val);
            if mod(Z,2) == 0
                plus = ones(1, (Z/2))*(-1);
                aLine = [plus, Val , plus];
            else
                plus = ones(1, fix(Z/2))*(-1);      % fix(), division returning integer
                aLine = [plus, Val , plus, -1];
            end
            fprintf( file_S , '%f\t', aLine );
            fprintf( file_S , '\n' );
            
        end %ff
    end %cc

    fclose(file_S) ;
    
end % Save .txt

end % MAIN fnc
