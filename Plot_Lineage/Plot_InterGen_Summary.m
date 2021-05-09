function Plot_InterGen_Summary(Path_txt, RGB_range)
% Here we Plot the whole genalogy tree for a single lineage
%
%
%
%
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


%% ----- Initialize -------------------------------------------------------
% Gather plotting parameters, the distribution of signal values and create
% a corresponding RGB colormap

% OLD_pole = Mask_PL_1
% NEW_pole = Mask_PL_2

global APP_opt ;

hg_pl = APP_opt.t3_height_PlotLine ;         % heights of the plot end lines
Dst_pl = APP_opt.t3_Dist_btw_PlotLines ;     % distances between plot end lines 
LinW =  APP_opt.t3_BorderLineWidth ;         % LineWidth of line connecting ancestor-to-daughters
Color_BorderLines = [.3 .3 .3]; 
LAST_Frm = 1;                                % highest value frame, use to find max x-axis length
fpmRate = 1 / APP_opt.t3_fpm_CH1 ;

for kk = 1 : size(Path_txt,2)
    for gg = 1 : size(Path_txt{kk},1)
        % Load the clone_List   
        ggData = importdata( Path_txt{1,kk}{1,gg} );
        

        

        

%% 
% ***** PLOT TRACKS ****************************************************************
    
% Make figure visible of invisible according to user choice
if APP_opt.t3_choose_Display_Plot == 0
    h1f = figure('Position', [100 300 1200 600], 'Visible','off'); 
elseif APP_opt.t3_choose_Display_Plot == 1
    h1f = figure('Position', [100 300 1200 600]);    
end
hold on;       
ax = gca;
% In independent clones plotting, there are as many as the number of tracks
y_H_Row = 1 ; 
                                       

% ----- Plot all Tracks independently
for cc = 1 : length(track_List)    
    % yrow = store the y position where the current cc-th cell is plotted
    %        (in independent clone-plot this is just cc)
    yrow = cc ;
    
    % ff = only one counter for "time", which is the same as array position
    for ff = 1 : length(track_List{cc})        
        if ff > LAST_Frm;     LAST_Frm = ff;     end

        % If choosen, plot HUD: cell poles orientation and cell_ID number
        if ff == 1      % at frame of birth;
           fnc_Plot_HUD( track_List{cc}{1}, ff, yrow, Dst_pl, hg_pl , LinW)
        end
        
        % Define the coordinates of the "box" square polygon and find the
        % signal value of the desired subcellular compartment to plot.
        [Xs, Ys, Value] = Find_Value_XYSquare( track_List{cc}{ff}, ff, yrow, Dst_pl, hg_pl, ChNum , SubLine, ...
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
            plot( ff -0.5, (yrow*Dst_pl +hg_pl) +hg_pl/4 , '.', 'MarkerSize', 8, 'LineWidth', 0.2, 'Color', El_Ec_Col);
        end
        
        % --- SAVE_step_1 ---- Collect plotted value of all clones
        if APP_opt.t3_choose_Save_txt == 1  
            for kk = SubLine
                % Save the value(s) plotted
                tree_Values{kk}{1, cc} = [tree_Values{kk}{1, cc} , Value{kk}'];                                    
                % [px] Cell axial-length       
                tree_Values{kk}{2, cc} = [tree_Values{kk}{2, cc}, track_List{cc}{ff}.geom.length] ;    
                % [px] Cell area
                tree_Values{kk}{3, cc} = [tree_Values{kk}{3, cc}, track_List{cc}{ff}.geom.area] ;
                % y-position in plot  
                tree_Values{kk}{4, cc} = yrow*Dst_pl ; 
            end
        end
                
    end % ff

    % Draw two horizontal lines that bound the plotted data for the clone-Line
    plot( [0, length(track_List{cc}) ], [yrow*Dst_pl,  yrow*Dst_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
    plot( [0, length(track_List{cc}) ], [(yrow*Dst_pl)+hg_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
    
    % Draw vertical line at the STRAT and END, to enclose all end_Lines
    plot( [0, 0], [yrow*Dst_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
    plot( [ff, ff], [yrow*Dst_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
    
    % If necessary draw a line to separate the two channels plot-lines
    if APP_opt.t3_choose_ChannelMode == 3
        plot( [ff, ff-length(track_List{cc})], [(yrow*Dst_pl)+(hg_pl*1/2),  (yrow*Dst_pl)+(hg_pl*1/2)], 'LineWidth', LinW/2, 'Color', Color_BorderLines);
    end
    
    % Plot a marker and track ID where the cc-th track generated new track(s)
    if ~isempty(info_Track{2,cc})  &&  APP_opt.t3_display_cID >= 1
      Dau_tracks = cellTrack{6,cc}( find(cellTrack{6,cc}~=0) );
      if ~isempty(Dau_tracks)
          for dd = 1 : length(info_Track{2,cc})
              fdiv = info_Track{2,cc}(dd);
              plot( (fdiv - info_Track{1,cc}) -0.5, (yrow*Dst_pl +hg_pl) +hg_pl/1.5 , ...
                  'v', 'MarkerSize', 6, 'MarkerFaceColor', Color_BorderLines, 'MarkerEdgeColor', Color_BorderLines);
              text( (fdiv - info_Track{1,cc}) +1, (yrow*Dst_pl +hg_pl) +hg_pl/1.5 , ...
                     num2str(Dau_tracks(dd)), 'FontSize', hg_pl);
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
% Create a little extra space at end of x axis, where we will place the
% colorbar legends (simply by extending of  MaxX+Delta_X)
XTicksPos = linspace(ax.XLim(1), ax.XLim(2), length(ax.XTick)) ./ fpmRate ;
XTicksPos = [XTicksPos , XTicksPos(end)+(XTicksPos(end)-XTicksPos(end-1))];
ax.XTickLabel  = XTicksPos ;
ax.XLim(2)     = XTicksPos(end) ;


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
        if ValuePos(vv) < 1             
            StrBar{vv} = num2str(ValuePos(vv), '%.3f') ;
        elseif ValuePos(vv) < 10000
            StrBar{vv} = num2str(round(ValuePos(vv)));
        else
            StrBar{vv} = num2str(ValuePos(vv), '%4.2e') ;
        end
    end          % StrBar = num2str(RGB_range{kk}(ValuePos,4), '%4.2e\n') ;
    
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
    elseif APP_opt.t3_choose_AutoSave_Plot == 2;       print(h1f, filename_plot ,'-r300', '-dtiffn');
    elseif APP_opt.t3_choose_AutoSave_Plot == 3;       print(h1f, filename_plot ,'-dsvg');
    end
end



% --- SAVE_step_2 ---- Save Plotted data in a .txt file -------------------
if APP_opt.t3_choose_Save_txt == 1  
    
  % Reorganize data to have all have same "length" and data at the correct "column" position.
    for kk = SubLine
        for cc = 1 : size(track_List,2)
          for rr = [1,2,3]
            nRows = size(tree_Values{kk}{rr,cc},1);        % number of rows == values saved in variable
            tree_Values{kk}{rr,cc} = [ tree_Values{kk}{rr,cc}, repmat(-1, nRows, last_timepoint - length( tree_Values{kk}{rr,cc})) ];
          end %rr
        end %cc

        % Sort data and save it as .txt in the same vertical order as it was plotted
        [~,indx] = sort( cell2mat(tree_Values{kk}(4,:)) , 'descend');

        % ---- Fluorescent Signal Data  ------------------------------------
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
    % ---- Cell axial-length Data  ------------------------------------
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

    % ---- Cell area Data  ------------------------------------
    file_A = fopen([Path_txt 'DataA_' terminator '.txt' ], 'w+');    
    fprintf( file_A , 'Clone ID\tTime Birth\tFrame-to-Min\t' );            % First row is the time point numbers
    fprintf( file_A , '%f\t', [1 : last_timepoint]' );
    fprintf( file_A , '\n' );
    for ii = indx
        fprintf( file_A , '%s\t', clone_List{ii}{1}.ID_clone );
        fprintf( file_A , '%f\t', clone_List{ii}{1}.frame_birth );
        fprintf( file_A , '%f\t', 1/ fpmRate );
        fprintf( file_A , '%f\t', tree_Values{1}{3,ii} );
        fprintf( file_A , '\n' );
    end
    fclose(file_A) ;
  
end % Save .txt


    end
end

end % MAIN fnc
