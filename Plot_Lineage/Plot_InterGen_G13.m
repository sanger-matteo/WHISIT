function Plot_InterGen_G13
% Here we Plot the whole genalogy tree for a single lineage
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

% Reorganize all the clone_list in the folder given by the user into
% gen_List, where each element stores all the clones of x-th generation
[gen_List, gen_Fold] = Create_GenList ;


% Using the entire distibution of values of all generations, create one RGB 
% colormap for each channel
if APP_opt.t3_choose_NormColormap == 1
    % Gather the distribution of all values (for each channel) needed to plot. 
    % Note if APP_opt.t3_choose_ChannelMode >= 4, Hist_Values_generationList
    % return the ratio between the two channels.
    for kk = ChNum
        [ Mu, Sigma, Distr_Vals ] = Hist_Values_generationList(0, gen_List);
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
end



% ----- INITIALIZE the summary file containing plotted data of all generations ----------------
% Find the latest time point for entire clone list, necessary to save data
% in .txt file by ensuring that each clone data has same length.
last_timepoint = -1 ;
for aa = 1 : size(gen_List,2)
    for bb = 1 : size(gen_List{1,aa},2)
        if last_timepoint <= size(gen_List{1,aa}{1,bb} ,2)
           last_timepoint  = size(gen_List{1,aa}{1,bb} ,2);
        end
    end
end

% Based on analysis performed, add a terminal identifier at the end of filename
terminator = [ '2_' num2str(APP_opt.t3_PlotOpt_B) '_' num2str(APP_opt.t3_PlotOpt_C)];
% Create a subfolder to store all the results, if it does not exsist, and
% all the subfolder for each generation, if is necessary.   
if exist(gen_Fold{2,1}) == 0
    mkdir(gen_Fold{2,1}) ;
end
if APP_opt.t3_choose_AutoSave_Plot ~= 0 
    for hh = 1 : size(gen_List,2)            
        if  exist( [gen_Fold{2,hh} '/' gen_Fold{3,hh}] ) == 0
            mkdir( [gen_Fold{2,hh} '/' gen_Fold{3,hh}] ) ;
        end            
    end
end

% Create handles for summary files, that carries the data for all generations
% [ with handles named File_A ]
for kk = SubLine                    % For 2 channels, create separate files
    file_A_S{kk} = fopen( [ gen_Fold{2,1} '/Summary_Sig_CH', num2str(ChNum(kk)) , '_', terminator, '.txt']  , 'w+');
    fprintf( file_A_S{kk} , 'Generation\tTrack ID\tTime Birth\tFrame-to-Min\t' );
    fprintf( file_A_S{kk} , '%f\t', [1 : last_timepoint]' );
    fprintf( file_A_S{kk} , '\n' );
end
   
file_A_L = fopen( [ gen_Fold{2,1} '/Summary_Len_' terminator '.txt']  , 'w+');
fprintf( file_A_L , 'Generation\tTrack ID\tTime Birth\tFrame-to-Min\t' );
fprintf( file_A_L , '%f\t', [1 : last_timepoint]' );
fprintf( file_A_L , '\n' );

file_A_A = fopen( [ gen_Fold{2,1} '/Summary_Area_' terminator '.txt'] , 'w+');
fprintf( file_A_A , 'Generation\tTrack ID\tTime Birth\tFrame-to-Min\t' );
fprintf( file_A_A , '%f\t', [1 : last_timepoint]' );
fprintf( file_A_A , '\n' );   
if APP_opt.t3_choose_Summary_Plot == 1
    Path_txt = {} ; 
end
    

% ***** PLOT ALL GENERATIONS ************************************************************************      
for gg = 1 : size(gen_List,2)
    
    % Take the gg-th generation
    clone_List = gen_List{1, gg};
    
    % --- SAVE_step_1 -----
    % At each generation initialize tree_Values to collect all values 
    % of data plotted [cell array in which each element is a channel; 
    % for ratios we collect their values and not individual channels]
    if APP_opt.t3_choose_Save_txt == 1
        for kk = SubLine       
            tree_Values{kk} = cell(4, size(clone_List,2));    
        end
    end
    
    % For each individual gg-th generation create a RGB colormap for each channel
    if APP_opt.t3_choose_NormColormap == 2
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
    end
    
    

    % ***** PLOT ALL CLONES of gg-th generation ***************************
    % Make figure visible of invisible according to user choice
    if APP_opt.t3_choose_Display_Plot == 0
        h1f = figure('Position', [100 300 1200 600], 'Visible','off'); 
    elseif APP_opt.t3_choose_Display_Plot == 1
        h1f = figure('Position', [100 300 1200 600]);    
    end
    hold on;       
    ax = gca;
    % In independent clones plotting, there are as many as the number of tracks
    y_H_Row = [1:size(clone_List,2)]; 
    

    for cc = 1 : size(clone_List, 2)   
    % yrow = store the y position where the current cc-th cell is plotted
    %        (in independent clone-plot this is just cc)
        yrow = cc ;

        % ff = only one counter for "time", since we plot each clone starting 
        % from x=0, array position and time (x-axis) coincides.
        for ff = 1 : size(clone_List{1,cc}, 2)
            if ff > LAST_Frm;     LAST_Frm = ff;     end

            % If choosen, plot HUD: cell poles orientation and cell_ID number
            if ff == 1      % at frame of birth;
               fnc_Plot_HUD( clone_List{1,cc}{1}, ff, yrow, Dst_pl, hg_pl , LinW)
            end

            % Define the coordinates of the "box" square polygon and find the
            % signal value of the desired subcellular compartment to plot.
            [Xs, Ys, Value] = Find_Value_XYSquare( clone_List{1,cc}{1,ff}, ff, yrow, Dst_pl, hg_pl, ChNum , SubLine, ...
                                                   APP_opt.t3_choose_ChannelMode, APP_opt.t3_PlotOpt_B, APP_opt.t3_PlotOpt_C, APP_opt.algorithm );

            % Now, evaluate the RGB_Color for each "Line" that we wish to plot for each clone
            for kk = SubLine            
                cRGB  = Assign_Value_RGB( Value{kk} , RGB_range{kk} );
                % Fill squares with evaluated RGB_Color
                fill(Xs{kk}, Ys{kk}, [cRGB(1), cRGB(2), cRGB(3)], 'LineStyle', 'none');
            end

            % --- Extra Info ---- Place dots over cell_line to represent cell's
            % eccentricity value [greys scale: black = 0/circle; white = 1/ellipse)]
            if APP_opt.t3_display_ExI == 1
                eXs = [clone_List{1,cc}{1,ff}.mesh(:,1); flipud(clone_List{1,cc}{1,ff}.mesh(:,3))];
                eYs = [clone_List{1,cc}{1,ff}.mesh(:,2); flipud(clone_List{1,cc}{1,ff}.mesh(:,4))];        
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

            % --- SAVE_step_2 ---- Collect plotted value of all clones
            if APP_opt.t3_choose_Save_txt == 1 
                for kk = SubLine            
                    tree_Values{kk}{1, cc} = [tree_Values{kk}{1, cc} , Value{kk}] ;                            
                    % [px] Cell axial-length       
                    tree_Values{kk}{2, cc} = [tree_Values{kk}{2, cc}, clone_List{1,cc}{1,ff}.geom.length] ;    
                    % [px] Cell area
                    tree_Values{kk}{3, cc} = [tree_Values{kk}{3, cc}, clone_List{1,cc}{1,ff}.geom.area] ;
                    % y-position in plot  
                    tree_Values{kk}{4, cc} = yrow*Dst_pl ;  
                end
            end

        end % ff

        % Draw two horizontal lines that bound the plotted data for the clone-Line
        plot( [0, size(clone_List{1,cc},2) ], [yrow*Dst_pl,  yrow*Dst_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
        plot( [0, size(clone_List{1,cc},2) ], [(yrow*Dst_pl)+hg_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);

        % Draw vertical line at the STRAT and END, to enclose all end_Lines
        plot( [0, 0], [yrow*Dst_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
        plot( [ff, ff], [yrow*Dst_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
        
        % If necessary draw a line to separate the two channels plot-lines
        if APP_opt.t3_choose_ChannelMode == 3
            plot( [ff, ff-size(clone_List{1,cc},2)], [(yrow*Dst_pl)+(hg_pl*1/2),  (yrow*Dst_pl)+(hg_pl*1/2)], 'LineWidth', LinW/2, 'Color', Color_BorderLines);
        end
    
    end % cc
    
    
    
    
    %% *** SAVE the FIGURE and DATA *******************************************

    % --- FINISH layout and SAVE Plot -----------------------------------------

    % --- FINISH layout and SAVE Plot -----------------------------------------
    h1f.Color = [1 1 1];
    ax.TickDir = 'out';                 ax.YTick = [] ;                 
    ax.XColor = Color_BorderLines ;     ax.YColor = 'none' ;
    ax.LineWidth = LinW ;               ax.FontSize = 15 ;
    % Create a little extra space at end of x axis, where we will place the
    % colorbar legends (simply by extending of  MaxX+Delta_X)
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
    
    
    if APP_opt.t3_choose_AutoSave_Plot ~= 0         
        filename_plot = [ gen_Fold{2,gg} '/' gen_Fold{3,gg} '/' 'G' num2str(gg) '_' terminator ];
        if     APP_opt.t3_choose_AutoSave_Plot == 1;       print(h1f, filename_plot ,'-dpdf');
        elseif APP_opt.t3_choose_AutoSave_Plot == 2;       print(h1f, filename_plot , '-r300', '-dtiffn');
        elseif APP_opt.t3_choose_AutoSave_Plot == 3;       print(h1f, filename_plot ,'-dsvg');
        end     
    end
    
    
        
    % --- SAVE_step_3 ---- Save in .txt file ----------------------------------
    if APP_opt.t3_choose_Save_txt == 1 
        
        % Reorganize data to have all have same "length" and data at the correct "column" position.
        for kk = SubLine
            for jj = 1 : size(clone_List,2)
              for rr = [1,2,3]
                  tree_Values{kk}{rr,jj} = [ tree_Values{kk}{rr,jj}, repmat(-1, 1, last_timepoint - length( tree_Values{kk}{rr,jj})) ];
              end
            end
            
            % Store path to the saved .txt file for each generation, so that we
            % can use them to create summary plots
            if APP_opt.t3_choose_Summary_Plot == 1
                Path_txt{kk}{gg} = [ gen_Fold{2,gg} '/' gen_Fold{3,gg} '/Signal_CH', num2str(ChNum(kk)) ,'_S_G' num2str(gg) '_' terminator '.txt'] ;
            end
            
            % ----- 3B ----- Save ONE GENERATION in a .txt file
            % ---- Fluorescent Signal Data ------------------------------------
            file_O_S = fopen([ gen_Fold{2,gg} '/' gen_Fold{3,gg} '/Signal_CH', num2str(ChNum(kk)) ,'_S_G' num2str(gg) '_' terminator '.txt'], 'w+');       
            % First row is intestation:
            fprintf( file_O_S , 'Generation\tTrack ID\tTime Birth\tFrame-to-Min\t' );
            fprintf( file_O_S , '%f\t', [1 : last_timepoint]' );
            fprintf( file_O_S , '\n' );
            for ii = size(tree_Values{kk},2) :-1: 1            % in reverse order organizes data in the same way as it was plotted
                fprintf( file_O_S , '%f\t', gg );
                fprintf( file_O_S , '%s\t', num2str(clone_List{1,ii}{1}.ID_ManualTrack) );
                fprintf( file_O_S , '%f\t', clone_List{1,ii}{1}.frame_birth );
                fprintf( file_O_S , '%f\t', 1/ fpmRate );
                fprintf( file_O_S , '%f\t', tree_Values{kk}{2,ii} );
                fprintf( file_O_S , '\n' );
            end
            fclose(file_O_S) ;
        end % kk
        

        % ---- Cell axial-length Data -------------------------------------
        file_O_L = fopen([ gen_Fold{2,gg} '/' gen_Fold{3,gg} '/' 'Data_L_G' num2str(gg) '_' terminator '.txt'], 'w+');       
        % First row is intestation:
        fprintf( file_O_L , 'Generation\tTrack ID\tTime Birth\tFrame-to-Min\t' );
        fprintf( file_O_L , '%f\t', [1 : last_timepoint]' );
        fprintf( file_O_L , '\n' );
        for ii = size(tree_Values{1},2) :-1: 1            % in reverse order organizes data in the same way as it was plotted
            fprintf( file_O_L , '%f\t', gg );
            fprintf( file_O_L , '%s\t', num2str(clone_List{1,ii}{1}.ID_ManualTrack) );
            fprintf( file_O_L , '%f\t', clone_List{1,ii}{1}.frame_birth );
            fprintf( file_O_L , '%f\t', 1/ fpmRate );
            fprintf( file_O_L , '%f\t', tree_Values{1}{3,ii} );
            fprintf( file_O_L , '\n' );
        end
        fclose(file_O_L) ;

        % ---- Cell area Data ---------------------------------------------
        file_O_A = fopen([ gen_Fold{2,gg} '/' gen_Fold{3,gg} '/' 'Data_A_G' num2str(gg) '_' terminator '.txt'], 'w+');       
        % First row is intestation:
        fprintf( file_O_A , 'Generation\tTrack ID\tTime Birth\tFrame-to-Min\t' );
        fprintf( file_O_A , '%f\t', [1 : last_timepoint]' );
        fprintf( file_O_A , '\n' );
        for ii = size(tree_Values{1},2) :-1: 1            % in reverse order organizes data in the same way as it was plotted
            fprintf( file_O_A , '%f\t', gg );
            fprintf( file_O_A , '%s\t', num2str(clone_List{1,ii}{1}.ID_ManualTrack) );
            fprintf( file_O_A , '%f\t', clone_List{1,ii}{1}.frame_birth );
            fprintf( file_O_A , '%f\t', 1/ fpmRate );
            fprintf( file_O_A , '%f\t', tree_Values{1}{1,ii} );
            fprintf( file_O_A , '\n' );
        end
        fclose(file_O_A) ;


        % ----- 3A ----- Save ALL GENERATIONS in one .txt file
        for ii = size(tree_Values{kk},2) :-1: 1             % in reverse order organizes data in the same way as it was plotted
            for kk = SubLine        % For 2 channels use separate files
                fprintf( file_A_S{kk} , '%f\t', gg );
                fprintf( file_A_S{kk} , '%s\t', num2str(clone_List{1,ii}{1}.ID_ManualTrack) );
                fprintf( file_A_S{kk} , '%f\t', clone_List{1,ii}{1}.frame_birth );
                fprintf( file_A_S{kk} , '%f\t', 1/ fpmRate );
                fprintf( file_A_S{kk} , '%f\t', tree_Values{kk}{1,ii} );
                fprintf( file_A_S{kk} , '\n' );
            end

            fprintf( file_A_L , '%f\t', gg );
            fprintf( file_A_L , '%s\t', num2str(clone_List{1,ii}{1}.ID_ManualTrack) );
            fprintf( file_A_L , '%f\t', clone_List{1,ii}{1}.frame_birth );
            fprintf( file_A_L , '%f\t', 1/ fpmRate );
            fprintf( file_A_L , '%f\t', tree_Values{1}{2,ii} );
            fprintf( file_A_L , '\n' );

            fprintf( file_A_A , '%f\t', gg );
            fprintf( file_A_A , '%s\t', num2str(clone_List{1,ii}{1}.ID_ManualTrack) );
            fprintf( file_A_A , '%f\t', clone_List{1,ii}{1}.frame_birth );
            fprintf( file_A_A , '%f\t', 1/ fpmRate );
            fprintf( file_A_A , '%f\t', tree_Values{1}{3,ii} );
            fprintf( file_A_A , '\n' );
        end 


    end % Save .txt

end % GG generation

% Close all summary files
for kk = SubLine
    fclose( file_A_S{kk} );
end
fclose(file_A_L) ; 
fclose(file_A_A) ; 

if APP_opt.t3_choose_Summary_Plot == 1
    Plot_InterGen_Summary(Path_txt, RGB_range);
end
                
end % MAIN fnc





