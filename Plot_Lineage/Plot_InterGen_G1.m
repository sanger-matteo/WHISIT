function Plot_InterGen_G1
% Here we Plot the whole genalogy tree for a single lineage
%
%
%
%
%

% OLD_pole = Mask_PL_1
% NEW_pole = Mask_PL_2

global APP_opt ;

% ----- INITIALIZE some general variables and the RGB_range -----------------------------------
hg_pl = APP_opt.t3_height_PlotLine ;         % heights of the plot end lines
Dst_pl = APP_opt.t3_Dist_btw_PlotLines ;     % distances between plot end lines 
LinW =  APP_opt.t3_BorderLineWidth ;         % LineWidth of line connecting ancestor-to-daughters
Color_BorderLines = [.3 .3 .3]; 
LAST_Frm = 1;                                % highest value frame, use to find max x-axis length
fpmRate_CH1 = APP_opt.t3_fpm_CH1 ;

% Load the color LUT. LUT folder should be in same folder of the scripts Plot_Lineage
pathParts_LUT = strsplit(mfilename('fullpath'), {'/','\'} ) ;        % mfilename = take path of currently running script.
path_LUT = fullfile(pathParts_LUT{1,1:end-1},'\') ;                  % fullfile  = build full filename from string parts
filename_LUT = [path_LUT '/LUT/' 'LUT_' APP_opt.t3_ColorMap_LUT '.txt'] ;
RGB_range = textread(filename_LUT) ;

% Reorganize all the clone_list in the folder given by the user into
% gen_List, where each element stores all the clones of x-th generation
[gen_List, gen_Fold] = Create_GenList ;


% ----- INITIALIZE the summary file containing plotted data of all generations ----------------
% Based on analysis performed, add a terminal identifier at the end of filename
terminator = [ '2_' num2str(APP_opt.t3_PlotOpt_B) '_' num2str(APP_opt.t3_PlotOpt_C)];
% Find the latest time point for entire clone list, necessary to save data
% in .txt file by ensuring that each clone data has same length.
max_fr_last = -1 ;
for aa = 1 : size(gen_List,2)
    for bb = 1 : size(gen_List{1,aa},2)
        if max_fr_last <= length(gen_List{1,aa}{1,bb})
           max_fr_last = length(gen_List{1,aa}{1,bb});
        end
    end
end

% First row intestation is necessary only once for File_A
file_A_S = fopen( [ gen_Fold{2,1} '/Summary_Sig_' terminator '.txt'] , 'w+');
fprintf( file_A_S , 'Generation\tTrack ID\tTime Birth\tFrame-to-Min\t' );
fprintf( file_A_S , '%f\t', [1 : max_fr_last]' );
fprintf( file_A_S , '\n' );

file_A_L = fopen( [ gen_Fold{2,1} '/Summary_Len_' terminator '.txt'] , 'w+');
fprintf( file_A_L , 'Generation\tTrack ID\tTime Birth\tFrame-to-Min\t' );
fprintf( file_A_L , '%f\t', [1 : max_fr_last]' );
fprintf( file_A_L , '\n' );

file_A_A = fopen( [ gen_Fold{2,1} '/Summary_Area_' terminator '.txt'] , 'w+');
fprintf( file_A_A , 'Generation\tTrack ID\tTime Birth\tFrame-to-Min\t' );
fprintf( file_A_A , '%f\t', [1 : max_fr_last]' );
fprintf( file_A_A , '\n' );   

    

% ***** PLOT ALL GENERATIONS ************************************************************************      
for gg = 1 : size(gen_List,2)
    
    clone_List = gen_List{1, gg};

    % Color Range Setup and determining the Extreme Values
    [IC_V , Rel_V, Distr_Vals] = Eval_RawValue_Extrema(clone_List);
    V_Min = Rel_V(1);     % 1 --> min;
    V_Max = Rel_V(2);     % 2 --> max;

    Value_Range = strsplit(APP_opt.t3_Value_Range, ':') ;
    [nmin, status] = str2num(Value_Range{1}) ;    % if status == 0, conversion was not successful  
    if status ~= 0 ;         MinRange = nmin ;
    else ;                   MinRange = Rel_V(1);
    end
    [nmax, status] = str2num(Value_Range{2}) ;    % if status == 0, conversion was not successful  
    if status ~= 0 ;         MaxRange = nmax ;
    else ;                   MaxRange = Rel_V(2);
    end
    Ext_Rng = 0.0 ;          % Extend Range value of X% in both directions
    t_range = linspace(MinRange, MaxRange, size(RGB_range,1)+ (size(RGB_range,1)*2*Ext_Rng) ) ;
    RGB_range(:,4) = t_range( 1+ (size(RGB_range,1)*Ext_Rng)  :  end- (size(RGB_range,1)*Ext_Rng) )';


    % ----- PLOT all clones in gg-th generation ---------------------------
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

    % --- SAVE_step_1 ---- Initialize
    if APP_opt.t3_choose_Save_txt == 1
        % Collect value of data plotted in temporary array
        tree_Values = cell(4,size(clone_List,2));       
    end

    % ----- Plot all clones independently
    for cc = 1 : length(clone_List)    
        % yrow = store the y position where the current cc-th cell is plotted
        %        (in independent clone-plot this is just cc)
        yrow = cc ;

        % ff = only one counter for "time", which is the same as array position
        for ff = 1 : length(clone_List{1,cc})        
            if ff > LAST_Frm;     LAST_Frm = ff;     end

            % If choosen, plot HUD: cell poles orientation and cell_ID number
            if ff == 1      % at frame of birth;
               fnc_Plot_HUD( clone_List{1,cc}{1}, ff, yrow, Dst_pl, hg_pl , LinW)
            end

            % Create the square for plotting the value at time point ff
            Xs = [ ff-1, ff, ff, ff-1, ff-1 ];
            Ys = [ yrow*Dst_pl,  yrow*Dst_pl, (yrow*Dst_pl)+hg_pl, (yrow*Dst_pl)+hg_pl, yrow*Dst_pl ];

            % Evaluate the Value and the RGB_Color to apply inside the square
            switch APP_opt.t3_PlotOpt_C
                case 1      % Whole Cell
                    Value = mean(clone_List{1,cc}{ff}.CH1.IC(  clone_List{1,cc}{ff}.Mask_wCell )) ; 
                    [Value, cRGB] =  Assign_Value_RGB( Value , V_Min, V_Max, RGB_range );
                case 2      % Membrane
                    Value = mean(clone_List{1,cc}{ff}.CH1.IC(  clone_List{1,cc}{ff}.Mask_Memb )) ; 
                    [Value, cRGB] =  Assign_Value_RGB( Value , V_Min, V_Max, RGB_range ); 
                case 3      % Old_PL                   
                    Value = mean(clone_List{1,cc}{ff}.CH1.IC(  clone_List{1,cc}{ff}.CH1.Mask_PL_1 )) ; 
                    [Value, cRGB] =  Assign_Value_RGB( Value , V_Min, V_Max, RGB_range ); 
                case 4      % New_PL                    
                    Value = mean(clone_List{1,cc}{ff}.CH1.IC(  clone_List{1,cc}{ff}.CH1.Mask_PL_2 )) ; 
                    [Value, cRGB] =  Assign_Value_RGB( Value , V_Min, V_Max, RGB_range ); 
            end

            % Fill squares with evaluated RGB_Color
            fill(Xs,Ys, [cRGB(1), cRGB(2), cRGB(3)], 'LineStyle', 'none');

            % --- Extra Info ---- Place dots over cell_line to represent cell's
            % eccentricity value [greys scale: black = 0/circle; white = 1/ellipse)]
            if APP_opt.t3_display_ExI == 1
                eXs = [clone_List{1,cc}{ff}.mesh(:,1); flipud(clone_List{1,cc}{ff}.mesh(:,3))];
                eYs = [clone_List{1,cc}{ff}.mesh(:,2); flipud(clone_List{1,cc}{ff}.mesh(:,4))];        
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
                tree_Values{1, cc} = [tree_Values{1, cc} , Value] ;                                            
                % [px] Cell axial-length       
                tree_Values{2, cc} = [tree_Values{2, cc}, clone_List{1,cc}{ff}.geom.length] ;    
                % [px] Cell area
                tree_Values{3, cc} = [tree_Values{3, cc}, clone_List{1,cc}{ff}.geom.area] ;
                % y-position in plot  
                tree_Values{4, cc} = yrow*Dst_pl ;
            end

        end % ff

        % Draw two horizontal lines that bound the plotted data for the clone-Line
        plot( [0, length(clone_List{1,cc}) ], [yrow*Dst_pl,  yrow*Dst_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
        plot( [0, length(clone_List{1,cc}) ], [(yrow*Dst_pl)+hg_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);

        % Draw vertical line at the STRAT and END, to enclose all end_Lines
        plot( [0, 0], [yrow*Dst_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
        plot( [ff, ff], [yrow*Dst_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);

    end % cc

    % --- FINISH layout of plot and RENDER the figure ------
    h1f.Color = [1 1 1];
    ax.TickDir = 'out';                 ax.YTick = [] ;                 
    ax.XColor = Color_BorderLines ;     ax.YColor = [1 1 1] ;
    ax.LineWidth = LinW ;               ax.FontSize = 15 ;
    ax.XTickLabel = linspace(ax.XLim(1),ax.XLim(2),length(ax.XTick)) .* fpmRate_CH1 ;
    h1f.Renderer = 'Painters';
    
    if APP_opt.t3_choose_AutoSave_Plot ~= 0 
      if     APP_opt.t3_choose_AutoSave_Plot == 1
          print(h1f, [ gen_Fold{2,gg} '/' gen_Fold{3,gg} '/' 'G' num2str(gg) '_' terminator ] ,'-dpdf');
      elseif APP_opt.t3_choose_AutoSave_Plot == 2
          print(h1f, [ gen_Fold{2,gg} '/' gen_Fold{3,gg} '/' 'G' num2str(gg) '_' terminator ] ,'-r300', '-dtiffn');
      end      
    end
    
    
    
    % --- SAVE_step_3 ---- Save in .txt file ----------------------------------
    if APP_opt.t3_choose_Save_txt == 1 
        
        % Reorganize data to have all have same "length" and data at the correct "column" position.
        for jj = 1 : size(clone_List,2)
            for rr = [1,2,3]
                tree_Values{rr,jj} = [ tree_Values{rr,jj}, repmat(0, 1, max_fr_last - length( tree_Values{rr,jj})) ];
            end
        end
        
        % ----- 3B ----- Save ONE GENERATION in a .txt file
        % ---- Fluorescent Signal Data -----
        file_O_S = fopen([ gen_Fold{2,gg} '/' gen_Fold{3,gg} '/' 'Data_S_G' num2str(gg) '_' terminator '.txt'], 'w+');       
        % First row is intestation:
        fprintf( file_O_S , 'Generation\tTrack ID\tTime Birth\tFrame-to-Min\t' );
        fprintf( file_O_S , '%f\t', [1 : max_fr_last]' );
        fprintf( file_O_S , '\n' );
        for ii = size(tree_Values,2) :-1: 1            % in reverse order organizes data in the same way as it was plotted
            fprintf( file_O_S , '%f\t', gg );
            fprintf( file_O_S , '%s\t', num2str(clone_List{1,ii}{1}.ID_ManualTrack) );
            fprintf( file_O_S , '%f\t', clone_List{1,ii}{1}.fr_birth );
            fprintf( file_O_S , '%f\t', fpmRate_CH1 );
            fprintf( file_O_S , '%f\t', tree_Values{2,ii} );
            fprintf( file_O_S , '\n' );
        end
        fclose(file_O_S) ;   
                
        % ---- Cell axial-length Data -----
        file_O_L = fopen([ gen_Fold{2,gg} '/' gen_Fold{3,gg} '/' 'Data_L_G' num2str(gg) '_' terminator '.txt'], 'w+');       
        % First row is intestation:
        fprintf( file_O_L , 'Generation\tTrack ID\tTime Birth\tFrame-to-Min\t' );
        fprintf( file_O_L , '%f\t', [1 : max_fr_last]' );
        fprintf( file_O_L , '\n' );
        for ii = size(tree_Values,2) :-1: 1            % in reverse order organizes data in the same way as it was plotted
            fprintf( file_O_L , '%f\t', gg );
            fprintf( file_O_L , '%s\t', num2str(clone_List{1,ii}{1}.ID_ManualTrack) );
            fprintf( file_O_L , '%f\t', clone_List{1,ii}{1}.fr_birth );
            fprintf( file_O_L , '%f\t', fpmRate_CH1 );
            fprintf( file_O_L , '%f\t', tree_Values{3,ii} );
            fprintf( file_O_L , '\n' );
        end
        fclose(file_O_L) ;
        
        % ---- Cell area Data -----
        file_O_A = fopen([ gen_Fold{2,gg} '/' gen_Fold{3,gg} '/' 'Data_A_G' num2str(gg) '_' terminator '.txt'], 'w+');       
        % First row is intestation:
        fprintf( file_O_A , 'Generation\tTrack ID\tTime Birth\tFrame-to-Min\t' );
        fprintf( file_O_A , '%f\t', [1 : max_fr_last]' );
        fprintf( file_O_A , '\n' );
        for ii = size(tree_Values,2) :-1: 1            % in reverse order organizes data in the same way as it was plotted
            fprintf( file_O_A , '%f\t', gg );
            fprintf( file_O_A , '%s\t', num2str(clone_List{1,ii}{1}.ID_ManualTrack) );
            fprintf( file_O_A , '%f\t', clone_List{1,ii}{1}.fr_birth );
            fprintf( file_O_A , '%f\t', fpmRate_CH1 );
            fprintf( file_O_A , '%f\t', tree_Values{1,ii} );
            fprintf( file_O_A , '\n' );
        end
        fclose(file_O_A) ;
        
                
        % ----- 3A ----- Save ALL GENERATIONS in one .txt file
        for ii = size(tree_Values,2) :-1: 1             % in reverse order organizes data in the same way as it was plotted
            fprintf( file_A_S , '%f\t', gg );
            fprintf( file_A_S , '%s\t', num2str(clone_List{1,ii}{1}.ID_ManualTrack) );
            fprintf( file_A_S , '%f\t', clone_List{1,ii}{1}.fr_birth );
            fprintf( file_A_S , '%f\t', fpmRate_CH1 );
            fprintf( file_A_S , '%f\t', tree_Values{1,ii} );
            fprintf( file_A_S , '\n' );
            
            fprintf( file_A_L , '%f\t', gg );
            fprintf( file_A_L , '%s\t', num2str(clone_List{1,ii}{1}.ID_ManualTrack) );
            fprintf( file_A_L , '%f\t', clone_List{1,ii}{1}.fr_birth );
            fprintf( file_A_L , '%f\t', fpmRate_CH1 );
            fprintf( file_A_L , '%f\t', tree_Values{2,ii} );
            fprintf( file_A_L , '\n' );
            
            fprintf( file_A_A , '%f\t', gg );
            fprintf( file_A_A , '%s\t', num2str(clone_List{1,ii}{1}.ID_ManualTrack) );
            fprintf( file_A_A , '%f\t', clone_List{1,ii}{1}.fr_birth );
            fprintf( file_A_A , '%f\t', fpmRate_CH1 );
            fprintf( file_A_A , '%f\t', tree_Values{3,ii} );
            fprintf( file_A_A , '\n' );
        end            
    end % Save .txt

end % GG generation

fclose(file_A_S) ; 
fclose(file_A_L) ; 
fclose(file_A_A) ; 

end % MAIN fnc







%% 
% *********************************************************************************************************
% -----> SCRIP-RELATED FUNCTIONS --------------------------------------------------------------------------
% *********************************************************************************************************

function [IC_V , Rel_V, Distr_Vals] = Eval_RawValue_Extrema(track_List)
% ----- Fluorescence VALUE SCALING ----------------------------------------
% ESTABLISH RELATIVE MIN MAX, ACCORDING TO...
% Find Min and Max in all fluorescence Frames
global APP_opt ;
I_max = -1 ;           I_min = 2*10^20 ;
R_avg_max = -1 ;       R_avg_min = 2*10^20 ;
Distr_Vals = [];
for cc = 1 : length(track_List)
    for ff = 1 : length(track_List{1,cc})
        Value = [];     
        
        if I_max < max( max(track_List{1,cc}{ff}.CH1.IC))
            I_max = max( max(track_List{1,cc}{ff}.CH1.IC) );
        end
        if I_min > min( min(track_List{1,cc}{ff}.CH1.IC))
            I_min = min( min(track_List{1,cc}{ff}.CH1.IC));
        end 

        % Evaluate the Value and the RGB_Color to apply inside the square
        switch APP_opt.t3_PlotOpt_C
            case 1      % Whole Cell
                Value = mean(track_List{1,cc}{ff}.CH1.IC(  track_List{1,cc}{ff}.Mask_wCell )) ; 
            case 2      % Membrane
                Value = mean(track_List{1,cc}{ff}.CH1.IC(  track_List{1,cc}{ff}.Mask_Memb )) ; 
            case 3      % Old_PL                   
                Value = mean(track_List{1,cc}{ff}.CH1.IC(  track_List{1,cc}{ff}.CH1.Mask_PL_1 )) ; 
            case 4      % New_PL                    
                Value = mean(track_List{1,cc}{ff}.CH1.IC(  track_List{1,cc}{ff}.CH1.Mask_PL_2 )) ; 
        end
            
        % Relative min-max to the avg values Cyto and Poles
        if R_avg_max < Value
            R_avg_max = Value;       end
        if R_avg_min > Value
            R_avg_min = Value;       end
        
        Distr_Vals = [Distr_Vals , Value] ;            

    end
end
IC_V = [I_min, I_max];
Rel_V = [R_avg_min,R_avg_max];

end



function [Value RGB] = Assign_Value_RGB( Value, V_min, V_max, RGB_range ) 
%     Value = Value - V_min ;
    % use abs() because half will be negative values: we search for the one
    % where the difference is closest to zero.
    diff = (abs(RGB_range(:,4) - Value));
    idx_min = find(diff == min(diff));
    if length(idx_min) >= 1           % in the rare exception when there are two ...    
        idx_min = idx_min(1);         % equal minima diff, take just the first
    end
    R = RGB_range(idx_min,1);      G = RGB_range(idx_min,2);    B = RGB_range(idx_min,3); 
    if Value > RGB_range(end,4) 
        R = RGB_range(end,1);      G = RGB_range(end,2);        B = RGB_range(end,3);        
    elseif Value < RGB_range(1,4)
        R = RGB_range(1,1);        G = RGB_range(1,2);          B = RGB_range(1,3);
    end
    RGB(1) = R/256;     RGB(2) = G/256;     RGB(3) = B/256;
end


