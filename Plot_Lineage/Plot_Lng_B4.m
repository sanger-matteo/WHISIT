function Plot_Lng_B4
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

% OLD_pole = Mask_PL_1
% NEW_pole = Mask_PL_2

global APP_opt ;

load([APP_opt.t3_path_cloneList, APP_opt.t3_file_cloneList]);   % load clone_List

hg_pl = APP_opt.t3_height_PlotLine ;         % heights of the plot end lines
Dst_pl = APP_opt.t3_Dist_btw_PlotLines ;     % distances between plot end lines 
LinW =  APP_opt.t3_BorderLineWidth ;         % LineWidth of line connecting ancestor-to-daughters
Color_BorderLines = [.3 .3 .3]; 
LAST_Frm = 1;                                % highest value frame, use to find max x-axis length
fpmRate_CH1 = APP_opt.t3_fpm_CH1 ;

% Load the color LUT. LUT folder should be in same folder of the scripts Plot_Lineage
pathParts_LUT = strsplit(mfilename('fullpath'), {'/','\'} ) ;   % mfilename = take path of currently running script.
path_LUT = fullfile(pathParts_LUT{1,1:end-1},'\') ;                  % fullfile  = build full filename from string parts
filename_LUT = [path_LUT '/LUT/' 'LUT_' APP_opt.t3_ColorMap_LUT '.txt'] ;
RGB_range = textread(filename_LUT) ;

% Create ordered list of the clones (cs_All) and the end-of-line clone (cS_EndLine) 
% This will allow to plot in the correct order, from the "latest" to the founder cell
[cS_EndLine, cN_EndLine, cS_All] = fnc_Organize_cloneList(clone_List) ;

% Color Range Setup and determining the Extreme Values
[IC_V , Abs_V, Distr_Vals, Seg_IDX, Seg_VAL,   Rel_cl_max, Rel_cl_min] = ...
    Eval_RawValue_Segment(clone_List, APP_opt.t3_B4_LenStep);
% V_Min = Abs_V(1);     % 1 --> min;
% V_Max = Abs_V(2);     % 2 --> max;

Value_Range = strsplit(APP_opt.t3_Value_Range, ':') ;
[nmin, status] = str2num(Value_Range{1}) ;    % if status == 0, conversion was not successful  
if status ~= 0 ;         MinRange = {nmin} ;
else ;                   MinRange = Rel_cl_min ;
end
[nmax, status] = str2num(Value_Range{2}) ;    % if status == 0, conversion was not successful  
if status ~= 0 ;         MaxRange = {nmax} ;
else ;                   MaxRange = Rel_cl_max ;
end

Ext_Rng = 0.0 ;     % Extend Range value of X in both directions
% Create absolute RGB scale that fit all clones
% t_range = linspace(MinRange, MaxRange, size(RGB_range,1)+ (size(RGB_range,1)*2*Ext_Rng) ) ;
% RGB_range(:,4) = t_range( 1+ (size(RGB_range,1)*Ext_Rng)  :  end- (size(RGB_range,1)*Ext_Rng) )';
% Create as many RGB scales as clones, each relative the specific clone
for cc = 1 : length(MaxRange)
   t_range = linspace(MinRange{cc}, MaxRange{cc}, size(RGB_range,1)+(size(RGB_range,1)*2*Ext_Rng) ) ;
   RGB_range(:,cc+3) = t_range( 1+ (size(RGB_range,1)*Ext_Rng)  :  end- (size(RGB_range,1)*Ext_Rng) )';
end

% find the Max_CellLen, the number of segments the longest cell has
Max_CellLen = 0;
for ii = 1:length(Seg_IDX)
    if Max_CellLen <  max(cellfun(@length, Seg_IDX{ii})) 
    Max_CellLen = max(cellfun(@length, Seg_IDX{ii}));      end
end
stp_Seg = Dst_pl / Max_CellLen ;


% ***** PLOT LINEAGE TREE ************************************************************************
if APP_opt.t3_PlotOpt_A == 1            % LINEAGE TREE
    
% Make figure visible of invisible according to user choice
if APP_opt.t3_choose_Display_Plot == 0
    h1f = figure('Position', [100 300 1200 600], 'Visible','off'); 
elseif APP_opt.t3_choose_Display_Plot == 1
    h1f = figure('Position', [100 300 1200 600]);    
end
hold on;       
ax = gca;
LAST_Frm = 1;                        % highest value frame, use to find max x-axis length
y_H_Row = [1:length(cN_EndLine)];   


%% --- STEP 1 --- Plot all End_Line clones ---------------------------------
% First we plot all end line clones, meaning all the clones that do not
% generate daughter cells.
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
    
    % % Plot all time points. We use two counters for "time":
    % - ff  = stores the position of each element (/frame) of the array
    %         clone_List{idx_EndL}: simply goes from 1 to last frame
    % - Frm = real time variable that plots the data extracted at clone_List{idx_EndL}{ff} 
    %         in the correct time point (x-axis) of the lineage tree
    for ff = 1 : length(clone_List{idx_EndL}) 
        seg_Xs = [];
        seg_Ys = [];
        Frm = ff + clone_List{idx_EndL}{1}.fr_birth;
        if Frm > LAST_Frm;     LAST_Frm = Frm;     end

        % Create the segments that define a "cell"
        for bb = 0 : length(Seg_IDX{idx_EndL}{ff})
            seg_Xs(bb+1,:) = [ Frm-1, Frm, Frm, Frm-1, Frm-1 ];
            seg_Ys(bb+1,:) = [(yrow*Dst_pl)+(stp_Seg*bb)    , (yrow*Dst_pl)+(stp_Seg*bb),    ... % last point = first point
                              (yrow*Dst_pl)+(stp_Seg*(bb+1)), (yrow*Dst_pl)+(stp_Seg*(bb+1)),    (yrow*Dst_pl)+(stp_Seg*bb) ];
        end
        
        % Fill ecah SEGMENT with corresponding evaluated RGB_Color
        [Value, seg_RGB] =  Assign_Clones_RGB_value( Seg_VAL{idx_EndL}{ff} , RGB_range, idx_EndL );
        for bb = 1 : length(Seg_IDX{idx_EndL}{ff})-1
            fill(seg_Xs(bb,:),seg_Ys(bb,:), [seg_RGB(bb,1), seg_RGB(bb,2), seg_RGB(bb,3)], 'LineStyle', 'none');
        end
        
        % Draw two horizontal lines that bound the plotted data for the clone-Line
        plot( [Frm-1, Frm], seg_Ys(1,1:2), 'LineWidth', LinW, 'Color', Color_BorderLines);
        plot( [Frm-1, Frm], seg_Ys(length(Seg_IDX{idx_EndL}{ff})-1, 3:4), 'LineWidth', LinW, 'Color', Color_BorderLines);

        % place dots on top of plotted cell_line to represent epsilon value
        % in scale of grays (black: circle/stand up; white:ellipse)
        if APP_opt.t3_display_ExI == 1
            epsl = clone_List{idx_EndL}{ff}.geom.Epsilon ; 
            FaceCol = [epsl epsl epsl] ;            
            plot( Frm -0.5, (yrow*Dst_pl +hg_pl) +hg_pl/3 , '.', 'MarkerSize', 8, 'LineWidth', 0.2, 'Color', FaceCol);
        end

        % If choosen, plot HUD: cell poles orientation and cell_ID number
        if ff == 1      % at frame of birth;
           fnc_Plot_HUD( clone_List{idx_EndL}{1}, Frm, yrow, Dst_pl, hg_pl , LinW)
           plot( [Frm-1, Frm-1], [seg_Ys(1,1),seg_Ys(length(Seg_IDX{idx_EndL}{ff})-1,4)], 'LineWidth', LinW, 'Color', Color_BorderLines);
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
            plot( Frm -0.5, (yrow*Dst_pl) -hg_pl/3 , '.', 'MarkerSize', 8, 'LineWidth', 0.2, 'Color', El_Ec_Col);
        end
        
    end % ff
    
    % draw vertical line at the end, to enclose all end_Lines
    plot( [Frm, Frm], [seg_Ys(1,1),seg_Ys(length(Seg_IDX{idx_EndL}{ff})-1,4)], 'LineWidth', LinW, 'Color', Color_BorderLines);

end % cc




%% --- STEP 2 --- Plot all remaining clones ---------------------------------
% This step will plot all the remaining clones. The algorithm find clone
% couples using cS_EndLine: clone couple that share same ID root (.XYZ ).
% From this it is possible to identify the common ancestror and plot the
% such cell (.XY ). The cS_EndLine list is updated and the algorithm
% proceed with the next couple, and so on until we end up to the progenitor
% cell itself. [This is why we must use while-loop for this step]

while length(cS_EndLine)>=2
    % cc is the counter of the cell position in arrays, such as cS_EndLine and clone_List 
    cc = 0 ; 
    
    while cc < length(cS_EndLine)
        cc = cc+1;
        
        % Find the two twins: ID name string as well as index positionin array
        TWstr_A = cS_EndLine{cc}(2:end);
        if      TWstr_A(end) == '1';   TWstr_B = [cS_EndLine{cc}(2:end-1) '2'];
        elseif  TWstr_A(end) == '2';   TWstr_B = [cS_EndLine{cc}(2:end-1) '1'];
        end
        idx_TA =  find(strcmp(cS_EndLine, ['.' TWstr_A] ) == 1) ;
        idx_TB =  find(strcmp(cS_EndLine, ['.' TWstr_B] ) == 1) ;
        
        % Plot the common ancestor
        if ~isempty(idx_TA) && ~isempty(idx_TB)  
            % Find ancestor ID name and postion in array
            str_Anc = [APP_opt.t3_PlotOpt_NLineage '.' TWstr_A(1:end-1)] ;
            idx_Anc = find(strcmp(cS_All , str_Anc ) == 1);
            % yrow = store the y position for plotting the current cc-th clone
            yrow = (y_H_Row(idx_TA) + y_H_Row(idx_TB)) /2 ;

            % The last time point of the ancestor is the birth of the two twins
            idx_birth = clone_List{idx_Anc}{1}.fr_last ;
            % Plot vertical line at idx_birth, connecting the two daughter cells (twins) in lineage tree
            seg_Xs = [ idx_birth+1 , idx_birth+1 ];
            if y_H_Row(idx_TA) < y_H_Row(idx_TB)
                Ys = [ y_H_Row(idx_TA)*Dst_pl,  (y_H_Row(idx_TB)*Dst_pl)+hg_pl];
            elseif y_H_Row(idx_TA) > y_H_Row(idx_TB)
                Ys = [ y_H_Row(idx_TB)*Dst_pl,  (y_H_Row(idx_TA)*Dst_pl)+hg_pl];
            end
            plot(seg_Xs, Ys, '-', 'LineWidth', LinW, 'Color', Color_BorderLines); 

            % Plot all time points. We use two counters for "time":
            % - ff  = stores the position of each element (/frame) of the array
            %         clone_List{idx_EndL}: simply goes from 1 to last frame
            % - Frm = real time variable that plots the data extracted at clone_List{idx_EndL}{ff} 
            %         in the correct time point (x-axis) of the lineage tree
            for ff = 1 : length(clone_List{idx_Anc})
                seg_Xs = [];
                seg_Ys = [];
                Frm = ff + clone_List{idx_Anc}{1}.fr_birth;
                if Frm > LAST_Frm;     LAST_Frm = Frm;     end

                % Create the segments that define a "cell"
                for bb = 0 : length(Seg_IDX{idx_Anc}{ff})
                    seg_Xs(bb+1,:) = [ Frm-1, Frm, Frm, Frm-1, Frm-1 ];
                    seg_Ys(bb+1,:) = [(yrow*Dst_pl)+(stp_Seg*bb)    , (yrow*Dst_pl)+(stp_Seg*bb),    ... % last point = first point
                                     (yrow*Dst_pl)+(stp_Seg*(bb+1)) , (yrow*Dst_pl)+(stp_Seg*(bb+1)),    (yrow*Dst_pl)+(stp_Seg*bb) ];
                end
                
                % Fill ecah SEGMENT with corresponding evaluated RGB_Color
                [Value, seg_RGB] =  Assign_Clones_RGB_value( Seg_VAL{idx_Anc}{ff} , RGB_range, idx_Anc );
                for bb = 1 : length(Seg_IDX{idx_Anc}{ff})-1
                    fill(seg_Xs(bb,:),seg_Ys(bb,:), [seg_RGB(bb,1), seg_RGB(bb,2), seg_RGB(bb,3)], 'LineStyle', 'none');
                end
                
                % Draw two horizontal lines that bound the plotted data for the clone-Line
                plot( [Frm-1, Frm], seg_Ys(1,1:2), 'LineWidth', LinW, 'Color', Color_BorderLines);
                plot( [Frm-1, Frm], seg_Ys(length(Seg_IDX{idx_Anc}{ff})-1,3:4), 'LineWidth', LinW, 'Color', Color_BorderLines);

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
                    plot( Frm -0.5, (yrow*Dst_pl) -hg_pl/3 , '.', 'MarkerSize', 8, 'LineWidth', 0.2, 'Color', El_Ec_Col);
                end                

                % If choosen, plot HUD: cell poles orientation and cell_ID number
                if ff == 1      % at frame of birth;
                   fnc_Plot_HUD( clone_List{idx_Anc}{1}, Frm, yrow, Dst_pl, hg_pl , LinW)
                   plot( [Frm-1, Frm-1], [seg_Ys(1,1),seg_Ys(length(Seg_IDX{idx_Anc}{ff})-1,4)], 'LineWidth', LinW, 'Color', Color_BorderLines);
                end  
                
            end % ff
            
            % Draw vertical line at the end, to enclose all end_Lines
            plot( [Frm, Frm], [seg_Ys(1,1),seg_Ys(length(Seg_IDX{idx_Anc}{ff})-1,4)], 'LineWidth', LinW, 'Color', Color_BorderLines);
            
            % UPDATe cS_EndLine ----------------------------------------
            % Remove the ID of the couple we just plotted and...
            % (we do not need to update cN_End_Line)
            cS_EndLine{idx_TA} = '';
            cS_EndLine{idx_TB} = '';
            y_H_Row(idx_TA) = yrow;            
            y_H_Row(idx_TB) = [];
            % ... put their common ancestor ID_anc in the list, since now it a "end_line"                  
            [String] = strsplit(str_Anc, '.');
            cS_EndLine{idx_TA} = [ '.' String{2} ];
            % ... and we remove empty elements from the list
            cS_EndLine = cS_EndLine(~cellfun(@isempty, cS_EndLine)) ;
            
        end
    end  % while cc
end % while length

% Last idx_Anc is lineage progenitor: draw vertical line at the start, to enclose Plot_Line
plot( [clone_List{idx_Anc}{1}.fr_birth, clone_List{idx_Anc}{1}.fr_birth], [yrow*Dst_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);

% --- FINISH layout of plot and RENDER the figure ------
h1f.Color = [1 1 1];
ax.TickDir = 'out';                 ax.YTick = [] ;                 
ax.XColor = Color_BorderLines ;     ax.YColor = [1 1 1] ;
ax.LineWidth = LinW ;               ax.FontSize = 15 ;
ax.XTickLabel = linspace(ax.XLim(1),ax.XLim(2),length(ax.XTick)) .* fpmRate_CH1 ;
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
    end
end




% ------------------------------------------------------------------------------------------------
% ***** PLOT CLONES INDEPENDENTLY ****************************************************************
elseif APP_opt.t3_PlotOpt_A == 2            % if GUI Plot Type: '2 - Individual Clones' 
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

% ----- Plot all clones
for cc = 1 : length(clone_List)
    % yrow = store the y position where the current cc-th cell is plotted
    %        (in independent clone-plot this is just cc)
    yrow = cc ;
    
    % ff = only one counter for "time", which is the same as array position
    for ff = 1 : length(clone_List{cc})        
        if ff > LAST_Frm;     LAST_Frm = ff;     end

        % If choosen, plot HUD: cell poles orientation and cell_ID number
        if ff == 1      % at frame of birth;
           fnc_Plot_HUD( clone_List{cc}{1}, ff, yrow, Dst_pl, hg_pl , LinW)
        end

        seg_Xs = [ ff-1, ff, ff, ff-1, ff-1 ];
        Ys = [ yrow*Dst_pl,  yrow*Dst_pl, (yrow*Dst_pl)+hg_pl, (yrow*Dst_pl)+hg_pl, yrow*Dst_pl ];

        % Evaluate the Value and the RGB_Color to apply inside the square
        % and fill squares with evaluated RGB_Color
        switch APP_opt.t3_PlotOpt_C
            case 1      % Whole Cell
                Value = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.Mask_wCell )) ; 
                [Value, seg_RGB] =  Assign_Clones_RGB_value( Value , RGB_range );
            case 2      % Membrane
                Value = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.Mask_Memb )) ; 
                [Value, seg_RGB] =  Assign_Clones_RGB_value( Value , RGB_range ); 
            case 3      % Old_PL                   
                Value = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.CH1.Mask_PL_1 )) ; 
                [Value, seg_RGB] =  Assign_Clones_RGB_value( Value , RGB_range ); 
            case 4      % New_PL                    
                Value = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.CH1.Mask_PL_2 )) ; 
                [Value, seg_RGB] =  Assign_Clones_RGB_value( Value , RGB_range ); 
        end

        % fill squares with evaluated RGB_Color
        fill(seg_Xs,Ys, [seg_RGB(1), seg_RGB(2), seg_RGB(3)], 'LineStyle', 'none');
        
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
            plot( Frm -0.5, (yrow*Dst_pl) -hg_pl/3 , '.', 'MarkerSize', 8, 'LineWidth', 0.2, 'Color', El_Ec_Col);
        end        
       
        % draw horizontal lines that bound the squares of the cell_Line
        plot( [ff, ff-1], [yrow*Dst_pl,  yrow*Dst_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
        plot( [ff, ff-1], [(yrow*Dst_pl)+hg_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
    
    end % ff
    
    % draw vertical line at the STRAT and END, to enclose all end_Lines
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
    end
end

end % if Plot

end % Main fnc







%% 
% *********************************************************************************************************
% -----> SCRIP-RELATED FUNCTIONS --------------------------------------------------------------------------
% *********************************************************************************************************

function [IC_V , Abs_V,   Distr_Vals,   Seg_IDX, Seg_VAL,   Rel_cl_max, Rel_cl_min] = ...
         Eval_RawValue_Segment(clone_List , L_stp)
% ----- Fluorescence VALUE SCALING ----------------------------------------
% ESTABLISH RELATIVE MIN MAX, ACCORDING TO...
% Find Min and Max in all fluorescence Frames
I_max = -1 ;           I_min = 2*10^20 ;
A_avg_max = -1 ;       A_avg_min = 2*10^20 ;
Rel_cl_max = {};       Rel_cl_min = {};
Distr_Vals = [];
Seg_IDX = {};
for cc = 1 : length(clone_List)
  Rel_cl_max{cc} = -1;       Rel_cl_min{cc} = 2*10^20 ;
  for ff = 1 : length(clone_List{cc}) 
      
    % Absolute Min and Max of the cell IC fluorescence
    if I_max < max( max(clone_List{cc}{ff}.CH1.IC))
        I_max = max( max(clone_List{cc}{ff}.CH1.IC) );        end
    if I_min > min( min(clone_List{cc}{ff}.CH1.IC))
        I_min = min( min(clone_List{cc}{ff}.CH1.IC));        end 
    % Left and Right set of coordinates of the cell mesh
    X_1 = clone_List{cc}{ff}.R_mesh(:,1) ;        Y_1 = clone_List{cc}{ff}.R_mesh(:,2) ;
    X_2 = clone_List{cc}{ff}.R_mesh(:,3) ;        Y_2 = clone_List{cc}{ff}.R_mesh(:,4) ;

    % According to the size of segments' width, we have find the
    % indexes that define the boundaries of each segment.         
    % (Segment width can only go from 2 to 5.)
    idx =  L_stp+2 : L_stp : length(clone_List{cc}{ff}.R_mesh(:,1)) -(L_stp);
    Rest = mod(length(clone_List{cc}{ff}.R_mesh(:,1)) ,L_stp);
    % There is a Rest if meshList length is not multiple of step_size. 
    % In the following switch-block we arrange all segments intervals to be 
    % the same width and increasing or decresing the size of the segments 
    % in the very center of the cell if there is extra length. Moreover,  
    % since the mesh(1) and mesh(end) are the same point for left and right 
    % set of coordinates points, we make the two pole_block one set of 
    % coordinate longer than the rest  (i.e. blocks width in a cell can be:
    %  5 4 4 4 3 4 4 4 5   or   4 3 3 3 3 4   or   4 3 3 2 3 4 )

    switch L_stp
        case 2     % if L_stp == 2 --------------------------------------------
            switch Rest
                case 0
                    mid = round(length(idx)/2) ;
                    idx(mid) = idx(mid) +1;
                    idx(mid+1:end) = idx(mid+1:end) -1;
                    idx(mid:end-1) = idx(mid+1:end) ;
                    idx = [1, idx];     idx(end) = length(clone_List{cc}{ff}.R_mesh(:,1));
                case 1
                    idx = [1, idx, length(clone_List{cc}{ff}.R_mesh(:,1))]; 
            end

        case 3     % if L_stp == 3 --------------------------------------------
            switch Rest
                case 0
                    idx = [1, idx, length(clone_List{cc}{ff}.R_mesh(:,1)) ];
                case 1
                    mid = round(length(idx)/2) ;
                    idx(mid:end) = idx(mid:end) +1 ;
                    idx = [1, idx, length(clone_List{cc}{ff}.R_mesh(:,1))];
                case 2
                    mid = round(length(idx)/2) ;
                    idx(mid) = idx(mid) +1 ;
                    idx = [1, idx, length(clone_List{cc}{ff}.R_mesh(:,1))];
            end

        case 4     % if L_stp == 4 --------------------------------------------
            switch Rest         
                case 0
                    mid = round(length(idx)/2) ;
                    idx(mid+1:end) = idx(mid+1:end) +1;
                    idx = [1, idx, length(clone_List{cc}{ff}.R_mesh(:,1)) ];
                case 1
                    mid = round(length(idx)/2) ;
                    idx(mid:end)   = idx(mid:end)   +1;
                    idx(mid+1:end) = idx(mid+1:end) +1;
                    idx = [1, idx, length(clone_List{cc}{ff}.R_mesh(:,1))]; 
                case 2
                    mid = round(length(idx)/2) ;
                    idx(mid+1:end) = idx(mid+1:end) -1;
                    idx = [1, idx, length(clone_List{cc}{ff}.R_mesh(:,1))]; 
                case 3
                    idx = [1, idx, length(clone_List{cc}{ff}.R_mesh(:,1))]; 
            end
    end

    jj = 1;         Avg_pxVal = [];        Area_Segm = [];
    for ii = 1 : 1 : length(idx)-1
        x1 = X_1(idx(ii):idx(ii+1));
        x2 = X_2(idx(ii):idx(ii+1));
        y1 = Y_1(idx(ii):idx(ii+1));
        y2 = Y_2(idx(ii):idx(ii+1));
        Mask_Seg = roipoly( clone_List{cc}{ff}.CH1.IC, [x1;flipud(x2)],[y1;flipud(y2)] );

        Avg_pxVal(jj) = mean(clone_List{cc}{ff}.CH1.IC(Mask_Seg));
        Area_Segm(jj) = sum(sum(Mask_Seg)) ;
        % We take the average value normalized by the area
        Avg_SegVal(jj) = Avg_pxVal(jj);  %/Area_Segm(jj);

        % Relative min-max to the avg values Cyto and Poles
        if A_avg_max < Avg_SegVal(jj)
            A_avg_max = Avg_SegVal(jj);       end
        if A_avg_min > Avg_SegVal(jj)
            A_avg_min = Avg_SegVal(jj);       end
        Distr_Vals = [Distr_Vals , Avg_SegVal] ;   

        jj = jj+1;
    end       

    if Rel_cl_max{cc} < max(Avg_SegVal)
        Rel_cl_max{cc} = max(Avg_SegVal) ;       end
    if Rel_cl_min{cc} > min(Avg_SegVal)
        Rel_cl_min{cc} = min(Avg_SegVal) ;       end
    
    Seg_IDX{cc}{ff} = idx;    
    Seg_VAL{cc}{ff} = Avg_SegVal;
    
  end
end

IC_V = [I_min, I_max];
Abs_V = [A_avg_min,A_avg_max];

end



function [Value RGB] = Assign_Value_RGB( Value, V_min, V_max, RGB_range ) 
	for ii = 1 : length(Value)
        % use abs() because half will be negative values: we search for the one
        % where the difference is closest to zero.
        diff = (abs(RGB_range(:,4) - Value(ii)));
        idx_min = find(diff == min(diff));
        R = RGB_range(idx_min,1);      G = RGB_range(idx_min,2);    B = RGB_range(idx_min,3); 
        if Value(ii) > RGB_range(end,4) 
            R = RGB_range(end,1);      G = RGB_range(end,2);        B = RGB_range(end,3);        
        elseif Value(ii) < RGB_range(1,4)
            R = RGB_range(1,1);        G = RGB_range(1,2);          B = RGB_range(1,3);
        end
        RGB(ii,1) = R/256;     RGB(ii,2) = G/256;     RGB(ii,3) = B/256;
    end
end


function [Value RGB] = Assign_Clones_RGB_value( Value, RGB_range, c ) 
%     cc = c + 3 ;
	for ii = 1 : length(Value)
        % use abs() because half will be negative values: we search for the one
        % where the difference is closest to zero.
        diff = (abs(RGB_range(:,4) - Value(ii)));
        idx_min = find(diff == min(diff));
        idx(ii) = idx_min;
        if Value(ii) > RGB_range(end,4) 
            idx(ii) = length(RGB_range(:,4));
        elseif Value(ii) < RGB_range(1,4)
            idx(ii) = 1;
        end
    end
    RGB(:,1) = RGB_range(idx,1)./256;
    RGB(:,2) = RGB_range(idx,2)./256;     
    RGB(:,3) = RGB_range(idx,3)./256;   
%     fprintf(' %2.1f ', Value)
%     fprintf('\n')
%     fprintf(' %2.0f ', idx)
%     fprintf('\n')
end


