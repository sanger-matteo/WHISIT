function Plot_Lng_B2
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
pathParts_LUT = strsplit(mfilename('fullpath'), {'/','\'} ) ;        % mfilename = take path of currently running script.
path_LUT = fullfile(pathParts_LUT{1,1:end-1},'\') ;                  % fullfile  = build full filename from string parts
filename_LUT = [path_LUT '/LUT/' 'LUT_' APP_opt.t3_ColorMap_LUT '.txt'] ;
RGB_range = textread(filename_LUT) ;

% Create ordered list of the clones (cs_All) and the end-of-line clone (cS_EndLine) 
% This will allow to plot in the correct order, from the "latest" to the founder cell
[cS_EndLine, cN_EndLine, cS_All] = fnc_Organize_cloneList(clone_List) ;

% Color Range Setup and determining the Extreme Values
[IC_V , Rel_V] = Eval_RawValue_Extrema(clone_List);
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

Ext_Rng = 0.0 ;     % Extend Range value of X in both directions
t_range = linspace(MinRange, MaxRange, size(RGB_range,1)+ (size(RGB_range,1)*2*Ext_Rng) ) ;
RGB_range(:,4) = t_range( 1+ (size(RGB_range,1)*Ext_Rng)  :  end- (size(RGB_range,1)*Ext_Rng) )';




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
LAST_Frm = 1;                           % highest value frame, use to find max x-axis length
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
        Frm = ff + clone_List{idx_EndL}{1}.fr_birth;
        if Frm > LAST_Frm;     LAST_Frm = Frm;     end

        % If choosen, plot HUD: cell poles orientation and cell_ID number
        if ff == 1      % at frame of birth;
           fnc_Plot_HUD( clone_List{idx_EndL}{1}, Frm, yrow, Dst_pl, hg_pl , LinW)
        end

        % Divide the of the square-"cell" into 2 or 3 sub-areas
        switch APP_opt.t3_PlotOpt_C
            case {1, 2}
                Xs = [ Frm-1, Frm, Frm, Frm-1, Frm-1 ];
                % upper half
                Ys_UP = [ (yrow*Dst_pl)+(hg_pl*1/2), (yrow*Dst_pl)+(hg_pl*1/2) , ... % last point = first point
                          (yrow*Dst_pl)+(hg_pl*1)  , (yrow*Dst_pl)+(hg_pl*1),    (yrow*Dst_pl)+(hg_pl*1/2) ];
                % lower half
                Ys_DW = [ (yrow*Dst_pl)+(hg_pl*0)  , (yrow*Dst_pl)+(hg_pl*0) ,   ... % last point = first point
                          (yrow*Dst_pl)+(hg_pl*1/2), (yrow*Dst_pl)+(hg_pl*1/2),  (yrow*Dst_pl)+(hg_pl*0) ];
            case 3
                Xs = [ Frm-1, Frm, Frm, Frm-1, Frm-1 ];
                % upper part
                Ys_UP = [ (yrow*Dst_pl)+(hg_pl*2/3), (yrow*Dst_pl)+(hg_pl*2/3) ,  ... % last point = first point
                          (yrow*Dst_pl)+(hg_pl*1)  , (yrow*Dst_pl)+(hg_pl*1),     (yrow*Dst_pl)+(hg_pl*2/3) ];
                % middle part            
                Ys_ML = [ (yrow*Dst_pl)+(hg_pl*1/3), (yrow*Dst_pl)+(hg_pl*1/3) ,  ... % last point = first point
                          (yrow*Dst_pl)+(hg_pl*2/3), (yrow*Dst_pl)+(hg_pl*2/3),   (yrow*Dst_pl)+(hg_pl*1/3) ];
                % lower part
                Ys_DW = [ (yrow*Dst_pl)+(hg_pl*0)  , (yrow*Dst_pl)+(hg_pl*0) ,    ... % last point = first point
                          (yrow*Dst_pl)+(hg_pl*1/3), (yrow*Dst_pl)+(hg_pl*1/3),   (yrow*Dst_pl)+(hg_pl*0) ];
        end

        % Evaluate the Value and the RGB_Color to apply inside the square
        % and fill squares with evaluated RGB_Color
        switch APP_opt.t3_PlotOpt_C
            case 1      % Old_PL / New_PL 
                Val_UP = mean(clone_List{idx_EndL}{ff}.CH1.IC(  clone_List{idx_EndL}{ff}.CH1.Mask_PL_1 )) ;
                Val_DW = mean(clone_List{idx_EndL}{ff}.CH1.IC(  clone_List{idx_EndL}{ff}.CH1.Mask_PL_2 )) ;  
                A_UP = sum(sum(clone_List{idx_EndL}{ff}.CH1.Mask_PL_1));
                A_DW = sum(sum(clone_List{idx_EndL}{ff}.CH1.Mask_PL_2));
%                     fprintf('%2.0f  %2.0f  %2.0f\n',Val_UP , Val_ML, Val_DW)
                Weighted_Denom = (Val_UP/A_UP +Val_DW/A_DW) ;
                Val_UP = Val_UP / Weighted_Denom;
                Val_DW = Val_DW / Weighted_Denom;
                [Val_UP, RGB_UP] =  Assign_Value_RGB( Val_UP , V_Min, V_Max, RGB_range ); 
                [Val_DW, RGB_DW] =  Assign_Value_RGB( Val_DW , V_Min, V_Max, RGB_range );                  
                fill(Xs, Ys_UP, [RGB_UP(1), RGB_UP(2), RGB_UP(3)], 'LineStyle', 'none');
                fill(Xs, Ys_DW, [RGB_DW(1), RGB_DW(2), RGB_DW(3)], 'LineStyle', 'none');
            case 2      % Cytosol / Membr
                Val_UP = mean(clone_List{idx_EndL}{ff}.CH1.IC( clone_List{idx_EndL}{ff}.Mask_mCyto )) ;
                Val_DW = mean(clone_List{idx_EndL}{ff}.CH1.IC( clone_List{idx_EndL}{ff}.Mask_Memb )) ;
                A_UP = sum(sum(clone_List{idx_EndL}{ff}.Mask_mCyto ));
                A_DW = sum(sum(clone_List{idx_EndL}{ff}.Mask_Memb ));
%                     fprintf('%2.0f  %2.0f  %2.0f\n',Val_UP , Val_ML, Val_DW)
                Weighted_Denom = (Val_UP/A_UP +Val_DW/A_DW) ;
                Val_UP = Val_UP / Weighted_Denom;
                Val_DW = Val_DW / Weighted_Denom;                    
                [Val_UP, RGB_UP] =  Assign_Value_RGB( Val_UP , V_Min, V_Max, RGB_range );
                [Val_DW, RGB_DW] =  Assign_Value_RGB( Val_DW , V_Min, V_Max, RGB_range );               
                fill(Xs, Ys_UP, [RGB_UP(1), RGB_UP(2), RGB_UP(3)], 'LineStyle', 'none');
                fill(Xs, Ys_DW, [RGB_DW(1), RGB_DW(2), RGB_DW(3)], 'LineStyle', 'none');
            case 3      % Old_PL / Cytosol / New_PL                    
                Val_UP = mean(clone_List{idx_EndL}{ff}.CH1.IC(  clone_List{idx_EndL}{ff}.CH1.Mask_PL_1 )) ;
                Val_ML = mean(clone_List{idx_EndL}{ff}.CH1.IC(  clone_List{idx_EndL}{ff}.CH1.Mask_pCyto )) ;
                Val_DW = mean(clone_List{idx_EndL}{ff}.CH1.IC(  clone_List{idx_EndL}{ff}.CH1.Mask_PL_2 )) ;
                A_UP = sum(sum(clone_List{idx_EndL}{ff}.CH1.Mask_PL_1));
                A_ML = sum(sum(clone_List{idx_EndL}{ff}.CH1.Mask_pCyto));
                A_DW = sum(sum(clone_List{idx_EndL}{ff}.CH1.Mask_PL_2));
%                     fprintf('%2.0f  %2.0f  %2.0f\n',Val_UP , Val_ML, Val_DW)
                Weighted_Denom = (Val_UP/A_UP +Val_ML/A_ML +Val_DW/A_DW) ;
                Val_UP = Val_UP / Weighted_Denom;
                Val_ML = Val_ML / Weighted_Denom;
                Val_DW = Val_DW / Weighted_Denom;
%                     fprintf('%2.0f  %2.0f  %2.0f\n',Val_UP , Val_ML, Val_DW)
                [Val_UP, RGB_UP] =  Assign_Value_RGB( Val_UP , V_Min, V_Max, RGB_range );
                [Val_ML, RGB_ML] =  Assign_Value_RGB( Val_ML , V_Min, V_Max, RGB_range );
                [Val_DW, RGB_DW] =  Assign_Value_RGB( Val_DW , V_Min, V_Max, RGB_range );
                fill(Xs, Ys_UP, [RGB_UP(1), RGB_UP(2), RGB_UP(3)], 'LineStyle', 'none');
                fill(Xs, Ys_ML, [RGB_ML(1), RGB_ML(2), RGB_ML(3)], 'LineStyle', 'none');
                fill(Xs, Ys_DW, [RGB_DW(1), RGB_DW(2), RGB_DW(3)], 'LineStyle', 'none');
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

    end % ff
    
    % Draw two horizontal lines that bound the plotted data for the clone-Line
    plot( [Frm, Frm-length(clone_List{idx_EndL})], [yrow*Dst_pl,  yrow*Dst_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
    plot( [Frm, Frm-length(clone_List{idx_EndL})], [(yrow*Dst_pl)+hg_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
 
    % Draw vertical line at the end, to enclose all clone-Line
    plot( [Frm, Frm], [yrow*Dst_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
    
end %cc



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
            Xs = [ idx_birth+1 , idx_birth+1 ];
            if y_H_Row(idx_TA) < y_H_Row(idx_TB)
                Ys = [ y_H_Row(idx_TA)*Dst_pl,  (y_H_Row(idx_TB)*Dst_pl)+hg_pl];
            elseif y_H_Row(idx_TA) > y_H_Row(idx_TB)
                Ys = [ y_H_Row(idx_TB)*Dst_pl,  (y_H_Row(idx_TA)*Dst_pl)+hg_pl];
            end
            plot(Xs, Ys, '-', 'LineWidth', LinW, 'Color', Color_BorderLines); 

            % Plot all time points. We use two counters for "time":
            % - ff  = stores the position of each element (/frame) of the array
            %         clone_List{idx_EndL}: simply goes from 1 to last frame
            % - Frm = real time variable that plots the data extracted at clone_List{idx_EndL}{ff} 
            %         in the correct time point (x-axis) of the lineage tree
            for ff = 1 : length(clone_List{idx_Anc})
                Frm = ff + clone_List{idx_Anc}{1}.fr_birth;
                if Frm > LAST_Frm;     LAST_Frm = Frm;     end

                % If choosen, plot HUD: cell poles orientation and cell_ID number
                if ff == 1      % at frame of birth;
                   fnc_Plot_HUD( clone_List{idx_Anc}{1}, Frm, yrow, Dst_pl, hg_pl , LinW)
                end

                % Divide the of the square-"cell" into 2 or 3 sub-parts
                switch APP_opt.t3_PlotOpt_C
                    case {1, 2}
                        Xs = [ Frm-1, Frm, Frm, Frm-1, Frm-1 ];
                        % upper half
                        Ys_UP = [ (yrow*Dst_pl)+(hg_pl*1/2), (yrow*Dst_pl)+(hg_pl*1/2) ,  ... % last point = first point
                                  (yrow*Dst_pl)+(hg_pl*1)  , (yrow*Dst_pl)+(hg_pl*1),     (yrow*Dst_pl)+(hg_pl*1/2) ];
                        % lower half
                        Ys_DW = [ (yrow*Dst_pl)+(hg_pl*0)  , (yrow*Dst_pl)+(hg_pl*0) ,    ... % last point = first point
                                  (yrow*Dst_pl)+(hg_pl*1/2), (yrow*Dst_pl)+(hg_pl*1/2),   (yrow*Dst_pl)+(hg_pl*0) ];
                    case 3
                        Xs = [ Frm-1, Frm, Frm, Frm-1, Frm-1 ];
                        % upper part
                        Ys_UP = [ (yrow*Dst_pl)+(hg_pl*2/3), (yrow*Dst_pl)+(hg_pl*2/3) ,   ... % last point = first point
                                  (yrow*Dst_pl)+(hg_pl*1)  , (yrow*Dst_pl)+(hg_pl*1),      (yrow*Dst_pl)+(hg_pl*2/3) ];
                        % middle part            
                        Ys_ML = [ (yrow*Dst_pl)+(hg_pl*1/3), (yrow*Dst_pl)+(hg_pl*1/3) ,   ... % last point = first point
                                  (yrow*Dst_pl)+(hg_pl*2/3), (yrow*Dst_pl)+(hg_pl*2/3),    (yrow*Dst_pl)+(hg_pl*1/3) ];
                        % lower part
                        Ys_DW = [ (yrow*Dst_pl)+(hg_pl*0)  , (yrow*Dst_pl)+(hg_pl*0) ,     ... % last point = first point
                                  (yrow*Dst_pl)+(hg_pl*1/3), (yrow*Dst_pl)+(hg_pl*1/3),    (yrow*Dst_pl)+(hg_pl*0) ];
                end

                % Evaluate the Value and the RGB_Color to apply inside the square
                % and fill squares with evaluated RGB_Color
                switch APP_opt.t3_PlotOpt_C
                    case 1      % Old_PL / New_PL
                        Val_UP = mean(clone_List{idx_Anc}{ff}.CH1.IC(  clone_List{idx_Anc}{ff}.CH1.Mask_PL_1 )) ; 
                        Val_DW = mean(clone_List{idx_Anc}{ff}.CH1.IC(  clone_List{idx_Anc}{ff}.CH1.Mask_PL_2 )) ;
                        A_UP = sum(sum(clone_List{idx_Anc}{ff}.CH1.Mask_PL_1));
                        A_DW = sum(sum(clone_List{idx_Anc}{ff}.CH1.Mask_PL_2));
    %                     fprintf('%2.0f  %2.0f  %2.0f\n',Val_UP , Val_ML, Val_DW)    
                        Weighted_Denom = (Val_UP/A_UP +Val_DW/A_DW) ;
                        
                        Val_UP = Val_UP / Weighted_Denom;
                        Val_DW = Val_DW / Weighted_Denom;
                        [Val_UP, RGB_UP] =  Assign_Value_RGB( Val_UP , V_Min, V_Max, RGB_range ); 
                        [Val_DW, RGB_DW] =  Assign_Value_RGB( Val_DW , V_Min, V_Max, RGB_range );                  
                        fill(Xs, Ys_UP, [RGB_UP(1), RGB_UP(2), RGB_UP(3)], 'LineStyle', 'none');
                        fill(Xs, Ys_DW, [RGB_DW(1), RGB_DW(2), RGB_DW(3)], 'LineStyle', 'none');
                    case 2      % Cytosol / Membr
                        Val_UP = mean(clone_List{idx_Anc}{ff}.CH1.IC( clone_List{idx_Anc}{ff}.Mask_mCyto )) ;
                        Val_DW = mean(clone_List{idx_Anc}{ff}.CH1.IC( clone_List{idx_Anc}{ff}.Mask_Memb )) ;
                        A_UP = sum(sum(clone_List{idx_Anc}{ff}.Mask_mCyto ));
                        A_DW = sum(sum(clone_List{idx_Anc}{ff}.Mask_Memb ));
    %                     fprintf('%2.0f  %2.0f  %2.0f\n',Val_UP , Val_ML, Val_DW)
                        Weighted_Denom = (Val_UP/A_UP +Val_DW/A_DW) ;
                        
                        Val_UP = Val_UP / Weighted_Denom;
                        Val_DW = Val_DW / Weighted_Denom;                    
                        [Val_UP, RGB_UP] =  Assign_Value_RGB( Val_UP , V_Min, V_Max, RGB_range );
                        [Val_DW, RGB_DW] =  Assign_Value_RGB( Val_DW , V_Min, V_Max, RGB_range );               
                        fill(Xs, Ys_UP, [RGB_UP(1), RGB_UP(2), RGB_UP(3)], 'LineStyle', 'none');
                        fill(Xs, Ys_DW, [RGB_DW(1), RGB_DW(2), RGB_DW(3)], 'LineStyle', 'none');
                    case 3      % Old_PL / Cytosol / New_PL                    
                        Val_UP = mean(clone_List{idx_Anc}{ff}.CH1.IC(  clone_List{idx_Anc}{ff}.CH1.Mask_PL_1 )) ;
                        Val_ML = mean(clone_List{idx_Anc}{ff}.CH1.IC(  clone_List{idx_Anc}{ff}.CH1.Mask_pCyto )) ;
                        Val_DW = mean(clone_List{idx_Anc}{ff}.CH1.IC(  clone_List{idx_Anc}{ff}.CH1.Mask_PL_2 )) ;
                        A_UP = sum(sum(clone_List{idx_Anc}{ff}.CH1.Mask_PL_1));
                        A_ML = sum(sum(clone_List{idx_Anc}{ff}.CH1.Mask_pCyto));
                        A_DW = sum(sum(clone_List{idx_Anc}{ff}.CH1.Mask_PL_2));
                        Weighted_Denom = (Val_UP/A_UP+Val_ML/A_ML+Val_DW/A_DW) ;
                        
    %                         fprintf('%2.0f  %2.0f  %2.0f\n',Val_UP , Val_ML, Val_DW)
                        Val_UP = Val_UP / Weighted_Denom ;
                        Val_ML = Val_ML / Weighted_Denom ;
                        Val_DW = Val_DW / Weighted_Denom ;
    %                         fprintf('%2.0f  %2.0f  %2.0f\n',Val_UP , Val_ML, Val_DW)
                        [Val_UP, RGB_UP] =  Assign_Value_RGB( Val_UP , V_Min, V_Max, RGB_range );
                        [Val_ML, RGB_ML] =  Assign_Value_RGB( Val_ML , V_Min, V_Max, RGB_range );
                        [Val_DW, RGB_DW] =  Assign_Value_RGB( Val_DW , V_Min, V_Max, RGB_range );
                        fill(Xs, Ys_UP, [RGB_UP(1), RGB_UP(2), RGB_UP(3)], 'LineStyle', 'none');
                        fill(Xs, Ys_ML, [RGB_ML(1), RGB_ML(2), RGB_ML(3)], 'LineStyle', 'none');
                        fill(Xs, Ys_DW, [RGB_DW(1), RGB_DW(2), RGB_DW(3)], 'LineStyle', 'none');
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
                
            end % ff
                           
            % Draw two horizontal lines that bound the plotted data for the clone-Line
            plot( [Frm, Frm-length(clone_List{idx_Anc})], [yrow*Dst_pl,  yrow*Dst_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
            plot( [Frm, Frm-length(clone_List{idx_Anc})], [(yrow*Dst_pl)+hg_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
                        
            % Draw vertical line at the end, to enclose all end_Lines
            plot( [Frm, Frm], [yrow*Dst_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);

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
    end % while cc
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

        % Divide the of the square-"cell" into 2 or 3 sub-parts for
        % plotting the value(s) at time point ff
        switch APP_opt.t3_PlotOpt_C
            case {1, 2}
                Xs = [ ff-1, ff, ff, ff-1, ff-1 ];
                % upper half
                Ys_UP = [ (yrow*Dst_pl)+(hg_pl*1/2), (yrow*Dst_pl)+(hg_pl*1/2) , ... % last point = first point
                          (yrow*Dst_pl)+(hg_pl*1)  , (yrow*Dst_pl)+(hg_pl*1),    (yrow*Dst_pl)+(hg_pl*1/2) ];
                % lower half
                Ys_DW = [ (yrow*Dst_pl)+(hg_pl*0)  , (yrow*Dst_pl)+(hg_pl*0) ,   ... % last point = first point
                          (yrow*Dst_pl)+(hg_pl*1/2), (yrow*Dst_pl)+(hg_pl*1/2),  (yrow*Dst_pl)+(hg_pl*0) ];
            case 3
                Xs = [ ff-1, ff, ff, ff-1, ff-1 ];
                % upper part
                Ys_UP = [ (yrow*Dst_pl)+(hg_pl*2/3), (yrow*Dst_pl)+(hg_pl*2/3) ,  ... % last point = first point
                          (yrow*Dst_pl)+(hg_pl*1)  , (yrow*Dst_pl)+(hg_pl*1),     (yrow*Dst_pl)+(hg_pl*2/3) ];
                % middle part            
                Ys_ML = [ (yrow*Dst_pl)+(hg_pl*1/3), (yrow*Dst_pl)+(hg_pl*1/3) ,  ... % last point = first point
                          (yrow*Dst_pl)+(hg_pl*2/3), (yrow*Dst_pl)+(hg_pl*2/3),   (yrow*Dst_pl)+(hg_pl*1/3) ];
                % lower part
                Ys_DW = [ (yrow*Dst_pl)+(hg_pl*0)  , (yrow*Dst_pl)+(hg_pl*0) ,    ... % last point = first point
                          (yrow*Dst_pl)+(hg_pl*1/3), (yrow*Dst_pl)+(hg_pl*1/3),   (yrow*Dst_pl)+(hg_pl*0) ];
        end

        % Evaluate the Value and the RGB_Color to apply inside the square
        % and fill squares with evaluated RGB_Color
        switch APP_opt.t3_PlotOpt_C
            case 1      % Old_PL / New_PL
                Val_UP = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.CH1.Mask_PL_1 )) ; 
                Val_DW = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.CH1.Mask_PL_2 )) ;
                A_UP = sum(sum(clone_List{cc}{ff}.CH1.Mask_PL_1));
                A_DW = sum(sum(clone_List{cc}{ff}.CH1.Mask_PL_2));
                Weighted_Denom = (Val_UP/A_UP +Val_DW/A_DW) ;
                Val_UP = Val_UP / Weighted_Denom;
                Val_DW = Val_DW / Weighted_Denom;
                [Val_UP, RGB_UP] =  Assign_Value_RGB( Val_UP , V_Min, V_Max, RGB_range ); 
                [Val_DW, RGB_DW] =  Assign_Value_RGB( Val_DW , V_Min, V_Max, RGB_range );                  
                fill(Xs, Ys_UP, [RGB_UP(1), RGB_UP(2), RGB_UP(3)], 'LineStyle', 'none');
                fill(Xs, Ys_DW, [RGB_DW(1), RGB_DW(2), RGB_DW(3)], 'LineStyle', 'none');
            case 2      % Cytosol / Membr
                Val_UP = mean(clone_List{cc}{ff}.CH1.IC( clone_List{cc}{ff}.Mask_mCyto )) ;
                Val_DW = mean(clone_List{cc}{ff}.CH1.IC( clone_List{cc}{ff}.Mask_Memb )) ;
                A_UP = sum(sum(clone_List{cc}{ff}.Mask_mCyto ));
                A_DW = sum(sum(clone_List{cc}{ff}.Mask_Memb ));
                Weighted_Denom = (Val_UP/A_UP +Val_DW/A_DW) ;
                Val_UP = Val_UP / Weighted_Denom;
                Val_DW = Val_DW / Weighted_Denom;                    
                [Val_UP, RGB_UP] =  Assign_Value_RGB( Val_UP , V_Min, V_Max, RGB_range );
                [Val_DW, RGB_DW] =  Assign_Value_RGB( Val_DW , V_Min, V_Max, RGB_range );               
                fill(Xs, Ys_UP, [RGB_UP(1), RGB_UP(2), RGB_UP(3)], 'LineStyle', 'none');
                fill(Xs, Ys_DW, [RGB_DW(1), RGB_DW(2), RGB_DW(3)], 'LineStyle', 'none');
            case 3      % Old_PL / Cytosol / New_PL                    
                Val_UP = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.CH1.Mask_PL_1 )) ;
                Val_ML = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.CH1.Mask_pCyto )) ;
                Val_DW = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.CH1.Mask_PL_2 )) ;
                A_UP = sum(sum(clone_List{cc}{ff}.CH1.Mask_PL_1));
                A_ML = sum(sum(clone_List{cc}{ff}.CH1.Mask_pCyto));
                A_DW = sum(sum(clone_List{cc}{ff}.CH1.Mask_PL_2));
                Weighted_Denom = (Val_UP/A_UP +Val_ML/A_ML +Val_DW/A_DW) ;
                Val_UP = Val_UP / Weighted_Denom;
                Val_ML = Val_ML / Weighted_Denom;
                Val_DW = Val_DW / Weighted_Denom;
                [Val_UP, RGB_UP] =  Assign_Value_RGB( Val_UP , V_Min, V_Max, RGB_range );
                [Val_ML, RGB_ML] =  Assign_Value_RGB( Val_ML , V_Min, V_Max, RGB_range );
                [Val_DW, RGB_DW] =  Assign_Value_RGB( Val_DW , V_Min, V_Max, RGB_range );
                fill(Xs, Ys_UP, [RGB_UP(1), RGB_UP(2), RGB_UP(3)], 'LineStyle', 'none');
                fill(Xs, Ys_ML, [RGB_ML(1), RGB_ML(2), RGB_ML(3)], 'LineStyle', 'none');
                fill(Xs, Ys_DW, [RGB_DW(1), RGB_DW(2), RGB_DW(3)], 'LineStyle', 'none');
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
            plot( Frm -0.5, (yrow*Dst_pl +hg_pl) +hg_pl/3 , '.', 'MarkerSize', 8, 'LineWidth', 0.2, 'Color', El_Ec_Col);
        end

    end % ff
    
    % Draw two horizontal lines that bound the plotted data for the clone-Line
    plot( [0, length(clone_List{cc})], [yrow*Dst_pl,  yrow*Dst_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
    plot( [0, length(clone_List{cc})], [(yrow*Dst_pl)+hg_pl, (yrow*Dst_pl)+hg_pl], 'LineWidth', LinW, 'Color', Color_BorderLines);
    
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

end % MAIN fnc







%% 
% *********************************************************************************************************
% -----> SCRIP-RELATED FUNCTIONS --------------------------------------------------------------------------
% *********************************************************************************************************

function [IC_V , Rel_V, Distr_Vals] = Eval_RawValue_Extrema(clone_List)
% ----- Fluorescence VALUE SCALING ----------------------------------------
% ESTABLISH RELATIVE MIN MAX, ACCORDING TO...
% Find Min and Max in all fluorescence Frames
global APP_opt ;
I_max = -1 ;           I_min = 2*10^20 ;
R_avg_max = -1 ;       R_avg_min = 2*10^20 ;
Distr_Vals = [];
for cc = 1 : length(clone_List)
    for ff = 1 : length(clone_List{cc})
        Val_1 = [];     Val_2 = [];     Val_3 = [];
        
        if I_max < max( max(clone_List{cc}{ff}.CH1.IC))
            I_max = max( max(clone_List{cc}{ff}.CH1.IC) );
        end
        if I_min > min( min(clone_List{cc}{ff}.CH1.IC))
            I_min = min( min(clone_List{cc}{ff}.CH1.IC));
        end 

        % Evaluate the Value and the RGB_Color to apply inside the square
        switch APP_opt.t3_PlotOpt_C
            case 1      % Old_PL / New_PL
                Val_1 = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.CH1.Mask_PL_1 )) ;
                Val_2 = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.CH1.Mask_PL_2 )) ;
                A1 = sum(sum(clone_List{cc}{ff}.CH1.Mask_PL_1));
                A2 = sum(sum(clone_List{cc}{ff}.CH1.Mask_PL_2));
                Weighted_Denom = (Val_1/A1+Val_2/A2) ;
                Val_1 = Val_1 / Weighted_Denom;
                Val_2 = Val_2 / Weighted_Denom;
            case 2      % Cytosol / Membr
                Val_1 = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.Mask_mCyto )) ;
                Val_2 = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.Mask_Memb )) ;
                A1 = sum(sum(clone_List{cc}{ff}.Mask_mCyto));
                A2 = sum(sum(clone_List{cc}{ff}.Mask_Memb));
                Weighted_Denom = (Val_1/A1+Val_2/A2) ;
                Val_1 = Val_1 / Weighted_Denom;
                Val_2 = Val_2 / Weighted_Denom;
            case 3      % Old_PL / Cytosol / New_PL
                Val_1 = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.CH1.Mask_PL_1 )) ;
                Val_2 = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.CH1.Mask_pCyto )) ;
                Val_3 = mean(clone_List{cc}{ff}.CH1.IC(  clone_List{cc}{ff}.CH1.Mask_PL_2 )) ;
                A1 = sum(sum(clone_List{cc}{ff}.CH1.Mask_PL_1));
                A2 = sum(sum(clone_List{cc}{ff}.CH1.Mask_pCyto));
                A3 = sum(sum(clone_List{cc}{ff}.CH1.Mask_PL_2));
                Weighted_Denom = (Val_1/A1+Val_2/A2+Val_3/A3) ;
                Val_1 = Val_1 / Weighted_Denom ;
                Val_2 = Val_2 / Weighted_Denom ;
                Val_3 = Val_3 / Weighted_Denom ;
        end
        % Relative min-max to the avg values Cyto and Poles
        if R_avg_max < Val_1
            R_avg_max = Val_1;       end
        if R_avg_min > Val_1
            R_avg_min = Val_1;       end  
        if R_avg_max < Val_2
            R_avg_max = Val_2;       end
        if R_avg_min > Val_2
            R_avg_min = Val_2;       end
        if ~isempty(Val_3)  &&  R_avg_max < Val_3
            R_avg_max = Val_3;       end
        if ~isempty(Val_3)  &&  R_avg_min > Val_3
            R_avg_min = Val_3;       end
        
        if ~isempty(Val_3)
            Distr_Vals = [Distr_Vals , Val_1, Val_2, Val_3] ;            
        elseif isempty(Val_3)
            Distr_Vals = [Distr_Vals , Val_1, Val_2] ;
        end
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


