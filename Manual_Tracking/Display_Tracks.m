function Display_Tracks
%%Display_TrackAid - Show the cell-track(s) according to the GUI selection
%
% APP_opt.t5_display_Track store the option chosen via the GUI panel 
% ShowTracks_Radio. In case that APP_opt.t5_display_Track ==
%
% 0 - none - Display no track at all
%
% 1-4 - these option work ONLY on the selected scc-th cell-track
% 1 - the previous frame XY point, if tracked
% 2 - only present frame XY point, if tracked
% 3 - all the XY point tracked up to present frame
% 4 - all the XY point tracked in the entire movie
%
% 5-8 - these option show the tracking done for ALL cell-tracks stored
% 5 - same as 1, but for all cells
% 6 - same as 2, but for all cells
% 7 - same as 3, but for all cells
% 8 - same as 4, but for all cells
%
% Options 1-4 color tracks with RGB value APP_opt.t5_display_Color decided
% by the user in the GUI. Options 5-8 will color each track with a unique
% color. Newly created cells are assigned a color when plotted for the
% first time. This RGB color is soter in CellTracks{3,:} and it is reused
% thereafter for the specific cell-track. This consistency will aid the
% user during tracking.
%  
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


% Quick explanation for repeatedly used variables and commands
% idx_scc ---- Find the column position for scc-th cell in CellTracks array  
% t_range ---- Store the time window to plot (as frame numbers)
% list_IDct -- Create a list of all ID-numbers stored in CellTracks 
% CellTracks - a global variable where that stores the information of
%              every cell track. 
% scc -------- stores the ID-number of the currently selected cell-track


global APP_opt;	    global CellTracks;     global scc;

clr = APP_opt.t5_display_Color ;        % Display track with color RGB  
mk_sz = APP_opt.t5_display_Marker  ;    % Chosen marker size   


% ------- DISPLAY ONLY scc-th CELL --------------------------------------------------------
if APP_opt.t5_display_Track == 1          % DISPLAY '# cell: previous frame'
    if APP_opt.t5_ff -1 >= 1
        idx_scc = find(cell2mat(CellTracks( 1,: )) == scc );
        t_range = APP_opt.t5_ff -1; 
        
        plot(CellTracks{2,idx_scc}(t_range , 1) , CellTracks{2,idx_scc}(t_range , 2), ...
             '.-', 'Color', clr, 'MarkerSize', mk_sz );
    end
    
elseif APP_opt.t5_display_Track == 2      % DISPLAY '# cell: only present frame'
    if APP_opt.t5_ff >= 1
        idx_scc = find(cell2mat(CellTracks( 1,: )) == scc );
        t_range = APP_opt.t5_ff ; 
        
        plot(CellTracks{2,idx_scc}(t_range , 1) , CellTracks{2,idx_scc}(t_range , 2), ...
             '.-', 'Color', clr, 'MarkerSize', mk_sz );
    end
    
elseif APP_opt.t5_display_Track == 3       % DISPLAY '# cell: up to present frame'
    idx_scc = find(cell2mat(CellTracks( 1,: )) == scc );
    
    t_range = find(CellTracks{2,idx_scc}(:,1)~=0);
    if ~isempty(t_range) && t_range(1)<=APP_opt.t5_ff && APP_opt.t5_ff<=t_range(end)
        t_range = [t_range(1) : APP_opt.t5_ff];
        Xs = CellTracks{2, idx_scc}(t_range , 1);
        Ys = CellTracks{2, idx_scc}(t_range , 2);
        % remove '0' coordinates - not-tracked points will not be considered
        Xs = Xs(Xs~=0);
        Ys = Ys(Ys~=0);
        plot( Xs , Ys , '.-', 'Color', clr, 'MarkerSize', mk_sz );
        plot( Xs(1), Ys(1), '*', 'Color', clr, 'MarkerSize', mk_sz-3 , 'LineWidth', mk_sz/10);
    end   
    
elseif APP_opt.t5_display_Track == 4     % DISPLAY '# cell: all frames'
    idx_scc = find(cell2mat(CellTracks( 1,: )) == scc );
    t_range = find(CellTracks{2,idx_scc}(:,1)~=0);
    
    if ~isempty(t_range)
        t_range = [t_range(1) : t_range(end)];
        Xs = CellTracks{2, idx_scc}(t_range , 1);
        Ys = CellTracks{2, idx_scc}(t_range , 2);
        % remove '0' coordinates - not-tracked points will not be considered
        Xs = Xs(Xs~=0);
        Ys = Ys(Ys~=0);
        plot( Xs , Ys , '.-', 'Color', clr, 'MarkerSize', mk_sz );
        plot( Xs(1), Ys(1), '*', 'Color', clr, 'MarkerSize', mk_sz-3 , 'LineWidth', mk_sz/10);
    end
 
    
    
    
% ------- DISPLAY ALL CELLS ---------------------------------------
elseif APP_opt.t5_display_Track == 5     % DISPLAY 'all cells: previous frame'
    if APP_opt.t5_ff-1 >= 1
        t_range = APP_opt.t5_ff -1;
        list_IDct = cell2mat(CellTracks( 1,: ));
        
        for kk = 1 : length(list_IDct)            
            colNum = list_IDct(kk);        % take array position for kk-th track
            idx_scc = find(cell2mat(CellTracks( 1,: )) == colNum);
            Xs = CellTracks{2,idx_scc}(t_range , 1);
            Ys = CellTracks{2,idx_scc}(t_range , 2);
            if ~isempty(Xs) && (Xs~=0 && Ys~=0)         
                % Plot and choose correct track color

                if isempty(CellTracks{3,idx_scc})       % it is a new cell-track
                    hh = plot( Xs , Ys , '.', 'MarkerSize', mk_sz );
                    c_clr = get(hh, 'Color');           % take auto-generated color from previous plot
                    CellTracks{3,idx_scc} = c_clr ;     % save color for this cell-track
                else
                    c_clr = CellTracks{3,idx_scc} ;     % take color specific to the cell-track
                    plot( Xs , Ys , '.', 'MarkerSize', mk_sz, 'Color', c_clr);
                end
                text(Xs,Ys, [' ',num2str(colNum)], 'Color', c_clr, 'FontSize',mk_sz ); 
                
            end
        end            
    end
    
    
elseif APP_opt.t5_display_Track == 6     % DISPLAY 'all cells: only present frame'
    if APP_opt.t5_ff >= 1
        t_range = APP_opt.t5_ff ;
        list_IDct = cell2mat(CellTracks( 1,: ));
        
        for kk = 1 : length(list_IDct)
            colNum = list_IDct(kk);        % take array position for kk-th track
            idx_scc = find(cell2mat(CellTracks( 1,: )) == colNum);
            Xs = CellTracks{2,idx_scc}(t_range , 1);
            Ys = CellTracks{2,idx_scc}(t_range , 2);
            if ~isempty(Xs) && (Xs~=0 && Ys~=0)         
                % Plot and choose correct track color

                if isempty(CellTracks{3,idx_scc})       % it is a new cell-track
                    hh = plot( Xs , Ys , '.', 'MarkerSize', mk_sz );
                    c_clr = get(hh, 'Color');           % take auto-generated color from previous plot
                    CellTracks{3,idx_scc} = c_clr ;     % save color for this cell-track
                else
                    c_clr = CellTracks{3,idx_scc} ;     % take color specific to the cell-track
                    plot( Xs , Ys , '.', 'MarkerSize', mk_sz, 'Color', c_clr);
                end
                text(Xs,Ys, [' ',num2str(colNum)], 'Color', c_clr, 'FontSize',mk_sz );
                
            end
        end            
    end
    
elseif APP_opt.t5_display_Track == 7     % DISPLAY 'all cells: up to present frame'
    list_IDct = cell2mat(CellTracks( 1,: ));
    
    for kk = 1 : length(list_IDct)
        colNum = list_IDct(kk);        % take array position for kk-th track
        idx_scc = find(cell2mat(CellTracks( 1,: )) == colNum);
        t_range = find(CellTracks{2,idx_scc}(:,1)~=0);
        
        if ~isempty(t_range) && t_range(1)<=APP_opt.t5_ff && APP_opt.t5_ff<=t_range(end)
            t_range = [t_range(1) : APP_opt.t5_ff];          % selected correct time window         
            Xs = CellTracks{2, idx_scc}(t_range , 1);
            Ys = CellTracks{2, idx_scc}(t_range , 2);
            % remove '0' coordinates - not-tracked points will not be considered
            Xs = Xs(Xs~=0);
            Ys = Ys(Ys~=0);
            if ~isempty(Xs)         
                % Plot and choose correct track color
                
                if isempty(CellTracks{3,idx_scc})
                    hh = plot( Xs , Ys , '.-', 'MarkerSize', mk_sz );
                    c_clr = get(hh, 'Color');           % take auto-generated color from previous plot
                    plot( Xs(1), Ys(1), '*', 'Color', c_clr, 'MarkerSize', mk_sz-3 , 'LineWidth', mk_sz/10);
                    CellTracks{3,idx_scc} = c_clr ;     % save color for this cell-track
                else
                    c_clr = CellTracks{3,idx_scc} ;     % take color specific to the cell-track
                    plot( Xs , Ys , '.-', 'MarkerSize', mk_sz, 'Color', c_clr);
                    plot( Xs(1), Ys(1), '*', 'Color', c_clr, 'MarkerSize', mk_sz-3 , 'LineWidth', mk_sz/10);
                end
                text(Xs(end),Ys(end), [' ',num2str(colNum)], 'Color', c_clr, 'FontSize',mk_sz );
                
            end
        end
    end
    
elseif APP_opt.t5_display_Track == 8     % DISPLAY 'all cells: all frames' 
    list_IDct = cell2mat(CellTracks( 1,: ));
    
    for kk = 1 : length(list_IDct)
        colNum = list_IDct(kk);        % take array position for kk-th track
        idx_scc = find(cell2mat(CellTracks( 1,: )) == colNum);
        t_range = find(CellTracks{2,idx_scc}(:,1)~=0);
        
        if ~isempty(t_range)           
            t_range = [t_range(1) : t_range(end)];          % selected correct time window 
            Xs = CellTracks{2, idx_scc}(t_range , 1);
            Ys = CellTracks{2, idx_scc}(t_range , 2);
            % remove '0' coordinates - not-tracked points will not be considered
            Xs = Xs(Xs~=0);
            Ys = Ys(Ys~=0);
            if ~isempty(Xs)         
                % Plot and choose correct track color            
                
                if isempty(CellTracks{3,idx_scc})
                    hh = plot( Xs , Ys , '.-', 'MarkerSize', mk_sz );
                    c_clr = get(hh, 'Color');           % take auto-generated color from previous plot
                    plot( Xs(1), Ys(1), '*', 'Color', c_clr, 'MarkerSize', mk_sz-3 , 'LineWidth', mk_sz/10);
                    CellTracks{3,idx_scc} = c_clr ;     % save color for this cell-track
                else
                    c_clr = CellTracks{3,idx_scc} ;     % take color specific to the cell-track
                    plot( Xs , Ys , '.-', 'MarkerSize', mk_sz, 'Color', c_clr);
                    plot( Xs(1), Ys(1), '*', 'Color', c_clr, 'MarkerSize', mk_sz-3 , 'LineWidth', mk_sz/10);
                end
                text(Xs(end),Ys(end), [' ',num2str(colNum)], 'Color', c_clr, 'FontSize',mk_sz );
                
            end
        end
    end
    
end %/if Display_Track
    
    

end



