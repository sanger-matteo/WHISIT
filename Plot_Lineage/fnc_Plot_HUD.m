function fnc_Plot_HUD ( clone_List, Frm, yrow, Dst_pl, hg_pl, LinW )
%
%fnc_Plot_HUD = according to parameters and option choosen in the interface
%   of WHISIT, the function change the general outline of the plot. Those
%   include format of axis, thickness of line, spacing between "cell lines"
%   size of font etc...
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------
global APP_opt;

Y_HUD = [];     % y-position where to plot the ID-name and pole info

if     APP_opt.t3_display_cID <= 1 
    str_label = {clone_List.ID_clone};
elseif APP_opt.t3_display_cID ==  2
    str_label = {clone_List.ID_ManualTrack};
end

%  PLOT LINEAGE TREE ******************************************************************************
if APP_opt.t3_PlotOpt_A == 1   
    
    if  APP_opt.t3_display_PL == 1                         % if option Pole is selected
        if  clone_List.ID_clone(end) == '.'  ...           % if is progenitor
            || clone_List.ID_clone(end-1) == '.'           % or after progenitor division
            % The polarity inheritance is unknown, therefore we place an 'x' instead of a dot
            Y_HUD = (yrow*Dst_pl +hg_pl) -hg_pl/2;
            plot( Frm -hg_pl*0.5 , Y_HUD, 'xk', 'MarkerSize', hg_pl, 'LineWidth', LinW);  

        elseif clone_List.ID_clone(end) ~= '.'
            % Plot a "O" filled according to Old/New pole
            if clone_List.ID_clone(end) == '1'             % OLD pole
              Y_HUD = (yrow*Dst_pl) +hg_pl/3;              MarkFaceCol = [0 0 0];
            elseif clone_List.ID_clone(end) == '2'         % NEW pole 
              Y_HUD = (yrow*Dst_pl +hg_pl) -hg_pl/3;       MarkFaceCol = [1 1 1];
            end
            plot( Frm -hg_pl*0.5 , Y_HUD, 'ok', 'MarkerSize', hg_pl, 'LineWidth', LinW, 'MarkerFaceColor', MarkFaceCol);
        end
    end

    if APP_opt.t3_display_cID >= 1                         % if option ID plot is not "none"
       if   clone_List.ID_clone(end) == '.'  ...           % if is progenitor
            || clone_List.ID_clone(end-1) == '.'           % or after progenitor division
            Y_HUD = (yrow*Dst_pl +hg_pl) -hg_pl ;
            h = text(Frm -hg_pl*0.5 , Y_HUD,  str_label ,'FontSize',hg_pl); 
            h.HorizontalAlignment = 'right' ;

       elseif clone_List.ID_clone(end) ~= '.'
            if clone_List.ID_clone(end) == '1'             % OLD pole
               Y_HUD = (yrow*Dst_pl) -hg_pl/3 ;
            elseif clone_List.ID_clone(end) == '2'         % NEW pole 
               Y_HUD = (yrow*Dst_pl +hg_pl) +hg_pl/3 ;
            end
            h = text(Frm -hg_pl*0.5 , Y_HUD,  str_label ,'FontSize',hg_pl);
            h.HorizontalAlignment = 'right' ;
       end
    end


% PLOT TRACKS or CLONES INDEPENDENTLY *************************************************************
elseif APP_opt.t3_PlotOpt_A == 2  ||  APP_opt.t3_PlotOpt_A == 3    
    
   if  APP_opt.t3_display_PL == 1                      % if option Pole is selected
       if  clone_List.ID_clone(end) == '.'  ...        % if is progenitor
           || clone_List.ID_clone(end-1) == '.'        % or after progenitor division
           % Plot a "X" because progenitor polarity is unknown
           Y_HUD = (yrow*Dst_pl +hg_pl) -hg_pl/2;
           plot( -5 , Y_HUD, 'xk', 'MarkerSize', hg_pl, 'LineWidth', LinW);
           
       elseif clone_List.ID_clone(end) ~= '.'
           % Plot a "O" filled according to Old/New pole
           Y_HUD = (yrow*Dst_pl) +hg_pl/2; 
           if clone_List.ID_clone(end) == '1';         MarkFaceCol = [0 0 0];  % OLD pole   
           elseif clone_List.ID_clone(end) == '2';     MarkFaceCol = [1 1 1];  % NEW pole 
           end
           plot( -5 , Y_HUD, 'ok', 'MarkerSize', hg_pl, 'LineWidth', LinW, 'MarkerFaceColor', MarkFaceCol);
       end
   end
   
   if  APP_opt.t3_display_cID >= 1                     % if option ID plot is not "none"
       if  clone_List.ID_clone(end) == '.'  ...        % if is progenitor
           ||  clone_List.ID_clone(end-1) == '.'       % or after progenitor division
       
           Y_HUD = (yrow*Dst_pl) +hg_pl/2 ;
           h = text(Frm -hg_pl*1.5 , Y_HUD,  str_label ,'FontSize',hg_pl);
           h.HorizontalAlignment = 'right' ;
           h1 = get(h, 'position') ;
           h.Position = [-10 h1(2) h1(3)] ;

       elseif clone_List.ID_clone(end) ~= '.'
           Y_HUD = (yrow*Dst_pl) +hg_pl/2 ;
           h = text(Frm -hg_pl*1.5 , Y_HUD,  str_label ,'FontSize',hg_pl);
           h.HorizontalAlignment = 'right' ;
           h1 = get(h, 'position') ;
           h.Position = [-10 h1(2) h1(3)] ;
       end
   end
   
end % if APP_opt.t3_PlotOpt_A == 1 || 2 || 3

end % MAIN fnc








