function Main_Fx(hObject, eventdata, handles)
% Main Function of the program
%
% ------ RETRIVE CELL'S CONTOUR -------------------------------------------
global GUI_opt ;    global my_DB ;

if GUI_opt.choice_Retrive == 1          % --- RETRIVE CELL'S CONTOUR
    if GUI_opt.choice_Analyse == 0 && ~isempty(GUI_opt.exp_name);    % if not selected, Retrive_Cells need name to save DB
        GUI_opt.my_DB = '1_DB' ;
    end
    Retrive_Cells ;
end  

% ------ ANALYSIS SIGNAL --------------------------------------------------
if GUI_opt.choice_Analyse == 1 ;
    temp = load([GUI_opt.path_DIR ,'/' ,GUI_opt.my_DB ]);
    my_DB = temp.my_DB ;
    clear temp ;
    
    tot_N_cells = 0;        % Count number of cells to analyse in my_DB.mat
    for f = 1 : length(my_DB)    
       tot_N_cells = tot_N_cells + length(my_DB(f).cell) ;
    end

    n_c = 1;
    for f = 1 : length(my_DB)                   % go through all frames
         % for each frame calculate the average background signal and the minimum
        [BkGr, min_px, J_M] = BckGrd_calc( f );
        my_DB(f).M_BkGr_Frame = ~J_M ;          % Background Signal mask for whole frame
       
        for c = 1 : length(my_DB(f).cell)           % go through each cell
            if GUI_opt.choice_Analyse == 1
                Signal_Analysis( f,c ) ;
                my_DB(f).cell(c).geom.BkGr_Frame = BkGr;
                my_DB(f).cell(c).geom.min_px_Frame = min_px;
            end

            if GUI_opt.STOP == 1                    % if STOP button is pressed, 
                return                              % interrupt function
            end 

            if GUI_opt.choice_plot == 1             % PLOT
                Plot_Cells_9( f,c ) ;         
            end
            text_OUT = [ 'Cell ', num2str(n_c) ,' of ', num2str(tot_N_cells) ];
            set(handles.text_OUT , 'String', text_OUT );
            pause(0.001);           % apparently without small pause it goes too fast to change GUI
            n_c = n_c + 1 ;
        end
    end
    save( [ GUI_opt.path_DIR ,'/' , GUI_opt.exp_name , '_2_DB'], 'my_DB');
end

% ------ SAVE DATA --------------------------------------------------------
if GUI_opt.choice_Save == 1                    
    Save_Data ;
end

end


%--------------------------------------------------------------------------
%    Copyright (c) 2016; WHISIT, Matteo Sangermani, All rights reserved
%--------------------------------------------------------------------------

