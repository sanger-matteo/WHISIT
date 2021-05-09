function Plot_Lng_MAIN_Fnc(app)

global APP_opt ;
% global clone_List ;

Str_value = strsplit(app.t3_ListBox_PlotType.Value, '-'); 
APP_opt.t3_PlotOpt_A = str2num(Str_value{1}(1));           

Str_value = strsplit(app.t3_ListBox_LineType.Value, '-');
APP_opt.t3_PlotOpt_B = str2num(Str_value{1}(1));

Str_value = strsplit(app.t3_DropDownValue.Value, '-');
APP_opt.t3_PlotOpt_C = str2num(Str_value{1}(1));

% *** PLOT single clone_List.mat ******************************************
if APP_opt.t3_Generation_VS_Single == 1

    switch APP_opt.t3_PlotOpt_B
        case 1 
            if APP_opt.t3_PlotOpt_A <= 2
                Plot_Lng_B1;    
            elseif APP_opt.t3_PlotOpt_A == 3
                Plot_Lng_M1;
            end
        case 2
            if APP_opt.t3_PlotOpt_A <= 2
                Plot_Lng_B2;    
            elseif APP_opt.t3_PlotOpt_A == 3
                Plot_Lng_M2;
            end
        case 3
            if APP_opt.t3_PlotOpt_A <= 2
                Plot_Lng_B3;    
            elseif APP_opt.t3_PlotOpt_A == 3
                Plot_Lng_M3;
            end
        case 4     
            Plot_Lng_B4;        
    end

% *** ANALYSE Intragenerational ******************************************* 
elseif APP_opt.t3_Generation_VS_Single == 0 
     switch APP_opt.t3_PlotOpt_B
        case 1 
            Plot_InterGen_G1 ;
        case 2
            app.TextOUT.Value = sprintf('%s\n%s', 'Intergenerational analysis can only be done for "Line Type" 1 and 3!', '... Idle ...');                 
        case 3
            Plot_InterGen_G3 ;
        case 4     
            app.TextOUT.Value = sprintf('%s\n%s', 'Intergenerational analysis can only be done for "Line Type" 1 and 3!', '... Idle ...');                 
     end    
end


end