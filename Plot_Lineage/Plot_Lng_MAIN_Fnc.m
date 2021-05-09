function Plot_Lng_MAIN_Fnc(app)
%
%Plot_Lng_MAIN_Fnc = decide which plot function to run according to the
%   option choosen by the user
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

global APP_opt ;

% *** PLOT single clone_List.mat ******************************************
if APP_opt.t3_Generation_VS_Single == 1

    switch APP_opt.t3_PlotOpt_B
        case {1, 2, 3}
            if APP_opt.t3_PlotOpt_A <= 2
                Plot_Lng_B123;    
            elseif APP_opt.t3_PlotOpt_A == 3
                Plot_Lng_M123;
            end
        case 4     
            if APP_opt.t3_PlotOpt_A <= 2
                Plot_Lng_B4;        
            elseif APP_opt.t3_PlotOpt_A == 3
                Plot_Lng_M4; 
            end
    end

% *** ANALYSE Intragenerational ******************************************* 
elseif APP_opt.t3_Generation_VS_Single == 0 
     switch APP_opt.t3_PlotOpt_B
        case 1 
            Plot_InterGen_G13 ;
        case 2
            app.TextOUT.Value = sprintf('%s\n%s', 'Intergenerational analysis can only be done for "Line Type" 1 and 3!', '... Idle ...');                 
        case 3
            Plot_InterGen_G13 ;
        case 4     
            app.TextOUT.Value = sprintf('%s\n%s', 'Intergenerational analysis can only be done for "Line Type" 1 and 3!', '... Idle ...');                 
     end    
end


end