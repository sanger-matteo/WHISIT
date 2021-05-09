function fr_backward( app, n )
%
%fr_backward( app, n ) - Move the currently selected frame BACKWARD of -n  
%   and show update the figure.
%
%   The function allow to select the frame APP_opt.t5_ff +n and update the
%   figure with the correct current frame as well as all other visual 
%   options (tracks, outllines, frame,...)
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

global APP_opt ;

if APP_opt.START_t5 == 1      % only if we are actively tracking
    
    % If number is above 1 (there is no frame 0 or negative ones)
    if APP_opt.t5_ff -n >= 1
        APP_opt.t5_ff = APP_opt.t5_ff -n ;
        app.t5_Edit_Frame.Value = APP_opt.t5_ff;         
        Display_BF_Frame;

    else         % remain fixed on the first frame 
        APP_opt.t5_ff = 1 ;
        app.t5_Edit_Frame.Value = APP_opt.t5_ff;         
        Display_BF_Frame;
        
    end

    ReFresh_Frame;           % REFRESH and update displayed frame
    
end


end




