function fr_forward( app, n )
%
%fr_forward( app, n ) - Move the currently selected frame FORWARD of +n  
%   and show update the figure.
%
%   The function allow to select the frame APP_opt.t5_ff -n and update the
%   figure with the correct current frame as well as all other visual 
%   options (tracks, outllines, frame,...)
%
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-

global APP_opt ;

if APP_opt.START_t5 == 1      % only if we are actively tracking
    
    % If number is belowe the last frame in the movie
    if APP_opt.t5_ff+n <= length(APP_opt.t5_srcFiles_BF)
        APP_opt.t5_ff = APP_opt.t5_ff +n;
        app.t5_Edit_Frame.Value = APP_opt.t5_ff;         
        Display_BF_Frame;

    else         % remain fixed on the last frame
        APP_opt.t5_ff = length(APP_opt.t5_srcFiles_BF);
        app.t5_Edit_Frame.Value = APP_opt.t5_ff;        
        Display_BF_Frame;
               
    end

    ReFresh_Frame;           % REFRESH and update displayed frame

end


end





