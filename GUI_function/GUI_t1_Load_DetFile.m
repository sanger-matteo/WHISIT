function GUI_t1_Load_DetFile(app)
%GUI_t1_Load_DetFile = upon bush button, the user provide a Det.mat file
%   generated via OUFTI. THe function extract the path and filenames and
%   store them
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

global APP_opt ;

[FileName, PathName] = uigetfile({'*.mat'},'Select Oufti Detection.mat file'); 
if FileName == 0 
    FileName = [] ;
    PathName = [] ;
    app.TextOUT.Value = sprintf('\n%s\n%s',  'No detection.mat file provided !!!');
    return
end                        
APP_opt.t1_path_Det = PathName ;
APP_opt.t1_file_Det = FileName ;
app.t1_Edit_Det_file.Value = [APP_opt.t1_path_Det , APP_opt.t1_file_Det];




