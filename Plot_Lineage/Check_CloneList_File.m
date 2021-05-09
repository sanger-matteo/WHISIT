function warning = Check_CloneList_File(app)
%
%Check_CloneList_File = Check that the file provided by the user is a
% 	correct clone_List.mat files.If it does not comply to the criteria a 
%   warning is raised.
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


global APP_opt ;
warning = 0;            % warning variable

% check if  clone_List.mat file is given
if isempty(APP_opt.t3_path_cloneList) | isempty(APP_opt.t3_file_cloneList)
    app.TextOUT.Value = sprintf('\n%s',  'No cloneList.mat file provided!');
    warning = 1;
    return
else
    % Load WHISIT_param from clone_List.mat
    load([APP_opt.t3_path_cloneList, APP_opt.t3_file_cloneList], 'WHISIT_param' );
end


% Check that we have all correct data to create clone_List
varlist = who('-file', [APP_opt.t3_path_cloneList, APP_opt.t3_file_cloneList]);
check_Track = any(strcmp(varlist,'cellTrack')) ;
check_LG = any(strcmp(varlist,'clone_List')) ;

if check_Track == 0 | check_LG == 0  |  WHISIT_param.choice_ManualTrack == 0
    % if anything is missing, then we cannot create clone list
    app.TextOUT.Value = sprintf('\n%s',  'File provided is invalid cloneList.mat!!!');       
    warning = 1;
    return;
end


% ----- Create variable APP_opt.algorithm ---------------------------------
% If execution arrive at this point, we have a "proper" cloneList. We can 
% now check WHISIT_param to find out which algorithm was used during
% analysis and stores it in APP_opt.algorithm
if WHISIT_param.choice_M2P == 1
    APP_opt.algorithm = 1 ;
end
if WHISIT_param.choice_PL == 1
    APP_opt.algorithm = 2 ;
end
if WHISIT_param.choice_AIS == 1
    APP_opt.algorithm = 3 ;
end


end % MAIN fnc








