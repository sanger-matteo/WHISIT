function warning = Check_CloneList_Folder(app)
%
% Check_CloneList_Folder = Check that the folder provided by the user
%   contain proper clone_List.mat files.
%   Those who do not comply to the criteria are removed and not considrered
%   for analysis. In this case a warning is raised. If no file is left,
%   (none is a clone_List) then an error is rised.
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


global APP_opt ;

warning = 0 ;       % warn if any .mat files is not in proper format
counter = [] ;      % keep the index of .mat files that are not clone_List format
W_check = [] ;      % store each WHISIT_param and check that they are all the same

% Reset founder DropDown List, because we do not use it for intergenerational analysis
app.t3_FounderCell_DropDown.Items = {} ;        % List all founders in DropDown menu
APP_opt.t3_PlotOpt_NLineage =  [];              % Specific cell Lineage to plot

if isempty(APP_opt.t3_fold_cloneFolder)  |  isempty(APP_opt.t3_path_cloneFolder)
    warning = 1;
    return
else

    
% Access folder and index all the file inside the folder that ends with .mat
srcFiles = dir( fullfile([APP_opt.t3_path_cloneFolder '/' APP_opt.t3_fold_cloneFolder],'*.mat'));
% If no .mat files is inside, rise an error and return
if isempty(srcFiles)
    app.TextOUT.Value = sprintf('\n%s',  'No .mat file found in the folder!');
    app.TextOUT.BackgroundColor = [0.75 0.3 0.3] ;
    % reset the variables
    app.t3_Edit_cloneFolder.Value = '' ;
    APP_opt.t3_intergen_srcFiles(:) = struct();
    warning = 1;
    return
end


jj = 0;
cell_srcFiles = struct2cell( srcFiles );
for ii = 1 : size( cell_srcFiles ,2)
    % Check if variables exsist inside .mat files without loading (this is the fastest solution)
    VarInfo = whos('-file', [cell_srcFiles{2,ii},'/',cell_srcFiles{1,ii}]) ;
    if ~isempty(VarInfo)  &&   size(VarInfo,1) == 3  && ...
            strcmp(VarInfo(1,1).name, 'WHISIT_param')  && ...
            strcmp(VarInfo(2,1).name, 'cellTrack')  && ...
            strcmp(VarInfo(3,1).name, 'clone_List')
        % The file is a clone_List.mat with possibly correct data inside
        load([cell_srcFiles{2,ii},'/',cell_srcFiles{1,ii}] , 'WHISIT_param')
        jj = jj+1;
        % store the WHISIT_param
        W_check(:,jj) = [ WHISIT_param.choice_AIS; ...
                          WHISIT_param.choice_M2P;  ...
                          WHISIT_param.choice_PL];                     
    else
        % Else, it is not a clone_List.mat file
        counter = [counter , ii] ;      % save index
        warning = 1 ;
    end
end


% Check that the WHISIT_param are all the same. We can do intergenerational 
% analysis only if all the clone_List where analysed with same algorithm.
% any() = check which rows have non-zero elements. If all clone_Lists where
% analysed the same, then only one row return 1.
if sum(any(W_check,2)) >= 2
    app.TextOUT.Value = sprintf('\n%s',  'clone_Lists have not all been analysed with same algorithm!');
    app.TextOUT.BackgroundColor = [0.75 0.3 0.3] ;
    % reset the variables
    app.t3_Edit_cloneFolder.Value = '' ;
    APP_opt.t3_intergen_srcFiles = struct([]);     % Initialize empty struct
    warning = 1;
    return
end    
   
srcFiles(counter,:) = [];               % all indexed .mat files are removed
if isempty(srcFiles)
    app.TextOUT.Value = sprintf('\n%s',  ['No .mat files inside folder is a clone_List!']);
    app.TextOUT.BackgroundColor = [0.75 0.3 0.3] ;
    warning = 1;
    return;
end

% Store a list of the .mat files that are clone_List and fit for intergenerational analysis 
APP_opt.t3_intergen_srcFiles = srcFiles ;


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








