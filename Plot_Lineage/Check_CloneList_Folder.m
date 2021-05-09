function Check_CloneList_Folder(app)
% Check_CloneList_Folder = Check that the folder provided by the user
%   contain proper clone_List.mat files.
%
%   Those who do not comply to the criteria are removed and not considrered
%   for analysis. In this case a warning is rised.
%   If no file is left, none is a clone_List, then an error is raised
%
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-

global APP_opt ;

Warning = 0 ;       % warn if any .mat files is not in proper format
counter = [] ;      % keep the index of .mat files that are not clone_List format
W_check = [] ;      % store each WHISIT_parameters and check that they are all the same
    
% Access folder and index all the file inside the folder that ends with .mat
srcFiles = dir( fullfile([APP_opt.t3_path_cloneFolder '/' APP_opt.t3_fold_cloneFolder],'*.mat'));
% If no .mat files is inside, rise an error and return
if isempty(srcFiles)
    app.TextOUT.Value = sprintf('\n%s\n%s',  'No .mat file found in the folder!');
    app.TextOUT.BackgroundColor = [0.75 0.3 0.3] ;
    % reset the variables
    app.t3_Edit_cloneFolder.Value = '' ;
    APP_opt.t3_intergen_srcFiles(:) = struct();
    return
end

jj = 0;
cell_srcFiles = struct2cell( srcFiles );
for ii = 1 : size( cell_srcFiles ,2)
    % Check if variables exsist inside .mat files without loading (this is the fastest solution)
    VarInfo = whos('-file', [cell_srcFiles{2,ii},'/',cell_srcFiles{1,ii}]) ;
    if ~isempty(VarInfo)  &&   size(VarInfo,1) == 3  && ...
            strcmp(VarInfo(1,1).name, 'WHISIT_parameters')  && ...
            strcmp(VarInfo(2,1).name, 'cellTrack')  && ...
            strcmp(VarInfo(3,1).name, 'clone_List')
        % The file is a clone_List.mat with possibly correct data inside
        load([cell_srcFiles{2,ii},'/',cell_srcFiles{1,ii}] , 'WHISIT_parameters')
        jj = jj+1;
        % store the WHISIT_parameters
        W_check(:,jj) = [ WHISIT_parameters.choice_AIS; ...
                          WHISIT_parameters.choice_M2C;  ...
                          WHISIT_parameters.choice_PL];                     
    else
        % Else, it is not a clone_List.mat file
        counter = [counter , ii] ;      % save index
        Warning = 1 ;
    end
end


% Check that the WHISIT_parameters are all the same. We can do intergenerational 
% analysis if all the clone_List where analysed the same way.
% any() = check which rows have non-zero elements. If all clone_Lists where
% analysed the same, then only one row return 1.
if sum(any(W_check,2)) >= 2
    app.TextOUT.Value = sprintf('\n%s\n%s',  'clone_Lists have not all been analysed with same algorithm!');
    app.TextOUT.BackgroundColor = [0.75 0.3 0.3] ;
    % reset the variables
    app.t3_Edit_cloneFolder.Value = '' ;
    APP_opt.t3_intergen_srcFiles(:) = struct();
    return
end    
    
srcFiles(counter,:) = [];               % all indexed .mat files are removed

if isempty(srcFiles)
    app.TextOUT.Value = sprintf('\n%s',  ['No .mat files inside folder is a clone_List!']);
    app.TextOUT.BackgroundColor = [0.75 0.3 0.3] ;
    return;
end

% Save a list of the .mat files that are clone_List and fit for intergenerational analysis 
APP_opt.t3_intergen_srcFiles = srcFiles ;

if Warning == 1
    app.TextOUT.Value = sprintf('\n%s',  ['Not all .mat files inside folder are clone_Lists!']);
else 
    app.TextOUT.Value = sprintf('\n%s',  ['... Idle ...']);
end


end % MAIN fnc








