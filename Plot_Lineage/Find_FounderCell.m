function warning = Find_FounderCell(app)
%
%Find_FounderCell = find the all the founder cells present in a clone_List
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


global APP_opt ;
warning = 0;            % warning variable

app.t3_FounderCell_DropDown.Items = {} ;        % List all founders in DropDown menu
APP_opt.t3_PlotOpt_NLineage =  [];              % Specific cell Lineage to plot

% check if clone_List.mat file is given
if isempty(APP_opt.t3_path_cloneList) | isempty(APP_opt.t3_file_cloneList)
    app.TextOUT.Value = sprintf('\n%s',  'No cloneList.mat file provided!');
    warning = 1;
    return
end

% Check that there is a clone_List inside the .mat file
load([APP_opt.t3_path_cloneList, APP_opt.t3_file_cloneList]);   % Load clone_List
if exist('clone_List') == 0
    app.TextOUT.Value = sprintf('\n%s',  'File .mat provided do not contain a cloneList !');
    warning = 1;
    return    
end

founders = {};
jj = 1 ;
for cc = 1 : size(clone_List ,2)    
    % split the ID_clone name string
    strID = strsplit( clone_List{cc}{1}.ID_clone , {'.'}) ;
    % If after the '.' there is nothing, cc-th clone is a founder cell
    if isempty(strID{2})
        founders{jj,1} = [repmat('0', 1, 3-length(num2str(cc))), num2str(cc)];
        jj = jj + 1 ;
    end    
end

if isempty(founders)
    % Warning, the file do not contain any founder
    APP_opt.t3_PlotOpt_NLineage =  [];                  % Cell Lineage to plot
    app.TextOUT.Value = sprintf('\n%s',  ['!!! Error: Clone list has no founder cell !!!']);
    app.TextOUT.BackgroundColor = [0.75 0.3 0.3] ;      % red bkgr
    warning = 1;
    return
else
    app.t3_FounderCell_DropDown.Items = founders;     
    APP_opt.t3_PlotOpt_NLineage = founders{1} ;
    warning = 0;
end



end