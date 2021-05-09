function Find_FounderCell(app)


global APP_opt ;

app.t3_FounderCell_DropDown.Items = {} ;        % List all founders in DropDown menu
APP_opt.t3_PlotOpt_NLineage =  [];              % Specific cell Lineage to plot

load([APP_opt.t3_path_cloneList, APP_opt.t3_file_cloneList]);   % Load clone_List

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
    return
else
    app.t3_FounderCell_DropDown.Items = founders;     
    APP_opt.t3_PlotOpt_NLineage = founders{1} ;
end



end