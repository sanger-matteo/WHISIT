function Dis_Enable_PlotOptions(app, warning, datatype )
%
%Dis_Enable_PlotOptions = enable or disable rather large sets of options 
%   in the GUI interface
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

global APP_opt ;

% Enable plotting options if there was no warning or error while loading

if datatype == 1         % --- Single Lineage / File selection ---
    if warning == 0
        act = 'on';
    elseif warning >= 1
        act = 'off';
    end
    app.PLOT_t3.Enable = act;               
    app.t3_PlotHystogram.Enable = act;
    % Enable Radio DropDown List
    app.t3_DropDownValue.Enable = act;
    % Enable Radio Plot
    app.LineageTreeButton.Enable = act;
    app.IndividualClonesButton.Enable = act;
    app.ManualTracksButton.Enable = act;
    % Enable Radio Line
    app.singlevalueButton.Enable = act;
    app.bivalueButton.Enable = act;
    app.singlevalueRatioButton.Enable = act;
    app.segmentationButton.Enable = act;    

    
elseif datatype == 2      % --- Intergenerational / Folder selection ---
    if warning == 0 
        act = 'on';
    elseif warning >= 1
        act = 'off';
    end
    app.PLOT_t3.Enable = act;               
    app.t3_PlotHystogram.Enable = act;
    % Enable Radio DropDown List
    app.t3_DropDownValue.Enable = act;
    % Enable Radio Plot -------------------- always disabled
    app.LineageTreeButton.Enable = 'off';
    app.IndividualClonesButton.Enable = 'off';
    app.ManualTracksButton.Enable = 'off';
    % Enable Radio Line
    app.singlevalueButton.Enable = act;
    app.bivalueButton.Enable = act;
    app.singlevalueRatioButton.Enable = act;
    app.segmentationButton.Enable = act;        
end


% --- For both Intergenerational and Single Lineage -----------------------
% If there was no warning, then we opened the file(s) and found what
% algorithm was used. We can we enable/disable some specific plot options
if warning == 0 
  switch APP_opt.algorithm
    case 1 | 2      % for algorithm M2P and PS
        app.singlevalueRatioButton.Enable = 'on';
        app.bivalueButton.Enable = 'on';
    case 3          % for algorithm AIS
        app.singlevalueRatioButton.Enable = 'off';
        app.bivalueButton.Enable = 'off';
        
        % Change selection if presently it is on an inactive radio button        
        value = app.t3_RadioBox_LineType.SelectedObject;
        N_value = strsplit(value.Text, '-'); 
        if str2num(N_value{1}) == 2 | str2num(N_value{1}) == 3
            app.t3_RadioBox_LineType.SelectedObject.Value = 0;
            % Update all plot options
            Update_PlotOptions('All_opt', app);
        end
  end
end


end % MAIN func






