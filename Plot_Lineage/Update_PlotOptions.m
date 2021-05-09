function Update_PlotOptions(oriopt, app)
%
% Update_PlotOptions = update the list of "Value to use for plotting"
%     (t3_DropDownValue) according to the value of the input oriopt. 
%
% oriopt = can carry a change for t3_PlotOpt_A or t3_PlotOpt_B
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

global APP_opt ;


% Refresh the values that determine what type of plot to create
if strcmp(oriopt, 'A_opt')  |  strcmp(oriopt, 'All_opt')
    value = app.t3_RadioBox_PlotType.SelectedObject;
    Str_value = strsplit(value.Text, '-'); 
    APP_opt.t3_PlotOpt_A = str2num(Str_value{1}(1));
end


if strcmp(oriopt, 'B_opt')  |  strcmp(oriopt, 'All_opt')
    value = app.t3_RadioBox_LineType.SelectedObject;
    Str_value = strsplit(value.Text, '-');
    APP_opt.t3_PlotOpt_B = str2num(Str_value{1}(1));  

    % Determine the DropDownLists according to the algorithm used
    switch APP_opt.algorithm
      case 1
        %(1) for algorithm M2P
        DDLists(:,1)= {{'1 - Whole Cell', '2 - Cytosol',...
                '3 - Old pole', '4 - New pole',...
                '5 - Lateral Membrane', '6 - Whole Membrane' } ,

               {'1 - Old pole | New pole', '2 - New Pole | Cytosol', ...
                '3 - Old Pole | Cytosol' , '4 - Poles | Cytosol',...
                '5 - All Membrane | Cytosol', '6 - Lateral Membrane | Cytosol',...
                '7 - Lateral Membrane | Poles' } ,

               {'1 - Cytosol / Poles'      , '2 - Poles / Cytosol' , ...
                '3 - Old pole / New pole'  , '4 - New pole / Old pole' , ...
                '5 - Cytosol / Membrane'   , '6 - Membrane / Cytosol' ,...
                '7 - Lateral Membrane / Poles'     , '8 - Poles / Lateral Membrane' } ,

               {'1 - Cell Segmentation', '2 - Axial Profile'}  };
      case 2
        % (2) for algorithm PL
        DDLists(:,2)= {{'1 - Whole Cell', '2 - Membrane', ...
                '3 - Old pole', '4 - New pole' } ,

               {'1 - Old pole | New_pole', '2 - New | Cytosol', ...
                '3 - Old pole | Cytosol' , '4 - Poles | Cytosol'} ,

               {'1 - Cytosol / Poles'      , '2 - Poles / Cytosol' , ...
                '3 - Old pole / New pole'  , '4 - New pole / Old pole' } ,                

               {'1 - Cell Segmentation', '2 - Axial Profile'}  };
      case 3
        % (3) for algorithm AIS
        DDLists(:,3)= {{'1 - Whole Cell' } ,
               {''} ,
               {''} ,
               {'1 - Segmentation (Align Bottom)'  , '2 - Segmentation (Align Middle)'  , '3 - Segmentation (Align Top)', ...
                '4 - Axial Profile (Align Bottom)' , '5 - Axial Profile (Align Middle)' , '6 - Axial Profile (Align Top)'} };
    end

    % Update the list in DropDownValue 
    app.t3_DropDownValue.Items = DDLists{APP_opt.t3_PlotOpt_B, APP_opt.algorithm};
    % Then rest it to the first element of the list
    app.t3_DropDownValue.Value = app.t3_DropDownValue.Items{1};
    Str_value = strsplit(app.t3_DropDownValue.Items{1}, '-');
    APP_opt.t3_PlotOpt_C = str2num(Str_value{1}(1)); 
    
end % if B_opt

end




