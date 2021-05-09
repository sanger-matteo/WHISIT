function Display_TrackAid
%
%Display_TrackAid - Show a colored frame around the figure to help
%   understand if the currently displayd frame was alread tracked or not
%   for the currently selected cell-track
%
%   frame is RED if it is not yet tracked
%   frame is BLU if it is already tracked
%
%   CellTracks = a global variable where that stores the information of
%   every cell track. 
%   scc = store the ID-number of the currently selected cell-track
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-

global APP_opt;	    global CellTracks;     global scc;

% Find the column position for the scc-th cell in CellTracks cell array       
idx_scc = find(cell2mat(CellTracks( 1,: )) == scc );     

if CellTracks{2,idx_scc}(APP_opt.t5_ff , 1:2) == [0,0] 
    clr_check = [.9 .2 .0];   % RED, frame yet not tracked
else
    clr_check = [.0 .8 .9];   % BLU, frame already tracked
end

hold on;
plot([1, APP_opt.t5_fr_W, APP_opt.t5_fr_W, 1, 1],...
     [1, 1, APP_opt.t5_fr_H, APP_opt.t5_fr_H, 1], ...
     'LineWidth', 1.25, 'Color', clr_check, 'Tag', 'RatioKeeper')
hold off;

end