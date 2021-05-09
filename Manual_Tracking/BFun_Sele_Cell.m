function BFun_Sele_Cell(app)
%
%BFun_Sele_Cell - Change the currently selected cell 
%   User change the current cell-track to work on using DropDown menu. The
%   function read the value of the current DropDown menu and update the 
%   value of scc
%
%   ---> Notes about variables used:
%   CellTracks = a global variable where that stores the information of
%   every cell track. 
%   scc = store the ID-number of the currently selected cell-track
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

global APP_opt;	    global CellTracks;     global scc;

% Take from GUI "Select Cell" the selected cell-track
value = str2double(app.t5_SelectCellDropDown.Value) ;
% Find the column position for the scc-th cell in CellTracks cell array 
idx_scc =  find(cell2mat(CellTracks(1,:)) == value) ;
scc = CellTracks{1 , idx_scc} ;         % update the scc value


if APP_opt.START_t5 == 1     % only if we started manual tracking   
    ReFresh_Frame;           % REFRESH and update displayed frame
end

end



