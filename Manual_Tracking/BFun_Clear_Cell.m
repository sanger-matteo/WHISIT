function  BFun_Clear_Cell(app)
%
%BFun_Clear_Cell = Reinitialize a specific cell-track to the empty state.
%
%   The function takes the currently selected scc-th cell-track and
%   reinitialize the tracked coordinates' list ( CellTracks{2,#) ) with 
%   zeros, as in the original state for newly created cell-tracks.   
%   The ID-number and RGB color are kept unchanged
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

% Find the column position for the scc-th cell in CellTracks cell array 
idx_scc = find(cell2mat(CellTracks( 1,: )) == scc);

% Initialize CellTracks at column scc-th, with ID- number, zeros and empty RGB color
CellTracks{ 1, idx_scc } = scc ;
CellTracks{ 2, idx_scc } = zeros( size(APP_opt.t5_srcFiles_BF,1) ,2);   
CellTracks{ 3, idx_scc } = [] ;

if APP_opt.START_t5 == 1     % only if we started manual tracking   
    ReFresh_Frame;           % REFRESH and update displayed frame
end

end





