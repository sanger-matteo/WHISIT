function  BFun_Remove_Point(app)
%
%BFun_Clear_Cell = Cancel all coordinates tracked in current cell-track in
% current frame and all following one untill end of movie.
%
% THe coordinates from time point APP_opt.t5_ff up to end of array are made
% into (0;0) coordinates.
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

% Total xy_coordinates' array length
ptsLen = size(CellTracks{2, idx_scc}, 1);      
% Zeros all coordinates from time point ff up to end of array 
CellTracks{2, idx_scc}(APP_opt.t5_ff : ptsLen ,:) = ...
            zeros( length(APP_opt.t5_ff : ptsLen) , 2);

if APP_opt.START_t5 == 1     % only if we started manual tracking   
    ReFresh_Frame;           % REFRESH and update displayed frame
end

end





