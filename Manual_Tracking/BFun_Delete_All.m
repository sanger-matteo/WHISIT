function  BFun_Delete_All(app)
%
%BFun_Delete_All - Delete all cell-tracks from the CellTracks array 
%
%   The function delete all the cell-tracks from the CellTracks array and
%   reinitialize a new empty first element to be ready to restart tracking
%   completely anew
%
%   CellTracks = a global variable where that stores the information of
%   every cell track. 
%   scc = stores the ID-number of the currently selected cell-track
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-

global APP_opt;	    global CellTracks;     global scc;

% Delete ALL the elements
CellTracks(:) = [] ;

% CellTracks array must always have at least one element
scc = 1 ;
% simply initialize a new, empty first element
CellTracks{ 1, scc } = scc ;
CellTracks{ 2, scc } = zeros( size(APP_opt.t5_srcFiles_BF,1) ,2);
CellTracks{ 3, scc } = [] ;

% update the select cell DropDown menu
app.t5_SelectCellDropDown.Items = {'1'};
app.t5_SelectCellDropDown.Value = {mat2str(scc)};

if APP_opt.START_t5 == 1     % only if we started manual tracking   
    ReFresh_Frame;           % REFRESH and update displayed frame
end

end




