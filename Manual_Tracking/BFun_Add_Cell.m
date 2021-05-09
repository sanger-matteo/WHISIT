function BFun_Add_Cell(app)
%
%BFun_Add_Cell - Add a new empty cell-track to our list.
%
%   The function add a new element at the end of the cell array CellTracks.
%   The ID-number of the new cell-track is defined as the last ID-number +1
%   which guarantee that cell-tracks are addedd incrementally and there
%   won't be cell-tracks sharing the same ID-number
%
%   The global variable scc is also updated to store the ID-number of the
%   newly created cell-track
%
%   CellTracks = a global variable where that stores the information of
%   every cell track. 
%   scc = store the ID-number of the currently selected cell-track
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-

global APP_opt;	    global CellTracks;     global scc;

% create a list of all ID-numbers stored in CellTracks 
list_IDct = cell2mat(CellTracks( 1,: ));
% selected the last ID-number on the list, add 1 to create a new cell-track
scc = list_IDct(end)+1 ;

% Initialize CellTracks at column scc-th, with ID- number, zeros and empty RGB color
CellTracks{ 1, scc } = scc ;
CellTracks{ 2, scc } = zeros( size(APP_opt.t5_srcFiles_BF,1) ,2);   
CellTracks{ 3, scc } = [] ;

% Add scc-th cell-track to GUI "Select Cell" DropDown list (as STRING cell array format)
new_list = [list_IDct, scc];
app.t5_SelectCellDropDown.Items =...
    cellfun(@mat2str, mat2cell( new_list, ...
    1 , ones( 1,length(new_list)) ),'UniformOutput',false);      
% Update GUI "Select Cell" to point to scc-th cell
app.t5_SelectCellDropDown.Value = {mat2str(scc)}; 

if APP_opt.START_t5 == 1     % only if we started manual tracking   
    ReFresh_Frame;           % REFRESH and update displayed frame
end

end





