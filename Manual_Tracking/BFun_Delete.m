function  BFun_Delete(app)
%
%BFun_Delete - Delete a specific cell-track from the CellTracks array 
%
%   The function takes the currently selected scc-th cell-track and
%   removes it from the cell array CellTracks
%
%   CellTracks = a global variable where that stores the information of
%   every cell track. 
%   scc = stores the ID-number of the currently selected cell-track
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-

global APP_opt;	    global CellTracks;     global scc;

% Find the column position for the scc-th cell in CellTracks cell array 
idx_scc = find(cell2mat(CellTracks( 1,: )) == scc);

% remove the specific CellTrack by deliting entire column
CellTracks( :, idx_scc ) = [] ;         
% this ID-number disappear. Since ID numbers are created incrementally,
% there is no possibility of having 2 equal IDs 

% create an updated list of all ID-numbers stored in CellTracks 
list_IDct = cell2mat(CellTracks( 1,: ));
% Update scc-th cell-track to GUI "Select Cell" DropDown list (as STRING cell array format)
new_list = list_IDct;
if ~isempty(new_list)
    app.t5_SelectCellDropDown.Items =...
         cellfun(@mat2str, mat2cell( new_list, ...
         1 , ones( 1,length(new_list)) ),'UniformOutput',false);

end

% --- UPDATE Selected Cell and DropDown menu ------------------------------
% Since we removed the scc-th track, we now update scc to point to the
% previous or next cell, depending where the deleted track was positioned
% in the array
if isempty(list_IDct)
    % if the only element left was deleted, the function simply reinitialize
    % a new first element with ID-number = 1
    scc = 1 ;
    % simply initialize a new first element
    CellTracks{1 , scc} = scc ; 
    CellTracks{2 , scc} = zeros(size(APP_opt.t5_srcFiles_BF,1) ,2);
    CellTracks{3 , scc} = [] ;

    app.t5_SelectCellDropDown.Value = {mat2str(scc)};       % Update GUI "Select Cell" to point to scc-th cell
       
elseif length(list_IDct) == 1
    % simply select the only element left available in CellTracks array
    scc = list_IDct(1) ;
    app.t5_SelectCellDropDown.Value = {mat2str(scc)}; 
    
elseif idx_scc == 1
    % if we delete first element in CellTracks array, take the new "first" element
    scc = list_IDct(1) ;
    app.t5_SelectCellDropDown.Value = {mat2str(scc)}; 
    
else
    % for any other position deleted, take the previous element in CellTracks array
    scc = list_IDct(idx_scc - 1) ;
    app.t5_SelectCellDropDown.Value = {mat2str(scc)};
end

if APP_opt.START_t5 == 1     % only if we started manual tracking   
    ReFresh_Frame;           % REFRESH and update displayed frame
end

end





