function BFun_SaveTracks(app)
%
%BFun_SaveTracks - Save CellTracks array into a tab separated .txt file
%
%   The function takes the list of tracked coordinates and converts them
%   into a matrix. Each cell-track occupy two columns, in order [x,y], and
%   the number of rows is equal to the number of BF images. The first row 
%   contains the ID-number of the track.
%   An example of the .txt file output is shown below:
%
%   1	1	2	2	3	3	4	4	5	5
%   0	0	29	28	0	0	38	46	0	0
%   0	0	25	36	49	37	28	69	0	0
%   20	47	47	54	26	53	63	44	0	0
%   47	47	41	29	53	68	0	0	39	34
%   62	37	0	0	60	43	0	0	41	57
%   50	30	0	0	14	63	0	0	35	50
%   0	0	0	0	32	77	0	0	42	45
%   0	0	0	0	0	0	0	0	0	0
%
%   Not tracked/empty coordinates have value 0.0
%   (N.B.: there cannot be a pixel position (0;0))
%
%   CellTracks = a global variable where that stores the information of
%   every cell track. 
%   scc = stores the ID-number of the currently selected cell-track
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-

global APP_opt;	    global CellTracks;     global scc;

cellID = [];
% Convert to a matrix storing CellTracks coordinates 
XY_tracks = cell2mat(CellTracks(2,:));
% Create a list of all ID-numbers stored in CellTracks 
list_IDct = cell2mat(CellTracks( 1,: ));

% we do not save RGB colors of tracks --> CellTracks(3,:)

% Delete any track that is empty (all points are (0,0)).
% We go in reverse order in for-loop to avoid skipping CellTracks elements
% and errors when trying to access positions that have been eliminated.
for ii = size(CellTracks, 2) : 1      
    if ~any(CellTracks{2,ii},1)       % if NOT all array elements are zeros
        CellTracks( : , ii ) = [] ;   % delete cell track
    end
end

% create a single row array with ID-numeber identifing cell-tracks
for ii = 1 : length(list_IDct)
   cellID = [cellID,  list_IDct(ii), list_IDct(ii)] ;
end

% --- SAVE cloneList.mat
app.TextOUT.Value = sprintf('\n%s',  'Saving clone_List file ... ');

MAT = [cellID; XY_tracks];
if isempty(APP_opt.t5_exp_name)
    file_name = [APP_opt.t5_path_BF ,'/', 'tracks.txt'];
else
    file_name = [APP_opt.t5_path_BF ,'/', APP_opt.t5_exp_name '_tracks.txt'];
end
dlmwrite(file_name , MAT, 'delimiter','\t') ;           % ,'precision',3 );

app.TextOUT.Value = sprintf('\n%s',  '... Idle ...');

end




