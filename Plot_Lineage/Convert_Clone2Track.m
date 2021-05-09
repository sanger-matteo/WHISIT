function [track_List , info_Track] = Convert_Clone2Track(clone_List)
%

% Gather general information of all clones in the List
TrackInfo = cell(3,length(clone_List));
for cc = 1 : size(clone_List,2)
    TrackInfo{1,cc} = clone_List{cc}{1}.ID_ManualTrack ;     % cell-track ID number of the clone
    TrackInfo{2,cc} = clone_List{cc}{1}.fr_birth ;           % when is the first frame of the clone
    TrackInfo{3,cc} = clone_List{cc}{1}.fr_last ;            % when is the last frame of the clone
end
% Find the "unique" cell-tracks in the clone_List
uniqueTracks = unique( cell2mat(TrackInfo(1,:)) ) ;

% Combine the different clones that have come from the same manual track,
% making sure to combine in the correct temporal order
track_List = cell( 1 , length(uniqueTracks) );
info_Track = cell( 3 , length(uniqueTracks) );
for ii = 1 : length( uniqueTracks )
    trk = uniqueTracks(ii) ; 
    % Find the position of the clones (idx_trk) as well as the time "order"
    % in which they have to be concatenated (idx_birth)
    idx_trk = find( cell2mat(TrackInfo(1,:)) == trk) ;
    [O,idx_birth] = sort( cell2mat(TrackInfo(2,idx_trk)) ) ;
    
    % Concatenate clones according to the sorted idx_birth order 
    track_List{1, ii} = [ clone_List{ idx_trk(idx_birth) } ] ;     
    if length(idx_trk) == 1
        info_Track{1, ii} = cell2mat( TrackInfo( 2, idx_trk(idx_birth(1)) )) ;       % first tracked frame
        info_Track{2, ii} = [] ;                                                     % divisions
        info_Track{3, ii} = cell2mat( TrackInfo( 3, idx_trk(idx_birth(1)) )) ;       % last tracked frame
        
    elseif length(idx_trk) > 1
        info_Track{1, ii} = cell2mat( TrackInfo( 2, idx_trk(idx_birth(1)) )) ;       % first tracked frame
        info_Track{2, ii} = cell2mat( TrackInfo( 2, idx_trk(idx_birth(2:end)) )) ;   % divisions
        info_Track{3, ii} = cell2mat( TrackInfo( 3, idx_trk(idx_birth(end)) )) ;     % last tracked frame
    
    end        
end


end
