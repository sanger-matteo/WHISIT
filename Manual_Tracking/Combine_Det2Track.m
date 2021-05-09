function Combine_Det2Track
%
%Combine_Det2Track - The function use the manual tracking data and combine
%  it with the detection done in Oufti. The tracking performed using WHISIT
%  is used to reorder the cellList to create a new time tracked cellsList.
%  The algorithm also track the polarity of the cell, from each new cell
%  division onward, in order to keep track which is the new and which is
%  the old pole
%
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-


% ----- INPUT -------------------------------------------------------------
% Track.txt = the function load the tab-separated track file. Each single
%      track is a pair columns (X,Y) and they alternates one after the
%      other crating a single matrix. Tracks stored in cTrack.
%          1   1   4   4   5   5   6   6
%          X   Y   X   Y   X   Y   X   Y
%      Untracked positions are filled with 0 (since no pixel position 0
%      can exist). The first row contains the track ID number.
%
% cTracks = [CellTracks in t5_Manual_Tracking, it is our main variable that
%       store all information needed to combine detection with manual track
%       It is a cell-array, where each column correspond to a single
%       specific cell-track. Each column is organized as follow:
%     row 1 - ID-number of the cell-track
%     row 2 - point coordinates, stored as [x,y] matrix of length LenStack
%     row 3 - ---- empty ---- (RGB color in manual tracking)
%
%     The following are all a colum array (1 x total_number_frames) 
%     row 4 - stores at each time point the cellList ID number that
%             correspond to the cell tracked (or the closest one to the
%             tracked point). Untracked points are filled with zeros.
%     row 5 - stores at the time points when the cell-track was created the
%             cell-track ID-number of ancestor (mother cell) from which it
%             was generated (see N.B.-DIVISION below).
%     row 6 - stores cell-track ID-number of daughter cell(s) that divided
%             , generated from current cell-track (see N.B.-DIVISION below)
%   
% cellList = the main variable that contain all informations about the
%           detection done in Oufti. We assume detection was done with
%           "independent frames" option selected. Therefore, each cell has
%           been ramdomly assigned a number and there is no correlation
%           between different frames.
%
% ----- OUTPUT ------------------------------------------------------------
% The entire det.mat provided with the updated cellList and cellListN, and
% the addition of cellTrack
%
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-

%%
% N.B.-DIVISION:
% The user manually track ensuring that mother and daugther share cell
% outline ONE FRAME BEFORE DIVISION event. This means that on the first
% tracked frame they are not physically separate yet and mother-daughter
% have the same cell outline.
% This allow to automatically identify cell divisions of daughter cell and
% the corresponding mother cell by simply looking with whom does the cell
% share cell outline. (cTrack{4, cc} and M_coID).
% However, this also means real "cell life" for daugther cells BEGIN at the 
% SECOND frame that is manually tracked.


global APP_opt ;


%--- Load DET.mat file ---------------------------------------------------%
% Load the DETection.mat file from Oufti
load([APP_opt.t5_path_Det , APP_opt.t5_fileName_Det ]);
 

%--- Load TRACK.txt file -------------------------------------------------%
% Load a list of cell-tracks, which are stored as tab-separated .txt file
temp_T = load([APP_opt.t5_path_Track , APP_opt.t5_fileName_Track ]);     
t_IDs = temp_T(1,2:2:end);          % row array with ID-numbers
t_trk = temp_T(2:end,:) ;           % matrix with tracked points
[Rtrk, Ctrk] = size(t_trk);         % extract row and column numbers    
% convet the matrixes in cell-array and place them correctly in
% columns each representing a single cell-track
cTrack(1,:) = mat2cell( t_IDs , 1    , ones( 1,length(t_IDs)) );
cTrack(2,:) = mat2cell( t_trk , Rtrk , 2.*ones(1, Ctrk/2) );
cTrack(3,:) = cell(1,length(t_IDs));
clearvars temp_T t_IDs t_trk Rtrk Ctrk;



%% --- STEP 1 --- Cell-outline ID -----------------------------------------------------
% First part the algorithm identify the cellList outline corresponding to
% each cell-track at every time point and store the respective ID number
% (position at cellList(frame_num)) in variable --- cTrack{4, :}
                    
for cc = 1 : size(cTrack, 2)  
    
  % Store detection outline cell_ID that correspond to the manual tracking
  cTrack{4, cc} = zeros( size(cTrack{2,cc},1) ,1) ;  
  % We can avoid testing where there is not tracking done (0;0)
  t_range = find(cTrack{2,cc}(:,1)~=0)';  
  
  if ~isempty(t_range)                                      % is cell-track is empty, ignore it
  for ff = t_range(1) :1: t_range(end)
    % Store XY_t point for all cell-tracks in ff-th frame
    ff_X = cTrack{2, cc}(ff , 1) ;
    ff_Y = cTrack{2, cc}(ff , 2) ;
    % remove '0' coordinates - not-tracked points will not be considered
    ff_X(ff_X==0) = NaN;     
    ff_Y(ff_Y==0) = NaN;        
    dist = [];      % store distances bewteen cell outline and track point at ff frame

    % Go through all cell's outline detected in ff-th frame. Calculate
    % distances to find closest outline to the tracked point of cc-th cell-track
    for dd = 1 : length( cellList.meshData{ff} )
        if  ~isempty(cellList.meshData{ff}{dd})  && ...
            ~isempty( cellList.meshData{ff}{dd}.mesh )  && ...
            size(cellList.meshData{ff}{dd}.mesh,2) == 4 
            
            % Store XY_o points for kk_th outline and plot it
            X_o = [cellList.meshData{ff}{dd}.mesh(:,1) ; flipud(  cellList.meshData{ff}{dd}.mesh(:,3)) ];
            Y_o = [cellList.meshData{ff}{dd}.mesh(:,2) ; flipud(  cellList.meshData{ff}{dd}.mesh(:,4)) ];
            % Calculate R2 between each XY_t and all XY_o points of kk-th cell 
            % ontour outline and store the minumum value
            dist(dd) = min(double( sqrt( abs(ff_X - X_o).^2 + abs(ff_Y - Y_o).^2 ) ));

            if inpolygon(ff_X,ff_Y, X_o,Y_o)                % and lot cell outline in selected_clr if XY_t is inside it             
                isinc = 1 ;
                cTrack{4, cc}(ff) = dd ;
            else
                isinc = 0 ; 
            end
        end
    end %/for dd  

    % If XY_t was not inside any cell outline, find the closest cell outline 
    dist(dist == 0) = NaN ;         % "empty" cells give a distance of zero which creates an error
    idx_closest =  find(dist == min(dist) );     
    if isinc == 0  &&  ~isempty(idx_closest)  &&  cTrack{4, cc}(ff) == 0     % if not inside cell ...      
        cTrack{4, cc}(ff) = idx_closest ;
    end   

  end % ff 
  end % if t_range notempty
end % cc


%% --- STEP 2 --- Identify DIVision events -----------------------------------------------------
% Second part of the algorithm find the frame of birth for each cell-track
% and who is the mother/ancestor cell.
% - cTrack{5, cc} --- ancestor
% - cTrack{6, cc} --- daughter(s)
% N.B.: cell-track cannot share same outline unless. Therefore, we can
% identify division events as when two cell-tracks share the same ID number
% at the same time point

M_coID = cell2mat(cTrack(4,:));                             % matrix storing all cell-outline ID of all cell-tracks

for cc = 1 : size(cTrack, 2) 
        
    cTrack{5, cc} = zeros( size(cTrack{2,cc},1) ,1) ;       % cell-track ID-number of ancestor    
    cTrack{6, cc} = zeros( size(cTrack{2,cc},1) ,1) ;       % cell-track ID-number of daughter(s) 
    
    % time frame where cc-th cell-track was actually tracked (avoid testing at (0;0))
    t_range = find( M_coID(:,cc) ~= 0 );        
    
    if ~isempty(t_range)                                    % is cell-track is empty, ignore it
        
        % find who share the same cell-outline with cc-th cell. This is the ancestor cell-track. 
        % Simply, look at rows in M_coID at each time point.
        co_share = find( M_coID(t_range(1),:) == M_coID(t_range(1),cc));
        Fr_Brt = t_range(2);                                % frame of birth of cell-track    
        if length(co_share) > 2                             % ERROR - there can be only 2 cells-track sharing same outline

        elseif length(co_share) == 1                        % it is a lineage founder cell

        elseif length(co_share) == 2 
            % N.B.: co_chare and cc are indexes to the column of cTrack.
            % We want to store the cell-track ID number (cTrack{1, x })
            co_share(co_share == cc) = [];                  % we exclude current cc-th cell-track we examine
            cTrack{5, cc}(Fr_Brt)       = cTrack{1, co_share};
            cTrack{6, co_share}(Fr_Brt) = cTrack{1, cc};
        end
    end % if t_range notempty    
end    


%% --- STEP 3 --- COMBINE CellTrack and cellList -----------------------------------------------------
% Now we have what is necessary to combine manual tracking and cell
% detection to create a new time-tracked cellList (new_cList).
% IMPORTANT: the algorithm is reordering the new_cellList as well:
% - any cell outline that is not used/tracked, will be discarded
% - at each frame cell outlines are places in the array position that
%   match their track ID number (this will be useful to create clone_List)

% --- Indexing (pcd) and counter (cc) -------------------------------------
% cc  = counter used to go through each cell-track in for-loops. 
% pcd = store where to position the cc-th data in the new_cList.
% tid = store the track ID number of interested.
% It is easier to go through cell-tracks column-by-column using cc as index.
% However, columns order is unrelated to the track ID number, which we need
% to correctly reorder the new_cellList.
% Hence, for simplicity, we create a separate variable to carry the index 
% where to reposition data in new_cList: 
% >> pcd = length(new_cList{ff}) +1 ;

new_cList = cell(1,size(cellList.meshData,2));      % reordered tracked cellList
old_cList = cellList.meshData;                      % we give a shorter practical name

% Store the position of each cell-det in each frame-array at time point ff
pcd_trk(1, :) = cTrack(1, :) ; 
for cc = 1 : size(cTrack, 2)   
    pcd_trk{2, cc} = zeros( size(cellList.meshData,2) ,1);
end

for cc = 1 : size(cTrack, 2)                        % go through each cell-track      
  % time frame where cc-th cell-track was actually tracked (avoid testing at (0;0))
  t_range = find(cTrack{2,cc}(:,1)~=0)';
  
  if ~isempty(t_range)                                % is cell-track is empty, ignore it 
    for ff = t_range(1) : t_range(end)                % <<<----- CATCH_1 (at line ~370)
      if ~isempty( old_cList{ff}{M_coID(ff,cc)})  &&  size( old_cList{ff}{M_coID(ff,cc)}.mesh, 2) == 4

        pcd = length(new_cList{ff}) +1 ;
        pcd_trk{2,cc}(ff) = pcd ;
        
        new_cList{ff}{pcd} = old_cList{ff}{M_coID(ff,cc)};
        
        % [see N.B.-DIVISION above] we must take the second frame tracked 
        % as the moment cells become "independent", unless it is a founder
        if isempty( cTrack{5, cc}(cTrack{5, cc}(:)~=0 ) )          % it is a founder cell   
            new_cList{ff}{pcd}.track_ancestor = 0 ;                                  
            new_cList{ff}{pcd}.fr_birth = t_range(1) ;             % take first frame tracked
        else
            new_cList{ff}{pcd}.track_ancestor = cTrack{5, cc}(cTrack{5, cc}(:)~=0 ) ;  
            new_cList{ff}{pcd}.fr_birth = t_range(2) ;             % take the second frame tracked
        end
        
        new_cList{ff}{pcd}.tracks_offspring  = cTrack{6, cc}(cTrack{6, cc}(:)~=0 ) ;                           
        new_cList{ff}{pcd}.tracks_fr_off_div = find(cTrack{6, cc}(:)~=0) ;
        new_cList{ff}{pcd}.fr_last           = t_range(end) ;
        new_cList{ff}{pcd}.ID_ManualTrack    = cTrack{1, cc} ;
        new_cList{ff}{pcd}.polarity          = 1 ;
        
      end
    end % ff 
  end % if t_range notempty  
end % cc


%% --- STEP 4 --- Determine POLE at DIVision ---------------------------------------------------
% For every division event the algorithm identifies which will be the new- 
% and the old-pole in the two daugthers.
% We have 2 cell-outlines: two daughter cells at division frame. They have
% 2 poles each makinf 4 possible pairings and distances to calculate. The
% pair with the shorter distance determine the two new poles.
% ---> OLD POLE coordinates are at first position in cellList.mesh(1,_)
%
% N.B.: This method should be no problem to correctly identifying as long
% cells do not move too much and/or time lapses are frequent, there


for cc = 1 : size(cTrack, 2)                  % go through each cell-track 
    
    % time frame where cc-th cell-track was actually tracked (avoid testing at (0;0))
    t_range = find(cTrack{2,cc}(:,1)~=0)';

    if ~isempty(t_range)                      % is cell-track is empty, ignore it 

        % Find positions of cc-th cell-track  and its ancestor in new_cList
        pcd = pcd_trk{2,cc}(t_range(1));                             % find the position of cc-th cell-track in new_cList
        tid_Anc = new_cList{1,t_range(1)}{pcd}.track_ancestor ;      % find the track ID of Ancestor   
            
        if tid_Anc == 0                       % if it is a founder cell,...
            Fr_Brt = t_range(1) ;             % take first frame tracked

        elseif tid_Anc ~= 0                   % if it is not a progenitor   
            Fr_Brt = t_range(2) ;             % take the second frame tracked
            colum = find( cell2mat(pcd_trk(1,:)) == tid_Anc) ;         % find the colum position of Anc in manual track
            pcd_Anc = pcd_trk{2,colum}(t_range(1)) ;                       % find the Anc position in new_cList

            %%%--->>> insert here CHECK SUBPLOT 1.1 <<<---%%%  

            % CELL_1 pole 1 and 2 x-y coordinates
            C1P1x = new_cList{Fr_Brt}{pcd}.mesh(end,1);       %end
            C1P1y = new_cList{Fr_Brt}{pcd}.mesh(end,2);
            C1P2x = new_cList{Fr_Brt}{pcd}.mesh(1,1);         % 1
            C1P2y = new_cList{Fr_Brt}{pcd}.mesh(1,2);
            % CELL_2 pole 1 and 2 x-y coordinates
            C2P1x = new_cList{Fr_Brt}{pcd_Anc}.mesh(end,1);       %end
            C2P1y = new_cList{Fr_Brt}{pcd_Anc}.mesh(end,2); 
            C2P2x = new_cList{Fr_Brt}{pcd_Anc}.mesh(1,1);         % 1
            C2P2y = new_cList{Fr_Brt}{pcd_Anc}.mesh(1,2);
            % four distances are calculated [ P11-P21 , P11-P22  , P12-P21 , P12-P22 ]
            dist = [ sqrt( abs(C1P1x - C2P1x).^2 + abs(C1P1y - C2P1y).^2 ), ...
                     sqrt( abs(C1P1x - C2P2x).^2 + abs(C1P1y - C2P2y).^2 ), ...
                     sqrt( abs(C1P2x - C2P1x).^2 + abs(C1P2y - C2P1y).^2 ), ...
                     sqrt( abs(C1P2x - C2P2x).^2 + abs(C1P2y - C2P2y).^2 ) ];
            [V, P_cls] = min(dist);

            % find two closest poles and reorient where necessary to have Old_Pole
            % positioned at .mesh(1,_) for both cells
            switch P_cls
              case 1   % P11-P21
                % nothing, because both have new pole correctly oriented

              case 2   % P11-P22
                new_cList{Fr_Brt}{pcd_Anc}.mesh(:,1) = flipud( new_cList{Fr_Brt}{pcd_Anc}.mesh(:,1)) ; 
                new_cList{Fr_Brt}{pcd_Anc}.mesh(:,2) = flipud( new_cList{Fr_Brt}{pcd_Anc}.mesh(:,2)) ; 
                new_cList{Fr_Brt}{pcd_Anc}.mesh(:,3) = flipud( new_cList{Fr_Brt}{pcd_Anc}.mesh(:,3)) ; 
                new_cList{Fr_Brt}{pcd_Anc}.mesh(:,4) = flipud( new_cList{Fr_Brt}{pcd_Anc}.mesh(:,4)) ; 

              case 3   % P12-P21
                new_cList{Fr_Brt}{pcd}.mesh(:,1) = flipud( new_cList{Fr_Brt}{pcd}.mesh(:,1)) ; 
                new_cList{Fr_Brt}{pcd}.mesh(:,2) = flipud( new_cList{Fr_Brt}{pcd}.mesh(:,2)) ; 
                new_cList{Fr_Brt}{pcd}.mesh(:,3) = flipud( new_cList{Fr_Brt}{pcd}.mesh(:,3)) ; 
                new_cList{Fr_Brt}{pcd}.mesh(:,4) = flipud( new_cList{Fr_Brt}{pcd}.mesh(:,4)) ; 

              case 4   % P12-P22
                new_cList{Fr_Brt}{pcd}.mesh(:,1) = flipud( new_cList{Fr_Brt}{pcd}.mesh(:,1)) ; 
                new_cList{Fr_Brt}{pcd}.mesh(:,2) = flipud( new_cList{Fr_Brt}{pcd}.mesh(:,2)) ; 
                new_cList{Fr_Brt}{pcd}.mesh(:,3) = flipud( new_cList{Fr_Brt}{pcd}.mesh(:,3)) ; 
                new_cList{Fr_Brt}{pcd}.mesh(:,4) = flipud( new_cList{Fr_Brt}{pcd}.mesh(:,4)) ; 

                new_cList{Fr_Brt}{pcd_Anc}.mesh(:,1) = flipud( new_cList{Fr_Brt}{pcd_Anc}.mesh(:,1)) ; 
                new_cList{Fr_Brt}{pcd_Anc}.mesh(:,2) = flipud( new_cList{Fr_Brt}{pcd_Anc}.mesh(:,2)) ; 
                new_cList{Fr_Brt}{pcd_Anc}.mesh(:,3) = flipud( new_cList{Fr_Brt}{pcd_Anc}.mesh(:,3)) ; 
                new_cList{Fr_Brt}{pcd_Anc}.mesh(:,4) = flipud( new_cList{Fr_Brt}{pcd_Anc}.mesh(:,4)) ; 
            end
            
            %%%--->>> insert here CHECK SUBPLOT 1.2 <<<---%%%
            
        end %if tid_Anc ~=0   
    end % if t_range notempty     
end % cc
                


%% --- STEP 5 --- TRACKING POLE over cell life -----------------------------------------------------
% Correctly assign the old and new poles throughout each cell-track time
% range. Based on the assignement at division event, the algorithm try to 
% keep the correct polarity untill next division event (or end of track)
% when a new assignment may be required
%
% At each time point ff a cell is compared with the previous one ff-1.
% There are 2 pairs of pole to compare and, as in STEP 4, we calculate
% distances between all possible combinations. Then cellList.mesh is
% rearrenged if necessary necessary in frame ff.
% ---> OLD POLE coordinates are at first position in cellList.mesh(1,_)
%
% N.B.: This method should be no problem to correctly identifying as long
% cells do not move too much and/or time lapses are frequent.

for cc = 1 : size(cTrack, 2)                    % go through each cell-track  
  % time frame where cc-th cell-track was actually tracked (avoid testing at (0;0))
  t_range = find(cTrack{2,cc}(:,1)~=0)';
  
  if ~isempty(t_range)                          % if cell-track is empty, ignores it 
      
    % Find positions of cc-th cell-track  and its ancestor in new_cList
    pcd = pcd_trk{2,cc}(t_range(1));                                        % find the position of cc-th cell-track in new_cList
    tid_Anc = new_cList{1,t_range(1)}{pcd}.track_ancestor ;                 % find the track ID of Ancestor      
    Fr_Brt = new_cList{1,t_range(1)}{pcd}.fr_birth ;                        % N.B.: we put +1 in next for-loop f_ii

    if tid_Anc ~= 0                                                         % if it is not a progenitor
        for fr_ii = Fr_Brt+1 : t_range(end)
          new_cList = Lng_track_polarity(new_cList, fr_ii, cc , pcd_trk);   
        end
    elseif tid_Anc == 0                                                     % if it is Founder cell
        Fr_OffDivs = find(cTrack{6, cc}(:)~=0) ;                            % frames where offspring separate from cell-track
        for gg_d = 1 : length(Fr_OffDivs)                                   % for all division events
          if gg_d < length(Fr_OffDivs)
              for fr_ii = Fr_OffDivs(gg_d) : Fr_OffDivs(gg_d+1)             % from division to division
                  new_cList = Lng_track_polarity(new_cList, fr_ii, cc , pcd_trk);  
              end

          elseif gg_d == length(Fr_OffDivs)
              for fr_ii = Fr_OffDivs(gg_d) : t_range(end)                   % from division to end of track
                  new_cList = Lng_track_polarity(new_cList, fr_ii, cc , pcd_trk);                     
              end 
          end
        end
    end % if Anc~=0
    
  end % if t_range notempty     
end % cc

%%%--->>> insert here CHECK PLOT 2 <<<---%%%


%% --- STEP 6 --- UPDATING and final corrections -------------------------------------
% ----->>> CATCH_1
% Above in STEP 3 for-ff-loop we went through t_range(1) : t_range(end).
% At t_range(1) two cell-track share same cell-outline and it is not the
% true "beginning" for the new track: t_range(2). Now we need to correct
% all daughter (ignore founder cells):
% - remove the first tracked point in cTrack{2,cc}
% - delete new_cList{t_range(1)}{cc} for every cc.

for cc = 1 : size(cTrack, 2)        % go through each cell-track 
    t_range = find(cTrack{2,cc}(:,1)~=0)';
    pcd = pcd_trk{2,cc}(t_range(1));                       % find the position of cc-th cell-track in new_cList
    
    if  new_cList{t_range(1)}{pcd}.track_ancestor ~= 0     % not a founder cell
        new_cList{t_range(1)}(pcd) = [];                   % remove cell detection at first tracked point
        cTrack{2,cc}(t_range(1),:) = [0,0];                % remove first tracked point
    end
end

% cTrack{4,...} so far carried the cell outline number that matched the now
% "old" cell_List at each time point. pcd_trk had the outline number for
% new_cList, but the last update step scrambled them one again. 
% Here, we update cTrack{4,...} to indicate the cellTracks ID number in the
% final new_cList at the every time where a cell was tracked.

for cc = 1 : size(cTrack, 2)                            % go through each cell-track      
  cTrack{4, cc} = zeros( size(cTrack{2,cc},1) ,1) ;     % reinitialize
  tid = cTrack{1,cc} ;                                  % Current track ID of interest

  % time frame where cc-th cell-track was actually tracked (avoid testing at (0;0))
  t_range = find(cTrack{2,cc}(:,1)~=0)';
  
  if ~isempty(t_range)                                  % is cell-track is empty, ignore it 
    for ff = t_range(1) : t_range(end)
        for tt = 1 : size(new_cList{ff},2)
            
            if tid == new_cList{ff}{tt}.ID_ManualTrack
                cTrack{4,cc}(ff) = tt ;
            end
            
        end % tt
    end % ff
  end
end % cc


%% --- STEP 7 --- SAVE Combined_Det2Track.mat file --------------------------------------------------------
% Saving in final variables that end up in .mat file
cellTrack = cTrack ;
cellList.meshData = new_cList(1:length(cTrack{4,1})) ;
cellListN = cellfun(@length, new_cList);
% Update:
% - cellList.detId,  which contains the cell detection ID in cellList for
%                    every timepoint
% - cellList.cellId, which contains the cell ID_ManualTrackID of every 
%                    detected cell in cellList at every timepoint
cellList.detId  = cellList.cellId(1:length(cTrack{4,1})) ;
cellList.cellId = cellList.cellId(1:length(cTrack{4,1})) ;
all_cc = cell2mat(cTrack(4,:));             % create temporary matrix
for ff = 1 : size(all_cc,1)
    cellList.detId{ff} = all_cc( ff, find(all_cc(ff,:)~=0));
    idx = [];
    for cc = cellList.detId{ff}
        idx = [ idx , new_cList{ff}{cc}.ID_ManualTrack] ;
    end
    cellList.cellId{ff} = idx ;
end


if isempty(APP_opt.t5_exp_name)     % if not given, use standard file name
    save( [ APP_opt.t5_path_Track ,'/', 'Combined_Det2Track.mat'] , ...
        'cellTrack', 'cellList', 'cellListN', 'coefPCA', 'mCell', 'p', 'paramString', 'rawPhaseFolder', 'shiftfluo', 'shiftframes', 'weights');
else
    save( [ APP_opt.t5_path_Track ,'/' , APP_opt.t5_exp_name , '_Det2Track.mat'], ...
        'cellTrack', 'cellList', 'cellListN', 'coefPCA', 'mCell', 'p', 'paramString', 'rawPhaseFolder', 'shiftfluo', 'shiftframes', 'weights');
end

      
end % MAIN fnc







%% ----------- CHECK PLOTs ------------------------------------------------
% insert at indicated points of the code to test and debug

%       %%%--->>> CHECK SUBPLOT 1.1 - Check polarity at cell division
%       clf;
%       subplot(1,2,1)
%       xs =  [new_cList{ff}{cc}.mesh(:,1) ; flipud( new_cList{ff}{cc}.mesh(:,3)) ];
%       ys =  [new_cList{ff}{cc}.mesh(:,2) ; flipud( new_cList{ff}{cc}.mesh(:,4)) ];
%       fill(xs,ys,[0 1/cc 1]);    alpha(0.5);    hold on;    axis equal;
%       plot(xs(1),ys(1),'ro', 'MarkerFaceColor', [0.6 0 0]);  
%       xs =  [new_cList{ff}{Anc}.mesh(:,1) ; flipud( new_cList{ff}{Anc}.mesh(:,3)) ];
%       ys =  [new_cList{ff}{Anc}.mesh(:,2) ; flipud( new_cList{ff}{Anc}.mesh(:,4)) ];


%      %%%--->>> CHECK SUBPLOT 1.2 - Check polarity is corrected if needed
%      subplot(1,2,2)
%      xs =  [new_cList{ff}{cc}.mesh(:,1) ; flipud( new_cList{ff}{cc}.mesh(:,3)) ];
%      ys =  [new_cList{ff}{cc}.mesh(:,2) ; flipud( new_cList{ff}{cc}.mesh(:,4)) ];
%      fill(xs,ys,[0 0 1/cc]);    alpha(0.5);    hold on;    axis equal;
%      plot(xs(1),ys(1),'ro', 'MarkerFaceColor', [0.6 0 0]);  
%      xs =  [new_cList{ff}{Anc}.mesh(:,1) ; flipud( new_cList{ff}{Anc}.mesh(:,3)) ];
%      ys =  [new_cList{ff}{Anc}.mesh(:,2) ; flipud( new_cList{ff}{Anc}.mesh(:,4)) ];
%      fill(xs,ys,[0 1/cc 1]);    alpha(0.5);    hold on;    axis equal;
%      plot(xs(1),ys(1),'ro', 'MarkerFaceColor', [0.6 0 0]);         
%      pause(2); 


%--------------------------------------------------------------------------

%   %%%--->>> CHECK PLOT 2 - check cell polarity after division untill end of movie
%   clf
%   for Fr = 15: length(new_cList)
%      for cc = 1:length(new_cList{Fr})
%      if ~isempty( new_cList{1, Fr}{1, cc})
%          xs =  [new_cList{1,Fr}{cc}.mesh(:,1) ; flipud( new_cList{1,Fr}{cc}.mesh(:,3)) ];
%          ys =  [new_cList{1,Fr}{cc}.mesh(:,2) ; flipud( new_cList{1,Fr}{cc}.mesh(:,4)) ];
%          fill(xs,ys,[0 1/cc 1]);    alpha(0.5);    hold on;    axis equal;
%          plot(xs(1),ys(1),'ro', 'MarkerFaceColor', [0.6 0 0]);  
%       end
%       end
%       pause(0.5);
%       clf;
%   end



