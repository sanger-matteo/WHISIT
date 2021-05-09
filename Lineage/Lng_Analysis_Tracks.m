function [cList, LG] = Lng_Analysis_Tracks( cList , cTrack )
% The function creats a a new reformatted database, called clone_List. It 
% add a clone ID number to each cell, according to the cell ancestor and 
% division(s) performed.
%
%
% ---> What is a cell-track?
% A cell-track is a single cell tracking, where the user tried to follow
% over time the same physical cell. Each cell-track can encompass multiple
% division events and at each one the user chose to keep following one of
% the two daughter cells, while creating a new cell-track for the other
% daughter cell.
%
% ---> What is a clone?
% Each cell track is the concatenation of multiple "cell-lifes". We use the
% term clone to define the time segment of a cell-track that:
% - starts with the division that generated the cell (or start of tracking)
% - ends with the next cell division where the cell generate two daughters
%   (or end of tracking)
% 
% ---> Aim of the script?
% The main steps of this function are:
% -  STEP 1 - Identify the clones that constitutes each cell-track and
%             their respective time windows
% -  STEP 2 - Find the founder clone(s) and define the lineage prefix ID
% -  STEP 3 - Determine which pole (old or new) clones inherits when born
% -  STEP 4 - Assigning an unique ID name to the clones
% -  STEP 5 - Split cellList to crate a list of clones, independent of the
%             detection file , track file and stack of images
%
% ---> Why are we doing such complex steps?
% In order to easily reach the one and only possible solution for the
% lineage, we need to split each cell-track in the component "clones". This
% is because each individual clone is a branch of the tree, the minimal
% "component" of the tree. Cell-tracks are related to the lineage tree, but
% done randomly. Thus, splitting in their minimal components, the clones,
% allow to find the unique solution for the lineage tree, irrespective of
% the order and number of tracks done. 
%
%
% ---> UNIQUE ID_clone name
% When a clone life begins, ends, the ancestor and ofspring it generates
% are relatively easy to get. The most important information is to give
% a unique ID name for each clone because this will allow to understand the
% relationship between clones and where to place them in the lineage tree:
% The ID_clone is in the form  XXX.ZZZZZ....Zn
% 
% prefix :  XXX.          indicates the lineage the clone belongs to. This
%                         usualy corresponds to the ID number of the ancestor  
%                         (usually the first tracked cell)
% postfix : .ZZZ....Zn    the postfix, is specific a set of numbers that
%                         are either 1 or 2.
%
% At every divisions event a clone inherit the ID_clone from the ancestor
% and add one digit (Z) at the end of it. The Z can only have one of two
% values:
% --> 1 = cell originated from Old_Pole of the "mother" cell
% --> 2 = cell originated from New_Pole of the "mother" cell
% In this just from the ID_clone we know the history and all the ancestors
% of a specific clone and who is the "twin" daughter clone
%
% i.e., 
% if we have a given clone with ID_clone:               -->  002.12112
% - We know that the immediats ancestor :               -->  002.1211
% - and the clone before that was :                     -->  002.121
% - so on untill we find the founder clone :            -->  002.
% 
% We also know that:                                    -->  002.1211  
% - must have generated two daughter cell,              -->  002.12112 
% and its "twin"                                        -->  002.12111
% 
%
% ----- INPUT -------------------------------------------------------------
% cellList = the main variable that contain all informations about the
%           detection done in Oufti. We assume detection was done with
%           "independent frames" option selected. Therefore, each cell has
%           been ramdomly assigned a number and there is no correlation
%           between different frames.
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
% ----- OUTPUT ------------------------------------------------------------
% clone_List = a cell array. Each element is a cell array to, containing
% the cell detection data, fluorescence analysis... of an individual clone
% for its entire lifespan (birth : next div).
% - cellList structure is organized as {Frame#}{cell#}
% - clone_List is reorganized in       {cell#}{Frame#}
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


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


%% --- Step 0 --- Initialize variable -------------------------------------
% The variable below are widely used through the code ant allow to write
% more concise and efficient code. Here we initialize and explain what
% would otherwise be obscure variable;

% Logical matrix that show the time window for each cell-track. This
% simplify future search and avoid using longer command, such as:
% find( cTrack{4, LG{cc}.Offspring(dd)}(:,1) ~=0 );    
tw_tracks = ( cell2mat(cTrack(4,:)) ~=0 ) ;                         

% Lineage variable LG = based on cTrack, we expand and create a struct
% variable that summarize all information for each specific cell-track,
% division events, clone segments, number ID etc...
LG = {} ;


% --- Indexing (tid) and counter (cc) -------------------------------------
% cc  = counter used to go through each cell-track in for-loops. 
% It is easier to go through cTrack and LG column-by-column using cc as
% index. However, columns order is unrelated to the track ID number, which
% we need to correctly access cList. cList elements are reorganized such
% that their position match the corresponding cell-track ID number.
% Hence, for simplicity, we create a separate variable to carry the track
% ID number and use it as index to reposition and access new_cList: 
% >> tid = cTrack{1,cc};

% % % % % % % % cc2tid = a  two row matrix that carries the cell-tracks ID number and the
% % % % % % % % respective column index position in cTrack, LG etc... We need two index
% % % % % % % % systems in parallel because the element position in cList is according to
% % % % % % % % cell-tracks ID number, but this is not so for cTrack, LG etc...
% % % % % % % % see above --- Indexing (tid) and counter (cc)
% % % % % % % cc2tid = [ 1 : length( cell2mat(cTrack(1,:)) ) ;  cell2mat(cTrack(1,:)) ];     
% % % % % % % % cTrack(1,:) = mat2cell( [1 : length(cell2mat(cTrack(1,:)))] , 1 , ones(1,length(cTrack(1,:))) );

cc_det = cell2mat(cTrack(4,:));        % the cell detection position in cList at each time point fore each cell-track (column)
cc2tid = zeros(size(cc_det));          % repeat the manual track ID number at each time point where cell-track (column) was tracked
for kk = 1 : size(cc_det,2)
    idx = find(cc_det(:,kk) ~= 0) ;
    cc2tid(idx,kk) = cTrack{1,kk};    
end
    



%% --- STEP 1 --- Time Windows Segments ----------------------------------%
% This step divides each cell-track in "segments" which identify the time
% window of individual "clone". A standard "clone" is a unique single cell
% whose "life" time-window goes from division event that generated it to
% division event where it generate two new clones. Each clone therefore has
% a start (division or start of manual tracking) and an end (either divide
% and generate two new clones or end of manual tracking).
% A cell-track can be composed of concatenated segments, each an individual
% clones. Objectives of the function is to identify and divide clones from
% cell-tracks.
   
for cc = 1 : size(tw_tracks,2)
    % Gather some general info for each cell-track
    LG{cc}.ID_cellTrack = cTrack{1,cc} ;                              % cell-Track ID, unique number (given during Manual Tracking)
    LG{cc}.tw_cellTrack = find(cTrack{2,cc}(:,1)~=0) ;                % cell-Track time windows when it was tracked
    LG{cc}.Fr_Birth = LG{cc}.tw_cellTrack(1) ;                        % begin and ...
    LG{cc}.Finish = LG{cc}.tw_cellTrack(end) ;                        % end of maual tracking for cell-track
    LG{cc}.cDetID_Birth = cc_det(LG{cc}.Fr_Birth,cc);                 % cell detection position in cList at time of birth
    
    LG{cc}.Off_Div  = [];                                             % when offspring was generated (division of new cell-track)    
    LG{cc}.Offspring = cTrack{6, cc}(cTrack{6, cc}(:)~=0 )' ;         % descendant/daughter cells that the given cell have generated
    
    if isempty(cTrack{5, cc}(cTrack{5, cc}(:)~=0 )')             	  % ancestor, mother cell that created cell track
        LG{cc}.Ancestor  = 0 ;                                        % it is a founder cell
    else
        LG{cc}.Ancestor  = cTrack{5, cc}(cTrack{5, cc}(:)~=0 )' ;   
    end    
    
    % Division Frame for Offspring = find frame number when offspring is
    % generate, simply time_window_offspring(1)
    for dd = 1:length(LG{cc}.Offspring)
        % t_cln = find the column where cell-track .Offspring(dd) is stored
        t_cln = find(cell2mat(cTrack(1,:)) == LG{cc}.Offspring(dd)) ;
        tw_offspring = find(tw_tracks(:, t_cln) == 1) ;
        LG{cc}.Off_Div(dd) = tw_offspring(1);
    end
    
    % --- Clones Time Windows Segments
    % tw_Seg_clones = stores two paired numbers: as string in form ---> start#:last#
    
    % ss = counter for clone segment. First and last segments are a bit "special" we need to use while-loop
    ss = 1;                                         
    if  isempty(LG{cc}.Off_Div)        
        % cell has no offspring then entire track is a single clone
        LG{cc}.Clones_Seg{ss} = [num2str(LG{cc}.Fr_Birth) ':' 'end'];
        
    elseif ~isempty(LG{cc}.Off_Div)     
        % Clone segment end at .Off_Div(1)-1 because the .Off_Div(1) is 
        % frame of birth of next generation, therefore, it must end one frame earlier
        LG{cc}.Clones_Seg{ss} = [num2str(LG{cc}.Fr_Birth) ':' num2str(LG{cc}.Off_Div(1)-1)];
        ss = ss +1;
        while ss <= length(LG{cc}.Off_Div)      % for all next clone's offspring     
            LG{cc}.Clones_Seg{ss} = [num2str(LG{cc}.Off_Div(ss-1)) ':' num2str(LG{cc}.Off_Div(ss)-1)];
            ss = ss +1;
        end
        % last clone segment is until the manual tracking is finished: 'end'
        % (i.e. cell was not anymore detected, swim/reach end movie...)
        LG{cc}.Clones_Seg{ss} = [num2str(LG{cc}.Off_Div(ss-1)) ':' 'end'];
    end
end % for cc




%% --- STEP 2 --- Founder Prefix -------------------------------------------------
% We find the founder clone(s) and define the lineage prefix the entire 
% time window(s). The founder clone define the Prefix of the lineage ID as
% 'XXX.' (three digits). Because we do not know polarity it received from
% ancestor, the ID has nothing following the '.' dot

Pfx_Lineage = 0;
for cc = 1:size(tw_tracks,2)                   % go thourgh each cell_line
    Pfx_Lineage = Pfx_Lineage+1 ;
    if LG{cc}.Ancestor == 0                    % if it is a founder clone  
       if isempty(LG{cc}.Offspring)            % in case founder has no offspring       
           sg_end = LG{cc}.Finish ;            % we go untill end of track

       elseif ~isempty(LG{cc}.Offspring)       % else ...        
           sg_end = LG{cc}.Off_Div(1) -1 ;     % end at .Off_Div(1)-1 because the .Off_Div(1) ...
       end                                     % is frame of birth of next generation

       if     Pfx_Lineage >= 100;    null = '';
       elseif Pfx_Lineage >= 10;     null = '0';
       elseif Pfx_Lineage < 10;      null = '00';
       end
       for ff = LG{cc}.Fr_Birth : sg_end          % set prefix for the entire clone's time window 
           ndet_cc = cc_det(ff ,cc) ;           % position of cc in cList                           
           cList{ff}{ndet_cc}.ID_clone = [null num2str(Pfx_Lineage) '.'];
      end            
    end
end % cc




%% --- STEP 3 --- INHERITED POLARITY at DIVISION (ID_div) ------------------
% As mentione above, the ID name give to cell will depend on which pole
% (old or new) a cell inherits when is born.
% This step takes every division event and determine which pole (old or
% new) the two daughter clones inherits when born. A matrix (ID_div) is
% created that store the information (value Z) for every division event of
% every celltrack. This will be the base from which we can univocally place
% an ID for each clone in for STEP3.
% ID_div is a matrix that has as many rows as cell_lines and as many
% columns as max number of divisions events in the lineage. Each row stores
% the pole inherited at each division (clone segment) of a cell-track. The
% value (Z) can be:
% --> 1 = cell originated from Old_Pole of the "mother" cell 
% --> 2 = cell originated from New_Pole of the "mother" cell
% --> 0 = no division event
%
% N.B.: In a progenitor with 5 segments, only 4 are considered for ID_div
% naming (.ZZZ...) because the first is the progenitor clone, with name XX.
% In non-progenitor cells they could have 4 descentants, but this means 5
% segments and all produce a clone whose naming must be (.ZZZ...).
%     we need to use ---> maxx = length(LG{cc}.Clones_Seg) 
%     rather than    ---> maxx = length(LG{cc}.Offspring)

maxd = -1 ; 
for cc = 1:size(tw_tracks,2)
    if length(LG{cc}.Offspring) > maxd    
        maxd = length(LG{cc}.Clones_Seg);     
    end
end

% initialize with zeros
ID_div = zeros(size(tw_tracks,2), maxd);         

for cc = 1:size(tw_tracks,2)               % go thourgh each cell_line    
  if LG{cc}.Ancestor ~= 0                  % skip progenitor cells (polarity is unknown in founder clone)
    
    ff = LG{cc}.Fr_Birth ;   
    ndet_cc = cc_det(ff ,cc) ;             % position of cc in cList  
    % column where ancestor is positioned (ndet_Anc carries cell-detection ID number)
    col_Anc = find( cell2mat(cTrack(1,:)) == LG{cc}.Ancestor) ;
    % ndet_Anc carries cell-detection position of the ancestor.
    ndet_Anc = cc_det(ff ,col_Anc) ;       % position of Ancestor in cList      

    % We have 3 cell-outlines and evaluate the old pole in all of them: 
    % - two new clones (daughter) cells at division frame, 
    % - the "united" ancestor at the frame before cell division 
    % Indexes: New clone 1 continues the "cell-track" (we use ndet_Anc),
    %          while the new clone 2 is a new track   (we use cc).
    OP1x  = cList{1,ff  }{1, ndet_Anc}.mesh(1,1);      OP1y  = cList{1,ff  }{1, ndet_Anc}.mesh(1,2);
    OP2x  = cList{1,ff  }{1, ndet_cc}.mesh(end,1);     OP2y  = cList{1,ff  }{1, ndet_cc}.mesh(end,2);
    % Old pole in "united" ancestor clone
    OPA_x = cList{1,ff-1}{1, ndet_Anc}.mesh(1,1);      OPA_y = cList{1,ff-1}{1, ndet_Anc}.mesh(1,2); 
    
    % distances between old_poles in ancestor-daughter pairs:   
    d1 = sqrt( abs(OP1x - OPA_x).^2 + abs(OP1y - OPA_y).^2 );
    d2 = sqrt( abs(OP2x - OPA_x).^2 + abs(OP2y - OPA_y).^2 );

    % emp : append the new division polarity in matrix ID_div, finding the 
    %       first empty (zero) element in the XY-th row.              
    % Update polarity of both clone 1 and 2.
    if d1 < d2                                                % New clone 1 (ndet_Anc) is nearer to Old_pole
        empt = min(find( ID_div(cc, 1:end)==0));              % find first empty place in row cc	
        ID_div(cc, empt) = 2;
        empt = min(find( ID_div(col_Anc, 1:end)==0));         % find first empty place in row ndet_Anc	
        ID_div(col_Anc , empt) = 1;
    elseif d1 > d2                                            % New clone 2 (cc) is nearer to Old_pole
        empt = min(find(ID_div(cc, 1:end)==0));             
        ID_div(cc, empt) = 1;
        empt = min(find( ID_div(col_Anc, 1:end)==0));     
        ID_div(col_Anc, empt) = 2;
    end
  end
end % for cc


DIV = {};
for cc = 1:size(tw_tracks,2)
    DIV{cc}(1,:) = LG{cc}.Clones_Seg ;
    lng = size(DIV{cc}(1,:),2) ;
    DIV{cc}(2,:) = mat2cell( [LG{cc}.Fr_Birth, LG{cc}.Off_Div] , 1, ones(1,lng) );
    if LG{cc}.Ancestor == 0        
        DIV{cc}(3,:) = mat2cell( [0, ID_div(cc, find(ID_div(cc,:)~=0))] , 1, ones(1,lng) );
    else
        DIV{cc}(3,:) = mat2cell( ID_div(cc, find(ID_div(cc,:)~=0)), 1, ones(1,lng) );
    end
        
end




%% --- STEP 4 --- ASSIGN c_ID (ID_clone) ----------------------------------------------
% Using ID_div this step is able to univocally determine the identity of
% each clone. The name is base on the ancestors and pole inherited.
% To create ID_clone we need to univocally identify:
% anc_cloneID : XX.ZZZ_     the ID of the ancestor is inherited from the ancestor
% pole_inher : Z            the pole inherited by the clone at the specific
%                           cell division determine the last digit of ...
% new_cloneID               = [anc_cloneID  pole_inher]      

for cc = 1 : size(tw_tracks,2)                   % go thourgh each cell-track
       
    if LG{cc}.Ancestor == 0                      % Founder clone do not have a postfix (unknown polarity)
        seg_start = 2;     seg_diff = 1;         % ... start from second segment, 
    elseif LG{cc}.Ancestor ~= 0    
        seg_start = 1;     seg_diff = 0;         % start from first segment
    end
    
    for seg = seg_start : length(LG{cc}.Clones_Seg)            % go through each segment of a cell_line   
        b = 0;
        lid = '' ;
        Vseg = seg;
        Vcc = cc;
        
        % --- Establish the cloneID -----------------------------------
        while b ~= 1 
            int_lid = int2str(cell2mat(DIV{Vcc}(3,1:Vseg)));
            lid = [ int_lid(~isspace(int_lid)), lid];
            if lid(1) == '0'
                b = 1;
                ndet_cc = cc_det(LG{Vcc}.Fr_Birth ,Vcc);
                anc_cloneID = cList{LG{Vcc}.Fr_Birth}{ndet_cc}.ID_clone ;
            else
                tid = LG{Vcc}.Ancestor;
                frdiv = LG{Vcc}.Fr_Birth;
                Vcc = find( cell2mat(cTrack(1,:)) == tid ) ;
                Vseg = find(cell2mat(DIV{Vcc}(2,:)) == frdiv) -1;   
                % [-1, we want the to go the back to "ancestor" segment, before division]
            end
        end
        lid = lid(2:end);
        new_cloneID = [anc_cloneID  lid];         

        % --- Establish start:last  save in all time points of the clone (segment) ------
        [fr] = strsplit(LG{cc}.Clones_Seg{seg},':');
        fr_str = str2num(fr{1}) ;    
        if strcmp(fr{end}, 'end');    fr_end = LG{cc}.Finish;
        else ;                        fr_end = str2num(fr{end}) ;
        end

        for ff = fr_str : fr_end
            ndet_cc = cc_det(ff ,cc) ;                % position of cc in cList 
            cList{ff}{ndet_cc}.ID_clone = new_cloneID;
        end
        
    end % for seg
    
end % for cc



end







