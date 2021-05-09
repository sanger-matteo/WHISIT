function Report_Track2Det
%
%Report_Track2Det - The function take the track.txt file and Create a .txt
%                   report file that contains information related to cell
%                   divisiona and daughter-mother connection for each
%                   specific cell-track.
%                   
% Errors such as having no mother cell, multiple one, sharing cell-outline
% with cells beyond the first "birth" frame, are all written down in the
% report.
% This output file can help the user to quickly the error in the tracking 
% and adress them accordingly
%
% The Error-free output for a specific cell-track consist of the following
% infos in the format shown below:
%
% ---------------------------------------------------
% CELL TRACK -----------> CC ------------------------
% 
% Mother cell:     X 
% Lifespan   :    XY  frames
% From frame :   XYZ 
% To frame   :   XYZ 
% Cell-tracks generated ->     XX    YY    ZZ 
% At frame -------------->     WW    KK    JJ 
%
% WARNING - ( here are indicated conflicting frames in which cell CC
%             share same cell-outline with more than one cell )
%
% ---------------------------------------------------
%
% Step (1) Cell-outline ID and (2) Identify DIVision events are the same as
% in function  Combine_Det2Track.m. This is necessary to identify cell
% track univocally with the detection done in Oufti. Using cell outlines
% rather than xy points we can find easily common errors:
%  - establish the relationship between mother and daughter cells, 
%  - when a cell outilne is used for more than one track
%  - who is a lineage founder cell
%
% INPUT data:
%   CellTracks = a global variable where that stores the information of
%   every cell track. 
%
%
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-


%% --- 1 --- Cell-outline ID -----------------------------------------------------
% First part of function's algorithm identify the ID-number of the cell-outline
% for each cell-track at each time point

% WE RECEIVE CellTracks as a .txt file
global APP_opt;


%--- Load DET.mat file ---------------------------------------------------%
% Load the DETection.mat file from Oufti. We this only to find the cell
% outline ID of tracks, in order to distinguuish cells according to the
% detection done in Oufti
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
    
    
% Initialize the Report.txt file, with appropriate filename
if isempty(APP_opt.t5_exp_name)
    filename_Report_txt = [ APP_opt.t5_path_Track ,'/', 'Report_Tracking.txt'];
else
    filename_Report_txt = [ APP_opt.t5_path_Track ,'/' , APP_opt.t5_exp_name , '_', 'Report_Tracking.txt'];
end
file_R = fopen(filename_Report_txt, 'w+');
fprintf( file_R, '--->>> MANUAL TRACKING REPORT <<<---\n\n' );




%% --- 1 --- Cell-outline ID -----------------------------------------------------
% First part of function's algorithm identify the ID-number of the cell-outline
% for each cell-track at each time point


dist = [];          % store the distances of each cell outline from the 
                    % track point of scc-th track in current frame
                    
for cc = 1 : size(cTrack, 2)        % go through each cell-track    
  % Store detection outline cell_ID that correspond to the manual tracking
  cTrack{4, cc} = zeros( size(cTrack{2,cc},1) ,1) ;  

  % We can avoid testing where there is not tracking done (0;0)
  t_range = find(cTrack{2,cc}(:,1)~=0)';
  if ~isempty(t_range)      % is cell-track is empty, ignore it 
      
    for ff = t_range(1) :1: t_range(end)
    % Store XY_t point for all cell-tracks in ff-th frame
    ff_X = cTrack{2, cc}(ff , 1) ;
    ff_Y = cTrack{2, cc}(ff , 2) ;
    % remove '0' coordinates - not-tracked points will not be considered
    ff_X(ff_X==0) = NaN;     
    ff_Y(ff_Y==0) = NaN;        
    dist = [];     % store the distances of each cell outline from the 
                   % track point of scc-th track in current frame

    % Go through all cell's outline detected in ff-th frame. Calculate
    % distances to find closest outline to the tracked point of cc-th cell-track
    for dd = 1 : length( cellList.meshData{ff} )
        if ~isempty( cellList.meshData{ff}{dd}.mesh  )
            % Store XY_o points for kk_th outline and plot it
            X_o = [cellList.meshData{ff}{dd}.mesh(:,1) ; flipud(  cellList.meshData{ff}{dd}.mesh(:,3)) ];
            Y_o = [cellList.meshData{ff}{dd}.mesh(:,2) ; flipud(  cellList.meshData{ff}{dd}.mesh(:,4)) ];
            % Calculate R2 between each XY_t and all XY_o points of kk-th cell 
            % ontour outline and store the minumum value
            dist(dd) = min(double( sqrt( abs(ff_X - X_o).^2 + abs(ff_Y - Y_o).^2 ) ));

            if inpolygon(ff_X,ff_Y, X_o,Y_o)       % and lot cell outline in selected_clr if XY_t is inside it             
                isinc = 1 ;
                cTrack{4, cc}(ff) = dd ;
            else
                isinc = 0 ; 
            end
        end
    end %/for dd  

    % If XY_t was not inside any cell outline, find the closest cell outline 
    idx_closest =  find(dist == min(dist) );     
    if isinc == 0  &&  ~isempty(idx_closest)  &&  cTrack{4, cc}(ff) == 0     % if not inside cell ...      
        cTrack{4, cc}(ff) = idx_closest ;
    end
    end % ff
  end %if t_range notempty
end % cc


%% --- 2 --- Identify DIVision events -----------------------------------------------------
% Second part of function's algorithm identify the frame of birth of each cell-track,
% who is the mother/ancestor cell and when those events occured

M_coID = cell2mat(cTrack(4,:));    % matrix storing all cell-outline ID of all cell-tracks
err = 0;
N_Founders = 0;                     % number of lineager founder cell

for cc = 1 : size(M_coID, 2)
    
    % Store cell-track ID-number of  ancestor/mother cell
    cTrack{5, cc} = zeros( size(cTrack{2,cc},1) ,1) ;
    % Store cell-track ID-number of daughter cell(s) 
    cTrack{6, cc} = zeros( size(cTrack{2,cc},1) ,1) ;
    
    % find time window and first number is frame of birth
    t_range = find( M_coID(:,cc) ~= 0 );         % Store time frame where cc-th cell was actually Tracked
    if ~isempty(t_range)      % is cell-track is empty, ignore it 
    
    Fr_Brt = t_range(2);                         % frame of birth of cell-track
    
    % Find who share the same cell-outline with cc-th cell. This is the mother cell
    co_share = find( M_coID(t_range(1),:) == M_coID(t_range(1),cc));
    if length(co_share) > 2                       % ERROR - there can be only 2 cells-tracks sharing same outline
        err = err +1;
        fprintf(file_R , 'ERROR - cell-track %i share cell outline at birth with too many cell-tracks: \n', cc );
        fprintf(file_R , '-------%s   at frame %3i\n\n',  strjoin(pad(split(num2str(co_share(co_share ~= cc))), 4,'left')) , t_range(1) );
        
    elseif length(co_share) == 1                  % it is a lineage founder cell
         N_Founders = N_Founders +1; 
         
    elseif length(co_share) == 2 
        co_share(co_share == cc) = [];            % we exclude current cc-th cell-track we examine
        cTrack{5, cc}(Fr_Brt)  = co_share;
        cTrack{6, co_share}(Fr_Brt)  = cc;
    end
    
    end %if t_range notempty
end % cc

if N_Founders == 1
    fprintf(file_R , 'There is only %i Founder cell-track %s \n', N_Founders );
elseif N_Founders >= 2
    fprintf(file_R , 'WARNING: the number of Founder cell-tracks is %3i \n', N_Founders );
end
fprintf(file_R , ' \n \n ') ;



%% --- 3 --- Create REPORT_Tracks.txt -----------------------------------------------------
% Here we create the major part of the Report file

M_Mot = cell2mat(cTrack(5,:));    % matrix contining cell-track ID-number of  ancestor/mother cell
M_Dau = cell2mat(cTrack(6,:));    % matrix contining cell-track ID-number of daughter cell(s) generated 
Warns = 0 ;                           % Count number of warning

for cc = 1 : size(M_coID, 2)          % go through each cell-track 
    
    % Heading for cc-th cell-track 
    fprintf(file_R , '--------------------------------------------------\n' );
    fprintf(file_R , 'CELL TRACK -----------> %i ------------------------\n', cc );
    
    Warns = 0 ;                          % Count number of warning for cc-th cell-track
    t_range = find(M_coID(:,cc)~=0);     % find time window for cell-track "lifespan"
    
    if ~isempty(t_range)
        % Info on mother-track from which cc-th was generated
        birth = find(M_Mot(:,cc)~=0);
        if ~isempty(birth)
            Mot = M_Mot(birth, cc); 
            fprintf(file_R , 'Mother cell: %4i \n', Mot);    
        elseif isempty(birth)
            fprintf(file_R , 'Mother cell: NONE -----> FOUNDER cell-track \n');
        end
        
        % Info on Lifespan of cell-track
        fprintf(file_R , 'Lifespan   : %4i  frames\n', t_range(end)-t_range(1) );  
        fprintf(file_R , 'From frame : %4i \n', t_range(1) );  
        fprintf(file_R , 'To frame   : %4i \n', t_range(end) );
    
        % Info on daughter-tracks generated by cc-th cell-track
        gen_Dau = find(M_Dau(:,cc)~=0);
        if ~isempty(gen_Dau)
            Dau = M_Dau(gen_Dau, cc) ;    
            fprintf(file_R , 'Cell-tracks generated -> %s \n',  strjoin(pad(split(num2str(Dau'))', 5,'left')) );
            fprintf(file_R , 'At frames -------------> %s \n',  strjoin(pad(split(num2str(gen_Dau'))', 5,'left')) );
        elseif isempty(gen_Dau)
            fprintf(file_R , 'Cell-tracks generated -> NONE \n' );
        end

        fprintf(file_R , '\n') ;
        
        % Check if multiple cell-track share same cell-outline ID behiond the
        % division frame. Although it should not generate errors when creating
        % a lineage tree, this could unintentionally alter the results.
        for ff = t_range(2) : t_range(end)
            co_share = find(M_coID(ff,:) == M_coID(ff,cc));            
            if length(co_share) >= 2
                if Warns == 0 
                    fprintf(file_R , 'WARNING: multiple cell-tracks sharing same cell-outline ID: \n' );
                    fprintf(file_R , '!!! at frame %4i  between tracks %s \n', ff, strjoin(pad(split(num2str(co_share)), 4,'left')) );
                else
                    fprintf(file_R , '!!! at frame %4i  between tracks %s \n', ff, strjoin(pad(split(num2str(co_share)), 4,'left')) );
                end
                 Warns = Warns +1;
            end
        end % ff
    
    else isempty(t_range)     % if cell-track is empty, it was not tracked at any frame
        fprintf(file_R , ' !!! ERROR: Cell-track %3i  has Lifespan of 0 frames !!! \n', cc );  
    end % if t_range isempty

    fprintf(file_R , ' \n \n' );
    
end % cc
    
    




