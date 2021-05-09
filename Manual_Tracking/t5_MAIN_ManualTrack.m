function t5_MAIN_ManualTrack(app)
%
%ManualTrack_MAIN - it is the Main Function of the tab "Manual Tracking"
%   that run all major events during the manual tracking
%
%   The function create a figure for the manual tracking that display, one
%   at a time, the frames of the selected stack of images. The user will 
%   be able to click using the mouse inside the displayed frame in order to 
%   track the position of a single microbial cell. 
%   The x and y pixel coordinates of the point where the cursor was during
%   clicking are stored. The function will ensure the point is stored in
%   the correct cell-track and at each click the frames will automatically
%   update and pregressively move forward, going one frame at a time 
%   through the entire movie 
%
%   ManualTrack_MAIN will pause at "grabinput" function and wait as long as
%   (1) there is a click inside the figure or (2) STOP button is pushed. 
%   At any moment the user can move between frames or change visualization
%   options, etc.... This does not affect the working of ManualTrack_MAIN
%   and since the most important variables are global they will all be
%   up-to-date when called.
%
%   --- MAJOR VARIABLE ----------------------------------------------------
%   Several variables are valid across all functions used during the
%   Manual_Tracking and will be discussed extensively only here below.
%   
%   CellTracks = a cell-track is the manual tracking for of a single 
%       microbial cell. CellTracks is a cell-array, where each column
%       correspond to a single specific cell-track. It stores a list of
%       coordinates of the position of the cell at each time point tracked.
%     Each column is organized as follow:
%       Row 1 - ID-number of the cell-track
%       Row 2 - point coordinates, stored as [x,y] matrix of length LenStack
%       Row 3 - the specific [R G B] color asigned to the cell-track
%
%     When a new CellTracks is initialize
%       Row 1 - the first element is ID-number "1"
%       Row 2 - (xy points) with as many zeros as total number of time points
%           (N.B.: there is no position (0;0), the "first" must be (1;1))
%       Row 3 - empty, no color assigned. A specific color will be determined 
%           for the cell-track the first time it is plotted
%
%   scc = store the ID-number of the currently selected cell-track that the
%       user is working on. When saving the xy points of clicking, it will be
%       only for this cell-track
%
%   oDet = will store the detection.mat file from Oufti, if provided. 
%       It will be usefull to visualize cell's outlines, in order to aid 
%      during trackinshow
%
%   APP_opt = Store most information for the GUI, choises and options, as
%       well as global variables that are necessary for the working of 
%       multiple functions.
%       For example APP_opt.t5_ff stores the number of the currently 
%       displayed frame
%
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-



%% ******* INITIALIZATION ******************************************************
% Load images, data and files necessary to perform manual tracking

global APP_opt ;        global oDet ;      
global CellTracks;      global scc ;


% Update the total stack number in the GUI
app.t5_Label_totFrames.Text = [ 'of ' num2str(size( APP_opt.t5_srcFiles_BF, 1))];

% Number the stack of images in the movie provided
LenStack = size(APP_opt.t5_srcFiles_BF,1);


%--- Find number of DIGITS in BF filename ------------------------------------------------% 
% Find .tif files inside BF stack and calcolate the total number of digits 
% at ending number (filename_xxxx.tif)

APP_opt.t5_T_dig = TotDigits_in_Filename(APP_opt.t5_srcFiles_BF, APP_opt.name_delimiters);

if APP_opt.t5_T_dig == -1       % then we have not found any tif correctly named
    app.TextOUT.Value = sprintf('\n%s\n%s',  '!!! FL stack contain uncorrect .tif filename !!!' ,...
                                'Format should be in form:   filename_xxxx.tif' );
    APP_opt.ERROR = 1;          
    return;                     % from the ManualTrack_MAIN
end


%--- Load DET.mat file -------------------------------------------------------------------%
% If provided load the DETection.mat file from Oufti
if APP_opt.t5_display_Det == 1
    app.TextOUT.Value = sprintf('\n%s',  'Loading Det.mat file' );
    oDet = load([APP_opt.t5_path_Det , APP_opt.t5_fileName_Det ]);
    app.TextOUT.Value = sprintf('\n%s',  'Processing ...' );
end  


%--- Load TRACK.txt file -----------------------------------------------------------------%
% Load a list of cell-tracks, which are stored as tab-separated .txt file
if app.t5_choose_Track.Value == 1 
    % Create a clear CellTracks
    CellTracks(:) = [] ;
    temp_T = load([APP_opt.t5_path_Track , APP_opt.t5_fileName_Track ]);    
    t_IDs = temp_T(1,2:2:end);          % row array with ID-numbers
    t_trk = temp_T(2:end,:) ;           % matrix with tracked points
    [Rtrk, Ctrk] = size(t_trk);         % extract row and column numbers    
    % convet the matrixes in cell-array and place them correctly in
    % columns each representing a single cell-track
    CellTracks(1,:) = mat2cell( t_IDs , 1    , ones( 1,length(t_IDs)) );
    CellTracks(2,:) = mat2cell( t_trk , Rtrk , 2.*ones(1, Ctrk/2) );
    CellTracks(3,:) = cell(1,length(t_IDs));
    clearvars temp_T t_IDs t_trk Rtrk Ctrk;
end


%--- Initialize CellTracks ---------------------------------------------------------------%
% If empty, initialize CellTracks with a single clear cell-track 
% If not empty, choose first cell in the array as the currently selected cell
% NOTE: if tracking is stopped, CellTracks variable is not deleted.
%       Therefore, it is not empty at restart and the tracking will resume
%       with the previously stored CellTracks (unless a .txt file is loaded)
%       In a way, the program "remember" (unsaved) data from previous tracking

if isempty(CellTracks)
    scc = 1 ;
    CellTracks{ 1, scc } = scc ;                    % set ID-number of track
    CellTracks{ 2, scc } = zeros( LenStack ,2) ;    % initialize coordinates
    CellTracks{ 3, scc } = [] ;                     % [R G G] color
    
    % Update GUI "Select Cell" to point to "1" cell
    app.t5_SelectCellDropDown.Items = {'1'};
    app.t5_SelectCellDropDown.Value = {mat2str(scc)}; 

elseif ~isempty(CellTracks)
    % create a list of all ID-numbers stored in CellTracks 
    list_IDct = cell2mat(CellTracks( 1,: ));
    scc = list_IDct(1) ;
    % Add scc-th cell-track to GUI "Select Cell" DropDown list (as STRING cell array format)
    app.t5_SelectCellDropDown.Items =...
        cellfun(@mat2str, mat2cell( list_IDct, ...
        1 , ones( 1,length(list_IDct)) ),'UniformOutput',false);
    % Update GUI "Select Cell" to point to scc-th cell
    app.t5_SelectCellDropDown.Value = {mat2str(scc)};
    clearvars list_IDct ;
end


% --- Initialize FIGURE for BF stack -----------------------------------------------------%
APP_opt.fig = figure(99);   % Create figure and make it currently handle
gca(APP_opt.fig);           % Ensure the figure has axes




%% ***** MANUAL TRACKING *****************************************************
% The main code contain a while-loop that run continuously to track the
% points clicked and move betwenn frames


% Every time START is pushed, ensure to begin at first frame
APP_opt.t5_ff = 1;           % the currently displayed frame
app.t5_Edit_Frame.Value = APP_opt.t5_ff ;
           
% Calculate the number of digits and zeros to create correct filename.
N_dig = length(num2str(APP_opt.t5_ff)); 
N_null = repmat('0', [1, APP_opt.t5_T_dig - N_dig]);
% Read the image and store the frame sizes
imgBF = imread([APP_opt.t5_path_BF, APP_opt.t5_foldName_BF ,'/' APP_opt.t5_Prefix_BF '_' N_null num2str(APP_opt.t5_ff) '.tif']);
[APP_opt.t5_fr_H, APP_opt.t5_fr_W] = size(imgBF);

app.TextOUT.Value = sprintf('\n%s',  'Manual Tracking mode ...' );
% WHILE-loop run untill the user push STOP button in GUI of the program
while APP_opt.BREAK == 0
    ReFresh_Frame;           % REFRESH and update displayed frame

    % return the coordinates of the point clicked using the mouse cursor,
    % only if it was inside the tracking figure
    [x,y,button] = grabinput(1);
    
    % x and y are returned empty when the figure is closed manually by the user
    % This would create an error when grabinput(1) try to take the handle of
    % the now closed figure. Thus, we break the while-loop
    if (isempty(x) && isempty(y))
       break 
    end
    
    % We save coordinates only if user clicks inside the "frame", clicking 
    % outideand in the surrounding grey area has no effects.
    if ~(x<0 || x>APP_opt.t5_fr_W || y<0 || y>APP_opt.t5_fr_H)        
        
      % The time points tracked so far, those that are non-zero
      tracked_range = find(CellTracks{2,scc}(:,1)~=0);
      
      % Check that the tracking is continuous and has no interruption,
      % meaning that it is continous over time (/frames) and no stretch of zeros 
      % break the tracking in multiple pieces. (it cannot disappear for 
      % some frames and reappear later). We accept the [x,y] IF:
      if isempty(tracked_range) ||...                 % - it is the first point addedd to an empty track
         any(tracked_range == APP_opt.t5_ff) ||...    % - we modify any of the previously tracked frames
         APP_opt.t5_ff - tracked_range(end) == 1      % - the frame tracked comes after the last one in the cell-track  

            CellTracks{2 , scc}(APP_opt.t5_ff , 1:2) = [x,y] ;     % store [x,y] 
            
            APP_opt.t5_ff = APP_opt.t5_ff+1;          % point to next frame
            if APP_opt.t5_ff > LenStack               % avoid exciding stack length
                APP_opt.t5_ff = APP_opt.t5_ff-1;      % go to previous frame
            end
            app.t5_Edit_Frame.Value = APP_opt.t5_ff;
            
      else      % issue a Warning and do not store the last [x,y] input
          app.TextOUT.Value = sprintf('%s\n%s',  'You can only modify previously tracked frames,', ...
                      'or add more after last time point in current track, ', num2str(tracked_range(end)) );
          
          % Move and display the frame after the last time point tracked      
          APP_opt.t5_ff = tracked_range(end)+1;          
          app.t5_Edit_Frame.Value = APP_opt.t5_ff;
          ReFresh_Frame;    
      end
      
    end
    
end %/while
app.TextOUT.Value = sprintf('\n%s',  'Idle ...' );

% If the track is not a single consecutive block .... (the cell-track is fragmented)
l_Block = double( CellTracks{2 , scc}(:,1)~=0);
ori1 = strfind([0,l_Block'==1],[0 1]);
end1 = strfind([l_Block'==1,0],[1 0]);
% Then we issue a Warning
if length(end1 - ori1 + 1) > 1
    app.TextOUT.Value = sprintf('%s\n%s',  '   WARNING !!!   ' , ...
                                'The track is fragmented in', num2str(length(end1)), ' pieces!!! ' );
end

% Reset to reinitiate new round if wished
APP_opt.t5_ff = 1; 
% Reset total stack number in the GUI to undetermined
app.t5_Label_totFrames.Text = [ '   ' ];

end








