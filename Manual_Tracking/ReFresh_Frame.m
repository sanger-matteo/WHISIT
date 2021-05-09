function ReFresh_Frame
%
%ReFresh_Frame - Refresh the image in the manual tracking figure, including
%   all selected option (tracks, frame,...)
%
%   The function will simply run all the major function for the
%   visualization of the manual tracking figure:
%   - Display_BF_Frame
%   - Display_TrackAid
%   - Display_Det_Outline
%   - Display_Tracks
%   In doing so update any change that was made or selection is now visible
%   and updated in the tracking figure
%
%   ---> Notes about variables used:
%   CellTracks = a global variable where that stores the information of
%   every cell track. 
%   scc = store the ID-number of the currently selected cell-track
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

global APP_opt;	    global CellTracks;     global scc;

% The manual tracking figure has a DUMMY border --- tagged 'RatioKeeper'
% Keeping a dummy plot(line) allows to manually resize the figure and
% during refreshing the figure size is maintained and does not revert to
% original size

% Deletes from the figure all elements except (dummy) frame border
children = get(gca, 'children');    % find all objects in the figure
for ii = 1 : length(children)
    if isempty(children(ii).Tag)    % delete objects with no Tag
        delete(children(ii));       % the one with Tag = 'RatioKeeper'
    end                             % allow to keep the frame resized
end

% displayed current bright field image
Display_BF_Frame;   

% If selected, show colored frame to aid tracking 
if APP_opt.t5_display_TrackAid ~= 0
    Display_TrackAid;
else                % plot DUMMY line, to maintained resized figure
    hold on;        % with the latest dimentions
    plot([0,0], [0,1], 'LineWidth', 0.1, 'Tag', 'RatioKeeper' );
    hold off;
end

% if checked, show cells outlines.
if APP_opt.t5_display_Det ~= 0      
    hold on ;
    Display_Det_Outline;
    hold off ;
end

% if checked, display track(s)
if APP_opt.t5_display_Track ~= 0
    hold on ;
    Display_Tracks;
    hold off ;
end

end



