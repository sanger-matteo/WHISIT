function  Lng_Create_clone_List(app)
% Lng_Create_clone_List = creats a new reorganized database: clone_List
% 
% The function load the life provided by the user:
%  >> [APP_opt.t1_path_Res_D2T , APP_opt.t1_file_Res_D2T]
%
% The algorithm main step will separate all "segments" that constiture each
% cell-track. A segment is the lifespan of a clone: going from birth/first
% appereance to next div/moving out of field of view. Each clone is given a 
% unique ID that allows to trace back all its lineage.
%
% The output clone_List.mat contains the cellTrack and clone_List, with the
% latter being an array of cellList, where each element is a single clone
% tracked at each time point.
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

%% 

global APP_opt ;
 
app.TextOUT.Value = sprintf('\n%s',  'Creating clone_List ... ');

load([APP_opt.t1_path_Res_D2T , APP_opt.t1_file_Res_D2T]);
cList = cellList.meshData;

% Check that we have all correct data to create clone_List
varlist = who('-file', [APP_opt.t1_path_Res_D2T , APP_opt.t1_file_Res_D2T]);
check_Track = any(strcmp(varlist,'cellTrack')) ;
check_LG = any(strcmp(varlist,'LG')) ;

if check_Track == 0 || check_LG == 0  ||  WHISIT_param.choice_ManualTrack == 0
    % if anything is missing, then we cannot create clone list
    app.TextOUT.Value = sprintf('\n%s',  'File provided is invalid for creating cloneList.mat!!!');
    APP_opt.ERROR = 1;          
    return;
end

% Logical matrix that show the time window for each cell-track. This
% simplify future search and avoid using longer command, such as:
% find( cellTrack{4, LG{cc}.Offspring(dd)}(:,1) ~=0 );    
tw_tracks = ( cell2mat(cellTrack(4,:)) ~=0 ) ;                         

% Lineage variable 
% LG = based on cellTrack, summarize all information for each specific
% cell-track, division events, clone segments, number ID etc...
cc_det = cell2mat(cellTrack(4,:));        % the cell detection position in cList at each time point fore each cell-track (column)



%% --- CREATE clone_List --------------------------------------------------
% This step separate all the cell-track segments that constitutes each
% cell-track and rearrange them in a simple struct array, the clone List.
% Each clone lifespan goes from birth : next div and carries a unique ID
% that allow to trace back all its lineage.

nn = 0;                                             % CLONE counter
for cc = 1:size(tw_tracks,2)                        % go thourgh each cell-track
    
    for seg = 1 : length(LG{cc}.Clones_Seg)         % go in each clone
        nn = nn+1;                                  % every new segment, it is a new clone.
        [fr] = strsplit(LG{cc}.Clones_Seg{seg},':');
        sg_str = str2num(fr{1}) ;    
        if strcmp(fr{end}, 'end');    sg_end = LG{cc}.Finish;
        else ;                        sg_end = str2num(fr{end}) ;          % f_end = str2num(fr{end})-1 ; 
        end
                   
        ndet_cc = cc_det(sg_str ,cc) ;             % position of cc in cList  
        
        if seg < length(LG{cc}.Clones_Seg)
            % Find position of the offspring in cList 
            cln = find( cell2mat(cellTrack(1,:)) == LG{cc}.Offspring(seg));
            ndet_offspring = cc_det(sg_end+1, cln) ;                         
            clone_descendants = { cList{sg_end+1}{ndet_offspring}.ID_clone ,...
                                  cList{sg_end+1}{ndet_cc}.ID_clone };
        elseif seg == length(LG{cc}.Clones_Seg)
            clone_descendants = [];
        end
            
        for ff = sg_str : sg_end  
            % position of cc in cList
            ndet_cc = cc_det(ff ,cc) ;             
            
            clone_List{nn}{ff-sg_str+1} = cList{ff}{ndet_cc} ;             
            % Below we update some variables:
            clone_List{nn}{ff-sg_str+1}.offspring_ID_clone = clone_descendants ;            
            clone_List{nn}{ff-sg_str+1}.frame_birth = sg_str ;
            clone_List{nn}{ff-sg_str+1}.frame_last  = sg_end ;       
        end  

    end % for seg
end % for cc


% Remove some fields that are related to manual tracking which are unimporant
% or misguiding in future use of clone_List
 remove_fields = {'divisions','descendants','ancestors','stage', 'birthframe',...
                  'frame_nextgen','tracks_offspring','track_ancestor', ...
                  'timelapse', 'algorithm'};
              
for cc = 1 : size(clone_List, 2)
    for ff = 1 : size(clone_List{cc}, 2)       
        clone_List{cc}{ff} = rmfield( clone_List{cc}{ff}, remove_fields) ;
    end
end


%% --- SAVE cloneList.mat -------------------------------------------------
app.TextOUT.Value = sprintf('\n%s',  'Saving clone_List file ... ');
if isempty(APP_opt.t1_exp_name)
    save([APP_opt.t1_path_Res_D2T, '/clone_List.mat'] ,...
        'clone_List', 'cellTrack', 'WHISIT_param');
else
    save([APP_opt.t1_path_Res_D2T, '/'  APP_opt.t1_exp_name '_clone_List.mat'] ,...
        'clone_List', 'cellTrack', 'WHISIT_param');
end

app.TextOUT.Value = sprintf('\n%s',  '...Idle...');
   

end % MAIN fnc