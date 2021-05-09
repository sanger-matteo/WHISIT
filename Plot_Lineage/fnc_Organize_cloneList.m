function [CL_s_EndLine, CL_n_EndLine, CL_s_All] = fnc_Organize_cloneList(clone_List)

% Take all clone_List.ID_clone and create several Lists that will be useful 
% for quickly searching clone_List, ordering, sorting them, convert to 
% numeric etc...

CL_s_All = {};         % all unique clone_List.c_ID
ID_pre = {};           % all prefix ID_clone that identify the lineag(s) 
ID_dot_post = {};      % postfix (with '.')
ID_post = {};          % postfix (without the '.')
jj = 1; 
for cc = 1 : size(clone_List,2)
    CL_s_All{cc} = clone_List{cc}{1}.ID_clone;
    [String] = strsplit(CL_s_All{cc}, '.');
    ID_pre{cc} = String{1};
    if ~isempty(String{2})              % if isempty, then is ancestor cell
        ID_post{jj} = String{2};
        ID_dot_post{jj} = ['.' String{2}];
        jj = jj +1;
    end
end

% Create a list of unique lineage prefix(-es)
Lin_Num = intersect(ID_pre, ID_pre);

%%% !!! INCOMPLETE !!!
%%%---> UPDATE GUI and let user choose which lineage he wants to plot

% Take all ID_post and fill with zeros at the end, such as they become 
% number-strings all haveing same number of digits (and as many as the 
% longest ID_clone). In such way we can compare and sort them numerically.
ID_Num_post = [];
N_Max = max(cellfun( @length, ID_post));
for cc = 1 : length(ID_post)
    len_str = length(ID_post{cc}) ;
    if len_str < N_Max
        ID_Num_post(cc) = str2num([ID_post{cc} num2str(zeros(1,N_Max-len_str),'%d')]);
        ID_post{cc} = [ID_post{cc} num2str(zeros(1,N_Max-len_str),'%d')];
    elseif len_str == N_Max
        ID_Num_post(cc) = str2num([ID_post{cc}]);
        ID_post{cc} = [ID_post{cc}];
    end
end
ID_Num_post = ID_Num_post';   


% For each ID_clone find how many share same root using ID_dot_post list.
% Roughly correspond to the number of all descendants a given ID_clone have. 
% Now we are looking for those who have 0, because those are the clones
% thare are terminal line (end_line) in the lineage tree
N_root = [];
for cc = 1 : length(ID_post)
    c1 = cellfun( @length, strfind(ID_dot_post, ID_dot_post{cc}));
    N_root(cc) = length(c1(c1 == 1)) -1;       % -1 because there is always one corresponding root, that of the cell itself
end


% Now we have end_lines, but we still have to decide where to position each
% within the lineage tree and have (twin) daughter cells together. 
% Solution is gather the end_line cells and sort them numerically. We then
% transform back in string form
idx_end_line = find(N_root==0);
CL_n_EndLine = ID_Num_post(idx_end_line);       % C_EndLine in numeric form
CL_n_EndLine = sort(CL_n_EndLine,'ascend') ;
CL_s_EndLine = {};                              % C_EndLine in string form
for cc = 1 : length(CL_n_EndLine)
    C_zeros = num2str(CL_n_EndLine(cc,:));              % we remove the zeros ...
    CL_s_EndLine{cc} = [ '.' C_zeros(C_zeros~='0') ];   % and add the '.'
end

%--------------------------------------------------------------------------
% FROM NOW ON we will use only cN_End_Line and cS_End_Line to navigate us
% in which order to plot the lineage tree. This will be updated and change
% every time we draw clone(s), common ancestors... progressing "back in
% time". They simply keep track of the "EndLine" clones and latest couples
% that need to be plotted.
% We also  keep using cS_All which  simply store ID_Clone names in same
% order as in clone_List, making searching far more easier.
%--------------------------------------------------------------------------


end