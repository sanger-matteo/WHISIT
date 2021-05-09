function [ g_List , g_Fold ] = Create_GenList
% Create_GenList = Take all the .mat files in the folder provided by the
%   user and create a single List of clones, organized according to the
%   generation
%
%   generation_List = each element stores all the clones of x-th generation
%
%   generation_Fold = store the Generation number and folder path where to
%                     place data (.txt file and figures)
%
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-

global APP_opt ;

srcFiles = APP_opt.t3_intergen_srcFiles ;
t_clone_List = [];

if size( srcFiles ,1) == 1                                      % if there is only one file    
    load([srcFiles(1).folder,'/',srcFiles(1).name]);            % load the only file present
    clone_List(2, 1:size(clone_List,2)) = {srcFiles(1).folder} ;
    clone_List(3, 1:size(clone_List,2)) = {srcFiles(1).name} ;
    t_clone_List = clone_List;
    
elseif size( srcFiles ,1) >= 2                                  % if there are more than two file    
    load([srcFiles(1).folder,'/',srcFiles(1).name]);            % load first file
    clone_List(2, 1:size(clone_List,2)) = {srcFiles(1).folder} ;
    clone_List(3, 1:size(clone_List,2)) = {srcFiles(1).name} ;
    t_clone_List = clone_List;
    % load all the remaning and concatenate the clones in t_clone_List, in
    % order to create a single comprehensive list
    for ii = 2 :  size( srcFiles ,1)
        load([srcFiles(ii).folder,'/',srcFiles(ii).name]);
        clone_List(2, 1:size(clone_List,2)) = {srcFiles(ii).folder} ;     % assign Folder name and ...
        clone_List(3, 1:size(clone_List,2)) = {srcFiles(ii).name} ;       % file name the clone belongs to the 
        t_clone_List = [t_clone_List , clone_List] ;
    end    
end

% Assign the generation number to each clone, based on the length of ID_clone
for tt = 1 : size(t_clone_List,2)
    IDstr = strsplit( t_clone_List{1,tt}{1}.ID_clone, '.' );
    t_clone_List(4, tt) = { (size(IDstr{2},2) +1) } ;     % assigne generation number 
end

% Find all unique generations presents in entire t_clone_List
uniqueGen =  unique( cell2mat( t_clone_List(4, :))) ;

% generation_List = each element stores all the clones of the uu-th generation
g_List = {};
for uu = 1 : length(uniqueGen)
    % Find all clones of the uu-th generation in t_clone_List and gather them together in new cell array
    idx_uu = find(cell2mat(t_clone_List(4,:)) == uniqueGen(uu) );
    g_List{uu} = t_clone_List(:, idx_uu);
end



% generation_Fold store the Generation number and folder path where to store data
g_Fold = {};
g_Fold(1, 1:size(g_List,2)) = mat2cell( [1:size(g_List,2)], 1, ones(1,size(g_List,2))); 
        
% Based on the size of generation_List, we create the folders for storing the results
for gg = 1 : size(g_List,2)
    if  isempty(APP_opt.t3_exp_name)
        pathname = [ g_List{gg}{2,1} '/Res_Generations' ];
    else
      if APP_opt.t3_choose_Save_txt == 1
          pathname = [ g_List{gg}{2,1} '/' APP_opt.t3_exp_name '_Res_Generations' ];
      else
          pathname = [ g_List{gg}{2,1} '/Res_Generations' ];
      end            
    end  
    mkdir([ pathname '/G_' num2str(g_List{gg}{4,1}) ]) ;
    g_Fold(2, gg) = {[ g_List{gg}{2,1} '/Res_Generations' ]};
    g_Fold(3, gg) = {[ 'G_' num2str(g_List{gg}{4,1}) ]};
end




end % Main fnc






