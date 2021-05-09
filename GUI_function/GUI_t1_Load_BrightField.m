function GUI_t1_Load_BrightField(app)
%GUI_t1_Load_BrightField = upon bush button, the user provide a folder with 
%   a stack of images whose path, size and filenames are checked and stored
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


global APP_opt ;

FullPath = uigetdir; 
if FullPath == 0
       FullPath = [];
       app.TextOUT.Value = sprintf('\n%s\n%s',  'No path provided for Bright Field stack !!!');
       return
end            
[PathName, FoldName]  = fileparts(FullPath) ;            
APP_opt.t1_path_BF = PathName ;
APP_opt.t1_foldName_BF = FoldName ;
app.t1_Edit_BF_Fold.Value = [APP_opt.t1_path_BF  filesep  APP_opt.t1_foldName_BF];

% Access stack folder and index all the file inside the folder that ends with .tif
APP_opt.t1_srcFiles_BF = dir( fullfile([APP_opt.t1_path_BF  filesep  APP_opt.t1_foldName_BF], '*.ti*f'));
% Ensure that the list of files is sorted
cell_srcFiles = struct2cell( APP_opt.t1_srcFiles_BF );   
[sort_name,idx] = sort( cell_srcFiles(1,:) );
APP_opt.t1_srcFiles_BF = APP_opt.t1_srcFiles_BF(idx);

% Take first .tif filename as prefix. 
% Filename format must be in form:  prefix_NNNNN.tif ; separation by '_' is necessary.
cnt = -1;            yy = 0;
while cnt == -1  &&  yy < size(APP_opt.t1_srcFiles_BF, 1)
    yy = yy+1;
    % split first filename that is a .tif
    strFile = strsplit( APP_opt.t1_srcFiles_BF(yy).name , APP_opt.name_delimiters);
    if length(strFile) >=3  &&  strcmp(strFile(end), 'tif')
        APP_opt.t1_Prefix_BF = strsplit(sort_name{yy} , '_');
        APP_opt.t1_Prefix_BF = APP_opt.t1_Prefix_BF(1:end-1);
        APP_opt.t1_Prefix_BF = strjoin(APP_opt.t1_Prefix_BF , '_');
        cnt = +2;
    end 
end
if cnt == -1        % then we have not found any tif correctly named
    app.TextOUT.Value = sprintf('\n%s\n%s',  ['!!! Bright Field stack provided contains uncorrect .tif filename !!!'],...
                               ['Format should be in form:   filename_xxxx.tif']);
    app.TextOUT.BackgroundColor = [0.75 0.3 0.3] ;
    return;
end