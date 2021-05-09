function GUI_t1_Load_Channel_2(app)
%GUI_t1_Load_Channel_2 = upon bush button, the user provide a folder with 
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
    app.TextOUT.Value = sprintf('\n%s\n%s',  'No path provided for Channel 2 stack !!!');
    return
end            
[PathName, FoldName]  = fileparts(FullPath) ;
APP_opt.t1_path_CH2 = PathName;
APP_opt.t1_foldName_CH2 = FoldName ;
app.t1_Edit_CH2_Fold.Value = [APP_opt.t1_path_CH2  filesep  APP_opt.t1_foldName_CH2];

% Access stack folder and index all the file inside the folder that ends with .tif
APP_opt.t1_srcFiles_CH2 = dir( fullfile([APP_opt.t1_path_CH2  filesep  APP_opt.t1_foldName_CH2], '*.ti*f'));
% Ensure that the list of files is sorted
cell_srcFiles = struct2cell( APP_opt.t1_srcFiles_CH2 );   
[sort_name,idx] = sort( cell_srcFiles(1,:) );
APP_opt.t1_srcFiles_CH2 = APP_opt.t1_srcFiles_CH2(idx);

% Take first .tif filename as prefix. Filename format must be in form:  
% prefix_NNNNN.tif  ---  The program takes the number between '_' an '.tif'
cnt = -1;            yy = 0;
while cnt == -1  &&  yy < size(APP_opt.t1_srcFiles_CH2, 1)
    yy = yy+1;
    % split first filename that is a .tif
    strFile = strsplit( APP_opt.t1_srcFiles_CH2(yy).name , APP_opt.name_delimiters);
    if length(strFile) >=3  &&  strcmp(strFile(end), 'tif')
        APP_opt.t1_Prefix_CH2 = strsplit(sort_name{yy} , '_');
        APP_opt.t1_Prefix_CH2 = APP_opt.t1_Prefix_CH2(1:end-1);
        APP_opt.t1_Prefix_CH2 = strjoin(APP_opt.t1_Prefix_CH2 , '_');
        cnt = +2;
    end 
end
if cnt == -1        % then we have not found any tif correctly named
    app.TextOUT.Value = sprintf('\n%s\n%s',  ['!!! Stack Ch. 2 contains uncorrect .tif filename !!!'],...
                               ['Format should be in form:   filename_xxxx.tif']);
    app.TextOUT.BackgroundColor = [0.75 0.3 0.3] ;
    return;
end 




