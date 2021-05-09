function digits = TotDigits_in_Filename(FileList, delimiters)
%
%TotDigits_in_Filename - examine filename and find how namy digits are 
%   present in filename between the '_' and '.tif' delimiter
%
%   The take a list of filenames from a folder and search the first .tiff
%   filename availble. then using delimiter select the section that contain
%   the stack number, count and return the number of digits in the name. 
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

digits = -1;            
ii = 0;
while digits == -1  &&  ii < size(FileList, 1)
    ii = ii+1;
    % split first filename that is a .tif
    strFile = strsplit( FileList(ii).name , delimiters);
    if length(strFile) >=3  &&  strcmp(strFile(end), 'tif')
        digits = length(strFile{end-1});                 
    end 
end
    
end


