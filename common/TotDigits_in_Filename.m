function digits = TotDigits_in_Filename(FileList, delimiters)
%
%TotDigits_in_Filename - examine filename and how namy digits are present
%in a .tiff filename
%
%   The take a list of filenames from a folder and search the first .tiff
%   filename availble. then using delimiter select the section that contain
%   the stack number, count and return the number of digits in the name. 
%
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-

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


