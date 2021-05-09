function Display_BF_Frame
%
%Display_BF_Frame - Show the "APP_opt.t5_ff" bright field image in the
%   manual tracking figure
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------
global APP_opt;

% Read the image an show it in the figure
imgBF  = double( imread( [ APP_opt.t5_path_BF  filesep   APP_opt.t5_foldName_BF  filesep   APP_opt.t5_srcFiles_BF(APP_opt.t5_ff).name]) ); 
imshow(imgBF, [min(min(imgBF)) max(max(imgBF))]);

end