function Display_BF_Frame
%
%Display_BF_Frame - Show the "APP_opt.t5_ff" bright field image in the
%   manual tracking figure
%
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-

global APP_opt;

% Calculate the number of digits and zeros to create filename.
N_dig = length(num2str(APP_opt.t5_ff)); 
N_null = repmat('0', [1, APP_opt.t5_T_dig - N_dig]);
% Read the image an show it in the figure
imgBF = imread([APP_opt.t5_path_BF, APP_opt.t5_foldName_BF ,'/' APP_opt.t5_Prefix_BF '_' N_null num2str(APP_opt.t5_ff) '.tif']);
imshow(imgBF, [min(min(imgBF)) max(max(imgBF))]);

end