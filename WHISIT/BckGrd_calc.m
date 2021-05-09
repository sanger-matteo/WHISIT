function [BkGr, min_px, J_M] = BckGrd_calc( f )
% Background calculator
%
% The script analize the backround signal of a given frame (f) and calculates:
% J_M  ---> create a mask of the background noise, excluding any signal above a
%          treshold which is considered noise
%          ( S_Lim is calculated for avery frame, by analysing the frequency 
%          distribution of the signal pixel value)
% BkGr ---> calculate the average signal in the backroung of the frame
% min_px -> calculate the minimum pixel value of the whole frame
%

%% --- Initialize ---------------------------------------------------------
%--------------------------------------------------------------------------
global GUI_opt;
global my_DB;

% Quick tips about functions used:
% --> mat2gray : transform a uintX format into a double and change greyscale range to 0 - 1
% --> imopen   : morpholoical erosion, allow to find the image background: shape
%                and area size of STREL allow to erode image elements and leave only backgr.
% --> bwlabel  : each independent connected region of 1's is given a label-number
%                like a mask, but the 1's are substituted with label-numbers
% --> im2bw    : set a threshold to detect elements and create a 0-1 mask
                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                        
%% --- Body of the Function -----------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if f <10                            null = '000';
elseif f >= 10 && f < 100           null = '00';
elseif f >= 100 && f < 1000         null = '0';
else                                null = '';
end
I = double(imread([ GUI_opt.path_FL, '/FL_' , null num2str(f), '.tif'])) ;    % handle to f-th frame

J = mat2gray(I) ;
J2 = imopen(J, strel('disk',40)) ;
J_nobg = J - J2 ;          % remove background

% --> Get histogram of frequency of pixels values and fit to Normal Gaussian 
% Find mu and delta and use to consider signal below S_Lim as noise      
[counts,x] = imhist(J_nobg,1000);       % histogram of intensity of image pixels
% stem(x,counts);
f = fit(x, counts, 'gauss1');
coeff =  coeffvalues(f);                % [a, mu, delta]
S_Lim = coeff(2) + 3*coeff(3);          % mu + 3*delta

t_J_M = im2bw(J_nobg, S_Lim) ;   
J_M = bwlabel(t_J_M) ;  
% Apply a not(Cell_Mask) to isolate all the background.
BkGr = sum(sum(I.* (~J_M))) / sum(sum(~J_M));    %mean(mean(I.* (~J_M))) ;
min_px = min(min(I)) ;

% --- PLOT ----------------------------------------------------------------
if GUI_opt.choice_plot == 1            
    % Display frame, background mask and average pixel value of background
    figure(3);  clf(3);  hold on;    
    annotation('textbox', [0 0.9 1 0.08], 'String', 'Background Subtraction', ...             
               'EdgeColor', 'none', 'HorizontalAlignment', 'center',...
               'Color',[0.4 0.4 0.4], 'FontSize', 14 );
    %title('Background Subtraction', 'Color',[0.5 0.5 0.5], 'FontSize', 14 );
    subplot(2,2,1);  imshow(J, [min(min(J)), max(max(J))] )
    subplot(2,2,2);  imshow(imadjust(J_nobg))
    subplot(2,2,3);  imshow(J.*(~J_M), [min(min(J)), max(max(J))]); 
    subplot(2,2,4);  hold on;   plot(f,x,counts);  
                     plot([S_Lim,S_Lim],[0,max(counts)], '--g', 'LineWidth',2);
                     text(S_Lim, max(counts)/2, ['Avg noise value :  ', num2str(BkGr, '%6.0f\n')] , 'Color',[0.5 0.5 0.5], 'FontSize', 14 );
                     ylabel('Counts', 'FontSize',12);     xlabel('Pixel value', 'FontSize',12);
                     %imshow(~J_M)    
    colormap jet ;
    set(gcf, 'Color', [1,1,1]);
    %pause(GUI_opt.plot_pause) ;       
      
end

end

%--------------------------------------------------------------------------
%       Copyright (c) 2016, Matteo Sangermani, All rights reserved
%--------------------------------------------------------------------------

