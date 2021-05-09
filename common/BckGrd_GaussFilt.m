function [BkGr, min_px, max_px] = BckGrd_GaussFilt( I, sigma, X_delta )
% BckGrd_calc = Calculate the average background signal of an image
%
% INPUT ------------------------------------------------------------------
% I     = the image to analyse, must be provided as type double().
%
% sigma = the standard deviation value for 2-D Gaussian smoothing kernel
%         function imgaussfilt( I, sigma) [2 or 3 are usually good values].
%         The higher the more "blurred" is the resulting image.
%
% X_delta = we will use a threshold to discriminate pixel that are "signal"
%         from noise (value == 0). This is done by fitting the histogram of
%         pixel values for I to a gaussian diastribution. The obtained mu,
%         average value, and sigma, the standard deviation, are use to find  
%         the threshold as >> mu + X_delta*sigma
%         X_delta is the multiplication factor for sigma: i.e. 2*sigma
%         cover 95% of the gaussian distribution, 3*sigma about 99%
%         [N.B.: in an image I the vast majority of pixel are "background",
%         unless it is super crowded by bacteria/signal. Therefore, we
%         usually the entire gaussian distribution is noise and we aim to
%         exclude it and consider only the exteme/highest pixel values]
%
%
% OUTPUT ------------------------------------------------------------------
% J_M  ---> create a mask of the background noise, excluding any signal above a
%          treshold which is considered noise
%          ( S_Lim is calculated for avery frame, by analysing the frequency 
%          distribution of the signal pixel value)
% BkGr ---> calculate the average signal in the backroung of the frame
% min_px -> calculate the minimum pixel value of the whole frame
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-


% Quick tips about functions used:
% --> mat2gray : transform a uintX format into a double and change 
%                greyscale range to 0 - 1
% --> imgaussfilt : filters image I with a 2-D Gaussian smoothing kernel
%                with standard deviation specified by sigma.
% --> imopen   : morpholoical erosion, allow to find the image background
%                shape and area size of STREL allow to erode image elements
%                and leave only backgr.
% --> bwlabel  : each independent connected region of 1's is given a
%                label-number like a mask, but the 1's are substituted with
%                label-numbers
% --> im2bw    : set a threshold to detect elements and create a 0-1 mask
 

%% --- Body of the Function -----------------------------------------------

J = mat2gray(I);
Jgauss = imgaussfilt(J, sigma);        % sigma: best values are 2 or 3

% --> Get histogram of frequency of pixels values and fit to Normal Gaussian 
% Find mu and delta and use to consider signal below S_Lim as noise      
[counts,x] = imhist(Jgauss,500);        % histogram of intensity of image pixels
f = fit(x, counts, 'gauss1');
coeff =  coeffvalues(f);                % [a, mu, delta]
% S_Lim is the signal threshold separating signal and noise. Anything above
% is valid real signal, while everything below is noise.
S_Lim = coeff(2) + X_delta*coeff(3);      % mu + 3*delta

% We use S_Lim as threshold to create a mask that define the signal area 
% (value == 1) and noise (value == 0).
t_J_M = im2bw(Jgauss, S_Lim) ;   
J_M = bwlabel(t_J_M) ;  
% Apply a not(Cell_Mask) to isolate all the background from original I.
BkGr = sum(sum(I.* (~J_M))) / sum(sum(~J_M));    %mean(mean(I.* (~J_M))) ;
min_px = min(min(I)) ;
max_px = max(max(I)) ;

% --- PLOT ----------------------------------------------------------------
% if APP_opt.choice_plot == 1            
%     % Display frame, background mask and average pixel value of background
%     figure(2);  clf(2);  hold on;    
%     annotation('textbox', [0 0.9 1 0.08], 'String', 'Background Subtraction', ...             
%                'EdgeColor', 'none', 'HorizontalAlignment', 'center',...
%                'Color',[0.4 0.4 0.4], 'FontSize', 14 );
%     %title('Background Subtraction', 'Color',[0.5 0.5 0.5], 'FontSize', 14 );
%     subplot(2,2,1);  imshow(J, [min(min(J)), max(max(J))] )
%     subplot(2,2,2);  imshow(imadjust(J_nobg))
%     subplot(2,2,3);  imshow(J.*(~J_M), [min(min(J)), max(max(J))]); 
%     subplot(2,2,4);  hold on;   plot(f,x,counts);  
%                      plot([S_Lim,S_Lim],[0,max(counts)], '--g', 'LineWidth',2);
%                      text(S_Lim, max(counts)/2, ['Avg noise value :  ', num2str(BkGr, '%6.0f\n')] , 'Color',[0.5 0.5 0.5], 'FontSize', 14 );
%                      ylabel('Counts', 'FontSize',12);     xlabel('Pixel value', 'FontSize',12);
%                      %imshow(~J_M)    
%     colormap jet ;
%     set(gcf, 'Color', [1,1,1]);
%     %pause(GUI_opt.plot_pause) ;            
% end

end





