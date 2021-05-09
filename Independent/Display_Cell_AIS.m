function Display_Cell_AIS( cData, Img_BF )
%Display_Cell_AIS = Display the results from analysis via algorithm AIS
% (Average Signal).
% Show all main results for the cc-th cell and help to visualize and check
% the work in progress
%
% INPUTS ------------------------------------------------------------------
% cData = the specific cc-th cell to plot ( normally in the provided as
%         cellList.meshData{1,ff}{1,cc})
% Img_BF = the Bright field image to plot
% -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|- 


global APP_opt ;                        % Variable storing WHISIT options
extra_border = APP_opt.BorderBox ;      % extra border area to add to a cell cropped image


%% -- FLUORO CHANNEL 1 -----------------------------------------------------------------------------------
%  -------------------------------------------------------------------------------------------------------
    if ~isempty(cData.model)  &&  size(cData.mesh,2) == 4
    hFig2 = figure(2);           % set(hFig2, 'Position', [50 200 600 500]);
    clf(2);
    set(gca, 'Box', 'off');
    set(gca, 'Color', [1,1,1]);
    set(gcf, 'Color', [1,1,1]);
    
%-1--Cell Raw Signal (false color)
    subplot(1,3, 1);     hold on;	axis equal;
    imshow(cData.CH1.IC , [min(min(cData.CH1.IC)), max(max(cData.CH1.IC))])
    title('Raw Signal (false color)', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    colormap jet ;
    freezeColors ;        %freeze this plot's colormap
 
%-2--Raw Signal and Mesh
    subplot(1,3, 2);     hold on;	axis equal;
    imshow(cData.CH1.IC , [min(min(cData.CH1.IC)), max(max(cData.CH1.IC))])
    title('Raw Signal and Mesh', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    plot( cData.R_mesh(:,1) , cData.R_mesh(:,2),'.', 'Color',[0.9, 0.2, 0]);
    plot( cData.R_mesh(:,3) , cData.R_mesh(:,4),'.', 'Color',[0.9, 0.8, 0]);
    plot( cData.geom.R_axis(:,1) , cData.geom.R_axis(:,2) ...
        , '.-', 'Color', [0 0.6 1]);     

%-3--Bright Field and Cell Contour
    min_x = min( cData.model(:, 1) );    max_x = max( cData.model(:, 1) );
    min_y = min( cData.model(:, 2) );    max_y = max( cData.model(:, 2) );    
    I_BF = imcrop( Img_BF , [(min_x-extra_border), (min_y-extra_border), ((max_x-min_x)+extra_border*2), ((max_y-min_y)+extra_border*2) ]);
    subplot(1,3, 3);     hold on;	axis equal;
    imshow(I_BF, [min(min(Img_BF)), max(max(Img_BF))]);  
    title('Bright Field and Cell Contour', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);    
    plot(cData.R_model(:, 1), cData.R_model(:, 2), '*y','MarkerSize' , 2); 

    
%% -- FLUORO CHANNEL 2 -----------------------------------------------------------------------------------
%  -------------------------------------------------------------------------------------------------------
    if APP_opt.t1_choose_Chan == 2        % ANALYSE Chan 2  
        hFig3 = figure(3);                % set(hFig3, 'Position', [700 200 600 500]);
        clf(3);
        set(gca, 'Box', 'off');
        set(gca, 'Color', [1,1,1]);
        set(gcf, 'Color', [1,1,1]);
%-1--Cell Raw Signal (false color)
        subplot(1,3, 1);     hold on;	axis equal;
        imshow(cData.CH2.IC , [min(min(cData.CH2.IC)), max(max(cData.CH2.IC))])
        title('Raw Signal (false color)', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
        colormap jet ;
        freezeColors ;        %freeze this plot's colormap

%-2--Raw Signal and Mesh
        subplot(1,3, 2);     hold on;	axis equal;
        imshow(cData.CH2.IC , [min(min(cData.CH2.IC)), max(max(cData.CH2.IC))])
        title('Raw Signal and Mesh', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
        plot( cData.R_mesh(:,1) , cData.R_mesh(:,2),'.', 'Color',[0.9, 0.2, 0]);
        plot( cData.R_mesh(:,3) , cData.R_mesh(:,4),'.', 'Color',[0.9, 0.8, 0]);
        plot( cData.geom.R_axis(:,1) , cData.geom.R_axis(:,2) ...
            , '.-', 'Color', [0 0.6 1]);     

%-3--Bright Field and Cell Contour
        min_x = min( cData.model(:, 1) );    max_x = max( cData.model(:, 1) );
        min_y = min( cData.model(:, 2) );    max_y = max( cData.model(:, 2) );   
        I_BF = imcrop( Img_BF , [(min_x-extra_border), (min_y-extra_border), ((max_x-min_x)+extra_border*2), ((max_y-min_y)+extra_border*2) ]);
        subplot(1,3, 3);     hold on;	axis equal;
        imshow(I_BF, [min(min(Img_BF)), max(max(Img_BF))]);  
        title('Bright Field and Cell Contour', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);    
        plot(cData.R_model(:, 1), cData.R_model(:, 2), '*y','MarkerSize' , 2); 
   
    end
    
%---WAIT to examin    
    set(gcf, 'Color', [1,1,1]);
    pause(APP_opt.plot_pause) ;
   
    end
    

end






