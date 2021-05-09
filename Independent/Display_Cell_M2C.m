function InF_Display_Cell_M2C( cData, Img_BF )
%Display_Cell_M2C = Display the results from analysis via algorithm M2C
% (Membrane 2 Cytosol).
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
    if ~isempty(cData.model)
    
    % Releative mesh coordinates of the cell
    m_Xs = [cData.mesh(:,1) ; flipud(cData.mesh(:,3))] ;
    m_Ys = [cData.mesh(:,2) ; flipud(cData.mesh(:,4))] ;    
    
    hFig2 = figure(2);           % set(hFig2, 'Position', [50 200 600 500]);
    clf(2);
    set(gca, 'Box', 'off');
    set(gca, 'Color', [1,1,1]);
    set(gcf, 'Color', [1,1,1]);
    
%-1--Cell Raw Signal (false color)
    subplot(6,6,[1 2 7 8]);     hold on;	axis equal;
    imshow(cData.CH1.IC , [min(min(cData.CH1.IC)), max(max(cData.CH1.IC))])
    title('Raw Signal (false color)', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    colormap jet ;
    freezeColors ;        %freeze this plot's colormap
 
%-2--Raw Signal and Mesh
    subplot(6,6,[3 4 9 10]);     hold on;	axis equal;
    imshow(cData.CH1.IC , [min(min(cData.CH1.IC)), max(max(cData.CH1.IC))])
    title('Raw Signal and Mesh', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    plot( cData.R_mesh(:,1) , cData.R_mesh(:,2),'.', 'Color',[0.9, 0.2, 0]);
    plot( cData.R_mesh(:,3) , cData.R_mesh(:,4),'.', 'Color',[0.9, 0.8, 0]);
    plot( cData.geom.R_axis(:,1) , cData.geom.R_axis(:,2) ...
        , '.-', 'Color', [0 0.6 1]);     

%-3--Bright Field and Cell Contour
    min_x = min( m_Xs );        max_x = max( m_Xs );
    min_y = min( m_Ys );        max_y = max( m_Ys );   
    I_BF = imcrop( Img_BF , [(min_x-extra_border), (min_y-extra_border), ((max_x-min_x)+extra_border*2), ((max_y-min_y)+extra_border*2) ]);
    subplot(6,6,[5 6 11 12]);     hold on;	axis equal;
    imshow(I_BF, [min(min(Img_BF)), max(max(Img_BF))]);  
    title('Bright Field and Cell Contour', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);    
    plot(cData.R_P_model(:,1) , cData.R_P_model(:,2),'*y','MarkerSize' , 2); 
 
%-4--Cell Membr. Area        (+ contour) 
    subplot(6,6,[13 14 19 20]);     hold on;	axis equal; 
    I_Mem =  cData.CH1.IC .* cData.Mask_Memb ;
    imshow(I_Mem, [min(min(cData.CH1.IC)), max(max(cData.CH1.IC))]) ;
    plot(cData.R_model(:,1) , cData.R_model(:,2),'*w','MarkerSize' , 2);        
    colormap jet ;
    freezeColors ;
    title('Cell Membr. Area', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
 
%-5--Cell Cytosol Area       (+ contour)
    subplot(6,6,[15 16 21 22]);     hold on;	axis equal; 
    I_Cyt =  cData.CH1.IC .* cData.Mask_mCyto ;
    imshow(I_Cyt, [min(min(cData.CH1.IC)), max(max(cData.CH1.IC))] ) ;  
    plot(cData.R_model(:,1) , cData.R_model(:,2),'*w','MarkerSize' , 2);       
    colormap jet ;
    freezeColors ;
    title('Cell Cytosol Area', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);

    
%% -- FLUORO CHANNEL 2 -----------------------------------------------------------------------------------
%  -------------------------------------------------------------------------------------------------------
    if APP_opt.t1_choose_Chan == 2        % ANALYSE Chan 2  
        hFig3 = figure(3);                % set(hFig3, 'Position', [700 200 600 500]);
        clf(3);
        set(gca, 'Box', 'off');
        set(gca, 'Color', [1,1,1]);
        set(gcf, 'Color', [1,1,1]);
%-1--Cell Raw Signal (false color)
        subplot(6,6,[1 2 7 8]);     hold on;	axis equal;
        imshow(cData.CH2.IC , [min(min(cData.CH2.IC)), max(max(cData.CH2.IC))])
        title('Raw Signal (false color)', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
        colormap jet ;
        freezeColors ;        %freeze this plot's colormap

%-2--Raw Signal and Mesh
        subplot(6,6,[3 4 9 10]);     hold on;	axis equal;
        imshow(cData.CH2.IC , [min(min(cData.CH2.IC)), max(max(cData.CH2.IC))])
        title('Raw Signal and Mesh', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
        plot( cData.R_mesh(:,1) , cData.R_mesh(:,2),'.', 'Color',[0.9, 0.2, 0]);
        plot( cData.R_mesh(:,3) , cData.R_mesh(:,4),'.', 'Color',[0.9, 0.8, 0]);
        plot( cData.geom.R_axis(:,1) , cData.geom.R_axis(:,2) ...
            , '.-', 'Color', [0 0.6 1]);     

%-3--Bright Field and Cell Contour
        min_x = min( m_Xs );        max_x = max( m_Xs );
        min_y = min( m_Ys );        max_y = max( m_Ys );   
        I_BF = imcrop( Img_BF , [(min_x-extra_border), (min_y-extra_border), ((max_x-min_x)+extra_border*2), ((max_y-min_y)+extra_border*2) ]);
        subplot(6,6,[5 6 11 12]);     hold on;	axis equal;
        imshow(I_BF, [min(min(Img_BF)), max(max(Img_BF))]);  
        title('Bright Field and Cell Contour', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);    
        plot(cData.R_P_model(:,1) , cData.R_P_model(:,2),'*y','MarkerSize' , 2); 


%-4--Cell Membr. Area        (+ contour)
        subplot(6,6,[13 14 19 20]);     hold on;	axis equal; 
        I_Mem =  cData.CH2.IC .* cData.Mask_Memb ;
        imshow(I_Mem, [min(min(cData.CH2.IC)), max(max(cData.CH2.IC))]) ;
        plot(cData.R_model(:,1) , cData.R_model(:,2),'*w','MarkerSize' , 2);        
        colormap jet ;
        freezeColors ;
        title('Cell Membr. Area', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);

%-5--Cell Cytosol Area        (+ contour)
        subplot(6,6,[15 16 21 22]);     hold on;	axis equal; 
        I_Cyt =  cData.CH2.IC .* cData.Mask_mCyto ;
        imshow(I_Cyt, [min(min(cData.CH2.IC)), max(max(cData.CH2.IC))] ) ;  
        plot(cData.R_model(:,1) , cData.R_model(:,2),'*w','MarkerSize' , 2);       
        colormap jet ;
        freezeColors ;
        title('Cell Cytosol Area', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
        
    end
    
%---WAIT to examin    
    set(gcf, 'Color', [1,1,1]);
    pause(APP_opt.plot_pause) ;
   
    end    

end






