function InF_Display_Cell_PL( cData, Img_BF )
%Display_Cell_PL = Display the results from analysis via algorithm PL
% (Polar Signal).
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
    
    % Renames mesh coordinates of the cell
    m_Xs = [cData.mesh(:,1) ; flipud(cData.mesh(:,3))] ;
    m_Ys = [cData.mesh(:,2) ; flipud(cData.mesh(:,4))] ;
    
    hFig2 = figure(2);           % set(hFig2, 'Position', [50 350 600 500]);   
    clf(2);
    set(gca, 'Box', 'off');
    set(gca, 'Color', [1,1,1]);
    set(gcf, 'Color', [1,1,1]);
    
%-1--Cell Raw Signal (false color)
    subplot(5,6,[1 2 7 8]);     hold on;	axis equal;
    imshow(cData.CH1.IC , [min(min(cData.CH1.IC)), max(max(cData.CH1.IC))])
    title('Raw Signal (false color)', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    %define the circle area where the to serch the poles for high signal
    temp = regionprops( bwperim( cData.CH1.M_sc_PL_1 ), 'PixelList');       idx1 = temp.PixelList;
    temp = regionprops( bwperim( cData.CH1.M_sc_PL_2 ), 'PixelList');       idx2 = temp.PixelList;
    temp = regionprops( bwperim( cData.CH1.M_sc_PL_1 ), 'Centroid');        C1 = temp.Centroid;
    temp = regionprops( bwperim( cData.CH1.M_sc_PL_2 ), 'Centroid');        C2 = temp.Centroid;
    plot( C1(1) , C1(2) , '.k');
    plot( idx1(:,1) , idx1(:,2), '*', 'MarkerSize', 2 ,'Color', [0 0 0] );
    plot( C2(1) , C2(2) , '.k');
    plot( idx2(:,1) , idx2(:,2), '*', 'MarkerSize', 2 ,'Color', [0 0 0] );
    colormap jet ;
    freezeColors ;        %freeze this plot's colormap

 %-2--Raw Signal and Mesh
    subplot(5,6,[3 4 9 10]);     hold on;	axis equal;
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
    subplot(5,6,[5 6 11 12]);     hold on;	axis equal;
    imshow(I_BF, [min(min(Img_BF)), max(max(Img_BF))]);  
    title('Bright Field and Cell Contour', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);    
    plot(cData.R_P_model(:,1) , cData.R_P_model(:,2),'.y');

    
%-4--Polar search area
    BW_PL_1 = cData.CH1.Mask_PL_1 ;
    BW_PL_2 = cData.CH1.Mask_PL_2 ;
    sc_1 = cData.CH1.IC .* cData.CH1.M_sc_PL_1 ;
    sc_2 = cData.CH1.IC .* cData.CH1.M_sc_PL_2 ;
    subplot(5,6,[13 14 19 20]);     hold on;	axis equal;
    imshow( sc_1 + sc_2 , [min(min(cData.CH1.IC)), max(max(cData.CH1.IC))] ); 
    title('Pole search area', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    colormap jet ;
    freezeColors ;        %freeze this plot's colormap

%-5--Pole areas
    subplot(5,6,[15 16 21 22]);     hold on;	axis equal;
    imshow( BW_PL_1 + BW_PL_2 );
    title('Pole area', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    temp_PL_1_Px = regionprops((BW_PL_1),'PixelList');
    temp_PL_2_Px = regionprops((BW_PL_2),'PixelList'); 
    plot(temp_PL_1_Px.PixelList(:,1),temp_PL_1_Px.PixelList(:,2),'.r');
    plot(temp_PL_2_Px.PixelList(:,1),temp_PL_2_Px.PixelList(:,2),'.m');
    plot(cData.R_P_model(:,1) , cData.R_P_model(:,2),'.y');

%-6--Cytosol area
    subplot(5,6,[17 18 23 24]);     hold on;	axis equal;
    imshow(cData.CH1.Mask_pCyto);
    title('Cytosol area', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);    
    plot(cData.R_P_model(:,1) , cData.R_P_model(:,2),'.r');
 
%- 8 & 9 -- Values
    subplot(5,6,[27 28]); hold on;	axis equal;  
    set(gca, 'Visible', 'off') 
    text( -1, 0.85, [ 'Pole 1  -  ', num2str(mean(cData.CH1.IC( BW_PL_1 )),'%6.f') ], 'FontSize', 16 , 'Color',[0.8, 0.2, 0] );
    text( -1, 0.45, [ 'Pole 2  -  ', num2str(mean(cData.CH1.IC( BW_PL_2 )),'%6.f') ], 'FontSize', 16 , 'Color',[0.8, 0.2, 0.8] );
    subplot(5,6,[29 30]); hold on;	axis equal;  
    set(gca, 'Visible', 'off')
    text( -1, 0.65, [ 'Cytosol  -  ', num2str(mean(cData.CH1.IC( cData.CH1.Mask_pCyto  )),'%6.f') ], 'FontSize', 16 , 'Color',[0, 0.2, 0.8] );


%% -- FLUORO CHANNEL 2 -----------------------------------------------------------------------------------
%  -------------------------------------------------------------------------------------------------------
    if APP_opt.t1_choose_Chan == 2         % ANALYSE Chan 2 
        hFig3 = figure(3);           % set(hFig3, 'Position', [700 350 600 500]);  
        clf(3);
        set(gca, 'Box', 'off');
        set(gca, 'Color', [1,1,1]);
        set(gcf, 'Color', [1,1,1]);

%-1--Cell Raw Signal (false color)
        subplot(5,6,[1 2 7 8]);     hold on;	axis equal;
        imshow(cData.CH2.IC , [min(min(cData.CH2.IC)), max(max(cData.CH2.IC))])
        title('Raw Signal (false color)', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
        %define the circle area where the to serch the poles for high signal
        temp = regionprops( bwperim( cData.CH2.M_sc_PL_1 ), 'PixelList');       idx1 = temp.PixelList;
        temp = regionprops( bwperim( cData.CH2.M_sc_PL_2 ), 'PixelList');       idx2 = temp.PixelList;
        temp = regionprops( bwperim( cData.CH2.M_sc_PL_1 ), 'Centroid');        C1 = temp.Centroid;
        temp = regionprops( bwperim( cData.CH2.M_sc_PL_2 ), 'Centroid');        C2 = temp.Centroid;
        plot( C1(1) , C1(2) , '.k');
        plot( idx1(:,1) , idx1(:,2), '*', 'MarkerSize', 2 ,'Color', [0 0 0] );
        plot( C2(1) , C2(2) , '.k');
        plot( idx2(:,1) , idx2(:,2), '*', 'MarkerSize', 2 ,'Color', [0 0 0] );
        colormap jet ;
        freezeColors ;        %freeze this plot's colormap


%-2--Raw Signal and Mesh
        subplot(5,6,[3 4 9 10]);     hold on;	axis equal;
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
        subplot(5,6,[5 6 11 12]);     hold on;	axis equal;
        imshow(I_BF, [min(min(Img_BF)), max(max(Img_BF))]);  
        title('Bright Field and Cell Contour', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);    
        plot(cData.R_P_model(:,1) , cData.R_P_model(:,2),'.y');
    
%-4--Polar search area
        BW_PL_1 = cData.CH2.Mask_PL_1 ;
        BW_PL_2 = cData.CH2.Mask_PL_2 ;
        sc_1 = cData.CH2.IC .* cData.CH2.M_sc_PL_1 ;
        sc_2 = cData.CH2.IC .* cData.CH2.M_sc_PL_2 ;
        subplot(5,6,[13 14 19 20]);     hold on;	axis equal;
        imshow( sc_1 + sc_2 , [min(min(cData.CH2.IC)), max(max(cData.CH2.IC))] ); 
        title('Pole search area', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
        colormap jet ;
        freezeColors ;        %freeze this plot's colormap

%-5--Pole areas
        subplot(5,6,[15 16 21 22]);     hold on;	axis equal;
        imshow( BW_PL_1 + BW_PL_2 );
        title('Pole area', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
        temp_PL_1_Px = regionprops((BW_PL_1),'PixelList');
        temp_PL_2_Px = regionprops((BW_PL_2),'PixelList'); 
        plot(temp_PL_1_Px.PixelList(:,1),temp_PL_1_Px.PixelList(:,2),'.r');
        plot(temp_PL_2_Px.PixelList(:,1),temp_PL_2_Px.PixelList(:,2),'.m');
        plot(cData.R_P_model(:,1) , cData.R_P_model(:,2),'.y');

%-6--Cytosol area
        subplot(5,6,[17 18 23 24]);     hold on;	axis equal;
        imshow(cData.CH2.Mask_pCyto);
        title('Cytosol area', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);    
        plot(cData.R_P_model(:,1) , cData.R_P_model(:,2),'.r');

    
%- 8 & 9 -- Values
        subplot(5,6,[27 28]); hold on;	axis equal;  
        set(gca, 'Visible', 'off') 
        text( -1, 0.85, [ 'Pole 1  -  ', num2str(mean(cData.CH2.IC( BW_PL_1 )),'%6.f') ], 'FontSize', 16 , 'Color',[0.8, 0.2, 0] );
        text( -1, 0.45, [ 'Pole 2  -  ', num2str(mean(cData.CH2.IC( BW_PL_2 )),'%6.f') ], 'FontSize', 16 , 'Color',[0.8, 0.2, 0.8] );
        subplot(5,6,[29 30]); hold on;	axis equal;  
        set(gca, 'Visible', 'off')
        text( -1, 0.65, [ 'Cytosol  -  ', num2str(mean(cData.CH2.IC( cData.CH2.Mask_pCyto  )),'%6.f') ], 'FontSize', 16 , 'Color',[0, 0.2, 0.8] );

    end
    
%---WAIT to examin 
    set(gcf, 'Color', [1,1,1]);
    pause(APP_opt.plot_pause) ;
       
    end

end






