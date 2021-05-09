function Display_Cell_PL( cellData )
%Display_Cell_PL = Display the masks and results found by  analysis with 
%   algorithm PL (Polar Signal).
%   Show all main results for the cc-th cell and help to visualize and 
%   check the work in progress
%
% INPUTS ------------------------------------------------------------------
% cData = the specific cc-th cell to plot ( normally in the provided as
%         cellList.meshData{1,ff}{1,cc})
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------

global APP_opt ;                            % Variable storing WHISIT options
displayBF = APP_opt.t1_choose_BrightField ;
sF = APP_opt.PL_ScaleFact;

if ~isempty(cellData.model)  &&  size(cellData.mesh,2) == 4

    Plot_PL(1, 11, cellData, displayBF, sF);             % Channel 1

    if APP_opt.t1_choose_Chan_2 == 1       
         Plot_PL(2, 22, cellData, displayBF, sF);        % Channel 2
    end
     if APP_opt.t1_choose_Chan_3 == 1  &&  APP_opt.t1_CH3_Marker ~= 1
         Plot_PL(3, 33, cellData, displayBF, sF);        % Channel 3
    end 

    %---WAIT to examin figure(s)
    pause(APP_opt.plot_pause) ;

end

end % main fnc




function Plot_PL( nChan, FigNum, cData, disp_BF , scaleF)

        
    % Calculate the varialbe necessary to plot the search circles, as in
    % Analysis_PL.m function
    cell_w = sqrt( (cData.R_mesh(:,1)-cData.R_mesh(:,3)).^2 + (cData.R_mesh(:,2)-cData.R_mesh(:,4)).^2 );
    % Calculate median and find which points have distance from if  above a specific value
    m_width = median(cell_w);
    % Create a circle coordinates with center (0;0) and radius that is half
    % the median cell width and scaled by scale_factor(sc_fc)
    ang = 0:0.1:2*pi;      
    x_cir = ((m_width/2)* scaleF) *cos(ang);      
    y_cir = ((m_width/2)* scaleF) *sin(ang);
    
    % We place the center of the cell cicle a bit off from the "pole" point
    % at about the second/3 point of the axis away from the pole
    x_Center_1 = cData.geom.R_axis( 3, 1);    % min(idx_Center)
    y_Center_1 = cData.geom.R_axis( 3, 2);
    x_Center_2 = cData.geom.R_axis( end-3, 1);
    y_Center_2 = cData.geom.R_axis( end-3, 2);
    
    % Renames mesh coordinates of the cell
    m_Xs = [cData.mesh(:,1) ; flipud(cData.mesh(:,3))] ;
    m_Ys = [cData.mesh(:,2) ; flipud(cData.mesh(:,4))] ;
    
    hFig = figure(FigNum);           % set(hFig2, 'Position', [50 350 600 500]);   
    clf(FigNum);
    set(gca, 'Box', 'off');
    set(gca, 'Color', [1,1,1]);
    set(gcf, 'Color', [1,1,1]);
    
%-1--Cell Raw Signal and range
    subplot(5,4,[1 2 5 6]);    
    imshow(cData.Fluor_Chan(nChan).IC , [min(min(cData.Fluor_Chan(nChan).IC)), max(max(cData.Fluor_Chan(nChan).IC))])
    hold on;	axis equal;
    title('Raw Signal and Mesh', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    plot( cData.R_mesh(:,1) , cData.R_mesh(:,2),'.', 'Color',[0.9, 0.2, 0]);
    plot( cData.R_mesh(:,3) , cData.R_mesh(:,4),'.', 'Color',[0.9, 0.8, 0]);
    plot( cData.geom.R_axis(:,1) , cData.geom.R_axis(:,2) ...
        , '.-', 'Color', [0 0.6 1]); 
    % display search circle
    plot(x_Center_1+x_cir , y_Center_1+y_cir, '.-', 'MarkerSize', 2, 'Color', [.2 .9 .0]);
    plot(x_Center_1 , y_Center_1 ,'*g', 'Color', [.0 1 .0]);
    plot(x_Center_2+x_cir , y_Center_2+y_cir, '.-', 'MarkerSize', 2, 'Color', [.2 .9 .0]); 
    plot(x_Center_2 , y_Center_2 ,'*g', 'Color', [.0 1 .0]);


 %-2--Polar Search Masks 
    subplot(5,4,[3 4 7 8]);     
    imshow( cData.Fluor_Chan(nChan).Mask.searchAPole_1 + ...
            cData.Fluor_Chan(nChan).Mask.searchAPole_2 );
    hold on;	axis equal;
    title('Search Area at Poles', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);  
    temp_PL_1_Px = regionprops((cData.Fluor_Chan(nChan).Mask.searchAPole_1),'PixelList');
    temp_PL_2_Px = regionprops((cData.Fluor_Chan(nChan).Mask.searchAPole_2),'PixelList'); 
    plot(temp_PL_1_Px.PixelList(:,1),temp_PL_1_Px.PixelList(:,2),'.c');
    plot(temp_PL_2_Px.PixelList(:,1),temp_PL_2_Px.PixelList(:,2),'.m');
    plot( [cData.R_mesh(:,1), flipud(cData.R_mesh(:,3))] , ...
          [cData.R_mesh(:,2), flipud(cData.R_mesh(:,4))] , '-y' );
    
    
%-3--Bright Field and Cell Contour
    if disp_BF == 1 
        subplot(5,4,[9 10 13 14]);        
        imshow( cData.Bright_Field.IC , [min(min(cData.Bright_Field.IC)), max(max(cData.Bright_Field.IC))]);  
        hold on;	axis equal;  
        title('Bright Field', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);    
        plot(cData.R_model(:, 1), cData.R_model(:, 2), '*y','MarkerSize' , 2); 
    else
        subplot(5,4,[9 10 13 14]);     hold on;	axis equal;
        set(gca, 'Visible', 'off') 
        text( -0.25, 0.5, [ 'Average Signal - ', num2str(mean(cData.Fluor_Chan(nChan).IC( cData.Mask.Cell_body )),'%6.f') ], 'FontSize', 14 , 'Color',[0.8, 0.2, 0] );
    end
    
%-4--Pole areas
    subplot(5,4,[11 12 15 16]);     
    BW_PL_1 = cData.Fluor_Chan(nChan).Mask.Pole_1 ;
    BW_PL_2 = cData.Fluor_Chan(nChan).Mask.Pole_2 ;
    imshow( BW_PL_1 + BW_PL_2 );
    hold on;	axis equal;
    title('Pole Masks', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    temp_PL_1_Px = regionprops((BW_PL_1),'PixelList');
    temp_PL_2_Px = regionprops((BW_PL_2),'PixelList'); 
    plot(temp_PL_1_Px.PixelList(:,1),temp_PL_1_Px.PixelList(:,2),'.c');
    plot(temp_PL_2_Px.PixelList(:,1),temp_PL_2_Px.PixelList(:,2),'.m');
    plot(cData.R_P_model(:,1) , cData.R_P_model(:,2),'.y');

 
%- 5 & 5 -- Values
    subplot(5,4,[17 18]);    hold on;	axis equal;  
    set(gca, 'Visible', 'off') 
    text( -1, 0.85, [ 'Pole 1  -  ', num2str(mean(cData.Fluor_Chan(nChan).IC( BW_PL_1 )),'%6.f') ], 'FontSize', 16 , 'Color',[0.2, 0.6, 0.8] );
    text( -1, 0.45, [ 'Pole 2  -  ', num2str(mean(cData.Fluor_Chan(nChan).IC( BW_PL_2 )),'%6.f') ], 'FontSize', 16 , 'Color',[0.8, 0.2, 0.8] );
    subplot(5,4,[19 20]);    hold on;	axis equal;  
    set(gca, 'Visible', 'off')
    text( -1, 0.65, [ 'Cytosol  -  ', num2str(mean(cData.Fluor_Chan(nChan).IC( cData.Fluor_Chan(nChan).Mask.Cytosol  )),'%6.f') ], 'FontSize', 16 , 'Color',[0, 0.2, 0.8] );
   
    set(gcf, 'Color', [1,1,1]);

end






