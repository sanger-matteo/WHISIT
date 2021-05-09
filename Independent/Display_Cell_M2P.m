function Display_Cell_M2P( cellData )
%Display_Cell_PL = Display the masks and results found by  analysis with 
%   algorithm M2P (Membrane 2 Polar).
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

if ~isempty(cellData.model)  &&  size(cellData.mesh,2) == 4

    Plot_M2P(1, 11, cellData );             % Channel 1

    if APP_opt.t1_choose_Chan_2 == 1       
         Plot_M2P(2, 22, cellData );        % Channel 2
    end
     if APP_opt.t1_choose_Chan_3 == 1  &&  APP_opt.t1_CH3_Marker ~= 1
         Plot_M2P(3, 33, cellData );        % Channel 3
    end 

    %---WAIT to examin figure(s)
    pause(APP_opt.plot_pause) ;

end

end % main fnc



function Plot_M2P( nChan, FigNum, cData )

  % Releative mesh coordinates of the cell
    m_Xs = [cData.mesh(:,1) ; flipud(cData.mesh(:,3))] ;
    m_Ys = [cData.mesh(:,2) ; flipud(cData.mesh(:,4))] ;    
    
    hFig = figure(FigNum);           % set(hFig2, 'Position', [50 200 600 500]);
    clf(FigNum);
    set(gca, 'Box', 'off');
    set(gca, 'Color', [1,1,1]);
    set(gcf, 'Color', [1,1,1]);
    
%-2-- Cell Raw Signal CH1 (false color)
    subplot(3,2, 2);   
    imshow(cData.Fluor_Chan(nChan).IC , [min(min(cData.Fluor_Chan(nChan).IC)), ...
                                         max(max(cData.Fluor_Chan(nChan).IC))])
    hold on;	axis equal;
    title('Raw Signal CH1', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    colormap jet ;
    freezeColors ;        %freeze this plot's colormap
    hold off;
        
%-4-- Cytosol Area
    subplot(3,2, 4);     
    imshow( cData.Mask.Cytosol )
    hold on;	axis equal; 
    title('Unified Membrane Area', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    hold off; 
    
%-1-- Unified Membrane Area
    subplot(3,2, 1);     
    imshow( cData.Mask.Memb_All )
    hold on;	axis equal; 
    title('Unified Membrane Area', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    hold off; 
 
%-3--Pole Membr. Areaw        (+ contour) 
    subplot(3,2, 3);      
    imshow( cData.Mask.MembPole_1 + cData.Mask.MembPole_2 );    
    hold on;	axis equal;   
    title('Pole Memb Areas', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    plot(cData.R_P_model(:,1) , cData.R_P_model(:,2),'.y');
    temp_PL_1_Px = regionprops((cData.Mask.MembPole_1),'PixelList');
    temp_PL_2_Px = regionprops((cData.Mask.MembPole_2),'PixelList'); 
    plot(temp_PL_1_Px.PixelList(:,1),temp_PL_1_Px.PixelList(:,2),'.r');
    plot(temp_PL_2_Px.PixelList(:,1),temp_PL_2_Px.PixelList(:,2),'.c');

%-5--Lateral Membr. Areas     (+ contour)
    subplot(3,2, 5);      
    imshow( cData.Mask.MembLateral_1 + cData.Mask.MembLateral_2 );
    hold on;	axis equal;
    title('Lateral Memb Areas', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    temp_PL_1_Px = regionprops((cData.Mask.MembLateral_1),'PixelList');
    temp_PL_2_Px = regionprops((cData.Mask.MembLateral_2),'PixelList'); 
    plot(cData.R_P_model(:,1) , cData.R_P_model(:,2),'.y');
    plot(temp_PL_1_Px.PixelList(:,1),temp_PL_1_Px.PixelList(:,2),'.r');
    plot(temp_PL_2_Px.PixelList(:,1),temp_PL_2_Px.PixelList(:,2),'.c');
    
     
%-6-- CH1 Values
    subplot(3,2, 6);        hold on;	axis equal;
    set(gca, 'Visible', 'off') 
    text( -0.5, 1.00, [ 'CH1 Pole 1  -  ', num2str(mean(cData.Fluor_Chan(nChan).IC( cData.Mask.MembPole_1 )),'%6.f') ], 'FontSize', 14 , 'Color',[0.8, 0.2, 0] );
    text( -0.5, 0.66, [ 'CH1 Pole 2  -  ', num2str(mean(cData.Fluor_Chan(nChan).IC( cData.Mask.MembPole_2 )),'%6.f') ], 'FontSize', 14 , 'Color',[0.0, 0.6, 0.8] );
    text( -0.5, 0.33, [ 'CH1 Side 1  -  ', num2str(mean(cData.Fluor_Chan(nChan).IC( cData.Mask.MembLateral_1 )),'%6.f') ], 'FontSize', 14 , 'Color',[0.8, 0.2, 0] );
    text( -0.5, 0.00, [ 'CH1 Side 2  -  ', num2str(mean(cData.Fluor_Chan(nChan).IC( cData.Mask.MembLateral_2 )),'%6.f') ], 'FontSize', 14 , 'Color',[0.0, 0.6, 0.8] );
 

%---WAIT to examin 
    set(gcf, 'Color', [1,1,1]);       

end






