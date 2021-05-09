function Display_Cell_M2C( cellData )
%Display_Cell_M2C = Display the masks and results found by  analysis with 
%   algorithm M2C (Membrane 2 Cytosol).
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

if ~isempty(cellData.model)  &&  size(cellData.mesh,2) == 4

    Plot_M2C(1, 11, cellData, displayBF );             % Channel 1

    if APP_opt.t1_choose_Chan_2 == 1       
         Plot_M2C(2, 22, cellData, displayBF );        % Channel 2
    end
     if APP_opt.t1_choose_Chan_3 == 1  &&  APP_opt.t1_CH3_Marker ~= 1
         Plot_M2C(3, 33, cellData, displayBF );        % Channel 3
    end 

    %---WAIT to examin figure(s)
    pause(APP_opt.plot_pause) ;

end

end % main fnc



function Plot_M2C( nChan, FigNum, cData, disp_BF )
    
    % Releative mesh coordinates of the cell
    m_Xs = [cData.mesh(:,1) ; flipud(cData.mesh(:,3))] ;
    m_Ys = [cData.mesh(:,2) ; flipud(cData.mesh(:,4))] ;    
    
    hFig = figure(FigNum);           % set(hFig2, 'Position', [50 200 600 500]);
    clf(FigNum);
    set(gca, 'Box', 'off');
    set(gca, 'Color', [1,1,1]);
    set(gcf, 'Color', [1,1,1]);
    
%-1--Cell Raw Signal (false color)
    subplot(6,6,[1 2 7 8]);     
    imshow(cData.Fluor_Chan(nChan).IC , [min(min(cData.Fluor_Chan(nChan).IC)), max(max(cData.Fluor_Chan(nChan).IC))])
    hold on;	axis equal;
    title('Raw Signal (false color)', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    colormap jet ;
    freezeColors ;        %freeze this plot's colormap
 
%-2--Raw Signal and Mesh
    subplot(6,6,[3 4 9 10]);     
    imshow(cData.Fluor_Chan(nChan).IC , [min(min(cData.Fluor_Chan(nChan).IC)), max(max(cData.Fluor_Chan(nChan).IC))])
    hold on;	axis equal;
    title('Raw Signal and Mesh', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    plot( cData.R_mesh(:,1) , cData.R_mesh(:,2),'.', 'Color',[0.9, 0.2, 0]);
    plot( cData.R_mesh(:,3) , cData.R_mesh(:,4),'.', 'Color',[0.9, 0.8, 0]);
    plot( cData.geom.R_axis(:,1) , cData.geom.R_axis(:,2) ...
        , '.-', 'Color', [0 0.6 1]);     

%-3--Bright Field and Cell Contour
    if disp_BF == 1 
        subplot(6,6,[5 6 11 12]);     
        imshow( cData.Bright_Field.IC , [min(min(cData.Bright_Field.IC)), max(max(cData.Bright_Field.IC))]);  
        hold on;	axis equal;
        title('Bright Field and Cell Contour', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);    
        plot(cData.R_P_model(:,1) , cData.R_P_model(:,2),'*y','MarkerSize' , 2); 
    else
        subplot(6,6,[5 6 11 12]);     hold on;	axis equal;
        set(gca, 'Visible', 'off') 
        text( -0.5, 0.5, [ 'Average Signal - ', num2str(mean(cData.Fluor_Chan(nChan).IC( cData.Mask.Cell_body )),'%6.f') ], 'FontSize', 14 , 'Color',[0.8, 0.2, 0] );
    end
 
%-4--Cell Membr. Area        (+ contour) 
    subplot(6,6,[13 14 19 20]);     
    I_Mem =  cData.Fluor_Chan(nChan).IC .* cData.Mask.Memb_All ;
    imshow(I_Mem, [min(min(cData.Fluor_Chan(nChan).IC)), max(max(cData.Fluor_Chan(nChan).IC))]) ;
    hold on;	axis equal; 
    title('Cell Membr. Area', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    plot(cData.R_model(:,1) , cData.R_model(:,2),'*w','MarkerSize' , 2);        
    colormap jet ;
    freezeColors ;
    
%-5--Cell Cytosol Area       (+ contour)
    subplot(6,6,[15 16 21 22]);     
    I_Cyt =  cData.Fluor_Chan(nChan).IC .* cData.Mask.Cytosol ;
    imshow(I_Cyt, [min(min(cData.Fluor_Chan(nChan).IC)), max(max(cData.Fluor_Chan(nChan).IC))] ) ;  
    hold on;	axis equal; 
    title('Cell Cytosol Area', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    plot(cData.R_model(:,1) , cData.R_model(:,2),'*w','MarkerSize' , 2);       
    colormap jet ;
    freezeColors ;
    

    
    set(gcf, 'Color', [1,1,1]);   

end






