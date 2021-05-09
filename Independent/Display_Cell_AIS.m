function Display_Cell_AIS( cellData )
%Display_Cell_AIS = Display the masks and results found by  analysis with 
%   algorithm AIS (Average Signal).
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
d_Prof = APP_opt.AIS_EvalAxisProfile ;
d_BF   = APP_opt.t1_choose_BrightField ;
d_Segm = APP_opt.AIS_EvalSegmentation ;
d_Peri = APP_opt.AIS_EvalPerim;

if ~isempty(cellData.model)  &&  size(cellData.mesh,2) == 4

    Plot_AIS(1, 11, cellData, d_Prof, d_Segm, d_Peri, d_BF );             % Channel 1

    if APP_opt.t1_choose_Chan_2 == 1       
         Plot_AIS(2, 22, cellData, d_Prof, d_Segm, d_Peri, d_BF );        % Channel 2
    end
     if APP_opt.t1_choose_Chan_3 == 1  &&  APP_opt.t1_CH3_Marker ~= 1
         Plot_AIS(3, 33, cellData, d_Prof, d_Segm, d_Peri, d_BF );        % Channel 3
    end 

    %---WAIT to examin figure(s)
    pause(APP_opt.plot_pause) ;

end

end % main fnc



function Plot_AIS( nChan, FigNum, cData, disp_Prof, disp_Segm, disp_Peri, disp_BF  )
      
    hFig = figure(FigNum);           % set(hFig2, 'Position', [50 200 600 500]);
    clf(FigNum);
    set(gca, 'Box', 'off');
    set(gca, 'Color', [1,1,1]);
    set(gcf, 'Color', [1,1,1]);
    
%-1--Cell Raw Signal (false color)
    subplot(2, 3, 1);    
    imshow(cData.Fluor_Chan(nChan).IC , [min(min(cData.Fluor_Chan(nChan).IC)), max(max(cData.Fluor_Chan(nChan).IC))])
    hold on;	axis equal;
    title('Raw Signal', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    colormap jet ;
    freezeColors ;        %freeze this plot's colormap
 
%-2--Raw Signal and Mesh
    subplot(2,3, 2);
    title('Signal and Mesh', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);
    hold on;
    imshow(cData.Fluor_Chan(nChan).IC , [min(min(cData.Fluor_Chan(nChan).IC)), max(max(cData.Fluor_Chan(nChan).IC))])
    hold on;	axis equal;
    plot( cData.R_mesh(:,1) , cData.R_mesh(:,2),'.', 'Color', [1 .9 .2]);
    plot( cData.R_mesh(:,3) , cData.R_mesh(:,4),'.', 'Color', [.9  1 .2]);
    % Plot entire axis
    plot( cData.geom.R_axis(:,1) , cData.geom.R_axis(:,2), '.-', 'Color', [.0 .6 .9] );   
    % Plot first axis point (should always coincide with POLE 1)
    plot(cData.geom.R_axis(1,1) , cData.geom.R_axis(1,2) , 'o' , 'MarkerFaceColor', [.9 .0 .6], 'MarkerSize', 9, 'Color',[1 .2 .6]);
    % Plot the POLE 1 (Old or marked)  
    plot(cData.R_mesh(1,1) , cData.R_mesh(1,2) , '*' , 'LineWidth', 0.9, 'Color',[.6  1 .3]);

%-3--Bright Field and Cell Contour
    if disp_BF == 1
        subplot(2, 3, 3);     
        imshow( cData.Bright_Field.IC , [min(min(cData.Bright_Field.IC)), max(max(cData.Bright_Field.IC))]);  
        hold on;	axis equal;
        title('Bright Field', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);    
        plot(cData.R_P_model(:,1) , cData.R_P_model(:,2),'*y','MarkerSize' , 2); 
    end
  
    
%-4--Cell bodyy axis signal profile
    if disp_Prof == 1
        subplot(2, 3, 4);          
        imshow(cData.Fluor_Chan(nChan).IC, [min(min(cData.Fluor_Chan(nChan).IC)), max(max(cData.Fluor_Chan(nChan).IC))] );
        hold on;        axis equal; 
        title('Profile line', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);    
        % Plot perimeter
        plot([cData.R_mesh(:,1) ; flipud(cData.R_mesh(:,3))] , ...
             [cData.R_mesh(:,2) ; flipud(cData.R_mesh(:,4))] , '.' , 'Color',[.9 .8 0]);
        % Prot the two profile axis
        plot(cData.geom.R_Naxis_1(:,1), cData.geom.R_Naxis_1(:,2), '-', 'Color', [.2 .6 .9], 'LineWidth', 1.25 ); 
        plot(cData.geom.R_Naxis_2(:,1), cData.geom.R_Naxis_2(:,2), '-', 'Color', [.2 .6 .9], 'LineWidth', 1.25 );
        % Plot the normal profile lines 
        plot( [cData.geom.R_Naxis_1(:,1), cData.geom.R_Naxis_2(:,1)]',  ...
              [cData.geom.R_Naxis_1(:,2), cData.geom.R_Naxis_2(:,2)]', '-', 'Color', [.3 .8 .0], 'LineWidth', 1.25);
    end
    
    
%-5--Cell Segmentation
    if disp_Segm == 1
        subplot(2, 3, 5); 
        imshow(cData.Fluor_Chan(nChan).IC, [min(min(cData.Fluor_Chan(nChan).IC)), max(max(cData.Fluor_Chan(nChan).IC))] );
        hold on;        axis equal; 
        title('Cell Segmentation', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]);    
        % Plot perimeter
        plot([cData.R_mesh(:,1) ; flipud(cData.R_mesh(:,3))] , ...
             [cData.R_mesh(:,2) ; flipud(cData.R_mesh(:,4))] , '.' , 'Color',[.9 .8 0]);
        % Plot the cell body axis
        plot(cData.geom.R_axis(:,1), cData.geom.R_axis(:,2), '-', 'Color', [.0 .6 .9] );    
        % Plot all segmentation lines   
        px = [cData.BodySegmentLines(2:end-1,1) , cData.BodySegmentLines(2:end-1,3)] ;
        py = [cData.BodySegmentLines(2:end-1,2) , cData.BodySegmentLines(2:end-1,4)] ; 
        plot(px' , py' , '.-' , 'Color', [.8 .2 .2], 'LineWidth', 1.25 )
    end
    
    
%-6--Perimeter Signal Profile
    if disp_Peri == 1
        subplot(2, 3, 6); 
        imshow(cData.Fluor_Chan(nChan).IC, [min(min(cData.Fluor_Chan(nChan).IC)), max(max(cData.Fluor_Chan(nChan).IC))] );
        hold on;        axis equal; 
        title('Perimeter profile', 'FontSize', 14 , 'Color',[0.4, 0.4, 0.4]); 
        outXY = cData.R_OutPeri ; 
        inXY  = cData.R_InnPeri ;
        % Plot inner perimeter
        plot(inXY(:,1), inXY(:,2), '.' , 'Color',[.0 .6 .9], 'MarkerSize', 6);
        plot(inXY(1,1), inXY(1,2), '*' , 'Color' ,[.0 .6 .9]);
        % Plot outer perimeter
        plot(outXY(:,1), outXY(:,2) , '.' , 'Color',[.9 .8 .0], 'MarkerSize', 6);
        plot(outXY(1,1), outXY(1,2) , '*' , 'Color',[.9 .8 .0]);
        % Prot Profile lines going from outer-to-inner perimeter
        plot( [inXY(:,1), outXY(:,1)]', [inXY(:,2), outXY(:,2)]' , '-' , 'Color',[1 .3 .3]);
    end    
    
    set(gcf, 'Color', [1,1,1]);

end






