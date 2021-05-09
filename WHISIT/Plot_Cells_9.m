function Plot_Cells_9( f, c )
% Plot the results for analysis by algorith Memb_2_Cyto
% create a 3x3 figure with all main passages of the analysis and help to
% visualize and check the work in progress
%

global GUI_opt ;
global my_DB ;

%% --- Visualize results (for testing) ------------------------------------
%  ------------------------------------------------------------------------
 
    if ~isempty(my_DB(f).cell(c).coord)
    figure(2)

%-1--Cell and Mesh
    subplot(3,3,1) ; 
    imshow(my_DB(f).cell(c).IC , [0, 2000])  ;  
    hold on;  
    for i = 1 : length(my_DB(f).cell(c).mesh)
        plot(my_DB(f).cell(c).R_mesh(i,1) , my_DB(f).cell(c).R_mesh(i,2),'.r');     end
    for i = 1 : length(my_DB(f).cell(c).mesh)
        plot(my_DB(f).cell(c).R_mesh(i,3) , my_DB(f).cell(c).R_mesh(i,4),'.g');     end  
    title('Cell and Mesh');      axis equal ;
    
%-2--Cell and Contour
    subplot(3,3,2) ; 
    % Inew = I.*M;   apply a mask on the image
    I_Cel =  my_DB(f).cell(c).IC .* my_DB(f).cell(c).M_Cel ;
    imshow(I_Cel, [0,2000])  ;  
    hold on;      plot(my_DB(f).cell(c).R_coord(:,1) , my_DB(f).cell(c).R_coord(:,2),'.r');       
    title('Cell and Contour');     axis equal ;
    
%-3--Contour, Central axis, Cell Center)
    subplot(3,3,3) ;
    imshow(my_DB(f).cell(c).M_Cel) ;    hold on; 
    plot(my_DB(f).cell(c).geom.R_axis(:,1), my_DB(f).cell(c).geom.R_axis(:,2), 'g.-')
    plot(my_DB(f).cell(c).R_coord(:,1) , my_DB(f).cell(c).R_coord(:,2),'.y')     
    plot(my_DB(f).cell(c).geom.C(1) , my_DB(f).cell(c).geom.C(2),'c.','MarkerSize' , 25);   
    title('Cell Area, Contour and Center');     axis equal ;
    
%-4--Cell Foci Area          (+ contour)
    subplot(3,3,4) ; 
    imshow(my_DB(f).cell(c).M_Foci  .* my_DB(f).cell(c).M_Cel)  ;  
    hold on;        plot(my_DB(f).cell(c).R_coord(:,1) , my_DB(f).cell(c).R_coord(:,2),'.g');    
  % plot(Foci(1).PixList(:,1) , Foci(1).PixList(:,2),'g.','MarkerSize' , 8);
    if  ~isempty( fieldnames( my_DB(1, f).cell(1, c).Foci ))  
    for i = 1 : length(my_DB(f).cell(c).Foci)
        plot(my_DB(f).cell(c).Foci(i).C(:,1) , my_DB(f).cell(c).Foci(i).C(:,2),'y.','MarkerSize' , 15);
    end
    end
    title('Cell Foci Area');        axis equal ;
 
%-5--Cell Membr. Area        (+ contour)
    subplot(3,3,5) ;
    I_Mem =  my_DB(f).cell(c).IC .* my_DB(f).cell(c).M_Mem ;
    imshow(I_Mem, [0,2000]) ;
    hold on;      plot(my_DB(f).cell(c).R_coord(:,1) , my_DB(f).cell(c).R_coord(:,2),'.g');       
    title('Cell Membr. Area');        axis equal ;

%-6--Cell Cytosol Area       (+ contour)
    subplot(3,3,6) ;
    I_Cyt =  my_DB(f).cell(c).IC .* my_DB(f).cell(c).M_Cyt ;
    imshow(I_Cyt, [0,2000] ) ;  
    hold on;       plot(my_DB(f).cell(c).R_coord(:,1) , my_DB(f).cell(c).R_coord(:,2),'.g');       
    title('Cell Cytosol Area');      axis equal ;
    
%-7---Profile Line Foci   (+ Peri contour)   
    subplot(3,3,7) ;
    imshow(my_DB(f).cell(c).M_Foci  .* my_DB(f).cell(c).M_Cyt) ;
    hold on;       plot(my_DB(f).cell(c).R_P_coord(:,1) , my_DB(f).cell(c).R_P_coord(:,2),'.y');       
    if  ~isempty( fieldnames( my_DB(1, f).cell(1, c).Foci ))
    for k = 1 : length( my_DB(1, f).cell(1, c).Foci )         
        plot( my_DB(f).cell(c).Foci(k).Prof_R_coord(:,1) , my_DB(f).cell(c).Foci(k).Prof_R_coord(:,2) , '-g', 'LineWidth', 2);
    end
    end
    title('Profile Line Foci ');     axis equal ;
    
%-8--Cell Mem Foci Area      (+ Peri contour)
    subplot(3,3,8) ;
    imshow(my_DB(f).cell(c).M_Foci  .* my_DB(f).cell(c).M_Mem) ;
    hold on;       plot(my_DB(f).cell(c).R_P_coord(:,1) , my_DB(f).cell(c).R_P_coord(:,2),'.y');       
    title('Cell Mem Foci Area');     axis equal ;    
    
%-9--Cell Cyto Foci Area     (+ Peri contour)
    subplot(3,3,9) ;
    imshow(my_DB(f).cell(c).M_Foci  .* my_DB(f).cell(c).M_Cyt) ;
    hold on;       plot(my_DB(f).cell(c).R_P_coord(:,1) , my_DB(f).cell(c).R_P_coord(:,2),'.y');       
    title('Cell Cytosol Foci Area');     axis equal ;

    colormap('jet');
    set(gcf, 'Color', [1,1,1]);
    
    pause(GUI_opt.plot_pause) ;
        
    clf(2)
    end


end


%--------------------------------------------------------------------------
%    Copyright (c) 2016; WHISIT, Matteo Sangermani, All rights reserved
%--------------------------------------------------------------------------

