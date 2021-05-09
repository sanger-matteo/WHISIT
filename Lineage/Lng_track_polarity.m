function cList = track_polarity(cList, ff, cc, pcd_track)

    % find the position of cc-th cell-track in new_cList
    now_cc  = pcd_track{2,cc}(ff);           % present (in ffh frame)
    past_cc = pcd_track{2,cc}(ff-1);         % past (in ff-1 frame)
    
    % coord pole 1 and 2 in new frame
    NP1x = cList{1,ff}{1,now_cc}.mesh(1,1);
    NP1y = cList{1,ff}{1,now_cc}.mesh(1,2);
    NP2x = cList{1,ff}{1,now_cc}.mesh(end,1);
    NP2y = cList{1,ff}{1,now_cc}.mesh(end,2);
    % coord old pole in previous frame
    OP_x = cList{1,ff-1}{1,past_cc}.mesh(1,1);
    OP_y = cList{1,ff-1}{1,past_cc}.mesh(1,2);
    % distances between poles in new frame vs old frame
    d1 = sqrt( abs(NP1x - OP_x).^2 + abs(NP1y - OP_y).^2 );
    d2 = sqrt( abs(NP2x - OP_x).^2 + abs(NP2y - OP_y).^2 );
    
    %%%--->>> insert here CHECK PLOT 1 <<<---%%%
    
    if d1 > d2
        cList{1,ff}{1,now_cc}.mesh(:,1) = flipud( cList{1,ff}{1,now_cc}.mesh(:,1)) ; 
        cList{1,ff}{1,now_cc}.mesh(:,2) = flipud( cList{1,ff}{1,now_cc}.mesh(:,2)) ; 
        cList{1,ff}{1,now_cc}.mesh(:,3) = flipud( cList{1,ff}{1,now_cc}.mesh(:,3)) ; 
        cList{1,ff}{1,now_cc}.mesh(:,4) = flipud( cList{1,ff}{1,now_cc}.mesh(:,4)) ; 
        
        %%%--->>> insert here CHECK PLOT 2 <<<---%%%
        
    end 
end


%% ----------- CHECK PLOTs ------------------------------------------------
% insert at different points of the code to test and debug

%     %%%--->>> CHECK PLOT part 1
%     clf;
%     subplot(1,2,1)
%     xs =  [cList{1,ff-1}{1,cc}.mesh(:,1) ; flipud( cList{1,ff-1}{1,cc}.mesh(:,3)) ];           
%     ys =  [cList{1,ff-1}{1,cc}.mesh(:,2) ; flipud( cList{1,ff-1}{1,cc}.mesh(:,4)) ];           
%     fill(xs,ys,[1 0 0]);    alpha(0.2);    hold on;    axis equal;
%     plot(OP_x, OP_y, 'ro', 'Markersize', 10, 'MarkerFaceColor', [0.6 0 0.6]); 
% 
%     xs =  [cList{1,ff}{1,cc}.mesh(:,1) ; flipud( cList{1,ff}{1,cc}.mesh(:,3)) ];
%     ys =  [cList{1,ff}{1,cc}.mesh(:,2) ; flipud( cList{1,ff}{1,cc}.mesh(:,4)) ];         
%     fill(xs,ys,[0 1 0]);    alpha(0.2);    hold on;    axis equal;
%     plot(NP1x, NP1y,'go', 'Markersize', 10, 'MarkerFaceColor', [0 0.4 0.7]); 
%     plot(NP2x, NP2y,'go', 'Markersize', 10, 'MarkerFaceColor', [0 0.6 0]); 


%         %%%--->>> CHECK PLOT part 2
%         subplot(1,2,2)
%         NP1x = cList{1,ff}{1,cc}.mesh(1,1);
%         NP1y = cList{1,ff}{1,cc}.mesh(1,2);
%         NP2x = cList{1,ff}{1,cc}.mesh(end,1);
%         NP2y = cList{1,ff}{1,cc}.mesh(end,2);            
%         fill(xs,ys,[0 1 0]);    alpha(0.2);    hold on;    axis equal;
%         plot(NP1x, NP1y,'go', 'Markersize', 10, 'MarkerFaceColor', [0 0.4 0.7]); 
%         plot(NP2x, NP2y,'go', 'Markersize', 10, 'MarkerFaceColor', [0 0.6 0]); 
%         pause(0.1)
