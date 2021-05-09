function cData = Adj_MarkedPole(cData, I_CH3)
%Adj_MarkedPole = The function receive the mesh of a cell and the
%   fluorescent frame from channel 3. If spots were detected (spotFinder in
%   Oufti), the algortihm find the spot with strongest intensity and to
%   which pole it is nearest. This pole will be tagged as "marked"
%   (cData.polarity = 1) and meshes reoriented correctly (first point/row 
%   in meshes is the marked pole)
%
% NOTE 1 : Normally, if polarity is set to 1, the meshes first row
% corresponds to the "old pole" (otherwise assignment is random/unknown).
% Now we set it as the "Marked" pole, and consider it "Marked" pole
%
% NOTE 2 : the axis will have same orientation as the meshes, because this
% is created at later point of analysis from the adjusted (or not) meshes 
%(row order btw axes and meshes is the same (first axes point == to first
% meshes' point)
%
% Called by t1_MAIN_Analysis.m
%
% INPUTS ------------------------------------------------------------------
% cData = Stores the data concerning a specific cell cc at frame ff
%         (normally in the provided as cellList.meshData{ff}{cc})
% I_CH3 = the image of the Channel 3 
%
% OUTPUT ------------------------------------------------------------------
% cData = update cData with set polarity and correctly oriented meshes
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


global APP_opt ;                        % Variable storing WHISIT options
extra_border = APP_opt.BorderBox ;      % extra border area to add to a cell cropped image

if ~isempty(cData.model)  &&  size(cData.mesh,2) == 4    
    ang = 0:0.1:2*pi;      
    x_cir = 3 *cos(ang);      
    y_cir = 3 *sin(ang);

    % Find Min and Max in X-Y axis for the c-th-cell modelinates
    Xs = [cData.mesh(:,1) ; flipud(cData.mesh(:,3))] ;
    Ys = [cData.mesh(:,2) ; flipud(cData.mesh(:,4))] ;   

    % Crop cell fluorescent signal of with additional (ex)-border
    IC = imcrop( I_CH3 , [(min(Xs)-extra_border), (min(Ys)-extra_border), ((max(Xs)-min(Xs))+extra_border*2), ((max(Ys)-min(Ys))+extra_border*2) ]);

    sf = 1;      % <<<----- correction factor for the shift of coordinates: origin at [0;0] or at [1;1]
    % All distances and coordinates must now be relative to the cropped image
    rel_X = (min(Xs)-sf - extra_border);
    rel_Y = (min(Ys)-sf - extra_border);

    R_mesh(:,1) = double( cData.mesh(:,1) - rel_X );
    R_mesh(:,2) = double( cData.mesh(:,2) - rel_Y );
    R_mesh(:,3) = double( cData.mesh(:,3) - rel_X );
    R_mesh(:,4) = double( cData.mesh(:,4) - rel_Y );
    
    %%%--->>> insert here CHECK PLOT 1 <<<---%%%

    if ~isempty(cData.spots.l)      % if there is no spot, we do nothing
        % Determine the spot to consider as Marker. If there is only one, it is
        % easy, otherwise we take the brightest one
        if length(cData.spots.l) == 1   
            x_spt = cData.spots.x - rel_X; 
            y_spt = cData.spots.y - rel_Y;    

        elseif length(cData.spots.l) >= 2
            for pp = 1 : length(cData.spots.l)
                x0 = cData.spots.x(pp) - rel_X;
                y0 = cData.spots.y(pp) - rel_Y;
                % create temp_mask of pole circle search area
                Mask_PL = roipoly( IC, x0+x_cir , y0+y_cir ) ;
                valInt(pp) = mean(IC(Mask_PL)) ;    
            end
            peak = find( valInt == max(valInt));
            x_spt = cData.spots.x(peak) - rel_X;
            y_spt = cData.spots.y(peak) - rel_Y;
        end

        % Polarity is set (polarity=1) or not set (polarity=0).
        cData.polarity = 1 ; 
        % Normally, if polarity is set, the meshes first row corresponds to the
        % "old pole" (otherwise the vertical orientation of the matrix is random).
        % Now we set it as the "Marked" pole, and still consider it "old pole"
        x_up = R_mesh(1,1) ;
        y_up = R_mesh(1,2) ;
        x_dw = R_mesh(end,2) ;
        y_dw = R_mesh(end,2) ;

        % two distances are calculated [ xy_up~~xy_spt , xy_dw~~xy_spt  ]
        dist = [ sqrt( (x_up - x_spt).^2 + (y_up - y_spt).^2 ), ...
                 sqrt( (x_dw - x_spt).^2 + (y_dw - y_spt).^2 ) ];
        [~, P_Mark] = min(dist);

        if P_Mark == 1
            % ---> do nothing, the meshes are already oriented correctly
        elseif P_Mark == 2
            % ---> Flip up-down meshe coordinates and change columns order
            %      to maintain columns order: [leftX, leftY, rightX, rightY]
            temp = flipud( [cData.mesh(:,3), cData.mesh(:,4), cData.mesh(:,1), cData.mesh(:,2)] );
            cData.mesh = temp ;  
        end
        
    else        %we found no spot, then we do not update meshes
        cData.polarity = 0 ; 
    end

%%%--->>> insert here CHECK PLOT 2 <<<---%%%

end %if ~isempty

end



% % % ---> Check subplot 1: BEFORE pole determination
%
%     imshow(IC, [min(min(IC)), max(max(IC))]);           hold on;
%     plot(R_mesh(:,1), R_mesh(:,2), '.-r');
%     plot(R_mesh(:,3), R_mesh(:,4), '.-b');
%     plot(R_mesh(1,1), R_mesh(1,2), '*g', 'MarkerSize', 10);
%     plot(R_mesh(end,1), R_mesh(end,2), '*y', 'MarkerSize', 10);
%     hold off;


% % % ---> Check subplot 2: AFTER pole determination
%
%     R_mesh(:,1) = double( cData.mesh(:,1) - rel_X );
%     R_mesh(:,2) = double( cData.mesh(:,2) - rel_Y );
%     R_mesh(:,3) = double( cData.mesh(:,3) - rel_X );
%     R_mesh(:,4) = double( cData.mesh(:,4) - rel_Y );
%     subplot(1,2,2)
%     imshow(IC, [min(min(IC)), max(max(IC))]);           hold on;
%     plot(R_mesh(:,1), R_mesh(:,2), '.-r');
%     plot(R_mesh(:,3), R_mesh(:,4), '.-b');
%     plot(R_mesh(1,1), R_mesh(1,2), '*g', 'MarkerSize', 10);
%     hold off;
% 
%     pause(0.5)
%     clf

