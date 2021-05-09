function Signal_Analysis( f, c )

% --- for WHISIT (v_1.0)
% Script encoding the algorithm Memb_2_Cyto
% Briefly, the algorithm is addressing the issue of the ratio of signal 
% between cytosol and membrane and how the signal profile along the cell width.
%
% Masks are generated to determine cell areas of cytosol of membrane and
% then the analysis of the signal in these two is performed and compared
%

%% --- Initialize ---------------------------------------------------------
%--------------------------------------------------------------------------
global GUI_opt;
global my_DB;

Thres_AreaFoci = GUI_opt.Thres_F_size ;     % [pixel], min Foci size to consider
Thres_SigFoci = GUI_opt.Thres_Signal;       % Threshold intensity signal that identify foci = '' ;
extra_border = GUI_opt.border_size ;    % extra border area to add to a cell cropped image
Thr_mem = GUI_opt.Thres_membrane ;      % number of circle to repeat to extend the membrane area: roughly each 
                                        % cicle add a pixel width from perimeter inward to the cell to membrane area
                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                        
%% --- Body of the Function -----------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if f <10                            null = '000';
elseif f >= 10 && f < 100           null = '00';
elseif f >= 100 && f < 1000         null = '0';
else                                null = '';
end

I = double(imread([ GUI_opt.path_FL, '/FL_' , null num2str(f), '.tif'])) ;    % handle to f-th frame
if ~ isempty(my_DB(f).cell(c).coord)            
    % Find Min and Max in X-Y axis for the c-th-cell coordinates
    min_x = min( my_DB(f).cell(c).coord(:, 1) );
    min_y = min( my_DB(f).cell(c).coord(:, 2) );
    max_x = max( my_DB(f).cell(c).coord(:, 1) );
    max_y = max( my_DB(f).cell(c).coord(:, 2) );     
    
    % Crop cell fluorescent signal of with additional (ex)-border
    my_DB(f).cell(c).IC = imcrop( I , [(min_x-extra_border), (min_y-extra_border), ((max_x-min_x)+extra_border*2), ((max_y-min_y)+extra_border*2) ]);
    
    sf = 1;      % <<<----- correction factor for the shift between BF and Fluor. due
                 % to difference in coordinates system: origin at [0;0] or at [1;1]
    % coordinates Relative to cropped image of the cell (.IC)
    my_DB(f).cell(c).R_mesh(:,1) = my_DB(f).cell(c).mesh(:,1) - (min_x-sf - extra_border) ;
    my_DB(f).cell(c).R_mesh(:,2) = my_DB(f).cell(c).mesh(:,2) - (min_y-sf - extra_border) ;
    my_DB(f).cell(c).R_mesh(:,3) = my_DB(f).cell(c).mesh(:,3) - (min_x-sf - extra_border) ;
    my_DB(f).cell(c).R_mesh(:,4) = my_DB(f).cell(c).mesh(:,4) - (min_y-sf - extra_border) ;
    
    my_DB(f).cell(c).R_coord(:,1) = my_DB(f).cell(c).coord(:,1) - (min_x-sf - extra_border) ;
    my_DB(f).cell(c).R_coord(:,2) = my_DB(f).cell(c).coord(:,2) - (min_y-sf - extra_border) ;
    
    my_DB(f).cell(c).geom.R_axis(:,1) = my_DB(f).cell(c).geom.axis(:,1) - (min_x-sf - extra_border) ;
    my_DB(f).cell(c).geom.R_axis(:,2) = my_DB(f).cell(c).geom.axis(:,2) - (min_y-sf - extra_border) ;

%% --- Extract cell's sectional signal ------------------------------------
%--------------------------------------------------------------------------
    sig = {} ;
    for p = 2 : length(my_DB(f).cell(c).mesh)      
        xs = [my_DB(f).cell(c).R_mesh(p,1) , my_DB(f).cell(c).R_mesh(p,3)] ;
        ys = [my_DB(f).cell(c).R_mesh(p,2) , my_DB(f).cell(c).R_mesh(p,4)] ;
        sig{p} = improfile(my_DB(f).cell(c).IC, xs ,ys) ;
    end    

    % Find longest section line that cut a cell
    Max = 0 ;
    for i = 1:length(sig)
        if Max < length(sig{i})
            Max = length(sig{i});   end
    end

    % Create a matrix-image containing all the section line s{p} intensity 
    % values, adding zeros where needed in order to keep all aligned around
    % center (which represent the cell main axis)
    temp_MSig = [];
    for i = 1 : length(sig)
        plus = [];
        Z = Max - length(sig{i});
        if mod(Z,2) == 0
            plus = zeros(1, (Z/2));
            temp_MSig(i,:) = [plus, (sig{i})' , plus];
        else
            plus = zeros(1, fix(Z/2));      % fix(), division returning integer
            temp_MSig(i,:) = [plus, (sig{i})' , plus, 0];
        end
    end
    my_DB(f).cell(c).Sig_Vec = sig ;
    my_DB(f).cell(c).Sig_Mat = temp_MSig ; 
    
    
%% --- Create the Masks that define cytosol and membrane area -------------
%--------------------------------------------------------------------------
    % For simplicity in coding we give shorter name to coordinates
    px = my_DB(f).cell(c).R_coord(:,1);
    py = my_DB(f).cell(c).R_coord(:,2);
    
    % ---> 1 Cell Mask
    BW = roipoly( my_DB(f).cell(c).IC, px, py );
    [R, C] = size(BW) ;
    
    %---------R_P_coord------Relative_Perimeter_coord----------------------
    % Using round(.mesh) to define .coord (see Retrive_BF_mesh.m), 
    % does not give a continuos line of pixel around the cell. But
    % with the mask we can can extract the perimeter and re-define cell 
    % coordinates so that they describe a continuous line.
    %
    tempPx = regionprops(bwperim(BW),'PixelList');
    my_DB(f).cell(c).R_P_coord = tempPx.PixelList ;   
    %
    % NOTE: coordinates are ordered in respect to x-axis. Is not a problem
    % for .coord, we never need them ordered (compared to .mesh)
    %----------------------------------------------------------------------
        
    % ---> 2 Temporary Masks
    % define contour perimeter
    t_BWm = bwperim(BW);
    for i = 1 : Thr_mem
        t_BWc = imfill(t_BWm,'holes') - t_BWm;  % create temporay cytosol mask:
        in_per =  bwperim(t_BWc) ;              % subtract 'holes' to perimeter to have new Cytosol
        t_BWm = t_BWm + in_per ;                % inner hole perimiter ...  
    end                                         % ... add it to Membrane area
    
    % ---> 3 Final Cytosol and Membrane Masks and Save them in my_DB     
    %BWm = t_BWm + in_per;
    BWm = logical( t_BWm );
    BWc = logical( imfill(t_BWm,'holes') - BWm );      
    my_DB(f).cell(c).M_Cel = BW;
    my_DB(f).cell(c).M_Cyt = BWc;
    my_DB(f).cell(c).M_Mem = BWm;

    
    
%% --- FOCI: detection and measurements -----------------------------------
%--------------------------------------------------------------------------
    temp = regionprops(my_DB(f).cell(c).M_Cel,'centroid');
    center = cat(1, temp.Centroid);         % Cell Centroid coorinates

    BW_Foci = (my_DB(f).cell(c).IC .*  my_DB(f).cell(c).M_Cel) > Thres_SigFoci ;
    temp = regionprops(BW_Foci,'centroid');      t_Foci.Centroid = cat(1, temp.Centroid);
    temp = regionprops(BW_Foci,'Area');          t_Foci.Area = cat(1, temp.Area); 
    temp = regionprops(BW_Foci,'Orientation');   t_Foci.Orient = cat(1, temp.Orientation) .* (pi/180) ;
    temp = regionprops(BW_Foci,'MinorAxisLength');  t_Foci.B_ax = cat(1, temp.MinorAxisLength) ;
    temp = regionprops(BW_Foci,'MajorAxisLength');  t_Foci.A_ax = cat(1, temp.MajorAxisLength) ;
    temp = regionprops(BW_Foci,'Eccentricity');  t_Foci.Epsilon = cat(1, temp.Eccentricity) ;
    t_tempPx = regionprops(BW_Foci,'PixelList');
    
    Foci = struct( );
    j = 1 ;
    for i = 1 : length(t_Foci.Area)
        if t_Foci.Area(i) >= Thres_AreaFoci
            Foci(j).C = t_Foci.Centroid(i,:) ;      % Centroid coorinates
            Foci(j).Area = t_Foci.Area(i) ;         % Area
            Foci(j).Beta = t_Foci.Orient(i) ;       % Angle orientation (rleative to x-axis) 
            Foci(j).Epsilon = t_Foci.Epsilon(i) ;   % Eccentricity of ellipse, (0:1)
            Foci(j).A_ax =  t_Foci.A_ax(i) ;        % major axis of fitted ellipse
            Foci(j).B_ax =  t_Foci.B_ax(i) ;        % minor axis of fitted ellipse
            % List of coordinates of all pixels in a given foci
            Foci(j).PixList = t_tempPx(i).PixelList;
            j = j+1 ;
        end
    end
    my_DB(f).cell(c).M_Foci = BW_Foci;
    my_DB(f).cell(c).Foci = Foci;
    my_DB(f).cell(c).geom.C = center;
    
%% --- Profile Linepassing through Foci -----------------------------------
%  ------------------------------------------------------------------------      
    % We have center point for each focii. What we do is to calculate the
    % mesh line that pass closest and do store the improfile of such line,
    % which is more significal than whole cell.
    if ~ isempty( fieldnames(my_DB(1,f).cell(1,c).Foci))      % if there are no fields, there are no Foci
    for jj = 1 : length(my_DB(f).cell(c).Foci)
        store_x = [];        store_y = [];
        store_d = 1000;      d = 1000 ;
        sig = [];
        Cx = my_DB(f).cell(c).Foci(jj).C(1) ;
        Cy = my_DB(f).cell(c).Foci(jj).C(2) ;
        for kk = 1 : length(my_DB(f).cell(c).R_mesh)
            x1 = my_DB(f).cell(c).R_mesh(kk,1) ;            
            y1 = my_DB(f).cell(c).R_mesh(kk,2) ; 
            x2 = my_DB(f).cell(c).R_mesh(kk,3) ;
            y2 = my_DB(f).cell(c).R_mesh(kk,4) ;
            % Distance from a point to a line using 3-points coordinates
            d = abs((y2-y1)*Cx - (x2-x1)*Cy + x2*y1 - y2*x1) / sqrt((y2-y1)^2 + (x2-x1)^2);
            if d < store_d
                store_x = [x1, x2] ;
                store_y = [y1, y2] ;
                store_d = d ;
            end
         end
         % Save Foci profile line coordinates and improfile
         my_DB(f).cell(c).Foci(jj).Prof_R_coord = [store_x' , store_y'] ;
         sig = improfile(my_DB(f).cell(c).IC, store_x ,store_y ) ;
         my_DB(f).cell(c).Foci(jj).Prof_line = sig' ;        
    end
    end

end

end

%--------------------------------------------------------------------------
%    Copyright (c) 2016; WHISIT, Matteo Sangermani, All rights reserved
%--------------------------------------------------------------------------


