function [mu, sigma, Distr_Vals] = Hist_Values_cloneList(btn, varargin)
%
%Hist_Values_cloneList = create a histogram of the values that would be 
%   plotted in a lineage tree, with the current options and parameters.
%
% The aim is to display in a pop-up figure the range of values, the average 
% value and standard deviation to the user. This will help the user to
% evaluate which min:max range to choose for the colormap in the plot
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


global APP_opt ;

% Assign the variable List according to inputs
if isempty(varargin)
    % Load clone_List, or if required, convert it to a track_List
    load([APP_opt.t3_path_cloneList, APP_opt.t3_file_cloneList]);
    List = clone_List;    
else
    List = varargin{1};
end

% Establish how many channels are used according to .t3_choose_ChannelMode
if APP_opt.t3_choose_ChannelMode == 1
    ChNum = 1;
elseif APP_opt.t3_choose_ChannelMode == 2
    ChNum = 2 ;
elseif APP_opt.t3_choose_ChannelMode >= 3
    ChNum = [1, 2] ;    
end

% Extract the signal value for the correct subcellular compartment from all
% the cells in cList. From it we will create a histogram
for kk = ChNum
    switch APP_opt.t3_PlotOpt_B
        case 1
            Distr_Vals{kk} = EvalDistr_1_SingleVals(   List, kk, APP_opt.t3_PlotOpt_C, APP_opt.algorithm );
        case 2
            Distr_Vals{kk} = EvalDistr_2_BiVals(       List, kk, APP_opt.t3_PlotOpt_C, APP_opt.algorithm );
        case 3        
            Distr_Vals{kk} = EvalDistr_3_SingleRatio(  List, kk, APP_opt.t3_PlotOpt_C, APP_opt.algorithm );
        case 4
            if APP_opt.algorithm == 3 & ...
                ( (APP_opt.t3_PlotOpt_C <= 3  &&  isfield( List{1,1}{1,1}.Fluor_Chan ,'SegmentSig')) | ...
                  (APP_opt.t3_PlotOpt_C >= 4  &&  isfield( List{1,1}{1,1}.Fluor_Chan ,'AxSig')) ) 
                Distr_Vals{kk} = EvalDistr_4_Segmentation( List, kk, APP_opt.t3_PlotOpt_C );
            end
    end
end


% If we have option "both 1 and 2", then we need to plot 2 histograms
if APP_opt.t3_choose_ChannelMode == 3
    nPlot = [1,2];
% If we want single channels or ratios, we plot only 1 histogram. 
elseif APP_opt.t3_choose_ChannelMode == 2
    nPlot = 2;
else
    nPlot = 1;
end

% If the option ratio is selected, calculate the ratio between the two channels
if APP_opt.t3_choose_ChannelMode == 4
    Distr_Vals = {Distr_Vals{1} ./ Distr_Vals{2}} ;
elseif APP_opt.t3_choose_ChannelMode == 5
    Distr_Vals = {Distr_Vals{2} ./ Distr_Vals{1}} ;
end


for kk = nPlot
    % If button is pressed plot distribution of value and make a Gaussian fit.
    % Moreover report the values for mu{kk} +- sigma{kk}
    if btn == 1    
        hf1 = figure(1111*kk);
        hf1.Position(1) = 50 +((kk-1)*hf1.Position(3));
        h1 = subplot(1,10,[1:8]);  % --- Histogram

        hnd = histogram(Distr_Vals{kk}, 20,  'Normalization', 'probability' );

        hnd.FaceColor = [ 1 .7 .3];   hnd.EdgeColor = [.2 .2 .2];    hnd.LineWidth = 1.00;    
        ax = gca;              
        ax.FontSize = 12;             ax.LineWidth = 1.25;           ax.Box = 'on';
        ax.XColor = [.2 .2 .2];       ax.YColor = [.2 .2 .2];
        xlabel('Values',  'Color',[.2 .2 .2], 'FontSize', 14 );
        ylabel('Probability',  'Color',[.2 .2 .2], 'FontSize', 14 );    

        yy = hnd.BinCounts;
        xx = (hnd.BinEdges(1:end-1) + hnd.BinEdges(2:end)) ./ 2;
               
        f = fit(xx.',yy.','gauss2');    
        mu{kk} = f.b1;    sigma{kk} = f.c1;

        ymax = get(gca,'ylim');
        hold on
        plot([(mu{kk}-sigma{kk}), (mu{kk}-sigma{kk})], ymax, '--', 'LineWidth', 2, 'Color', [.9 .3 .0]);
        plot([(mu{kk}+sigma{kk}), (mu{kk}+sigma{kk})], ymax, '--', 'LineWidth', 2, 'Color', [.9 .3 .0]);
        hold off


        h2 = subplot(1,10,[9,10]);  % --- Text
        cla(h2);
        hold on;	axis equal;
        set(gca, 'Visible', 'off') 
        text( 0, 2.20, 'Average Value:' , 'FontSize', 12 , 'Color', [.4 .2 .0] );
        text( 0, 2.00, num2str(mu{kk})      , 'FontSize', 12 , 'Color', [.4 .2 .0] );
        text( 0, 1.60, 'St. Dev.:'      , 'FontSize', 12 , 'Color', [.4 .2 .0] );
        text( 0, 1.40, num2str(sigma{kk})   , 'FontSize', 12 , 'Color', [.4 .2 .0] );

        text( 0, 0.80, '\mu{kk} - \delta ='  , 'FontSize', 12 , 'Color', [.9 .3 .0] );
        text( 0, 0.55, num2str(mu{kk}-sigma{kk}) , 'FontSize', 12 , 'Color', [.9 .3 .0] );
        text( 0, 0.25, '\mu{kk} + \delta ='  , 'FontSize', 12 , 'Color', [.9 .3 .0] );
        text( 0, 0.00, num2str(mu{kk}+sigma{kk}) , 'FontSize', 12 , 'Color', [.9 .3 .0] );

        text( 0, -0.50, '\mu{kk} - 2*\delta ='  , 'FontSize', 12 , 'Color', [.8 .4 .3] );
        text( 0, -0.75, num2str(mu{kk}-2*sigma{kk}) , 'FontSize', 12 , 'Color', [.8 .4 .3] );
        text( 0, -1.05, '\mu{kk} + 2*\delta ='  , 'FontSize', 12 , 'Color', [.8 .4 .3] );
        text( 0, -1.30, num2str(mu{kk}+2*sigma{kk}) , 'FontSize', 12 , 'Color', [.8 .4 .3] ); 

    else
        hf1 = figure(1111*kk);
        hf1.Position(1) = 50 +((kk-1)*hf1.Position(3));
        hnd = histogram(Distr_Vals{kk}, 20,  'Normalization', 'probability' ,'Visible', 'off');      
        yy = hnd.BinCounts;
        xx = (hnd.BinEdges(1:end-1) + hnd.BinEdges(2:end)) ./ 2;        
        close(hf1);
        
        f = fit(xx.',yy.','gauss2');    
        mu{kk} = f.b1;    sigma{kk} = f.c1;
        
        % If only channel 2 is selected, remove first empty variable
        if APP_opt.t3_choose_ChannelMode == 2
            Distr_Vals = {Distr_Vals{2}} ;
            mu = {mu{2}} ;
            sigma = {sigma{2}} ;
        end
        
    end
    
end %kk
    

end %Main func


% *********************************************************************************************************
% -----> SCRIP-RELATED FUNCTIONS --------------------------------------------------------------------------
% *********************************************************************************************************

% --1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1

function [Distr_Vals] = EvalDistr_1_SingleVals(clone_List, ChN, opt_C, algorithm)
Value = [];
Distr_Vals = [];        % Vars to store all the values found in clone_List
for cc = 1 : size(clone_List, 2)
  for ff = 1 : size(clone_List{1,cc}, 2)     
    
    % Evaluate the Value to apply inside the square    
    [Value] = ValueFinder_B1(clone_List{1,cc}{1,ff}, ChN, opt_C, algorithm);   
    % Store all values found
    Distr_Vals = [Distr_Vals , Value] ;
    
  end %-ff
end %-cc
end



% --2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2

function [Distr_Vals] = EvalDistr_2_BiVals(clone_List, ChN, opt_C, algorithm)
Val_1 = [];             % upper value (in y-axis of lineage plot)
Val_2 = [];             % lower value (in y-axis of lineage plot)
Distr_Vals = [];        % Vars to store all the values found in clone_List
for cc = 1 : size(clone_List, 2)
  for ff = 1 : size(clone_List{1,cc}, 2)     
    
    % Evaluate the Value to apply inside the square    
    [ Val_1, Val_2 ] = ValueFinder_B2(clone_List{1,cc}{1,ff}, ChN, opt_C, algorithm);    
    % Store all values found
    Distr_Vals = [Distr_Vals , Val_1, Val_2] ;
    
  end %-ff
end %-cc
end


% --3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3

function [Distr_Vals] = EvalDistr_3_SingleRatio(clone_List, ChN, opt_C, algorithm)
Val_A = [];             % Numerator (of the Ratio)
Val_B = [];             % Denominator (of the Ratio)
Distr_Vals = [];        % Vars to store all the values found in clone_List
for cc = 1 : size(clone_List, 2)
  for ff = 1 : size(clone_List{1,cc}, 2)     

    % Evaluate the Value to apply inside the square
    [Val_A , Val_B] = ValueFinder_B3(clone_List{1,cc}{1,ff}, ChN, opt_C, algorithm);       
    % Calculate the ratio and store all values found
    vRatio = (Val_A / Val_B) ;    
    Distr_Vals = [Distr_Vals , vRatio] ;
    
  end %ff
end %cc
end


% --4--4--4--4--4--4--4--4--4--4--4--4--4--4--4--4--4--4--4--4--4--4--4--4

function [Distr_Vals] = EvalDistr_4_Segmentation(clone_List, ChN, opt_C)
% Can only run for algorithm 3 (AIS), if we performed axial of segmentation
% analysis
Val = [];              % Denominator (of the Ratio)
Distr_Vals = [];       % Vars to store all the values found in clone_List
for cc = 1 : size(clone_List, 2)
  for ff = 1 : size(clone_List{1,cc}, 2)     

    % Evaluate the Value to apply inside the square
    [Val] = ValueFinder_B4(clone_List{1,cc}{1,ff}, ChN, opt_C);       
    % Store all values found   
    Distr_Vals = [Distr_Vals , Val] ;
    
  end %ff
end %cc
end

