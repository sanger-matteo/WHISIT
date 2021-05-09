function [Xs, Ys, Value] = Find_Value_XYSquare( cList, fXnum, yrow, Dst_pl, hg_stp,  ChNum , SubLine, ChannelMode, PlotOpt_B, PlotOpt_C, algorithm , NumSegments)
% Find_Value_XYSquare = Define the coordinates of the "box" square polygon
% and find the signal value of the subcellular compartment.
% This is achieve in the following steps:
% - 1 - Discriminate which PlotOpt_B was selected in the GUI 
% - 2 - Use Create_1_squareCell to define the coordinates of the box polygon
%       (using "position" fXnum, yrow and line parameters Dst_pl, hg_pl)
% - 3 - Use ValueFinder_BX to extract the signal Value for each channel,
%       using the original image and the cell Mask of the desired
%       subcellular compartment(s)

switch PlotOpt_B                                                                   % --- 1
  case 1      
      [Xs, Ys] = Create_1_SquareBox( SubLine, fXnum, yrow, Dst_pl, hg_stp );       % --- 2
      
      for kk = ChNum            
          Value{kk} = ValueFinder_B1(cList, kk, PlotOpt_C, algorithm);             % --- 3
      end

  case 3      
      [Xs, Ys] = Create_1_SquareBox( SubLine, fXnum, yrow, Dst_pl, hg_stp );       % --- 2       
      
      for kk = ChNum 
          [v_Num, vDen] = ValueFinder_B3( cList, kk, PlotOpt_C, algorithm );       % --- 3
          Value{kk} = v_Num / vDen ;
      end
      
  case 2      
      [Xs, Ys] = Create_2_SquareBox( SubLine, fXnum, yrow, Dst_pl, hg_stp );       % --- 2
              
      for kk = ChNum 
          [Val_UP, Val_DW] = ValueFinder_B2( cList, kk, PlotOpt_C, algorithm );    % --- 3
          Value{kk} = [Val_UP, Val_DW] ;
      end   
      
  case 4
      [Xs, Ys] = Create_Segment_SquareBox( NumSegments, fXnum, yrow, Dst_pl, hg_stp, PlotOpt_C);
      for kk = ChNum 
          Value{kk} = ValueFinder_B4( cList, kk, PlotOpt_C );
      end      
      
end %algorithm



% If a ratio between the channels is selected, we calculate it here
switch ChannelMode
case 2
  Value = {Value{2}} ;
case 4
  Value = {Value{1} ./ Value{2}} ;
case 5
  Value = {Value{2} ./ Value{1}} ;
end

    
end% main Func





