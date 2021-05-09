function [ Seg_Len, Max_Len ] = Find_NumberSegments(c_List, chn, compartment)
% Can only run for algorithm 3 (AIS), if we performed axial of segmentation
% analysis

Max_Len = 0;        % Denominator (of the Ratio)
Seg_Len = {};       % Vars to store the total number of segments in each cell of clone_List

for cc = 1 : size(c_List, 2)
  for ff = 1 : size(c_List{1,cc}, 2)     
      
      if compartment == 1
           Seg_Len{1,cc}{1,ff} = length(c_List{1,cc}{1,ff}.Fluor_Chan( chn ).SegmentSig) ;
           
      elseif compartment == 2
           Seg_Len{1,cc}{1,ff} = length(c_List{1,cc}{1,ff}.Fluor_Chan( chn ).AxSig) ;           
      end
      
      if Max_Len < Seg_Len{1,cc}{1,ff}
           Max_Len = Seg_Len{1,cc}{1,ff} ;
      end
  
  end %ff
end %cc
end