function [V_Min, V_Max, RGB_range] = Find_ColorValueRange( mu, sigma, Distr_Vals , RGB_range, vals_range)
% Find_ColorValueRange = 
%
%
%
%
%
% -------------------------------------------------------------------------
% Author: Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


V_Min = min(Distr_Vals);     % 1 --> min;
V_Max = max(Distr_Vals);     % 2 --> max;

Value_Range = strsplit( vals_range, ':' );

[nmin, status] = str2num(Value_Range{1}) ;    % if status == 0, conversion was not successful  
if status ~= 0 ;         MinRange = nmin ;
else ;                   MinRange = min(Distr_Vals);
end

[nmax, status] = str2num(Value_Range{2}) ;    % if status == 0, conversion was not successful  
if status ~= 0 ;         MaxRange = nmax ;
else ;                   MaxRange = max(Distr_Vals);
end

Ext_Rng = 0.0 ;          % Extend Range value of X% in both directions
t_range = linspace(MinRange, MaxRange, size(RGB_range,1)+ (size(RGB_range,1)*2*Ext_Rng) ) ;
RGB_range(:,4) = t_range( 1+ (size(RGB_range,1)*Ext_Rng)  :  end- (size(RGB_range,1)*Ext_Rng) )';


end