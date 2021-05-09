function qXY = evenspaced( XY, Npts )
% CURVSPACE = Create a set of points evenly spaced along a XY curve in 2D 
%   or 3D.
%   CURVSPACE(XY,Npts) generates Npts points that interpolates a curve
%   (represented by a set of points) with an equal spacing. Each
%   row of XY defines a point, which means that XY should be a n x 2
%   (2D) or a n x 3 (3D) matrix.
%
%   (Example)
%   x = -2*pi:0.5:2*pi;
%   y = 10*sin(x);
%   z = linspace(0,10,length(x));
%
%   Npts = 50;
%   points = [x',y',z'];
%   qXY = curvspace(XY, Npts );
%
%   plot3(p(:,1),p(:,2),p(:,3),'*b',q(:,1),q(:,2),q(:,3),'.r');
%
%   legend('Original Points','Interpolated Points');
%
%
% --> NOTE: Adapted from:   curvspace.m 
%             created by:   Yo Fukushima (22 Mar 2005)
%
% -------------------------------------------------------------------------
% Author (this version): Matteo Sangermani
% e-mail: matteo.sangermani@unibas.ch
% Release: 1.0
% Release date: 2019
% -------------------------------------------------------------------------


% Set Initial varialbes
currentpt = XY(1,:);    % current point
indfirst = 2;           % index of the next point in XY
len = size(XY,1);
qXY = currentpt;        % NOTE: Location of first point remain unchanged
ii = 0;

% Define the distance between points in XY
for k0 = 1:len-1
   dist_btw_pts(k0) = distance(XY(k0,:),XY(k0+1,:));
end
totaldist = sum(dist_btw_pts);
intv = totaldist./(Npts-1);     % Define interval between points


for ii = 1:Npts-1   
   newpt = []; distsum = 0;
   ptnow = currentpt;
   kk = 0;
   pttarget = XY(indfirst,:);
   
   % remainder of distance that should be accumulated
   remainder = intv;  
   
   while isempty(newpt)
      % calculate the distance from active point to the closest point in XY
      disttmp = distance(ptnow,pttarget);
      distsum = distsum + disttmp;
      % if distance is enough, generate newpt. else, accumulate
      % distance
      
      if distsum >= intv
         newpt = interpintv(ptnow,pttarget,remainder);
      else
         remainder = remainder - disttmp;
         ptnow = pttarget;
         kk = kk + 1;
         if indfirst+kk > len
            newpt = XY(len,:);
         else
            pttarget = XY(indfirst+kk,:);
         end
      end
      
   end %while
   
   % add to the output points
   qXY = [qXY; newpt];
   
   % update currentpt and indfirst
   currentpt = newpt;
   indfirst = indfirst + kk;
   
end %ii

end %main Func


% ---------------------------------------------------------------------------
% --- Script related Functions ----------------------------------------------  

function dist = distance(A,B)
% DISTANCE Calculate the distance.
% A is either 1 by 2 (2D) or 1 by 3 (3D) vector. 
% B is either n by 2 (2D) or n by 3 (3D) vector, where n is number of points.
% When n > 1, It is evaluated the distance between point A and all points 
% in vector B
%
% --- Adaptation from curvspace.m created by: Yo Fukushima (11 Mar 2005) --- 

if size(A,2) == 2
   dist = sqrt((A(1)-B(:,1)).^2+(A(2)-B(:,2)).^2);
elseif size(A,2) == 3
   dist = sqrt((A(1)-B(:,1)).^2+(A(2)-B(:,2)).^2+(A(3)-B(:,3)).^2);
end

end


function newpt = interpintv( pt1 ,pt2, intv)
% Generate a point between pt1 and pt2 in such a way that the distance
% between pt1 and new point is intv.
% pt1 and pt2 must be either 1 by 2 or 1 by 3 arrays (2D and 3D respectively)
dirvec = pt2 - pt1;
dirvec = dirvec./norm(dirvec);
l = dirvec(1); m = dirvec(2);
newpt = [intv*l+pt1(1),intv*m+pt1(2)];

if length(pt1) == 3
   n = dirvec(3);
   newpt = [newpt,intv*n+pt1(3)];
end

end



