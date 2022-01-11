function retImage = ShiftImage(inTarget,inShift)
%function retImage = ShiftImage(inTarget,inShift)
%
% Apply the shift parameters to a given image.
%
% Input:
%   inTarget -- the 2-dimensional image needs to be shifted;
%   inShift -- [x,y], the shift parameters;
%
% Output:
%   retImage -- the shifted image;
%
% Dabao Zhang     December 3, 2007
%
% Revised by:
%
% Copyright (c) 2007 by Zhang & Zhang.

if( nargin<2 ) error('Need two input arguments...'); end

x = inShift(1); y = inShift(2);
[rowI, colI] = size(inTarget);
retImage = zeros(rowI,colI);
if( (x>=0)&(y>=0) )
    retImage((1+y):rowI,(1+x):colI) = inTarget(1:(rowI-y),1:(colI-x));
elseif( (x>=0)&(y<0) )
    retImage(1:(rowI+y),(1+x):colI) = inTarget((1-y):rowI,1:(colI-x));
elseif( (x<0)&(y>=0) )
    retImage((1+y):rowI,1:(colI+x)) = inTarget(1:(rowI-y),(1-x):colI);
else % x<0, y<0
    retImage(1:(rowI+y),1:(colI+x)) = inTarget((1-y):rowI,(1-x):colI);
end
