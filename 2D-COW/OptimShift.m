function [retImage,retShift] = OptimShift(inTarget, inRef, xLimit, yLimit)
%function retImage = OptimShift(inTarget, inRef, xLimit, yLimit)
%
% Calculate the optimal shift to the target image in order to match the
% reference image.
%
% Input:
%   inTarget -- the 2-dimensional image needs to be shifted;
%   inRef -- the 2-dimensional reference image, rowI x colI;
%   xLimit -- [xL, xU], define the limit of shift along x
%   yLimit -- [yL, yU], define the limit of shift along y
%
% Output:
%   retImage -- the shifted image matching inRef;
%   retShift -- [xShift, yShift], shifting coefficients, which are applied
%               to inTarget to get retImage.
%
% Note: inTarget and inRef are assumed to be sampled at same time points.
%       The sampling interval is taken as the time unit.
%
% Dabao Zhang     Nov. 16, 2006
%
% Revised by:
%
% Copyright (c) 2007 by Zhang & Zhang.

if( nargin<2 ) error('Need more input arguments...'); end
[rowI, colI] = size(inRef);

if( (nargin<4)|isempty(yLimit) ) yLimit = [-1 1]*round(rowI/12); end
if( (nargin<3)|isempty(xLimit) ) xLimit = [-1 1]*round(colI/15); end

tmpCC = zeros(yLimit(2)-yLimit(1)+1,xLimit(2)-xLimit(1)+1);
j=0;
nL = rowI*colI;
for x=xLimit(1):xLimit(2)
    j = j+1; k = 0;
    disp(['column ' num2str(j) '...']);
    for y=yLimit(1):yLimit(2)
        k = k+1;
        
        tmpT = zeros(rowI,colI);
        if( (x>=0)&(y>=0) )
            tmpT((1+y):rowI,(1+x):colI) = inTarget(1:(rowI-y),1:(colI-x));
        elseif( (x>=0)&(y<0) )
            tmpT(1:(rowI+y),(1+x):colI) = inTarget((1-y):rowI,1:(colI-x));
        elseif( (x<0)&(y>=0) )
            tmpT((1+y):rowI,1:(colI+x)) = inTarget(1:(rowI-y),(1-x):colI);
        else % x<0, y<0
            tmpT(1:(rowI+y),1:(colI+x)) = inTarget((1-y):rowI,(1-x):colI);
        end
        %tmpCC(j,k) = corrcoef();
        tmpCC(j,k) = xcorr(reshape(tmpT,nL,1),reshape(inRef,nL,1), 0, 'coeff');
    end
end

[xID,yID] = find(tmpCC==max(max(tmpCC)));
x = xLimit(1)+xID-1;
y = yLimit(1)+yID-1;
retShift = [x, y];

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
