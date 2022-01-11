function retX = WarpImage(inX,inWarp)
% function retX = WarpImage(inX,inWarp)
%
% The two-dimensional correlation optimized warping (2DCOW) algorithm with
% linear interpolation, described in
%
% Zhang, D., Huang, X., Regnier, F. E. and Zhang, M. (2008) Two-dimensional
% correlation optimized warping algorithm for aligning GCxGC-MS data.
% Analytical Chemistry, 80: 2664-2671.
%
% inX -- the 2-dimensional image needs to be warped, rowX x colX;
% inWarp -- warping coefficients, which are applied to inX to get retX.
%
% retX -- the warped image.
%
% Note: inX and inWarp.xWarp/inWarp.yWarp are assumed to be with the same
%       size.
%
% Dabao Zhang     March 11, 2007
%
% Revised by:
%
% Copyright (c) 2007 by Zhang & Zhang.


if( nargin<2 )
    help WarpImage;
    return;
end

[rowX,colX] = size(inX);
if( (rowX~=size(inWarp.xWarp,1))|(colX~=size(inWarp.xWarp,2)) )
    % Note: we should use a better approach handle this situation...
    %warning('Two images are of different sizes:');
    rowX = size(inWarp.xWarp,1);
    colX = size(inWarp.xWarp,2);
    
    inX = inX(1:rowX,1:colX);
end

r = 0;  c = 0;
for j=1:rowX
    for k=1:colX
        r = inWarp.yWarp(j,k);
        c = inWarp.xWarp(j,k);
        rL = min(max(floor(r),1),rowX);
        rU = min(max(ceil(r),1),rowX);
        cL = min(max(floor(c),1),colX);  
        cU = min(max(ceil(c),1),colX);
        x = c-cL; y = r-rL;
        retX(j,k) = inX(rL,cL)*(1-x)*(1-y)+inX(rL,cU)*x*(1-y)+...
                        inX(rU,cL)*(1-x)*y+inX(rU,cU)*x*y;
    end
end
