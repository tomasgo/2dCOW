function [retX,retWarp] = TwoDCOW(inX,inRef,inSLen,inMaxW)
% function [retX,retWarp] = TwoDCOW(inX,inRef,inSLen,inMaxW)
%
% The two-dimensional correlation optimized warping (2DCOW) algorithm with
% linear interpolation, described in
%
% Zhang, D., Huang, X., Regnier, F. E. and Zhang, M. (2008) Two-dimensional
% correlation optimized warping algorithm for aligning GCxGC-MS data.
% Analytical Chemistry, 80: 2664-2671.
%
% Input:
%   inX -- the 2-dimensional image needs to be warped, rowX x colX;
%   inRef -- the 2-dimensional reference image, rowX x colX;
%   inSLen -- (1 x 2) segment lengths of all dimensions;
%   inMaxW -- (1 x 2) maximum degrees of warping in each segment of given
%             length, it is also called "slack" in N.-P.V. Nielsen, 
%             J.M. Carstensen and J. Smedegaard (1998) "Aligning of single
%             and multiple wavelength chromatographic profiles for
%             chemometric data analysis using correlation optimised
%             warping", Journal of Chromatography A, 805: 17-35.
%
% Output:
%   retX -- the corrected image which is warped to mach inRef;
%   retWarp -- warping coefficients, which are applied to inX to get retX.
%
% Note: inX and inRef are assumed to be sampled at same time points. The
%       sampling interval is taken as the time unit.
%
% Dabao Zhang     Sept. 16, 2006
%
% Revised by:
%       Dabao Zhang, May 16, 2007 -- Use the package WarpTB.
%
% Copyright (c) 2006-2007 by Zhang & Zhang.


if( nargin<4 )
    help TwoDCOW;
    return;
end

[rowX,colX] = size(inX);
if( (size(inRef,1)~=rowX) ||(size(inRef,2)~=colX) )
    % Note: we should use a better approach handle this situation...
    warning('Two images are of different sizes:');
    disp(['Reference image: ' num2str(size(inRef,1)) 'x' num2str(size(inRef,2))]);
    disp(['Work image: ' num2str(size(inX,1)) 'x' num2str(size(inX,2))]);
    
    rowX = min(rowX,size(inRef,1));
    colX = min(colX,size(inRef,2));
    inX = inX(1:rowX,1:colX);
    inRef = inRef(1:rowX,1:colX);
end

nRowSegs = floor(rowX/inSLen(1));   % Number of Segments
if( ~nRowSegs )
    inSLen(1) = floor(nRowSegs/2);
    nRowSegs = 2;
    warning('The segment length is forced to be smaller than warping-vector length!');
    disp(['inSLen(1) = ' num2str(inSLen(1))]);
end
if( inMaxW(1)>=(inSLen(1)*0.8) )
    inMaxW(1) = round(inSLen(1)*0.3);
    warning(['The slack (inMaxW(1)) is changed to: ', num2str(inMaxW(1))]);
end
T(1).Idx = 1:inSLen(1):(floor(rowX/inSLen(1))*inSLen(1));
if( T(1).Idx(end)~=rowX )
    T(1).Idx = [T(1).Idx rowX];
end

nColSegs = floor(colX/inSLen(2));   % Number of Segments
if( ~nColSegs )
    inSLen(2) = floor(nColSegs/2);
    nColSegs = 2;
    warning('The segment length is forced to be smaller than warping-vector length!');
    disp(['inSLen(2) = ' num2str(inSLen(2))]);
end
if( inMaxW(2)>=(inSLen(2)*0.8) )
    inMaxW(2) = round(inSLen(2)*0.3);
    warning(['The slack (inMaxW(2)) is changed to: ', num2str(inMaxW(2))]);
end
T(2).Idx = 1:inSLen(2):(floor(colX/inSLen(2))*inSLen(2));
if( T(2).Idx(end)~=colX )
    T(2).Idx = [T(2).Idx colX];
end

wDim = [length(T(1).Idx), length(T(2).Idx)];
retWarp.x = zeros(wDim(1),wDim(2));
retWarp.y = zeros(wDim(1),wDim(2));

for j=1:2
    tmpSeg = [T(3-j).Idx; T(3-j).Idx];
    for k=1:length(T(j).Idx)
        t = T(j).Idx(k);
        tmpRef = ProfileImage(inRef,j,t,inSLen(j));
        tmpX = ProfileImage(inX,j,t,inSLen(j));
        
        if( j==1 )
            [tmpW,~]       = COW_orig(tmpX',tmpRef',tmpSeg,inMaxW(2));
            retWarp.x(k,:) = tmpW(:,:,2);
        else
            %retWarp.y(:,k) = OneDCOW(tmpX,tmpRef,inSLen(1),inMaxW(1));
            [tmpW,~]       = COW_orig(tmpX',tmpRef',tmpSeg,inMaxW(1));
            retWarp.y(:,k) = tmpW(:,:,2)';
        end
    end
end

xWarp = zeros(rowX,colX);
for j=1:(length(T(1).Idx)-1)
    for k=T(1).Idx(j):(T(1).Idx(j+1)-1)
        % for each row of the image, get the warped one.
        tmpXW = retWarp.x(j,:)+(retWarp.x(j+1,:)-retWarp.x(j,:))...
                    *(k-T(1).Idx(j))/(T(1).Idx(j+1)-T(1).Idx(j));
        for jj=1:(length(T(2).Idx)-1)
            tmpIdxR = T(2).Idx(jj):(T(2).Idx(jj+1)-1);
            tmpSLen = T(2).Idx(jj+1)-T(2).Idx(jj);
            xWarp(k,tmpIdxR) = tmpXW(jj)+...
                   (0:(tmpSLen-1))*(tmpXW(jj+1)-tmpXW(jj)-1)/(tmpSLen-1);
        end
        % for the time at colX
        xWarp(k,colX) = tmpXW(nColSegs+1);
    end
end
% for k=rowX, get the warped one.
tmpXW = retWarp.x(nRowSegs+1,:);
for jj=1:(length(T(2).Idx)-1)
    tmpIdxR = T(2).Idx(jj):(T(2).Idx(jj+1)-1);
    tmpSLen = T(2).Idx(jj+1)-T(2).Idx(jj);
    xWarp(rowX,tmpIdxR) = tmpXW(jj)+...
           (0:(tmpSLen-1))*(tmpXW(jj+1)-tmpXW(jj)-1)/(tmpSLen-1);
end
% for the time at colX
xWarp(rowX,colX) = tmpXW(nColSegs+1);

yWarp = zeros(rowX,colX);
for j=1:(length(T(2).Idx)-1)
    for k=T(2).Idx(j):(T(2).Idx(j+1)-1)
        % for each column of the image, get the warped one.
        tmpXW = retWarp.y(:,j)+(retWarp.y(:,j+1)-retWarp.y(:,j))...
                    *(k-T(2).Idx(j))/(T(2).Idx(j+1)-T(2).Idx(j));
        for jj=1:(length(T(1).Idx)-1)
            tmpIdxR = T(1).Idx(jj):(T(1).Idx(jj+1)-1);
            tmpSLen = T(1).Idx(jj+1)-T(1).Idx(jj);
            yWarp(tmpIdxR,k) = tmpXW(jj)+...
                   (0:(tmpSLen-1))'*(tmpXW(jj+1)-tmpXW(jj)-1)/(tmpSLen-1);
        end
        % for the time at rowX
        yWarp(rowX,k) = tmpXW(nRowSegs+1);
    end
end
% for k=colX, get the warped one.
tmpXW = retWarp.y(:,nColSegs+1);
for jj=1:(length(T(1).Idx)-1)
    tmpIdxR = T(1).Idx(jj):(T(1).Idx(jj+1)-1);
    tmpSLen = T(1).Idx(jj+1)-T(1).Idx(jj);
    yWarp(tmpIdxR,colX) = tmpXW(jj)+...
           (0:(tmpSLen-1))'*(tmpXW(jj+1)-tmpXW(jj)-1)/(tmpSLen-1);
end
% for the time at rowX
yWarp(rowX,colX) = tmpXW(nRowSegs+1);

r = 0;  c = 0;
for j=1:rowX
    for k=1:colX
        r = yWarp(j,k);
        c = xWarp(j,k);
        rL = min(max(floor(r),1),rowX);
        rU = min(max(ceil(r),1),rowX);
        cL = min(max(floor(c),1),colX);  
        cU = min(max(ceil(c),1),colX);
        x = c-cL; y = r-rL;
        retX(j,k) = inX(rL,cL)*(1-x)*(1-y)+inX(rL,cU)*x*(1-y)+...
                        inX(rU,cL)*(1-x)*y+inX(rU,cU)*x*y;
    end
end

retWarp.xWarp = xWarp;
retWarp.yWarp = yWarp;