function [xAligned,warpPath] = TwoDCOW(X,T,segLen,slack,band)
% The two-dimensional correlation optimized warping (2DCOW) algorithm for multiple channel signals 
%
% [xAligned,warpPath] = TwoDCOW(X,T,segLen,slack)
%
% INPUT
% X     : sample (rt1 x rt2 x channels);
% T     : target (rt1 x rt2 x channels);
% segLen: segment lengths for mode(1) and (2)
% slack : slack parameter for mode(1) and (2)
% band  : max. absolute correction for mode (1) and (2)
%
% OUTPUT
% xAligned: aligned sample
% warpPath: warping coefficients
%
% References: 
% Zhang, D., et al. (2008) Analytical Chemistry, 80: 2664-2671.
% Nielsen,N.P.V. et al. (1998) Journal of Chromatography A, 805: 17-35.
%
% Note: inX and inRef are assumed to be sampled at same time points. The
%       sampling interval is taken as the time unit.
%
% Dabao Zhang     Sept. 16, 2006
%
% Revised by:
%       Dabao Zhang, May 16, 2007       -- Use the package WarpTB.
%       Giorgio Tomasi, March 14, 2019  -- Multichannel support and streamlining

narginchk(4,5)
dimX = size(X);
if (nargin < 5 || isempty(band)), band = [Inf,Inf]; end
dimT = size(T);
N    = ndims(X);
ind  = repmat({':'},1,ndims(X) - 2);
if(~isequal(dimX(1:2),dimT(1:2)))
    warning('TwoDCOW:nonMatchingSize','Reference and sample are of different sizes. Reference: [%i x %i] and sample: [%i x %i]',dimT(1:2),dimX(1:2));
    dimX(1:2) = min(dimX(1:2),dimT(1:2));
    X         = X(1:dimX(1),1:dimX(2),ind{:});
    T         = T(1:dimX(1),1:dimX(2),ind{:});
end
if (~isequal(dimX(3:end),dimX(3:end))),error('TwoDCOW:nonMatchingChannels','Number of channels does not match in reference and sample'); end
B                            = cell(2,1);
warpPath                     = cell(2,2);
nSeg2                        = zeros(1,2);
[segLen,slack,nSeg2(1),B{1}] = checkWarpParameters(dimX,segLen,slack,band,1);
[segLen,slack,nSeg2(2),B{2}] = checkWarpParameters(dimX,segLen,slack,band,2);
nSeg                         = cat(2,length(B{1}),length(B{2}));
[warpPath{1:4}]              = deal(zeros(nSeg));
[A,AL,AU]                    = deal(zeros(dimX(1),dimX(2),2));
for (j = 1:2)
    
    tmpSeg = B{3 - j}([1 1],:);
    ind    = {':',':'};
    rev    = [j,3 - j];
    ord    = [rev,3:N];
    for (k = 1:length(B{j}))
        
        ind{j}                = k;
        tmpRef                = permute(profImage(T,j,B{j}(k),segLen(j)),ord);
        tmpX                  = permute(profImage(X,j,B{j}(k),segLen(j)),ord);
        [tmpW,~]              = COW_orig_mw(tmpX,tmpRef,tmpSeg,slack(3 - j),[0 1 0 band(j) 0]); % Inverted because of the implementation of COW
        warpPath{j,1}(ind{:}) = permute(tmpW(:,:,2),rev);
        warpPath{j,2}(ind{:}) = permute(tmpW(:,:,1),rev);
        
    end
    [A(:,:,j),AL(:,:,j),AU(:,:,j)] = aWarp(permute(warpPath{j},rev),dimX(rev),nSeg(rev),nSeg2(rev),B(rev),rev);
    
end
xAligned = applyInterp(X,A(:,:,1),A(:,:,2),AL(:,:,1),AL(:,:,2),AU(:,:,1),AU(:,:,2),dimX);
warpPath = struct('x',warpPath{1,1},'y',warpPath{2,1},'xWarp',A(:,:,1),'yWarp',A(:,:,2),'xR',warpPath{1,2},'yR',warpPath{2,2});

end

function xAligned = applyInterp(X,Ax,Ay,AxL,AyL,AxU,AyU,dimX)
xAligned = NaN(dimX);
if (ismatrix(X)) % Non-negligible speedup
    
    for (j = 1:dimX(2))
        
        for (i = 1:dimX(1))
            
            cL = AxL(i,j);
            cU = AxU(i,j);
            rL = AyL(i,j);
            rU = AyU(i,j);
            x  = Ax(i,j) - cL;
            y  = Ay(i,j) - rL;
            dx = 1 - x;
            dy = 1 - y;
            xAligned(i,j) = (X(rL,cL) * dx + X(rL,cU) * x) * dy + (X(rU,cL) * dx + X(rU,cU) * x) * y;
            
        end
        
    end
    
else
    
    for (j = 1:dimX(2))
        
        for (i = 1:dimX(1))
            
            cL = AxL(i,j);
            cU = AxU(i,j);
            rL = AyL(i,j);
            rU = AyU(i,j);
            x  = Ax(i,j) - cL;
            y  = Ay(i,j) - rL;
            dx = 1 - x;
            dy = 1 - y;
            for (k = 1:dimX(3)), xAligned(i,j,k) = (X(rL,cL,k) * dx + X(rL,cU,k) * x) * dy + (X(rU,cL,k) * dx + X(rU,cU,k) * x) * y; end
            
        end
        
    end
    
end

end

function [A,AL,AU] = aWarp(rW,dimX,nSeg,nSeg2,B,rev)
A = zeros(dimX);
for (j = 1:(nSeg(1) - 1))
    a = B{1}(j+1)- B{1}(j);
    t = rW(j+1,:) - rW(j,:);
    for (k = B{1}(j):(B{1}(j + 1) - 1)), A(k,:) = bWarp(A(k,:),rW(j,:) + t * (k - B{1}(j))/a,B{2}); end
end
A(dimX(1),:) = bWarp(A(dimX(1),:),rW(nSeg2(1),:),B{2});
A            = permute(A,rev);
AL           = min(max(floor(A),1),dimX(2));
AU           = min(max(ceil(A),1),dimX(2));

    function c = bWarp(c,xW,b)
        for (m = 1:nSeg(2) - 1)
            ind    = b(m):(b(m + 1) - 1);
            s      = b(m + 1) - (b(m) + 1);
            c(ind) = xW(m) + (0:s) * (xW(m + 1) - xW(m) - 1) / s;
        end
        c(dimX(2)) = xW(nSeg2(2) + 1);
        
    end

end

function [sL,sP,nSeg,b] = checkWarpParameters(dimX,sL,sP,band,d)
nSeg = floor(dimX(d)/sL(d));   % Number of Segments
if (~nSeg)
    sL(d) = floor(nSeg/2);
    nSeg  = 2;
    id = sprintf('TwoDCOW:largeSegLen%i',d);
    warning(id,'The segment length in mode %i forced to %i',d,sL(d));
end
if(sP(d) >= (sL(d)*0.8) )
    sP(d) = round(sL(d)*0.3);
    id = sprintf('TwoDCOW:largeSlack%i',d);
    warning(id,'The slack parameter in mode %i forced to: %i',d,sP(d));
end
if (sP(d) > band(d))
    sP(d) = band(d);
    id = sprintf('TwoDCOW:smallBand%i',d);
    warning(id,'The slack parameter in mode %i reduced to: %i',d,sP(d));
end
s = sL(d) - 1;
b = 1:s:(floor((dimX(d) - 1)/s) * s); % The old version uses segments that are 1 point too long
if (b(end) ~= dimX(1)), b(end + 1) = dimX(d); end % FIXME: COW has an option to deal with this case

end

function prof = profImage(X,dim,pos,bW)
% Profile 2D chromatograms along dim at pos using the Epanechnikov kernel.
%
% prof = profImage(X,dim,pos,bW)
%
% Input:
% X  : data
% dim: dimension along which the profiles are obtained;
% pos: position of the profile
% bW : band width

dimX = size(X);
ind  = repmat({':'},1,ndims(X));
ind{dim} = max(1,pos - bW):min(dimX(dim),pos + bW); % CHECK: the width can be smaller than bW at boundaries
if (dim == 1), ind{dim} = ind{dim}'; end
if (verLessThan('matlab','9.3.0')), prof = bsxfun(@times,X(ind{:}),0.75 * (1-((ind{dim} - pos)/bW).^2));
else prof = X(ind{:}) .* (0.75 * (1-((ind{dim} - pos)/bW).^2));
end
prof = sum(prof,dim,'omitnan');
end
