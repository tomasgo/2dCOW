function [xAligned,warpPath] = TwoDCOW(X,T,segLen,slack,band,options)
% The two-dimensional correlation optimized warping (2DCOW) algorithm for multiple channel signals
%
% [xAligned,warpPath] = TwoDCOW(X,T,segLen,slack)
%
% INPUT
% X      : sample (rt1 x rt2 x channels x samples);
% T      : target (rt1 x rt2 x channels x 1);
% segLen : segment lengths for mode(1) and (2)
% slack  : slack parameter for mode(1) and (2)
% band   : max. absolute correction for mode (1) and (2)
% options: (1) extract all the slacks
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
if (N > 4), error('TwoDCOW:multiwayChannels','Multiway signals are currently not handled'); end
ind           = repmat({':'},1,ndims(X) - 2);
dimX(N + 1:4)   = ones(1,4 - N);
dimT(end + 1:3) = ones(1,3 - length(dimT)); % Target has no sample dimension
if(~isequal(dimX(1:2),dimT(1:2)))
    warning('TwoDCOW:nonMatchingSize','Reference and sample are of different sizes. Reference: [%i x %i] and sample: [%i x %i]',dimT(1:2),dimX(1:2));
    dimX(1:2) = min(dimX(1:2),dimT(1:2));
    X         = X(1:dimX(1),1:dimX(2),ind{:});
    T         = T(1:dimX(1),1:dimX(2),ind{:});
end
if (~isequal(dimT(3),dimX(3))),error('TwoDCOW:nonMatchingChannels','Number of channels does not match in reference and sample'); end
B                            = cell(2,1);
warpPath                     = cell(2,2);
nSeg                        = zeros(1,2);
for (d = 1:2), [segLen(d),slack(d),nSeg(d),B{d}] = checkWarpParameters(dimX(d),segLen(d),slack(d),band(d),d); end
nBound                       = cat(2,length(B{1}),length(B{2}));
[warpPath{1:4}]              = deal(zeros([nBound,dimX(4)]));
for (j = 1:2)
    
    tmpSeg = B{3 - j}([1 1],:);
    ind    = {':',':'};
    rev    = [3 - j,j];
    ord    = [4 rev,3];
    for (k = 1:length(B{j}))
        
        ind{j}                  = k;
        tmpRef                  = permute(profImage(T,j,B{j}(k),segLen(j)),ord);
        tmpX                    = permute(profImage(X,j,B{j}(k),segLen(j)),ord); % ord needs to be modified
        [tmpW,~]                = COW(tmpRef,tmpX,tmpSeg,slack(3 - j),[0 1 0 band(j) 0]); % Inverted because of the implementation of COW
        warpPath{j,1}(ind{:},:) = permute(tmpW(:,:,1),[rev + 1, 1]);
        warpPath{j,2}(ind{:},:) = permute(tmpW(:,:,2),[rev + 1, 1]);
        
    end
    
end
warpPath            = squeeze(struct('x',num2cell(warpPath{1,1},[1 2]),'y',num2cell(warpPath{2,1},[1 2]),'xR',num2cell(warpPath{1,2},[1 2]),'yR',num2cell(warpPath{2,2},[1 2])));
[xAligned,warpPath] = apply2DWarp(X,warpPath);

end

function [sL,sP,nSeg,b] = checkWarpParameters(dimX,sL,sP,band,d)
nSeg = floor((dimX - 1)/(sL - 1));   % Number of Segments
if (~nSeg)
    sL = floor(nSeg/2);
    nSeg  = 2;
    id = sprintf('TwoDCOW:largeSegLen%i',d);
    warning(id,'The segment length in mode %i forced to %i',d,sL);
end
if(sP >= (sL*0.8) )
    sP = round(sL*0.3);
    id = sprintf('TwoDCOW:largeSlack%i',d);
    warning(id,'The slack parameter in mode %i forced to: %i',d,sP);
end
if (sP > band)
    sP = band;
    id = sprintf('TwoDCOW:smallBand%i',d);
    warning(id,'The slack parameter in mode %i reduced to: %i',d,sP);
end
s = sL - 1;
b = 1:s:(floor((dimX - 1)/s) * s); % The old version uses segments that are 1 point too long
if (b(end) ~= dimX), b(end + 1) = dimX; end % FIXME: COW has an option to deal with this case

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
ind{dim} = max(1,pos - bW):min(dimX(dim),pos + bW); % CHECK [!!!!]: should not bW be half the segment?
% CHECK [!]: the width can be smaller than bW at boundaries
if (dim == 1), ind{dim} = ind{dim}'; end
if (verLessThan('matlab','9.3.0')), prof = bsxfun(@times,X(ind{:}),0.75 * (1-((ind{dim} - pos)/bW).^2));
else prof = X(ind{:}) .* (0.75 * (1-((ind{dim} - pos)/bW).^2));
end
prof = sum(prof,dim,'omitnan');
end
