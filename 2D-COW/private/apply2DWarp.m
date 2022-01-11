function [Xw,warpPath] = apply2DWarp(X, warpPath, Modes, Type)
% Apply warping to a signal
%
% [Xw] = Apply2DWarp(X, WarpPath, Modes, Type)
% 
% INPUT
% X       : data tensor/matrix
% warpPath: warping path structure
% Modes   : mode along which the warping is applied
% Type    : 'COW' warping generated with COW
% 
% OUTPUT
% Xw:
% 
% Author: Giorgio Tomasi
%         Giorgio.Tomasi@gmail.com
% 
% Created      : 25 March, 2007
% Last modified: 06 April, 2009; 11:03

% HISTORY
% 1.00.01 06 Apr 09 -> Added help

narginchk(2,4)
if (nargin < 3), Modes = 1:2;
elseif (numel(Modes) ~= 2 || any(fix(Modes) ~= Modes)), error('apply2DWarp:badModes','Modes must have two integers')
elseif (any(Modes > ndims(X))), error('apply2DWarp:badModes','One of the specified modes exceeds the array''s order') 
end
if (nargin < 4 || isempty(Type)), Type = 'cow'; end
m      = min(Modes);
M      = max(Modes);
N      = ndims(X);
ord    = [Modes,2:m - 1,m + 1:M - 1,M + 1:N];
doPerm = m ~= 1 || M ~= 2;
X      = permute(X,ord);
dimX   = size(X);
if (N > 4), error('apply2DCOW:multiwayChannels','Multiway signals are currently not handled'); end
dimX(N + 1:4) = ones(1,4 - N);
L             = numel(warpPath);
if (~isstruct(warpPath)), error('apply2DWarp:structWarping','Warping path must be provided as a structure')
elseif (~all(isfield(warpPath,{'x','y','xR','yR'}))), error('apply2DWarp:badFields','Invalid structure for warping')
elseif (~(L == 1 || L == dimX(4))), error('apply2DWarp:invalidWarping','Invalid warping')
elseif (L == 1), indWarp = ones(dimX(4),1);
else indWarp = 1:L;
end
predef = all(isfield(warpPath,{'xWarp','yWarp'}));
for (i = 1:dimX(4))
    if (warpPath(i).x(1,end) ~= dimX(2) || warpPath(i).y(end,1) ~= dimX(1)), error('apply2DWarp:incompatibleWarp','Warping path not compatible with X size'); end
end
Xw = NaN(dimX);
ind = {':',':',':'}; % Not only works for 4-way arrays
switch upper(Type)
   case 'COW'
       if (~predef)
           
           [A1,A2,AL1,AL2,AU1,AU2] = deal(zeros(dimX(1),dimX(2),1));  % CHECK [!!] is preallocation slower or faster? All samples are the same
           for (k = 1:dimX(4))
               
               w                           = indWarp(k);
               b                           = {warpPath(w).y(:,1),warpPath(w).x(1,:)}; % for legibility
               [A1(:,:),AL1(:,:),AU1(:,:)] = aWarp(warpPath(w).xR,b,false);  % Matrices are faster than 3 way (for loop)
               [A2(:,:),AL2(:,:),AU2(:,:)] = aWarp(warpPath(w).yR,b,true);
               Xw(ind{:},w)                = applyInterp(X,A1,A2,AL1,AL2,AU1,AU2,dimX(1:3));
               if (nargout > 1)
                   warpPath(w).xWarp = cat(3,A1,AL1,AU1);
                   warpPath(w).yWarp = cat(3,A2,AL2,AU2);
               end
               
           end
           
       else
           
           for (k = 1:dimX(4))
               w            = indWarp(k);
               Xw(ind{:},w) = applyInterp(X,warpPath(w).xWarp(:,:,1),warpPath(w).yWarp(:,:,1),...
                                            warpPath(w).xWarp(:,:,2),warpPath(w).yWarp(:,:,2),...
                                            warpPath(w).xWarp(:,:,3),warpPath(w).yWarp(:,:,3),dimX(1:3));
           end
           
       end

   otherwise
      error('Warping method not implemented')

end
if (doPerm), Xw = ipermute(reshape(Xw,dimX),ord); end

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

function [A,AL,AU] = aWarp(rW,B,rev)
if (rev)
    rW = rW';
    B  = B([2 1]);
end
nSeg   = size(rW);
dimX   = [B{1}(end),B{2}(end)];
A      = zeros(dimX);
nSeg2  = nSeg - 1;
segLen = diff(B{1},1);
T      = diff(rW,1,1);
for (j = 1:nSeg2(1))
    for (k = B{1}(j):(B{1}(j + 1) - 1)), A(k,:) = bWarp(A(k,:),rW(j,:) + T(j,:) * (k - B{1}(j))/segLen(j),B{2}); end
end
A(dimX(1),:) = bWarp(A(dimX(1),:),rW(nSeg2(1),:),B{2});
if (rev), A = A'; end
AL = min(max(floor(A),1),dimX(2));
AU = min(max(ceil(A),1),dimX(2));

    function c = bWarp(c,xW,b)
        for (m = 1:length(b) - 1)
            ind    = b(m):(b(m + 1) - 1);
            s      = b(m + 1) - (b(m) + 1);
            c(ind) = xW(m) + (0:s) * (xW(m + 1) - xW(m) - 1) / s;
        end
        c(dimX(2)) = xW(nSeg2(2) + 1);
        
    end

end
