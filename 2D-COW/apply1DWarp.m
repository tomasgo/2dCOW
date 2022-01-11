function [Xw] = apply1DWarp(X, WarpPath, Mode, Type)
% Apply warping to a signal
%
% [Xw] = ApplyWarp(X, WarpPath, Mode, Type)
% 
% INPUT
% X       : data tensor/matrix
% WarpPath: warping path
% Mode    : mode along which the warping is applied
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
if (nargin < 3), Mode = 2; 
elseif (Mode > ndims(X)), error('apply1DWarp:exceedingOrder','The specified mode exceeds the array''s order') 
end
if (nargin < 4 || isempty(Type)), Type = 'cow'; end
dimX = size(X);
ord     = [Mode,2:Mode - 1,Mode + 1:ndims(X),1];
X       = permute(X,ord);
dimXper = dimX(ord);
switch size(WarpPath,1) 
    case 1,       X  = reshape(X,dimXper(1),prod(dimXper(2:end)));
    case dimX(1), X  = reshape(X,[dimXper(1),prod(dimXper(2:end - 1)),dimXper(end)]);
    otherwise,    error('apply1DWarp:badWarping','Invalid warping path')
end
Xw = NaN(WarpPath(1,end,2),size(X,2),size(X,3));
switch upper(Type)
   case 'COW'
      DecFlag = ~any(rem(WarpPath(:),1)); % Allow fractional positions for boundaries
      for i_sam = 1:size(X,3)

         for i_seg = 1:size(WarpPath,2) - 1

            if DecFlag
               indX = WarpPath(i_sam,i_seg,1):WarpPath(i_sam,i_seg + 1,1);
               indT = WarpPath(i_sam,i_seg,2):WarpPath(i_sam,i_seg + 1,2);
            else
               indT = linspace(WarpPath(i_sam,i_seg,2),WarpPath(i_sam,i_seg + 1,2),WarpPath(i_sam,i_seg + 1,2) - WarpPath(i_sam,i_seg,2) + 1);
               indX = linspace(WarpPath(i_sam,i_seg,1),WarpPath(i_sam,i_seg + 1,1),WarpPath(i_sam,i_seg + 1,1) - WarpPath(i_sam,i_seg,1) + 1);
            end
            lenX             = WarpPath(i_sam,i_seg + 1,1) - WarpPath(i_sam,i_seg,1);
            lenT             = WarpPath(i_sam,i_seg + 1,2) - WarpPath(i_sam,i_seg,2);
            Xw(indT,:,i_sam) = interp1q(indX' - WarpPath(i_sam,i_seg,1) + 1,X(indX,:,i_sam),(0:lenT)'/lenT * lenX + 1);

         end

      end

   otherwise
      error('Warping method not implemented')

end
dimXper(1) = size(Xw,1);
Xw         = reshape(Xw,dimXper);
ord(ord)   = 1:length(dimXper);
Xw         = permute(Xw,ord);
