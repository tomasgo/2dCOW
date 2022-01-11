function retProfile = ProfileImage(inImage,inPDim,inT,inBWidth)
% function retProfile = ProfileImage(inImage,inPDim,inT,inBWidth)
%
% Profile a two-dimensional image at a specific location of a specific
% dimension, using the Epanechnikov kernel.
%
% Input:
%   inImage -- the 2-dimensional image;
%   inPDim -- the dimension in which the image is profiled;
%   inT -- the location at which the image is profiled;
%   inBWidth -- profiling band width;
%
% Output:
%   retProfile -- the profiled vector.
%
% Dabao Zhang     Sept. 16, 2006
%
% Revised by:
%
% Copyright (c) 2006 by Zhang & Zhang.


[row, col] = size(inImage);
if( inPDim==1 )
    % profiled out the row
    tmpIdx = [max(1,inT-inBWidth):min(row,inT+inBWidth)]';
    tmpX = inImage(tmpIdx,:);
    
    %Epanechnikov kernel
    tmpW = 0.75*(1-((tmpIdx-inT)/inBWidth).^2); % CHECK: Potential error here (as per row 26, the width can be smaller than inBWidth)
    wTensor = kron(ones(1,col),tmpW);           % FIXME: unnecessary
    retProfile = sum(tmpX.*wTensor,1)';         % retProfile = sum(inImage(tmpIdx,:) .* (1-((tmpIdx-inT)/size(inImage,1)).^2);<
else
    % profiled out the column
    tmpIdx = [max(1,inT-inBWidth):min(col,inT+inBWidth)];
    tmpX = inImage(:,tmpIdx);
    
    %Epanechnikov kernel
    tmpW = 0.75*(1-((tmpIdx-inT)/inBWidth).^2);
    wTensor = kron(ones(row,1),tmpW);
    retProfile = sum(tmpX.*wTensor,2);
end
