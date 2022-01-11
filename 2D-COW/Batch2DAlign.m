function [] = Batch2DAlign(inFRef,inSLen,inMaxW,inFTIC,inFSIC)
%function [] = Batch2DAlign(inFRef,inSLen,inMaxW,inFTIC,inFSIC)
%
% Align 2-dimensional TIC/SIC chromatographs using the 2D-COW algorithm
% described in
%
% Zhang, D., Huang, X., Regnier, F. E. and Zhang, M. (2008) Two-dimensional
% correlation optimized warping algorithm for aligning GCxGC-MS data.
% Analytical Chemistry, 80: 2664-2671.
%
% Input:
%   inFRef -- name of  the file including the reference image;
%   inSLen -- 1x2, which specifies the section lengths (row & column);
%   inMaxW -- 1x2, which specifies the maximum warpings (row & column);
%   inFTIC -- name pattern of all TIC files;
%   inFSIC -- name pattern of all SIC files.
%
% Note: the file lists specified by inFTIC and inFSIC should be matched if
%       inFSIC is not empty.
%
% Dabao Zhang     May 15, 2007
%
% Revised by:
%
% Copyright (c) 2007 by Zhang & Zhang.

if( nargin<5 ) inFSIC = []; end
if( (nargin<4)|isempty(inFTIC) ) inFTIC = '*.txt'; end
if( (nargin<3)|isempty(inMaxW) ) inMaxW = [2,2]; end
if( (nargin<2)|isempty(inSLen) ) inSLen = [20,20]; end
if( nargin<1 ) error('The reference file must be specified!'); end

tic;

imgRef = load(inFRef);
ticList = dir(inFTIC);
nFTICs = length(ticList);
if( isempty(inFSIC) )
    nFSICs = nFTICs;
else
    sicList = dir(inFSIC);
    nFSICs = length(sicList);
end

if( nFTICs~=nFSICs )
    error('# of SIC files is different from # of TIC files!');
end

for j=1:nFTICs
    tmpFTIC = ticList(j).name;
    disp(['Current Image: ', tmpFTIC, '...']);

    if( isempty(inFSIC) )
        % Warping homogeneous TIC images:
        if( ~strcmp(upper(tmpFTIC),upper(inFRef)) )
            imgTmp = load(tmpFTIC);
            [imgAlgn,wCoeff] = TwoDCOW(imgTmp,imgRef,inSLen,inMaxW);
            
            fSave = ['A' tmpFTIC];
            save(fSave,'imgAlgn','-ascii');
        else
            fSave = ['A' tmpFTIC];
            save(fSave,'imgRef','-ascii');
        end
    else
        % Warping non-homogeneous images using SIC images:
        tmpFSIC = sicList(j).name;
        if( ~strcmp(upper(tmpFSIC),upper(inFRef)) )
            % warp SICs:
            imgTmp = load(tmpFSIC);
            [imgAlgn,wCoeff] = TwoDCOW(imgTmp,imgRef,inSLen,inMaxW);

            fSave = ['A' tmpFSIC];
            save(fSave,'imgAlgn','-ascii');
        
            % warp TICs:
            imgTmp = load(tmpFTIC);
            imgAlgn = WarpImage(imgTmp,wCoeff);

            fSave = ['A', tmpFTIC];
            save(fSave,'imgAlgn','-ascii');
        else
            % warp SICs:
            fSave = ['A', tmpFSIC];
            save(fSave,'imgRef','-ascii');
        
            % warp TICs:
            imgTmp = load(tmpFTIC);
            fSave = ['A',tmpFTIC];
            save(fSave,'imgTmp','-ascii');
        end
    end
end

toc;
