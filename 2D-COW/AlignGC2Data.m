function [] = AlignGC2Data(inFRef,inFTrg,inSLen,inMaxW,inDPks)
%function [] = AlignGC2Data(inFRef,inFTrg,inSLen,inMaxW)
%
% Align the GCxGC chromatographs using the 2D-COW algorithm
% described in
%
% Zhang, D., Huang, X., Regnier, F. E. and Zhang, M. (2008) Two-dimensional
% correlation optimized warping algorithm for aligning GCxGC-MS data.
% Analytical Chemistry, 80: 2664-2671.
%
% Input:
%   inFRef -- name of  the file including the reference image;
%   inFTrg -- names of  the target files.
%   inSLen -- 1x2, which specifies the section lengths (row & column);
%   inMaxW -- 1x2, which specifies the maximum warpings (row & column);
%
% Dabao Zhang     March 25, 2008
%
% Revised by:
%
% Copyright (c) 2008 by Zhang & Zhang.

if( nargin<5 ) inDPks = []; end
if( (nargin<4)|isempty(inMaxW) ) inMaxW = [2,2]; end
if( (nargin<3)|isempty(inSLen) ) inSLen = [20,20]; end
if( nargin<2 ) error('Need to specify target files ...'); end
if( nargin<1 ) error('Need to specify reference file ...'); end

tic;

nDPks = length(inDPks);

imgRef = load(inFRef);
trgList = dir(inFTrg);
nFTrgs = length(trgList);

for j=1:nFTrgs
    tmpFTrg = trgList(j).name;
    disp(['Current Image: ', tmpFTrg, '...']);

    if( ~strcmp(upper(tmpFTrg),upper(inFRef)) )
        % align and warp CSICs:
        imgTmp = load(tmpFTrg);
        [imgAlgn,wCoeff] = TwoDCOW(imgTmp,imgRef,inSLen,inMaxW);
        fSave = ['A', tmpFTrg];
        save(fSave,'imgAlgn','-ascii');

        if( nDPks==0 )
            fSave = ['W', strrep(tmpFTrg,'.txt','.mat')];
            save(fSave,'wCoeff');
        end
        
        % warp SICs
        for k=1:nDPks
            % warp D[id].txt and save into AD[id].txt
            tmpDFN = strrep(tmpFTrg,'C',['D' num2str(inDPks(k,:))]);
            tmpDF = WarpImage(load(tmpDFN),wCoeff);
            tmpOFN = ['A' tmpDFN];
            save(tmpOFN,'tmpDF','-ascii');
        end
    else
        % warp SICs
        fSave = ['A', tmpFTrg];
        save(fSave,'imgRef','-ascii');
        
        if( nDPks==0 )
            wCoeff = [];
            fSave = ['W', strrep(tmpFTrg,'.txt','.mat')];
            save(fSave,'wCoeff');
        end
        
        for k=1:nDPks
            % save D[id].txt into AD[id].txt
            tmpDFN = strrep(tmpFTrg,'C',['D' num2str(inDPks(k,:))]);
            tmpDF = load(tmpDFN);
            tmpOFN = ['A' tmpDFN];
            save(tmpOFN,'tmpDF','-ascii');
        end
    end
end

toc;