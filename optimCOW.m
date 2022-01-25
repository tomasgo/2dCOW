function [optimPars,OS,diagnos] = optimCOW(X,optimSpace,varargin)
% [optim_pars,OS,diagnos] = optimCOW(X,optim_space,options,ref)
%
% Automatically optimizes the segment lengths and slack parameters for (2D)COW alignment.
% Optimality is found in terms of warping effect
%
% “Warping Effect = "Simplicity" + “Peak Factor”
%
% A grid search on the two dimensions ensued by "discrete-coordinates simplex" optimization
% Segment lengths and slack parameters for the COW alignment with so-called
% The number of dimensions to align depends on the number of rows (r) in optimSpace:
%     r = 1 -> Ordinary COW with alignment along mode 2 of X (i.e. columns)
%     r = 2 -> TwoDCow with alignment along modes 2 and 3 (i.e. columns and frontal slabs)
%
% NB r will matter only when X is at least of order (N) 3.
%    E.g. For N = 3, and r = 1, COW will only align along dimension 2 and will interpret dimension 3 as channels.
%         For N = 3, and r = 2, 2DCOW will align along dimension 2 and 3 and horizontal slabs in X are interpreted as single channel signals
%
%
% INPUT
% X           (s x t1 [x t2] x ch): s samples x t1 rt in mode 1 x t2 rt in mode 2 x ch channels
% optimSpace            (1|2 x 10): optimization space. row 1 -> mode 1, row 2 -> mode 2
%                                   columns 1:6: [min segLen, max segLen, min Slack, maxSlack, gridSizeSegLen, gridSizeSlack]
%                                   Columns 7:10 are relevant only for if simplex search is active
%                                           [segLen LB, segLen UB, Slack LB, Slack UB] (default: [4 lenX 1 lenX - 3])
%                                   where LB, UB and lenX are the lower bound, the upper bound and the length of the signal in the relevant mode
%                                   only columns 1 to 4 are required. Default grid size is min(5,b - a + 1) where a and
%                                   b denote lower and upper end for a parameter. Grid size may also be adjusted based
%                                   on options(7).
%     options                  (1): trigger plot and progress text
%                              (2): number of optimizations from grid maxima
%                              (3): maximum number of optimization steps
%                            (4:5): Maximum absolute correction as a fraction of signal length in modes 1 and 2 respectively
%                                   Default [0 3 50 0.15 0.15] (no plot; 3 starts from 3 maxima in grid search; maximum 50 steps; 15% max absolute correction)
%                              (6)  Use parallel toolbox for grid search
%                              (7)  Activate DP slack acceleration: the warping paths for all slacks within search space are calculated within COW
%         ref (1 x t1 [x t2] x ch): reference object used for alignment; if omitted first sample is used)
%
% OUTPUT
% optim_pars (1 x 4) optimal segment length and slack size
% OS         (7 x N) optimization sequence
%                    (1)  mode 1 segment length
%                    (2)  mode 1 slack parameter
%                    (3)  mode 2 segment length  (NaN in case of 1D alignment)
%                    (4)  mode 2 slack parameter (NaN in case of 1D alignment)
%                   (5:7) "Warping Effect", fourth "Simplicity", Fifth "Peak Factor")
%                    (8)  computation time
% diagnos (struct):  simplicity raw data, total run time, start points for optimization (columns in OS),
%                        "optim_space" and warping results for optimum (path + aligned matrix)
%
% uses ref_select.m, cow.m
%
% Authors:
% Giorgio Tomasi & Guilherme L. Alexandrino
%
% References
% OptimCOW...
% COW
% 2DCOW
% Code adapted from optimCOW by T. Skov and F.v.d.Berg

narginchk(2,4)
dimX = size(X);
N    = ndims(X);
ind  = repmat({':'},1,N - 1);

% Check optional input parameters
options    = checkOptions(varargin{:});
show       = options(1); % For legibility
[ref,refN] = checkReference(X,varargin{2:end});
if (~isempty(refN) && show), fprintf(1,'Object %i selected as reference',refN); end

% Grid search setup
[R,band,optimSpace] = checkOptimSpace(optimSpace,options,dimX); % Check grid parameters
OS                  = initOS(optimSpace,options(7));            % Form actual grid

% Initialise output
diagnos  = struct('baseSimplicity',[],'optimWarpFactor',[],'time_min',[],'optimStarts',[],'optimSteps',[],'optimSpace',optimSpace,...
    'reference',ref,'referenceSample',refN,'warping',[],'Xw',[]);

% Some initialisations
nGridRun = size(OS,2);
t00      = tic;
normX    = NaN(dimX(1) + 1,1); % Norms for peak factor
for (i = 1:dimX(1)), normX(i) = norm(X(i,:),2); end
normX(end) = norm(ref(:),2);
if (R == 2) % TwoDCOW requires samples to be the last mode
    X   = permute(X,[2:N,1]);
    ref = permute(ref,[2:N,1]);
end

% Actual grid search
if (show), fprintf(1,'Starting grid search\n'); end
nTopSlacks = prod(optimSpace(:,5),1);
nSlacks    = prod(optimSpace(:,6),1);
if (~options(7)), startGridWith = 1; else, startGridWith = nGridRun - nTopSlacks + 1; end
warpingPaths = cell(1,nGridRun);
if (options(6) && ~isempty(ver('parallel')))
    
    OSpar = num2cell(OS,1);
    parfor (iGrid = startGridWith:nGridRun)
        [OSpar{iGrid},warpingPaths{iGrid}] = optim_eval(X,OSpar{iGrid},ref,normX,band);
        if (show), showIter(0,iGrid,nGridRun,OSpar{iGrid}); end
    end
    OS = cat(2,OSpar{:});
    
else
    for (iGrid = startGridWith:nGridRun)
        [OS(:,iGrid),warpingPaths{iGrid}] = optim_eval(X,OS(:,iGrid),OS,ref,normX,band);
        if (show), showIter(0,iGrid,nGridRun,OS(:,iGrid)); end
        options(7) = false; % Just temporary fix
        if (options(7))
            indGrid = OS(:,1:R) == OS(1:R,iGrid);
            for (iGrid = 1:nGridRun - TopSlacks)
                OS(:,iGrid) = [];
            end
        end
        
    end
end
disp(warpingPaths)
% Simplex search
if (options(2))
    
    % Find best options(2) combinations
    if (verLessThan('matlab','9.3.0'))
        [~,starts] = sort(OS(5,:),'descend');
        tmp        = find(~isnan(starts),1,'first');
        starts     = starts(tmp + 1:tmp + options(2));
    else, [~,starts] = maxk(OS(5,:),options(2));
    end
    N        = nGridRun + 1;
    steps    = NaN(1,options(2));
    nVert    = 2 * R; % N. vertices polytope - 1
    upperLim = optimSpace(1,[8 10])';
    lowerLim = optimSpace(1,[7,9])';
    if (R > 1)
        upperLim = cat(1,upperLim,optimSpace(2,[8 10])'); 
        lowerLim = cat(1,lowerLim,optimSpace(2,[7  9])'); 
    end
    for a = 1:length(starts)
        
        if show
            showIter(3,a,options(2),OS(:,starts(a)))
            t0 = tic;
        end
        s        = 0;
        ps       = [starts(a) NaN(1,nVert)];
        polytope = OS(1:nVert,ps(1)) + 2 * eye(nVert) .* (double(OS(1:nVert,ps(1)) ~= upperLim) - 0.5); % Ensures it is within bounds
        for (i = 1:nVert)
            [OS,osPos] = calcVert(OS,polytope(:,i));
            if (osPos), ps(i + 1) = osPos; else, error('this should not happen!'); end % Error because all vertices in the polytope should be feasible
        end
        term                       = false;
        k                          = 1;
        [currentPolytope,ps,worst] = sortPolytope(OS,ps);
        while (~term)

            n          = newPoint(OS(5,ps),currentPolytope,k); % update the k-th worst result
            [OS,osPos] = calcVert(OS,n); % NOTE: if parallelised ensure that the same combination is not calculated twice
            if (osPos && OS(5,osPos) > worst) % Polytope is updated
                [currentPolytope,ps] = sortPolytope(OS,[osPos,ps(2:end)]);
                k                    = 1;
            else, k = k + 1; % The combination is invalid/worse than current worst result
            end
            term = s > options(3) || k > nVert + 1; % Max. iter or went through all vertices in polytope without updates
            
        end
        if (s > options(3) && k > nVert + 1), showIter(4,steps(a)); end
        steps(a) = s - 1;
        if (show)
            [~,aOptim] = max(OS(5,ps));
            showIter(5,a,options(2),OS(:,ps(aOptim)),steps(a),toc(t0)/60);
        end
                    
    end
    [diagnos.optim_warpFactor,optim] = max(OS(5,:));
    optimPars                        = OS(1:nVert,optim);
    diagnos.time_min                 = toc(t00)/60;       % to exclude the calculation of base simplicity
    if (show || nargout > 2)
        
        if (nargout > 2)
            
            diagnos   = struct('baseSimplicity',[],'optimWarpFactor',[],'time_min',[],'optimStarts',[],'optimSteps',[],'optimSpace',optimSpace,...
                'reference',ref,'referenceSample',refN,'warping',[],'Xw',[]);
            if (options(2))
                diagnos.optimStarts = starts;
                diagnos.optimSteps  = steps;
            end
            fprintf('Applying optimal warp\n');
            if (R == 1) % COW
                
            else % TwoDCOW
               
                try
                    
                    for (k = 1:dimX(1))
                        if (k == 1),[diagnos.Xw,diagnos.warping]              = TwoDCOW(X(ind{:},k),ref,optimPars(1:2),optimPars(3:4),band);
                        else,       [diagnos.Xw(ind{:},k),diagnos.warping(k)] = TwoDCOW(X(ind{:},k),ref,optimPars(1:2),optimPars(3:4),band);
                        end
                    end
                    diagnos.Xw = permute(diagnos.Xw,[N,1:(N-1)]);
                    
                catch, fprintf('[\bThe data set may be too big for the current 2DCOW implementation; final result not included in "diagnos".]\b\n');
                end
                
            end
            
        end
        diagnos.baseSimplicity = sum(svd(X/sqrt(sum(X(:).^2))).^4);         % Base simplicity
        if (show)
            
            plotOptim(OS,optim,optimSpace,diagnos.time_min,diagnos.baseSimplicity);    % Plotting
            % if (~isempty(diagnos.Xw))
            %
            %     figure
            %     subplot(2,1,1);
            %     plot(1:size(X,2),X);
            %     title('Data raw');
            %     grid;
            %     subplot(2,1,2); plot(1:size(diagnos.Xw,2),diagnos.Xw);
            %     title(['Data from optimal correction (segment ' num2str(optimPars(1)) ', slack ' num2str(optimPars(2)) ')']);
            %     grid;
            %
            % end
            
        end
        
    end
    
    
    
end

    function [OS,osPos] = calcVert(OS,n)
        [flag,osPos] = checkPoint(OS,n,lowerLim,upperLim);
        if (flag == 0)
            OS(:,N) = optim_eval(X,n,OS,ref,normX,band); 
            osPos   = N;
            N       = N + 1;
            s       = s + 1;
        elseif (flag > 1), osPos = 0; 
        end % if osPos == 0 -> The point has already been evaluated
        if (show && flag ~= 1), showIter(flag,s,options(3),OS(:,osPos)); end

    end

end

%%%
function [currentPolytope,ps,worst] = sortPolytope(OS,ps)
[~,ord]         = sort(OS(5,ps));
ps              = ps(ord);
currentPolytope = OS(:,ps);
worst           = currentPolytope(5,1); % Smallest warping factor

end


function [exitFlag,p] = checkPoint(OS,p,band,lowerBound,upperBound)
R        = length(band); 
index    = ~(isnan(OS(5,:)) | OS(5,:) == 0) & OS == p;
exitFlag = any(index);
if (exitFlag), p = find(index,1,'first'); 
else
    exitFlag = 2 * checkNonFeasibility(p(1),p(3),band(1),lowerBound(1:2),upperBound(1:2));
    if (R > 1 && exitFlag == 0)
        exitFlag = 2 * checkNonFeasibility(p(2),p(4),band(2),lowerBound(3:4),upperBound(3:4));
    end
end
    function flag = checkNonFeasibility(l,s,b,LB,UB)
       flag = l <= s + 3 || ...          % slack too large compared to segment length
           l < LB(1) || l > UB(1) || ... % segment length is out of bounds
           s < LB(2) || s > UB(2) || ... % slack is out of bounds
           s > b;                        % slack is larger than band-width
    end

end

function y = optim_eval(X,p,ref,normX,band)
y         = p;
y(5:8)    = NaN;
t0        = tic;
N         = ndims(X);
K         = size(X,max(3,N));
ind       = repmat({':'},max(N - 1,2),1);
normXcorr = NaN(K,1);
for (k = 1:K)
    xTmp         = TwoDCOW(X(ind{:},k),ref,p(1:2),p(3:4),band);
    normXcorr(k) = norm(xTmp(:),2);
    X(ind{:},k)  = xTmp;
    if (K == 1), s = sum(svd([X(:),ref(:)]/sqrt(normXcorr(k)^2 + normX(end)^2)).^4); end % Simplicity with one sample is calculated using the reference
end
if (K > 1), y(5:6) = sum(svd(reshape(X,numel(X)/K,K)'/sqrt(sum(normXcorr.^2))).^4); else, y(5:6) = s; end
y(7)   = mean((1 - min(abs((normXcorr - normX(1:K))./normX(1:K)),1)).^2);
y(5)   = y(7) + y(5);
y(8)   = toc(t0)/60;

end

function showIter(id,iter,total,varargin)
switch id
    case 0, fprintf(1,'   run %i/%i %s\n',iter,total,itStr(varargin{1}));
        %case 1, fprintf(1,'   run %i/%i: - min (segment/slack combination was already computed)\n',iter,total);
    case 2, fprintf(1,'   [\brun %i/%i: - illegal segment/slack combination)]\b\n',iter,total);
    case 3, fprintf(1,'Starting optimization %i/%i %s\n',iter,total,itStr(varargin{1}(1:7)));
    case 4, fprintf(1,'   Early termination after %i steps!',iter);
    case 5, fprintf(1,'   Optimization %i/%i terminated in %i steps %s\n',iter,total,varargin{2},itStr(cat(1,varargin{1}(1:7),varargin{3})));
    case 7
        fprintf(1,'Finished optimization\n');
        fprintf(1,'    Optimal values: %s\n',itStr(varargin{1}(1:7)));
        fprintf(1,'    Total time    : %3.2f min\n',varargin{2});
end

    function str = itStr(val)
        if (length(val) == 8), str = sprintf('seg. = {%i,%i}, sl. = {%i,%i} - %3.2f min',val(1:4),val(8));
        else, str = sprintf('seg. = {%i,%i}, sl. = {%i,%i}',val(1:4));
        end
    end

end

function showWarn(id,a,b,varargin)
switch id
    case 'smallOptimSegLen', msg = sprintf('Segment length grid reduced to %i values for mode %i',a,b);
    case 'badOptimSegLen',   msg = sprintf('Segment length max. in mode %i changed to %i to accommodate for meaningful grid',b,a);
    case 'smallOptimSlack',  msg = sprintf('Slack parameter grid reduced to %i values for mode %i',a,b);
    case 'badOptimSlack',    msg = sprintf('Slack parameter max. in mode %i changed to %i to accommodate for meaningful grid',b,a);
    case 'smallOptimBand',   msg = sprintf('Max. absolute correction (band constraint) in mode %i adjusted to accommodate for slack',b,a);
        
end
warning(sprintf('optimCOW:%s',id),msg)

end

function n = newPoint(L,X,s) % Assumes L is sorted from smallest to largest
if (any(diff(L,1) < 0)), error('Polytope not sorted'); end
D  = size(X,1);
cL = calcUpdate(sum(X < X(:,s),2));
cU = calcUpdate(sum(X > X(:,s),2));
n  = max(1,X(:,s) + cU - cL);

    function t = calcUpdate(v)
        t = zeros(D,1);
        t(v < D & v > 0) = 1;
        t(v == D)        = 2;
    end

end

function plotOptim(OS,optim,optimSpace,eTime,S)
P = cell(2,2);
for (i = 1:2)
    P{i,1} = fix(linspace(optimSpace(i,1),optimSpace(i,2),optimSpace(i,5)));
    P{i,2} = fix(linspace(optimSpace(i,3),optimSpace(i,4),optimSpace(i,6)));
end
Titles            = {'Warping factor','Simplicity','Peak Factor'};
[maxVal,posVal]  = max(OS(5:7,:),[],2);
for (i = 1:3)
    hFig = figure('Name',Titles{i});
    localPlot([1 3],4 + i,subplot(1,2,1,'NextPlot','add','XGrid','on','YGrid','on','ZGrid','on'),1)
    localPlot([2 4],4 + i,subplot(1,2,2,'NextPlot','add','XGrid','on','YGrid','on','ZGrid','on'),2)
end

    function localPlot(pred,resp,sh,dim)
        view(3)
        sh.Title.String = sprintf('Dimension %i',dim);
        av = accumarray(OS(pred,:)',OS(resp,:),[],@(a) mean(a,'omitnan'),NaN);
        m  = accumarray(OS(pred,:)',OS(resp,:),[],@(a) min(a,[],1,'omitnan'),NaN);
        M  = accumarray(OS(pred,:)',OS(resp,:),[],@(a) max(a,[],1,'omitnan'),NaN);
        av = av(P{dim,:});
        m  = m(P{dim,:});
        M  = M(P{dim,:});
        [A,B] = meshgrid(P{dim,2},P{dim,1});
        a = A(:)';
        b = B(:)';
        ph1 = plot3(sh,A,B,av,'Color','b','Marker','o','MarkerSize',6,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[.6 0 0],'DisplayName','Average');
        ph2 = plot3(sh,a([1 1],:),b([1 1],:),[m(:) M(:)]','Color','b','Marker','o','MarkerSize',2,'MarkerEdgeColor',[.6 0 0],'MarkerFaceColor',[.6 0 0],'LineStyle','--','DisplayName','min\rightarrowmax');
        xMax  = OS(pred(2),posVal(resp - 4));
        yMax  = OS(pred(1),posVal(resp - 4));
        zMax  = maxVal(resp - 4);
        xtMax = OS(7 - pred(2),posVal(resp - 4));
        ytMax = OS(3 - pred(1),posVal(resp - 4));
        ph3 = plot3(sh,xMax,yMax,zMax,'Marker','d','MarkerEdgeColor',[0 .6 0],'MarkerFaceColor',[0 .6 0],'MarkerSize',8,'LineStyle','none','DisplayName',sprintf('Max.:%1.2f',zMax));
        text(xMax,yMax,zMax,sprintf('  (%i,%i)',ytMax,xtMax))
        if (resp > 5)
            
            xMax  = OS(pred(2),optim);
            yMax  = OS(pred(1),optim);
            zMax  = OS(resp,optim);
            xtMax = OS(7 - pred(2),optim);
            ytMax = OS(3 - pred(1),optim);
            ph4 = plot3(sh,xMax,yMax,zMax,'Marker','^','MarkerEdgeColor',[.6 .6 0],'MarkerFaceColor',[.6 .6 0],'MarkerSize',8,'LineStyle','none','DisplayName','Optimum');
            text(xMax,yMax,zMax,sprintf('  (%i,%i)',ytMax,xtMax))
            
        else, ph4 = gobjects(0);
        end
        xlabel('Segment length'); ylabel('Slack size'); zlabel(hFig.Name);
        legend([ph1(1),ph2(1),ph3,ph4])
        axis tight
        offLim = [1.05 -0.05;-0.05 1.05];
        sh.XLim = sh.XLim * offLim;
        sh.YLim = sh.YLim * offLim;
        sh.ZLim = sh.ZLim * offLim;
        sh.XTick = P{dim,2};
        sh.YTick = P{dim,1};
        
    end

end

function [R,band,optimSpace] = checkOptimSpace(optimSpace,options,dimX)
% columns 1:6: [min segLen, max segLen, min Slack, maxSlack, gridSizeSegLen, gridSizeSlack]
%              [segLen LB, segLen UB, Slack LB, Slack UB] (default: [4 lenX 1 lenX - 3])
nanOS = ~isnan(optimSpace);
[R,c] = size(optimSpace);
d     = dimX(2:R + 1)';
if (~(size(optimSpace,2) >= 4 && all(nanOS(1:R * 4)) && isequal(fix(optimSpace(nanOS)),optimSpace(nanOS)))), error('optimCOW:badOptions','"optim_space" must have at least four columns and must contain only integers'); end
band                   = fix(options(4:R + 3) .* dimX(2:R + 1));
optimSpace(:,c + 1:10) = NaN;
[optimSpace,pNaN]      = fillOS(optimSpace);

% To prevent initial search is set as bounds if NaNs are present

% Compliance with boundaries checked before in case it is already partly defined
if (c > 6), optimSpace(:,1) = max(optimSpace(:,1),optimSpace(:, 7),'omitnan'); end % segLen lower bound
if (c > 7), optimSpace(:,2) = min(optimSpace(:,2),optimSpace(:, 8),'omitnan'); end % segLen upper bound
if (c > 8), optimSpace(:,3) = max(optimSpace(:,3),optimSpace(:, 9),'omitnan'); end % slack lower bound
if (c > 9), optimSpace(:,4) = min(optimSpace(:,4),optimSpace(:,10),'omitnan'); end % slack upper bound
[a,t]      = deal(optimSpace(:,2) - optimSpace(:,1) + 1);
[b,u]      = deal(optimSpace(:,4) - optimSpace(:,3) + 1);
a(a > 8)   = 5; % Default grid size = 5 if interval contains at least 9 values - i.e. 1:2:length
b(b > 8)   = 5;
optimSpace = fillOS(optimSpace,[a,b],pNaN);
for (i = 1:R)
    
    % Segment length
    if (t(i) < optimSpace(i,5))            % Ensure that enough grid points are available for the segment length search space
        optimSpace(i,5) = t + 1;
        showWarn('smallOptimSegLen',optimSpace(i,5),i)
    elseif (t(i) > optimSpace(i,5) && t(i) < (optimSpace(i,5) - 1) * 2)    % Ensures that all increments are at least 2 for tight grid
        optimSpace(i,2) = optimSpace(i,1) + (optimSpace(i,5) - 1) * 2 ;
        showWarn('badOptimSegLen',optimSpace(i,2),i)
    end
    
    % Slack
    if (options(7))
        optimSpace(i,3) = optimSpace(i,9); % All warping paths from the lower boundary to the maximum can be calculated
        optimSpace(i,6) = optimSpace(i,4) - optimSpace(i,3) + 1;
    else
        
        if (u(i) < optimSpace(i,6))            % Ensure that enough grid points are available for the slack parameter search space
            optimSpace(i,6) = u(i) + 1;
            showWarn('smallOptimSlack',optimSpace(i,6),i)
        elseif (u(i) > optimSpace(i,6) && u(i) < (optimSpace(i,6) - 1) * 2)    % Ensures that all increments are at least 2 for tight grid
            optimSpace(i,4) = optimSpace(i,3) + (optimSpace(i,6) - 1) * 2;
            showWarn('badOptimSlack',i,optimSpace(i,4))
        end
        
    end
    if (band(i) < optimSpace(i,4))
        band(i) = optimSpace(i,4);
        showWarn('smallOptimBand',i,band(i))
    end
    
end
if (any(optimSpace(:) <= 0)), error('optimCOW:badOptimSpace','Invalid choice of optimisation space'); end

    function [oS,pNaN] = fillOS(oS,a,pNaN)
        if (nargin < 2), a = NaN(R,2); end
        if (nargin < 3), pNaN = isnan(oS); end
        defValues = [NaN(R,4),a,4 * ones(R,1),d,ones(R,1),d - 3];
        oS(pNaN)  = defValues(pNaN);
        
    end

end

function OS = initOS(optimSpace,flag)
R = size(optimSpace,1);
P = cell(R,2);
for (i = 1:R)
    P{i,1} = fix(linspace(optimSpace(i,1),optimSpace(i,2),optimSpace(i,5)));
    P{i,2} = fix(linspace(optimSpace(i,3),optimSpace(i,4),optimSpace(i,6)));
end
[P{:}]          = ndgrid(P{:});
nGridRun        = numel(P{1});
OS              = NaN(8,nGridRun);
OS(1:(R * 2),:) = reshape(cat(2 * R + 1,P{:}),nGridRun,2 * R)';

end

function options = checkOptions(options)
optionsDef = [0 3 50 0.15 0.15 0 1];
L          = length(optionsDef); % for legibility
if (~nargin || isempty(options)), options = optionsDef; % Default options
elseif (numel(options) ~= length(options)), error('optimCOW:badOptions','Options must be a vector')
else, options(end + 1:L) = optionsDef(length(options) + 1:L);
end

end

function [ref,refN] = checkReference(X,ref)
refN = [];
N    = ndims(X);
ind  = repmat({':'},1,N - 1);
dimX = size(X);
if (nargin < 2 || isempty(ref))
    
    if (exist('ref_select','file')), [ref,~,refN] = ref_select(X,[],[5 0]);
    else
        refN = 1;
        ref  = X(refN,ind{:});
    end
    
elseif (numel(ref) == 1 && fix(ref) == ref && ref > 0 && ref <= dimX(1))
    refN = ref;
    ref  = X(ref,ind{:});
end

end
