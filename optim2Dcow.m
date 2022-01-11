function [optimPars,OS,diagnos] = optim2Dcow(X,optimSpace,options,ref)
% [optim_pars,OS,diagnos] = optim_cow(X,optim_space,options,ref)
%
% Automatically optimizes the segment lengths and slack parameters for 2DCOW alignment.
% Optimality is found in terms of warping effect
%
% “Warping Effect = "Simplicity" + “Peak Factor”
%
% A grid search on the two dimensions ensued by "discrete-coordinates simplex" optimization
% Segment lengths and slack parameters for the COW alignment with so-called
%
% INPUT
% X           (s x t1 x t2 x ch): s samples x t1 rt in mode 1 x t2 rt in mode 2 x ch channels
% optim_space            (2 x 6): optimization space. row 1 -> mode 1, row 2 -> mode 2
%                                 columns: [min segLen, max segLen, min Slack, maxSlack, gridSizeSegLen, gridSizeSlack]
%     options                (1): trigger plot and progress text
%                            (2): number of optimizations from grid maxima
%                            (3): maximum number of optimization steps
%                          (4:5): Maximum absolute correction as a fraction of signal length in modes 1 and 2 respectively
%                                 Default [0 3 50 0.15 0.15] (no plot; 3 starts from 3 maxima in grid search; maximum 50 steps; 15% max absolute correction)
%         ref (1 x t1 x t2 x ch): reference object used in 2DCOW alignment; if omitted first sample is used (i.e. X(1,:,:,:))
%
% OUTPUT
% optim_pars (1 x 4) optimal segment length and slack size
% OS         (7 x N) optimization sequence
%                    (1) mode 1 segment length
%                    (2) mode 1 slack parameter
%                    (3) mode 2 segment length
%                    (4) mode 2 slack parameter
%                    (5:7) "Warping Effect", fourth "Simplicity", Fifth "Peak Factor")
%                    (8) computation time
%      diagnos (struct): simplicity raw data, total run time, start points for optimization (columns in OS),
%                        "optim_space" and warping results for optimum (path + aligned matrix)
%
% uses ref_select.m, cow.m
%
% Authors:
% Giorgio Tomasi & Guilherme L. Alexandrino
%
% References
% OptimCOW...
% 2DCOW
% Code adapted from optimCOW by T. Skov and F.v.d.Berg

narginchk(2,5)
refN       = 1;
optionsDef = [0 3 50 0.15 0.15];
L          = length(optionsDef); % for legibility
if (nargin < 3 || isempty(options)), options = optionsDef; % Default options
elseif (numel(options) ~= length(options)), error('optim2DCOW:badOptions','Options must be a vector')
else options(end + 1:L) = optionsDef(length(options) + 1:L);
end
N                  = ndims(X);
ind                = repmat({':'},1,N - 1);
dimX               = size(X);
if (nargin < 4 || isempty(ref))
    if (exist('ref_select','file')), [ref,~,refN] = ref_select(X,[],[5 0]); else ref  = X(refN,ind{:}); end
    if (options(1)), fprintf(1,'Object %i selected as reference',refN); end
elseif (numel(ref) == 1 && fix(ref) == ref && ref > 0 && ref <= dimX(1)), ref = X(ref,ind{:});  % FIXME: Ignore the ref sample in the optimisation!
end

% 4D grid search setup
% Check grid parameters
band = fix(options(4:5) .* dimX(2:3));
if (~(isequal(size(optimSpace),[2 6]) && isequal(fix(optimSpace),optimSpace))), error('optim2DCOW:badOptions','"optim_space" must be of size (2 x 6) and must contain only integers');
else
    
    for (i = 1:2)
        
        t = optimSpace(i,2) - optimSpace(i,1) + 1;
        if (t < optimSpace(i,5))            % Ensure that enough grid points are available for the segment length search space
            optimSpace(i,5) = t + 1;
            showWarn('smallOptimSegLen',optimSpace(i,5),i)
        elseif (t > optimSpace(i,5) && t < (optimSpace(i,5) - 1) * 2)    % Ensures that all increments are at least 2 for tight grid
            optimSpace(i,2) = optimSpace(i,1) + (optimSpace(i,5) - 1) * 2 ;
            showWarn('badOptimSegLen',optimSpace(i,2),i)
        end
        t = optimSpace(i,4) - optimSpace(i,3) + 1;
        if (t < optimSpace(i,6))            % Ensure that enough grid points are available for the slack parameter search space
            optimSpace(i,6) = t + 1;
            showWarn('smallOptimSlack',optimSpace(i,6),i)
        elseif (t > optimSpace(i,6) && t < (optimSpace(i,6) - 1) * 2)    % Ensures that all increments are at least 2 for tight grid
            optimSpace(i,4) = optimSpace(i,3) + (optimSpace(i,6) - 1) * 2;
            showWarn('badOptimSlack',i,optimSpace(i,4))
        end
        if (band(i) < optimSpace(i,4))
            band(i) = optimSpace(i,4);
            showWarn('smallOptimBand',i,band(i))
        end
        
    end
    
end

% Form actual grid
P = cell(2,2);
for (i = 1:2)
    P{i,1} = fix(linspace(optimSpace(i,1),optimSpace(i,2),optimSpace(i,5)));
    P{i,2} = fix(linspace(optimSpace(i,3),optimSpace(i,4),optimSpace(i,6)));
end
[P{:}]    = ndgrid(P{:});
nGridRun  = numel(P{1});
OS        = NaN(8,nGridRun);
OS(1:4,:) = reshape(cat(5,P{:}),nGridRun,4)';
diagnos   = struct('baseSimplicity',[],'optimWarpFactor',[],'time_min',[],'optimStarts',[],'optimSteps',[],'optimSpace',optimSpace,...
    'reference',ref,'referenceSample',refN,'warping',[],'Xw',[]);

t00   = tic;
normX = NaN(dimX(1) + 1,1);
for (i = 1:dimX(1)), normX(i) = norm(X(i,:),2); end
normX(end) = norm(ref(:),2);
X          = permute(X,[2:N,1]);
ref        = permute(ref,[2:N,1]);
if (options(1)), fprintf(1,'Starting grid search\n'); end
for (iGrid = 1:nGridRun)
    [OS(:,iGrid),exitFlag] = optim_eval(X,OS(:,iGrid),OS,ref,normX,band);
    if (options(1)), showIter(exitFlag,iGrid,nGridRun,OS(:,iGrid)); end
end
if (options(2))
    
    if (verLessThan('matlab','9.3.0'))
        [~,starts] = sort(OS(5,:),'descend');
        tmp        = find(~isnan(starts),1,'first');
        starts     = starts(tmp + 1:tmp + options(2));
    else [~,starts] = maxk(OS(5,:),options(2));
    end
    N     = nGridRun + 1;
    steps = NaN(1,options(2));
    for a = 1:length(starts)
        
        Na = N;
        if options(1)
            showIter(3,a,options(2),OS(:,starts(a)))
            t0 = tic;
        end
        s  = 0;
        ps = [starts(a) NaN(1,4)];
        T  = bsxfun(@plus,eye(4),OS(1:4,ps(1)));
        for (i = 1:4)
            
            [res,exitFlag] = optim_eval(X,T(:,i),OS,ref,normX,band);
            if (~exitFlag)
                OS(:,N)   = res;
                ps(i + 1) = N;
                N         = N + 1;
                s         = s + 1;
            elseif (exitFlag == 1), ps(i + 1) = find(all(OS(1:4,:) == res(1:4)));
            end
            if (options(1)), showIter(exitFlag,s,options(3),res); end
            
        end
        term = false;
        while (~term)
            
            [n,worst,c]    = newPoint(OS(5,ps),OS(1:4,ps),1);
            [res,exitFlag] = optim_eval(X,n,OS,ref,normX,band);
            conv = res(5) > worst;
            if (~exitFlag)
                
                OS(:,N) = res;
                if (conv), ps(c) = N; end
                N = N + 1;
                s = s + 1;
                if (options(1)), showIter(exitFlag,s,options(3),res); end
                
            elseif (conv && exitFlag == 1), ps(c) = find(all(OS(1:4,:) == res(1:4)));
            elseif (~conv)
                
                k = 2;
                while (k < 5 && ~conv)
                    
                    [n,~,c]        = newPoint(OS(5,ps),OS(1:4,ps),k);
                    [res,exitFlag] = optim_eval(X,n,OS,ref,normX,band);
                    conv           = OS(5,N) > worst;
                    if (~exitFlag)
                        
                        OS(:,N) = res;
                        if (conv), ps(c) = N; end
                        N    = N + 1;
                        s    = s + 1;
                        if (options(1)), showIter(exitFlag,s,options(3),res); end
                        
                    elseif (conv && exitFlag == 1), ps(c) = find(all(OS(1:4,:) == res(1:4)));
                    end
                    k = k + 1;
                    
                end
                
            end
            steps(a) = s;
            term     = ~conv || s >= options(3);
            
        end
        if (~conv && steps(a) >= options(3)), showIter(4,steps(a)); end
        if (options(1))
            [~,aOptim] = max(OS(5,Na:end));
            showIter(5,a,options(2),OS(:,Na + aOptim - 1),steps(a),toc(t0)/60);
        end
        
    end
    
end
[diagnos.optim_warpFactor,optim] = max(OS(5,:));
optimPars                        = OS(1:4,optim);
diagnos.time_min                 = toc(t00)/60;       % to exclude the calculation of base simplicity
if (options(1) || nargout > 2)
    
    if (nargout > 2)
        
        diagnos   = struct('baseSimplicity',[],'optimWarpFactor',[],'time_min',[],'optimStarts',[],'optimSteps',[],'optimSpace',optimSpace,...
            'reference',ref,'referenceSample',refN,'warping',[],'Xw',[]);
        if (options(2))
            diagnos.optimStarts = starts;
            diagnos.optimSteps  = steps;
        end
        fprintf('Applying optimal warp\n');
        try
            
            for (k = 1:dimX(1))
                if (k == 1),[diagnos.Xw,diagnos.warping]              = TwoDCOW(X(ind{:},k),ref,optimPars(1:2),optimPars(3:4),band);
                else        [diagnos.Xw(ind{:},k),diagnos.warping(k)] = TwoDCOW(X(ind{:},k),ref,optimPars(1:2),optimPars(3:4),band);
                end
            end
            diagnos.Xw = permute(diagnos.Xw,[N,1:(N-1)]);
            
        catch, fprintf('[\bThe data set may be too big for the current 2DCOW implementation; final result not included in "diagnos".]\b\n');
        end
        
    end
    diagnos.baseSimplicity = sum(svd(X/sqrt(sum(X(:).^2))).^4);         % Base simplicity
    if (options(1))
        
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

%%%
function [y,exitFlag] = optim_eval(X,p,OS,ref,normX,band)
index3 = ~(isnan(OS(5,:)) | OS(5,:) == 0);
y      = p;
y(5:8) = NaN;
for (i = 1:4), index3 = index3 & OS(i,:) == p(i); end
exitFlag = any(index3);
if (exitFlag), y = OS(:,find(index3,1,'first'));
elseif ((p(1) <= p(3) + 3) || (p(3) < 1) || (p(2) <= p(4) + 3) || (p(4) < 1)), exitFlag = 2; % segment > slack OR slack < 1
else
    
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
    if (K > 1), y(5:6) = sum(svd(reshape(X,numel(X)/K,K)'/sqrt(sum(normXcorr.^2))).^4); else y(5:6) = s; end
    y(7)   = mean((1 - min(abs((normXcorr - normX(1:K))./normX(1:K)),1)).^2);
    y(5)   = y(7) + y(5);
    y(8)   = toc(t0)/60;
    
end

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
        else str = sprintf('seg. = {%i,%i}, sl. = {%i,%i}',val(1:4));
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
warning(sprintf('optim2DCOW:%s',id),msg)

end

function [n,M,c] = newPoint(L,X,s)
[M,c] = sort(L);
M     = M(s);
c     = c(s);
n     = X(:,c);
D     = size(X,1);
cL    = calcUpdate(sum(X < n,2));
cU    = calcUpdate(sum(X > n,2));
n     = max(1,n + cU - cL);

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
        
        else ph4 = gobjects(0);
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