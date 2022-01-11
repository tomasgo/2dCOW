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
% optim_pars (1 x 2) optimal segment length and slack size
% OS         (7 x N) optimization sequence
%                    (1) mode 1 segment length
%                    (2) mode 1 slack parameter
%                    (3) mode 2 segment length
%                    (4) mode 2 slack parameter
%                    (5:7) "Warping Effect", fourth "Simplicity", Fifth "Peak Factor")
%                    (8) computation time
%      diagnos (struct): simplicity raw data, total run time, start points for optimization (columns in OS),
%                        "optim_space" and warping results for optimum (path + aligned matrix + diagnostics)
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
if (nargin < 3 || isempty(options)), options = [0 3 50 0.15]; % Default options
elseif (numel(options) ~= length(options)), error('optim2DCOW:badOptions','Options must be a vector')
end
L                  = length(optionsDef); % for legibility
options(end + 1:L) = optionsDef(length(options) + 1:L);

if (nargin < 4 || isempty(ref))
    if (exist('ref_select','file')), [ref,~,refN] = ref_select(X,[],[5 0]); else ref  = X(refN,:,:,:); end
    if (options(1)), fprintf(1,'Object %i selected as reference',refN); end
end

% 4D grid search setup
% Check grid parameters
if (~(isequal(size(optimSpace),[2 6]) && isequal(fix(optimSpace),optimSpace))), error('optim2DCOW:badOptions','"optim_space" must be of size (2 x 6) and must contain only integers');
else
    
    for (i = 1:2) 
        
        t = optimSpace(i,2) - optimSpace(i,1);
        if (t < optimSpace(i,5))            % Ensure that enough grid points are available for the segment length search space
            optimSpace(i,5) = t + 1;
            showWarn('smallOptimSegLen',optimSpace(i,5),i)
        elseif (t < optimSpace(i,5) * 2)    % Ensures that all increments are at least 2 for tight grid
            optimSpace(i,2) = optimSpace(i,1) + optimSpace(i,5) * 2;
            showWarn('badOptimSegLen',optimSpace(i,2),i)
        end
        t = optimSpace(i,4) - optimSpace(i,3);
        if (t < optimSpace(i,6))            % Ensure that enough grid points are available for the slack parameter search space
            optimSpace(i,6) = t + 1;
            showWarn('smallOptimSlack',optimSpace(i,6),i)
        elseif (t < optimSpace(i,6) * 2)    % Ensures that all increments are at least 2 for tight grid
            optimSpace(i,4) = optimSpace(i,3) + optimSpace(i,6) * 2;
            showWarn('badOptimSlack',i,optimSpace(i,4))
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
nGridRun  = length(P{1});
OS        = zeros(8,nGridRun);
OS(1:4,:) = reshape(cat(5,P{:}),numel(Pp)/4,4)'; 
diagnos   = struct('base_simplicity',[],'time_min',[],'optim_starts_in_OS',[],'optim_steps_in_OS',[],'optim_space',optimSpace,'reference',ref,'reference_sample',refN,'warping',[],'Xw',[],'warping_diagnos',[]);

t00 = tic;
if (options(1)), fprintf(1,'Starting grid search\n'); end
for (iGrid = 1:nGridRun)
    [OS(:,iGrid),exitFlag] = optim_eval(X,OS(1:2,iGrid),OS,ref,round(length(ref)*options(4)));
    if options(1), showIter(exitFlag,iGrid,nGridRun,OS(8,iGrid)); end
end
if (verLessThan('matlab','9.3.0')), [~,starts] = maxk(OS(5,:),options(2));
else
    [~,starts] = sort(OS(5,:),'descend');
    starts     = starts(1:options(2));
end
N          = nGridRun + 1;
steps      = NaN(1,options(3));
for a = 1:length(starts)
    
    if options(1) 
        showIter(3,a,length(starts),OS(:,starts(a)))
        t0 = clock;
    end
    Na              = N - 1;
    ps              = [starts(a) N:N + 3];
    OS(1:4,N:N + 3) = eye(4) + OS(1:4,ps(1));
    N               = N + 3;
    for (i = 1:4)
        [OS(:,Na + i),exitFlag] = optim_eval(X,OS(1:4,Na + i),OS,ref,round(length(ref)*options(4)));
        showIter(exitFlag,N - Na,options(3),OS(:,Na + i));
    end
    term = false;
    while (~term)
        
        N           = N + 1;
        [n,worst,c] = newPoint(OS(5,ps),OS(1:4,ps),1);
        OS(:,N)     = optim_eval(X,n,OS,ref,round(length(ref) * options(4)));
        conv        = OS(5,N) > worst;
        if (conv),ps(c) = N;
        else
            
            k = 2;
            while (k < 5)
                
                N         = N + 1;
                [n,~,c]   = newPoint(OS(5,ps),OS(1:4,ps),k);
                [OS(:,N)] = optim_eval(X,n,OS,ref,round(length(ref)*options(4)));
                conv      = OS(5,N) > worst;
                k         = k + 1;
                if (conv)
                    ps(c) = N;
                    k     = 5;
                end
                
            end
            
        end
        steps(a) = N - Na - 1;
        term   = ~conv || steps(a) >= options(3);
        
    end
    if (~conv && steps(a) >= options(3)), showIter(4,steps(a)); end
    if (options(1)), showIter(5,a,options(2),OS(:,N),steps(a),etime(clock,t0)/60); end

end

[~,optim]        = max(OS(3,:));
optimPars        = OS(1:4,optim);
diagnos.time_min = toc(t00)/60;       % to exclude the calculation of base simplicity
if (options(1) || nargout > 2)
    
    if (nargout > 2)
        
        diagnos.optim_starts_in_OS = starts;
        diagnos.optim_steps_in_OS  = steps;
        fprintf('Computing diagnostics\n'); 
        try    [diagnos.warping,diagnos.Xw,diagnos.warping_diagnos] = cow(ref,X,optimPars(1),optimPars(2),[0 1 0 round(length(ref)*options(4)) 0]);
        catch, fprintf('[\bThe data set is too big for the current 2DCOW implementation; final result not included in "diagnos".]\b\n');
        end
        
    end
    diagnos.base_simplicity = sum(svd(X/sqrt(sum(X(:).^2))).^4);         % Base simplicity
    if (options(1))
        
        plotOptim(OS,optim,diagnos.time_min,diagnos.base_simplicity);    % Plotting
        if (~isempty(diagnos.Xw))
            
            figure
            subplot(2,1,1);
            plot(1:size(X,2),X);
            title('Data raw');
            grid;
            subplot(2,1,2); plot(1:size(diagnos.Xw,2),diagnos.Xw);
            title(['Data from optimal correction (segment ' num2str(optimPars(1)) ', slack ' num2str(optimPars(2)) ')']);
            grid;
            
        end 
        
    end
    
end

end

%%%
function [y,exitFlag] = optim_eval(X,p,OS,ref,losange)
t0     = clock;
index3 = true(1,size(OS,2) - 1);
y      = NaN(8,1);
for (i = 1:4), index3 = index3 && OS(i,1:end - 1) == p(i); end
exitFlag = any(index3);
if (exitFlag), y = OS(:,find(index3)); %#ok<FNDSB>
elseif ((p(1) <= p(2)+3) || (p(2) < 1) || (p(3) <= p(4)+3) || (p(4) < 1))
    exitFlag = 2;
    y(1:4)   = p;% segment > slack OR slack < 1
else
    
    K                 = size(X,1);
    [normX,normXcorr] = deal(NaN(K,1));
    for (k = 1:K), normX(k) = norm(X(k,:),2); end
    try
        [~,X,diagnos] = cow(ref,X,p(1),p(2),[0 1 0 losange 0]);
        for (k = 1:K), normXcorr(k) = norm(X(k,:),2); end
    catch
        for (k = 1:K)
            [~,X(k,:),diagnos] = cow(ref,X(k,:),p(1),p(2),[0 1 0 losange 0]);
            normXcorr(k)       = norm(X(k,:),2);
        end
    end
    y(1)   = diagnos.segment_length(1,1)+1;
    y(2)   = diagnos.slack;
    y(5:6) = sum(svd(X'/sqrt(sum(normXcorr.^2))),K).^4;
    y(7)   = mean((1 - min(abs((normXcorr - normX)./normX),1))^2);
    y(5)   = y(3) + y(5);
    y(8)   = etime(clock,t0);
    
end

end

function showIter(id,iter,total,varargin)
switch id
    case 0, fprintf(1,'run %i/%i %s\n',iter,total,itStr(varargin{1}));
    case 1, fprintf(1,'run %i/%i: - min (segment/slack combination was already computed)\n',iter,total);
    case 2, fprintf(1,'[\brun %i/%i: - illegal segment/slack combination)]\b\n',iter,total);
    case 3, fprintf(1,'Starting optimization %i/%i %s\n',iter,total,itStr(varargin{1}(1:7)));
    case 4, fprintf(1,'Early termination after %i steps!',iter);
    case 5, fprintf(1,'   Optimization %i/%i terminated in %i steps %s',iter,total,varargin{2},itStr(cat(1,varargin{1}(1:7),varargin{3})));
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
end
warning(sprintf('optim2DCOW:%s',id),msg)

end

function [n,M,c] = newPoint(L,X,s)
[M,c]      = sort(L);
M          = M(s);
c          = c(s);
n          = X(:,c);
D          = size(X,1);
cL         = sum(X < n,2);
if (cL < D), cL = 1; end
cU         = sum(X > n,2);
if (cU < D), cU = 1; end
n          = n + cU - cL;

end

function plotOptim(OS,optim,eTime,S)
    optim_pars = OS(1:4,optim); 
    f = figure;
    stem3(OS(1,:),OS(2,:),OS(3,:),'filled');
    xlabel('Segment length'); ylabel('Slack size'); zlabel('Warping Effect');
    s = ['Warping Effect(' num2str(optim_pars(1)) ',' num2str(optim_pars(2)) ')'];
    text(OS(1,optim),OS(2,optim),OS(3,optim),s);
    
    figure;
    stem3(OS(1,:),OS(2,:),OS(4,:),'filled');
    hold on;
    
    stem3(OS(1,:),OS(2,:),ones(size((OS(4,:))))*S);
    hold off;
    title('(Solid = Simplicity, Open = Simplicity Raw Data)');
    xlabel('Segment length'); ylabel('Slack size'); zlabel('Simplicity');
    [~,b] = max(OS(4,:));
    s = ['Simplicity(' num2str(OS(1,b)) ',' num2str(OS(2,b)) ')'];
    text(OS(1,b),OS(2,b),OS(4,b),s);
    s = ['Warping Effect(' num2str(optim_pars(1)) ',' num2str(optim_pars(2)) ')'];
    text(OS(1,optim),OS(2,optim),OS(4,optim),s);
    
    figure;
    stem3(OS(1,:),OS(2,:),OS(5,:),'filled');
    xlabel('Segment length'); ylabel('Slack size'); zlabel('Peak Factor');
    [~,b] = max(OS(5,:));
    s = ['Peak Factor(' num2str(OS(1,b)) ',' num2str(OS(2,b)) ')'];
    text(OS(1,b),OS(2,b),OS(5,b),s);
    s = ['Warping Effect(' num2str(optim_pars(1)) ',' num2str(optim_pars(2)) ')'];
    text(OS(1,optim),OS(2,optim),OS(5,optim),s);
    figure(f);
    showIter(7,NaN,NaN,OS(:,optim),eTime)
    
end