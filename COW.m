function [Warping,XWarped,Diagnos] = COW(T,X,Seg,Slack,Options)
% function [Warp,XWarped,Diagnos] = COW(T,X,Seg,Slack,Options);
% Correlation Optimized Warping function with linear interpolation
% Giorgio Tomasi / Frans van den Berg 090411 (GT)
%
% INPUT
% T (1 x nt x ...) : target. "..." denotes the channels
% X (mP x nP x ...): data to alignm for mP row vectors of length nP to be warped/corrected
% Seg (1 x 1)      : segment length (consequently: number of segments N = floor(nP/Seg))
%  or (2 x N+1)      matrix with segment (pre-determined) boundary-points
%                    1st row: indices in "T", must start with 1 and end with "nt"
%                    2nd row: indices in "xP", must start with 1 and end with "nP"
% Slack (1 x 1)      'slack' - maximum range or degree of warping in segment length "m"
% Options (1 x 6)  : options. NaN = Default is used
%                    1 : triggers plot and progress-text (default: 0). NB. Only last row/object in "xP" is plotted)
%                    2 : correlation power (default: 1). It must be comprised between 1 and 4.
%                    3 : force equal segment lengths in "T" and "X" instead of filling up "T" with N boundary-points
%                        (notice that different number of boundaries in "T" and "X" will generate an error)
%                    4 : fix maximum correction to + or - options(4) points from the diagonal
%                    5 : save in "diagnos" the table with the optimal values of loss function and predecessor (default: 0)
%                        (memory consuming for large problems)
%                        Row 1: boundary position
%                        Row 2: optimal cost function up to the current node given the constraints (i.e., Seg, Slack and Options(4))
%                        Row 3: pointer to the predecessor leading to the optimal cost function up to the current node.
%                    6 : 1 -> calculate all optimal warping paths for all slacks between 1 and Slack (default: 0)
% 
% OUTPUT
% Warp (mP x N x 2): interpolation segment starting points (in "nP"  units) after warping (first slab) and before warping (second slab)
%                    (difference of the two = alignment by repositioning segment boundaries; useful for comparing correction in different/new objects/samples)
% XWarped (mP x nt): corrected signals (from "xP" warped to mach "xt")
% Diagnos (struct) : warping diagnostics: options, segment, slack, 
%                    index in target ("xt", "warping" is shift compared to this) and sample ("xP"), search range in "xP", computation time
%                    (note: diagnostics are only saved for one - the last - signal in "xP")
%
% based on: Niels-Peter Vest Nielsen, Jens Micheal Carstensen and Jørn Smedegaard 'Aligning of singel and multiple
%           wavelength chromatographic profiles for chemometric data analysis using correlation optimised warping' 
%           J. Chrom. A 805(1998)17-35
%
% Reference: Correlation optimized warping and dynamic time warping as preprocessing methods for chromatographic Data
%            Giorgio Tomasi, Frans van den Berg and Claus Andersson, Journal of Chemometrics 18(2004)231-241
%
% Authors: 
% Giorgio Tomasi / Frans van den Berg
% Royal Agricultural and Veterinary University - Department of Food Science
% Quality and Technology - Spectroscopy and Chemometrics group - Denmark
% email: gto@life.ku.dk / fb@life.ku.dk - www.models.life.ku.dk

%% ERRORS
CORRPOWER  = '"Options(2)" (correlation power) must be in the range 1:4';
MISSING    = 'COW cannot handle missing values';
SLACKGTSEG = 'Slack is larger than the length of the segments';

%% Check Input values
if (nargin < 4), help cow; return; end
% X and T
if any(isnan(T(:))) || any(isnan(X(:))), error('COW:Missing',MISSING); end


defOptions = [0 1 0 0 0 0];
if (nargin < 5 || isempty(Options)), Options = defOptions;
else
    Options(end + 1:length(defOptions)) = defOptions(length(Options) + 1:end);
    Options(isnan(Options))             = defOptions(isnan(Options));
    if (Options(2) < 1) || (Options(2) > 4), error('COW:corrPower',CORRPOWER); end
end
doShow    = Options(1);
allSlacks = Options(6);

%% Initialise
dimX    = size(X);
nX      = dimX(1);                      % nX     : number of signals that are to be aligned
pX      = dimX(2);                      % pX     : number of data points in each signal
pT      = size(T,2);                    % pT     : number of data points in the target
XWarped = zeros([nX,pT,dimX(3:end)]);   % XWarped: initialise matrix of warped signals

%% Initialise segments
[LenSeg,nSeg,predefB] = defSegments(Seg,pT,pX,Options(3),doShow);
if any(LenSeg(:) <= Slack + 1), error('COW:slackGTSeg',SLACKGTSEG); end % Two points are the minimum required for linear interpolation
bT      = cumsum([1,LenSeg(1,:)]);
bP      = cumsum([1,LenSeg(2,:)]);

%% Check slack
if (length(Slack) > 1) % Different slacks for the segment boundaries will be implemented
   if (size(Slack,2) <= nSeg), error('The number of slack parameters is not equal to the number of optimised segments'); end
   fprintf('\n Multiple slacks have not been implemented yet')
   return
end
slackVec = -Slack:Slack;                     % All possible slacks for a segment boundary

%% Set feasible points for boundaries
% Slope Constraints
Bounds = defBounds(bT,bP,Slack,nSeg,Options(4),allSlacks);

%% Calculate first derivatives for interpolation
Tdiff = diff(T,1,2);

%% Calculate coefficients and indexes for interpolation
Int_Coeff = cell(nSeg,1);
Int_Index = Int_Coeff;
if ~predefB
   [A,B]                             = InterpCoeff(LenSeg(2,1) + 1,LenSeg(1,1) + slackVec + 1,slackVec); 
   [Int_Coeff{1:nSeg - 1}]           = deal(A);
   [Int_Index{1:nSeg - 1}]           = deal(B);
   [Int_Coeff{nSeg},Int_Index{nSeg}] = InterpCoeff(LenSeg(2,nSeg) + 1,LenSeg(1,nSeg) + slackVec + 1,slackVec);   
else
   for (i_seg = 1:nSeg), [Int_Coeff{i_seg},Int_Index{i_seg}] = InterpCoeff(LenSeg(2,i_seg) + 1,LenSeg(1,i_seg) + slackVec + 1,slackVec); end
end

%% Dynamic Programming Section
Table_Index    = cumsum([0,diff(Bounds(:,:,1)) + 1],2);     % Indexes for the first node (boundary point) of each segment in Table
nS             = size(Bounds,3);
Table          = zeros(1 + 2 * nS,Table_Index(nSeg + 2),nX);% Table: each column refer to a node
                                                            %        (1,i) position of the boundary point in the signal
                                                            %        (2,i) optimal
                                                            %        value of the loss function up to node (i)
                                                            %        (3,i) pointer to optimal preceding node (in Table)
Table(2,2:end,1:nX) = -Inf;                                 % All loss function values apart from node (1) are set to -Inf

if (nS > 1), tableFlag = false(nS - 1,size(Table,2)); end
for i_seg = 1:nSeg + 1                                      % Initialise Table
   v              = (Bounds(1,i_seg,1):Bounds(2,i_seg,1))';
   ind            = Table_Index(i_seg) + 1:Table_Index(i_seg + 1);
   Table(1,ind,:) = v(:,ones(nX,1));
   for (s = 2:nS), tableFlag(s - 1,ind) = v >= Bounds(1,i_seg,s) & v <= Bounds(2,i_seg,s); end
end
if (nS > 1), tableFlag = sum(tableFlag,1) + 1; else, tableFlag = ones(1,size(Table,2)); end
   
tic
% Reshape to allow for treatment of second order signals
% If X is a matrix and T a row vector, this is equivalent to a
% transposition for T.
T     = reshape(permute(T,[2 1 3:length(dimX)]),dimX(2),prod(dimX(3:end)));
Tdiff = reshape(permute(Tdiff,[2 1 3:length(dimX)]),dimX(2) - 1,prod(dimX(3:end)));
X     = reshape(permute(X,[1 3:length(dimX) 2]),[dimX(1),prod(dimX(3:end)),dimX(2)]);
d     = prod(dimX(3:end));
D     = true(2 * Slack + 1,nS);
for (i = 2:nS), D([1:(i - 1),end - i + 2:end],i) = false; end
D = D(:,[2:end,1]);

% Forward phase
for i_seg = 1:nSeg                             % Loop over segments

   a             = slackVec + LenSeg(1,i_seg);               % a,b,c: auxiliary values that depend only on segment number and not node
   b             = Table_Index(i_seg) + 1 - Bounds(1,i_seg,:);
   c             = LenSeg(2,i_seg) + 1;
   Count         = 1;                                          % Counter for local table for segment i_seg
   Node_Z        = Table_Index(i_seg + 2);                     % Last node for segment i_seg
   Node_A        = Table_Index(i_seg + 1) + 1;                 % First node for segment i_seg
   Bound_k_Table = zeros(nS * 2,Node_Z - Node_A + 1,nX);       % Initialise local table for boundary
   Int_Index_Seg = Int_Index{i_seg}' - (LenSeg(1,i_seg) + 1);  % Indexes for interpolation of segment i_seg
   Int_Coeff_Seg = Int_Coeff{i_seg}';                          % Coefficients for interpolation of segment i_seg
   dimX_Seg      = [dimX(1),(bP(i_seg + 1) - bP(i_seg) + 1) * d];
   XSeg          = X(:,:,bP(i_seg):bP(i_seg + 1));                  % Segments i_seg of signals X
   if isequal(XSeg(:,:,ones(bP(i_seg + 1) - bP(i_seg) + 1,1)),XSeg) % Constants
      XSeg_centred  = zeros(dimX_Seg);
   else
      XSeg          = reshape(XSeg,dimX_Seg);                     % Segments i_seg of signals X
      XSeg_mean     = sum(XSeg,2)/size(XSeg,2);                   % Centred XSeg (for correlation coefficients)
      XSeg_centred  = XSeg - XSeg_mean(:,ones(size(XSeg,2),1));
   end
   Norm_XSeg_cen = sqrt(sum(XSeg_centred.^2,2));               % (n - 1) * standard deviation of TSeg
   for i_node = Node_A:Node_Z                                  % Loop over nodes (i.e. possible boundary positions) for segment i_seg

      Prec_Nodes         = Table(1,i_node) - a;                                               % Possible predecessors given the allowed segment lengths
      Allowed_Arcs       = Prec_Nodes >= Bounds(1,i_seg,1) & Prec_Nodes <= Bounds(2,i_seg,1); % Arcs allowed by local and global constraints
      iAA                = Allowed_Arcs;
      Nodes_TablePointer = b(1) + Prec_Nodes(Allowed_Arcs);                                   % Pointer to predecessors in Table
      N_AA               = sum(Allowed_Arcs);                                                 % Number of allowed arcs
      Dloc               = ':';
      if N_AA > 0 % Sometimes boundaries are ineffective and few nodes are allowed that cannot be reached
              % It has to be further investigated

         Index_Node                    = Table(1,i_node) + Int_Index_Seg(:,Allowed_Arcs);                % Interpolation signal indexes for all the allowed arcs for node i_node
         Coeff_b                       = Int_Coeff_Seg(:,Allowed_Arcs);                                  % Interpolation coefficients for all the allowed arcs for node i_node
         Coeff_b                       = Coeff_b(:);
         Ti_Seg                        = T(Index_Node,:);
         Coeff_b                       = Coeff_b(:,ones(size(Ti_Seg,2),1));
         Ti_diff                       = Tdiff(Index_Node,:);
         Ti_Seg                        = reshape((Ti_Seg + Coeff_b .* Ti_diff)',c * d,N_AA);                   % Interpolate for all allowed predecessors
         Ti_Seg_mean                   = sum(Ti_Seg)/size(Ti_Seg,1);                                     % Means of the interpolated segments
         Norm_Ti_Seg_cen               = sqrt(sum(Ti_Seg.^2) - size(Ti_Seg,1) * Ti_Seg_mean.^2);         % Fast method for calculating the covariance of T and X (no centering of X is needed)
         CCs_Node                      = ((XSeg_centred * Ti_Seg)./(Norm_XSeg_cen * Norm_Ti_Seg_cen))';  % Correlation coefficients relative to all possible predecessors
         CCs_Node(~isfinite(CCs_Node)) = 0;                                                              % If standard deviation is zero, update is not chosen
         if (Options(2) == 1), CCs_Node = CCs_Node.^Options(2); end
         for (k = 1:tableFlag(i_node))
             
             f                            = k * 2;
             Cost_Fun                     = reshape(Table(f,Nodes_TablePointer,:),N_AA,nX) + CCs_Node(Dloc,:); % Optimal value of loss function from all predecessors
             [ind,pos]                    = max(Cost_Fun,[],1);
             Bound_k_Table(f - 1,Count,:) = ind;
             Bound_k_Table(f,Count,:)     = Nodes_TablePointer(pos);                                           % Pointer to optimal predecessor
             if (k < tableFlag(i_node))
                 
                 Prec_Nodes         = Prec_Nodes(2:end - 1);
                 Allowed_Arcs       = Prec_Nodes >= Bounds(1,i_seg,k + 1) & Prec_Nodes <= Bounds(2,i_seg,k + 1);
                 N_AA               = nnz(Allowed_Arcs);
                 if (N_AA == 0), break; end
                 Dloc               = D(:,k);
                 Dloc(Dloc)         = Allowed_Arcs;
                 Dloc               = Dloc(iAA);
                 Nodes_TablePointer = b(1) + Prec_Nodes(Allowed_Arcs);
             
             end
             
         end
         Count = Count + 1;

      end

   end % i_node
   Table(2:end,Node_A:Node_Z,:) = Bound_k_Table;       % Update general table (it turned out to be faster than using Table directly in the loop over nodes

end % i_seg
Time = toc;

%% Backward phase
Warping               = zeros(nX,nSeg + 1,nS);
Warping(:,nSeg + 1,:) = bT(end);
for (k = 1:nS)
    
    for (i = 1:nX)                                  % Loop over samples/signals
        
        Pointer = size(Table,2);           % Backtrace optimal boundaries using the pointers in Table
        for (j = nSeg:-1:1)
            Pointer        = Table(1 + 2 * k,Pointer,i);
            Warping(i,j,k) = Table(1,Pointer,i);
        end
        
    end
    
end

Warping = cat(3,bP(ones(nX,1),:),Warping);
% Warping          = Warping(:,:,[2 1]);
% Warping(:,end,2) = pT;

%% Output
if (nargout > 1)  % Reconstruct aligned signals

    X       = reshape(permute(X,[1 3 2]),dimX);
    XWarped = apply1DWarp(X,Warping(:,:,[1 end]),2);
    if nargout > 2    % Save some diagnostics if requested
        Diagnos = struct('indexP',bP,'indexT',bT,'Nsegments',nSeg,'options',Options,'rangeP',permute(Bounds,[2 1,3]),...
            'segment_length',LenSeg,'slack',Slack,'table',[],'time',Time);
        if (Options(5)), Diagnos.table = Table; end
    end
    
end

%% Plot
if (doShow)
    
    LEGPROP = {'Location','SouthOutside','Orientation','horizontal'}; 
    AXPROP  = {'NextPlot','add','XGrid','on','YGrid','on','Box','on'};
    figure('Name',sprintf('COW - sample %i',nX))
    yTarget = reshape(T,1,pT,[]);
    ySample = reshape(X(nX,:),1,pX,[]);
    if (length(dimX) == 2)
        
        YL         = [min([T(:)',X(nX,:)]) max([T(:)' X(nX,:)])] * [1.005 -0.005;-0.005 1.005];
        minmaxaxis = [1 max([pT pX]) YL];
        ah         = subplot(2,1,1,AXPROP{:});
        ph         = gobjects(4,1);
        ph(1)      = plot(1:pT,T,'Color','b','DisplayName','Target');
        ph(2)      = plot(bT,T(bT),'Color','b','Marker','+','LineStyle','none','DisplayName','Target bound.');
        ph(3)      = plot(1:pX,X(nX,:),'Color',[0.1 0.5 0.1],'DisplayName','Sample');
        ph(4)      = plot(bP,X(nX,bP),'Color',[0.1 0.5 0.1],'Marker','+','LineStyle','none','DisplayName','Sample bound.');
        for a = 2:length(Warping(nX,:,1))
            mh = plot([bT(a) Warping(nX,a,1)],[T(Warping(nX,a,2)) T(Warping(nX,a,2))],'Color','r','DisplayName','Correction');
            if (Warping(nX,a,2) > Warping(nX,a,1)), plot(Warping(nX,a,2),T(Warping(nX,a,2)),'>r');
            else plot(Warping(nX,a,2),T(Warping(nX,a,2)),'<r');
            end
        end
        axis(minmaxaxis)
        ah.Title.String = sprintf('Raw data (Sample #%i/%i)',nX,nX);
        
    end
    lh = legend([ph;mh],LEGPROP{:});
    lh.FontSize = 7;
    ah = subplot(2,1,2,AXPROP{:});
    ph = plot(1:pT,T,'b',1:pT,XWarped(nX,:),'g');
    set(ph,{'DisplayName'},{'Target','Aligned sample'}')
    set(legend(ph,LEGPROP{:}),'FontSize',7)
    axis(minmaxaxis);
    ah.Title.String = 'Aligned sample';

end
