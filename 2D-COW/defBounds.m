function Bounds = defBounds(bT,bP,Slack,nSeg,Band,allSlacks)
if (allSlacks)
    s       = (-Slack:Slack)';
    [~,ord] = sort(abs(s),'descend');
    s       = s(ord(1:end - 1));
else s = [-Slack;Slack];
end
nSlacks = length(s)/2;
Bounds  = ones(2,nSeg + 1,nSlacks);
offs    = s * (0:nSeg);
if (verLessThan('matlab','9.3'))
    Bounds_a    = bsxfun(@plus,bT(:,1:nSeg + 1),offs);
    Bounds_b    = bsxfun(@plus,bT(:,1:nSeg + 1),offs(:,nSeg + 1:-1:1));
else
    Bounds_a = bT(:,1:nSeg + 1) + offs;
    Bounds_b = bT(:,1:nSeg + 1) + offs(:,nSeg + 1:-1:1);
end
for (i = 1:nSlacks)
    k             = (i - 1) * 2 + 1;
    Bounds(1,:,i) = max(Bounds_a(k,:),    Bounds_b(k,:));
    Bounds(2,:,i) = min(Bounds_a(k + 1,:),Bounds_b(k + 1,:));
end
pT = bT(end);
pX = bP(end);

% Band Constraints
if (Band)
    
    if (abs(pT - pX) > (Band)), error('The band is too narrow and proper correction is not possible');end
    lowBound  = max(0,pT/pX * bP - Band);
    highBound = min(pT,pT/pX * bP + Band);
    if (verLessThan('matlab','9.3'))
        Bounds(1,:,:) = bsxfun(@max,Bounds(1,:,:),lowBound);
        Bounds(2,:,:) = bsxfun(@min,Bounds(2,:,:),highBound);
    else
        Bounds(1,:,:) = max(Bounds(1,:,:),lowBound);
        Bounds(2,:,:) = min(Bounds(2,:,:),highBound);
    end
    if (any(any(diff(Bounds < 0,1,1)))), error('The band is incompatible with the fixed boundaries'); end
    
end

