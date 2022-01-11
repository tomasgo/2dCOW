function [LenSeg,nSeg,predefB] = defSegments(Seg,pT,pX,splitRem,doShow)
BADENDPTS  = 'End points must be equal to 1 and to the length of the pattern/target';
TWOPTS     = 'Segments must contain at least two points';
OFFBLEN    = 'Segment length is %i than for signals of length %i';
DIFFNSEG   = 'Target and signal do not have the same number of segments';

predefB = length(Seg) > 1; % True if segment boundaries are predefined
Seg     = round(Seg);      % Only integers are currently allowed as segment boundaries
if (predefB) 
   if (~isequal(Seg(1,1),Seg(2,1),1) || Seg(1,end) ~= pT || Seg(2,end) ~= pX), error('COW:badEndpoints',BADENDPTS); end
   LenSeg = diff(Seg,1,2);    % LenSeg(1,:): Length of the segments in the - 1
   if (any(LenSeg < 2)), error('COW:TwoPts',TWOPTS); end
   nSeg   = size(LenSeg,2);   % nSeg: number of segments
else
   
    tmp = min(pX,pT);
    if (Seg > tmp), error('COW:OffBLen',OFFBLEN,Seg,tmp); end
    if (splitRem) % Segments in the signals can have different length from those in the target
        nSeg             = floor((pT - 1)/Seg);
        LenSeg(1,1:nSeg) = floor((pT - 1)/nSeg);
        LenSeg(2,1:nSeg) = floor((pX - 1)/nSeg);
        fprintf('\n Segment length adjusted to best cover the remainders')
    else
        nSeg               = floor((pT - 1) / (Seg - 1));
        LenSeg(1:2,1:nSeg) = Seg - 1;
        if (floor((pX - 1) / (Seg - 1)) ~= nSeg), error('COW:diffNSeg',DIFFNSEG); end
    end
    temp = rem(pT - 1,LenSeg(1,1)); % The remainders are attached to the last segment in the target and in the reference
    if (temp > 0)
        LenSeg(1,nSeg) = LenSeg(1,nSeg) + temp;
        if (doShow), fprintf('\n Segments: %i points x %i segments + %i (target)',LenSeg(1,1) + 1,nSeg - 1,LenSeg(1,end) + 1); end
    elseif (doShow), fprintf('\n Segments: %i points x %i segments (target)',LenSeg(2,1) + 1,nSeg);
    end
    temp = rem(pX - 1,LenSeg(2,1));
    if (temp > 0) 
        LenSeg(2,nSeg) = LenSeg(2,nSeg) + temp;
        if (doShow), fprintf('\n           %i points x %i segments + %i (signals)\n',LenSeg(2,1) + 1,nSeg - 1,LenSeg(2,end) + 1); end
    elseif (doShow), fprintf('\n           %i points x %i segments (signals)\n',LenSeg(2,1) + 1,nSeg);
    end

end
