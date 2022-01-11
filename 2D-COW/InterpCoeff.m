function [Coeff,Index] = InterpCoeff(n,nprime,offs)
% Function to calculate coefficients for interpolation
p     = length(nprime);
q     = n - 1;
Coeff = zeros(p,n);
Index = zeros(p,n);
for i_p = 1:p
   
   pp                  = 1:nprime(i_p);
   p                   = (0:q) * (nprime(i_p) - 1)/q + 1;
   [~,k]               = histc(p,pp);
   k(p < 1)            = 1;
   k(p >= nprime(i_p)) = nprime(i_p) - 1;
   Coeff(i_p,:)        = (p - pp(k));
   Index(i_p,:)        = k - offs(i_p);

end

