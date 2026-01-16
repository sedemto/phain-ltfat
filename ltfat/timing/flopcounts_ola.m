function [fcfacola,fcola]=flopcounts(a,M,L,Lg,Lb)
%
%   Compute the flop count for the factorization algorithm.
%
%   Url: http://ltfat.github.io/doc/timing/flopcounts_ola.html

% Copyright (C) 2005-2023 Peter L. Soendergaard <peter@sonderport.dk> and others.
% This file is part of LTFAT version 2.6.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
N=L/a;

[c,h_a,h_m]=gcd(a,M);
h_a=-h_a;
p=a/c;
q=M/c;
d=N/q;

Nb=L/Lb;
rho=(Lb+Lg)/Lb;
  
fcfacola  = rho*(8*L*q+4*L*(1+q/p)*log2(d*rho)+4*M*N*log2(M));

fcola     = rho*(8*L*M+4*L*log2(rho*Lb)+4*M*N*log2(rho*N));



