%TIME_INTOVERSAMP
%
%   The purpose of this test is to determine whether special code for
%   handling the integer oversampling case in the factorization routines
%   pays off in any way.
%
%   Url: http://ltfat.github.io/doc/timing/time_intoversamp.html

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


Lr=[1048576*2,2400];
ar=[      512,   2];
Mr=[     1024, 800];

for ii=1:length(Lr)

  L=Lr(ii);
  
  M=Mr(ii);
  a=ar(ii); 
  
  [L, a, M]
  
  N=L/a;
  c=gcd(a,M);
  p=a/c;
  q=M/c;
  d=N/q;

  f=rand(L,1);
  gf=rand(p*q,c*d);  
  c1=mex_dgt_fac_1(f,gf,a,M);

  f=rand(L,1);
  gf=rand(p*q,c*d);
  c2=mex_dgt_fac_2(f,gf,a,M);

end;


