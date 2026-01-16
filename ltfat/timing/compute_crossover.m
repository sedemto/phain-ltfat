if 1
  % This is the setup used in the paper
  a=40;
  M=60; 
  L=a*M;
  W=4;
  nrep=50;
else
  a=50;
  M=200; 
  L=a*M; 
  W=4;
  nrep=50;  
  
end;

system('rm crossover.log');
%for gl=M:M:20*M
%
%   Url: http://ltfat.github.io/doc/timing/compute_crossover.html

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
for gl=M:M:16*M
  s=sprintf('./time_dgt_fb %i %i %i %i %i %i >> crossover.log\n',a,M,L,W,gl,nrep);
  
  disp(s);
  system(s);
end;

s=sprintf('./time_dgt_fac %i %i %i %i %i > crossover.ref\n',a,M,L,W,nrep);

disp(s);
system(s);

system('rm crossover_real.log');
%for gl=M:M:20*M
for gl=M:M:16*M
  s=sprintf('./time_dgtreal_fb %i %i %i %i %i %i >> crossover_real.log\n',a,M,L,W,gl,nrep);
  
  disp(s);
  system(s);
end;

s=sprintf('./time_dgtreal_fac %i %i %i %i %i > crossover_real.ref\n',a,M,L,W,nrep);

disp(s);
system(s);


