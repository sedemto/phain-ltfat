%COMPUTE_LONGER  Vary the length of the transform
%
%  This script computes the running time for longer and longer
%  transforms. Use the script plot_longer to visualize the result.
%
%  All other parameters except L remain fixed. The window length for the
%  filter bank algorithm is also kept fixed.
%
%
%   Url: http://ltfat.github.io/doc/timing/compute_much_longer.html

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
a=40;
M=60; 
W=4;
nrep=20;
gl=2400;
bl=24000;

% Get test sizes. Use only test sizes from nextfastfft, as the others
% cause to large a variability.
[dummy,table]=nextfastfft(1);

system('rm much_longer_fb_real.log');
system('rm much_longer_fac_real.log');
system('rm much_longer_ola_real.log');

testrange=bl*table(find(table<=30));

for ii=1:length(testrange)
  L=testrange(ii);
    
  s=sprintf('./time_dgtreal_fb %i %i %i %i %i %i >> much_longer_fb_real.log\n',a,M,L,W, ...
            gl,nrep);
  disp(s);
  system(s);

  s=sprintf('./time_dgtreal_fac %i %i %i %i %i >> much_longer_fac_real.log\n',a,M,L,W,nrep);  
  disp(s);
  system(s);

  % Extend L to multiple of bl for the OLA algorithm.
  Lola = ceil(L/bl)*bl;

  s=sprintf('./time_dgtreal_ola %i %i %i %i %i %i %i >> much_longer_ola_real.log\n',a,M,Lola,W,gl,bl,nrep);  
  disp(s);
  system(s);

end;





