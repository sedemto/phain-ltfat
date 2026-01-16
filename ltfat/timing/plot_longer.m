lw=2;
fz=20;

data_fb_real =load('longer_fb_real.log');
data_fac_real=load('longer_fac_real.log');
data_ola_real=load('longer_ola_real.log');

% Columns in data, fb : a M L W gl time
% Columns in data, fac: a M L W time
%
%   Url: http://ltfat.github.io/doc/timing/plot_longer.html

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
Ls=data_fac_real(:,3);
t_fb_real =data_fb_real(:,6);
t_fac_real=data_fac_real(:,5);
t_ola_real=data_ola_real(:,7);

if 0
  % Color legend
  l1='b';
  l2='b--';
  l3='r';
  l4='r--';
else
  % bw legend
  l1='b';
  l2='b--';
  l3='b-.';
  l4='b:';
end;

figure(1);

a=data_fac_real(1,1);
M=data_fac_real(1,2);
L=data_fac_real(:,3);

plot(Ls,t_fb_real,l1,...  
     Ls,t_fac_real,l2,...
     Ls,t_ola_real,l3,'LineWidth',lw);
set(gca,'Fontsize',fz);

h=legend('Portnoff, real','Fac, real','Fac-OLA, real',...
       'Location','NorthWest');
% Grow the box a little, otherwise the export to .eps is messed up.
q=get(h,'Position');
set(h,'Position',[q(1)*0.90 q(2)*.95 q(3)*1.8 q(4)]);

xlabel('Signal length / samples','Fontsize',fz);
ylabel('Running time / seconds','Fontsize',fz);

print -deps plot_longer_1.eps




