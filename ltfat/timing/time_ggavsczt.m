clear all;

L = [300,1024,2000];
K = 1:5:100;
rep = 1000;


commands = {'time_ggareal','time_cztreal','time_cztreal_fact','time_ggareal_noplan','time_cztreal_nextpow2','time_cztreal_fact_nextpow2'};
abbrv = {'gga','czt_{nextfastfft}','czt(fact)_{nextfastfft}','gga(np)','czt_{nextpow2}','czt(fact)_{nextpow2}'};
times = cell(numel(commands),1);
for c=1:numel(commands)
   times{c} = zeros(numel(K),numel(L));
end

for l=1:numel(L)
   
 for k=1:numel(K)  
    for c =1:numel(commands)
       cmd = sprintf('%s %i %i %i',commands{c},L(l),K(k),rep);
       [status,out] = system(cmd);
       if ~isempty(out(strfind(out,'is not recognized')))  || status~=0 || isempty(out) || isempty(out(strfind(out,', ')))
          error('%s: Call to\n%s\nfailed with the following error message:\n%s',upper(mfilename),commands{c},out);
       end
       times{c}(k,l) = str2double(out(strfind(out,', ')+2:end));
    end
 end
 
figure(l); 
timesplot = zeros(numel(K),numel(commands));
for c = 1:numel(commands)
   timesplot(:,c) = times{c}(:,l);
end


h=plot(K,timesplot);
legend(abbrv);
xlabel('K\rightarrow');
ylabel('t[ms]');
title(sprintf('L=%i',L(l)));
set(gcf, 'Color', 'w');
set(h,'LineWidth',2);
%export_fig(sprintf('ggavscztL%i.png',L(l)))
%
%   Url: http://ltfat.github.io/doc/timing/time_ggavsczt.html

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
end


