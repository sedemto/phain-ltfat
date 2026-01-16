load time_fftw.log

x=time_fftw(:,1);
r=zeros(65536,1);
r(1:49999)=time_fftw(:,2);

L=size(r,1);

x=(1:L).';

x2=2.^nextpow2(x);
[x3,table]=nextfastfft(x);
r2=r(x2);
r3=r(x3);

plot(x,r,x,r2,x,r3)


% Repetitions for an accurate result
%
%   Url: http://ltfat.github.io/doc/timing/plot_nextfastfft.html

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
nrep=50;

% Time the FFT.
N=numel(r);
res=zeros(N,1);
for ii=1:N
  f=randn(ii,nrep);
  tic;    
  output=fft(f);
  res(ii)=toc/nrep;
end;

% Remap the data, so we only consider the next higher sizes.
res2=res(r2-range_min+1);
res3=res(r3-range_min+1);

% Plot the running time for all sizes.
figure(2);
plot(r,res,r,res2,r,res3);
xlabel('Problem size (samples)');
ylabel('Running time (seconds)');


figure(3);
plot(r,res2,r,res3)



