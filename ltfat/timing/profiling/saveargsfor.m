function saveargsfor(matfile,args,varargin)
%SAVEARGSFOR Save input arguments for mexExecuter
%
%   Url: http://ltfat.github.io/doc/timing/profiling/saveargsfor.html

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

if nargin<2
   error('%s: Too few input arguments.',upper(mfilename) );
end

if isempty(args)
   error('%s: No variables to be saved.',upper(mfilename) );
end

definput.keyvals.res=[];
[flags,kv]=ltfatarghelper({},definput,varargin);


varNo = numel(args);
varNames = cell(varNo,1);


for ii=1:numel(args)
   varNames{ii} = sprintf('arg%i',ii-1);
   eval([varNames{ii},' = args{ii};']);   
end

if ~isempty(kv.res)
   varNames{end+1} = 'res';
   eval([varNames{end},' = kv.res;']);
end

save(matfile,varNames{:});

