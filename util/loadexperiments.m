% Copyright 2013 Clemens Kirisits and Lukas Lang
%
% This file is part of OFD.
%
%    OFD is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    OFD is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with OFD.  If not, see <http://www.gnu.org/licenses/>.
function E = loadexperiments(path)
%LOADEXPERIMENTS Takes a path and loads all mat files into a cell array.

files = getFiles(path);

% Load experiments.
E = cell(length(files), 1);
for k=1:length(files)
    filename = files(k).name;
    fprintf('Loading run %d/%d: %s\n', k, length(files), filename);
    E{k} = load(fullfile(path, filename));
end

end