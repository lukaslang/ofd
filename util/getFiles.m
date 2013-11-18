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
function files = getFiles(path)
%GETFILES Returns a struct of .mat files in a dir.
%
%   files = GETFILES(path) takes a path and returns a struct of files.

% Get list of files and match for files matching '.mat'.
files = dir(path);
files = files([files.isdir] == 0);
files = files(cellfun(@(x) ~isempty(regexp(x, '\.mat$', 'once')), {files.name}));

end