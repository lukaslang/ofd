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

% This script sets up the paths of the libraries and adds all subfolders.

% Set library path.
libraryPath = 'Z:\libraries\';

% Export Figure is required for saving figures.
addpath(genpath(fullfile(libraryPath, 'export_fig\')));

% Add all subfolders.
y = dir('.');
y = y([y.isdir]);
y = y(~cellfun(@(x) strcmp(x, '.git') || strcmp(x, '.') || strcmp(x, '..') || strcmp(x, 'results') || strcmp(x, 'renderings'), {y.name}));
% Add to path.
cellfun(@(x) addpath(genpath(fullfile(pwd, x))), {y.name});

% Clean up.
clear y;
clear libraryPath;