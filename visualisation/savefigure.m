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
function savefigure(h, file, format, density, quality)
%SAVEFIGURE Saves a figure.
%
%   SAVEFIGURE(h, file) takes a figure handle h and a file and saves the
%   figure.

if(nargin == 2)
    export_fig(file, '-png', '-r300', '-zbuffer', '-transparent', h);    
elseif(nargin == 4)
    export_fig(file, format, density, '-zbuffer', '-transparent', h);
else
    export_fig(file, format, quality, density, '-zbuffer', '-transparent', h);
end