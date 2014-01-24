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
function f = createFigure(cmap, xmin, xmax, ymin, ymax)
%CREATEFIGURE Creates a 2D figure.
%
% f = CREATEFIGURE opens a figure with default settings.
%
% f = CREATEFIGURE(cmap) creates a figure with the given colormap cmap.
%
% f = CREATEFIGURE(cmap, xmin, xmax, ymin, ymax) sets axis limits
% accordingly.
%
% Note that CREATEFIGURE uses zbuffer as renderer!
f = figure('Renderer', 'zbuffer');
axis square;
daspect([1, 1, 1]);
hold on;
% Set colormap.
if(nargin > 0)
   colormap(cmap);
end
if(nargin == 5)
    set(gca, 'XLim', [xmin xmax]);
    set(gca, 'YLim', [ymin ymax]);
end
end