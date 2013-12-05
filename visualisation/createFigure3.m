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
function f = createFigure3(cmap, xmin, xmax, ymin, ymax, zmin, zmax)
%CREATEFIGURE3 Creates a 3D figure.
%
% f = CREATEFIGURE3 opens a figure with default settings.
%
% f = CREATEFIGURE3(cmap) creates a figure with the given colormap cmap.
%
% f = CREATEFIGURE3(cmap, xmin, xmax, ymin, ymax, zmin, zmax) sets axis
% limits accordingly.
%
% Note that CREATEFIGURE3 uses zbuffer as renderer!
    f = figure('Renderer', 'zbuffer');
    axis square;
    daspect([1, 1, 1]);
    hold on;
    view([60, 30]);
    % Set colormap.
    if(nargin > 0)
       colormap(cmap);
    end
    if(nargin == 7)
        set(gca, 'ZLim', [zmin zmax]);
        set(gca, 'XLim', [xmin xmax]);
        set(gca, 'YLim', [ymin ymax]);
    end
end