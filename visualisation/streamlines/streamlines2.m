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
function lh = streamlines2(P, v, S, h, maxit, cmap)
%STREAMLINES2 Computes integral curves of a static velocity field given on 
%a hemisphere and plots them in the plane. The curves are coloured 
%according to cmap (optional).
%
%   lh = STREAMLINES2(P, v, S, h, maxit, cmap) takes points P and a vector 
%   field v on a hemisphere and plots colored integral curves for an 
%   artificial time. The number of steps is given by maxit. The a cell 
%   array lh contains handles to the line segments lh{t} for each time t.
%   S are starting points of the integral curves in R^2. cmap is used to
%   interpolate colors.
%
%   Note that line segments are plotted with increasing time so that they
%   may overlap previously drawn segments.

assert(size(P, 2) == 3);
assert(size(v, 2) == 3);
assert(size(P, 1) == size(v, 1));
assert(h > 0);
assert(maxit > 0);

% Pull back vector field to plane.
w = pullback(P, v);

% Compute integral curves.
verts = integrateflow2(P(:, 1:2), w, S, h, maxit);

% Interpolate colours.
c = 1:maxit;
cmap = colormap(cmap);
y = linspace(max(c),min(c),size(cmap,1));
cm = spline(y, cmap', c);

for l=1:maxit-1
    X = [];
    Y = [];
    for k=1:length(verts)
        v = verts{k};
        if(l >= size(v, 1))
            continue;
        end        
        X = [X; v(l, 1), v(l+1, 1)];
        Y = [Y; v(l, 2), v(l+1, 2)];
    end
    if(isempty(X) || isempty(Y))
        break;
    end
    lh{l} = line(X', Y', 'color', cm(:, l), 'LineWidth', 2);
end

end