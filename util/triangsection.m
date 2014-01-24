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
function [F, V, Fidx, Vidx] = triangsection(F, V, lim)
%TRIANGSECTION Returns a section of a triangulation.
%
%   [F, V, Fidx, Vidx] = TRIANGSECTION(F, V, lim) returns faces F and 
%   vertices V of a triangulation F, V such that all vertices are in lim.
%   Fidx and Vidx are the indices of the faces and vertices, respectively,
%   which are in the area specified by lim.
%
%   Note that lim is a vector [xmin, xmax, ymin, ymax, zmin, zmax].

assert(isvector(lim));
assert(length(lim) == 6);

% Get faces on the northern hemisphere.
xmin = V(F(:, 1), 1) >= lim(1) & V(F(:, 2), 1) >= lim(1) & V(F(:, 3), 1) >= lim(1);
xmax = V(F(:, 1), 1) <= lim(2) & V(F(:, 2), 1) <= lim(2) & V(F(:, 3), 1) <= lim(2);
ymin = V(F(:, 1), 2) >= lim(3) & V(F(:, 2), 2) >= lim(3) & V(F(:, 3), 2) >= lim(3);
ymax = V(F(:, 1), 2) <= lim(4) & V(F(:, 2), 2) <= lim(4) & V(F(:, 3), 2) <= lim(4);
zmin = V(F(:, 1), 3) >= lim(5) & V(F(:, 2), 3) >= lim(5) & V(F(:, 3), 3) >= lim(5);
zmax = V(F(:, 1), 3) <= lim(6) & V(F(:, 2), 3) <= lim(6) & V(F(:, 3), 3) <= lim(6);
Fidx = xmin & xmax & ymin & ymax & zmin & zmax;
F = F(Fidx, :);

% Find vertices on the northern hemisphere.
Vidx = find(V(:, 1) >= lim(1) & V(:, 1) <= lim(2) & V(:, 2) >= lim(3) & V(:, 2) <= lim(4) & V(:, 3) >= lim(5) & V(:, 3) <= lim(6));
V = V(Vidx, :);

% Compute new vertex indices.
nids(Vidx) = 1:length(Vidx);
F = nids(F);

end