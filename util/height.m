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
function H = height(F, V)
%HEIGHT Computes the heights of the triangles for a given triangulation.
%
%   H = height(F, V) takes a triangulation F, V with n faces and returns a
%   n-by-3-by-3 matrix H with height vector H(k, i, :) for vertex i of
%   triangle k.

% Get coordinates of vertices.
V1 = V(F(:, 1), :);
V2 = V(F(:, 2), :);
V3 = V(F(:, 3), :);

% Compute triangle sides.
V12 = V2 - V1;
V13 = V3 - V1;
V23 = V3 - V2;

% Compute projection of vertex to opposite side.
P1 = repmat((dot(V1 - V2, V23, 2) ./ dot(V23, V23, 2)), 1, 3) .* V23 + V2;
P2 = repmat((dot(V2 - V1, V13, 2) ./ dot(V13, V13, 2)), 1, 3) .* V13 + V1;
P3 = repmat((dot(V3 - V1, V12, 2) ./ dot(V12, V12, 2)), 1, 3) .* V12 + V1;

% Compute height vectors.
H(:, 1, :) = P1 - V1;
H(:, 2, :) = P2 - V2;
H(:, 3, :) = P3 - V3;

end